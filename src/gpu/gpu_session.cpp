// Persistent interactive session. See include/gpu/gpu_session.h.
//
// The descriptor/pipeline boilerplate mirrors render_megakernel() in
// gpu_renderer.cpp - it is factored here into create_pipeline() so both
// passes (path trace and display modulate) share one implementation.

#include <gpu/gpu_session.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace tiny_renderer::gpu {

namespace {

// Push-constant blocks, mirroring the Slang structs they feed byte-for-byte.
#pragma pack(push, 1)
struct RtPathTracePushConstants {
    uint32_t frame_index;
    uint32_t max_depth;
    uint32_t rr_depth;
    uint32_t num_lights;
    uint32_t environment_light_id;
    uint32_t accumulate;
    float radiance_clamp;
    uint32_t spp_per_frame;
    uint32_t demodulate;
};

struct SvgfTemporalPushConstants {
    uint32_t width;
    uint32_t height;
    uint32_t reset;
    float alpha_color;
    float alpha_moments;
    uint32_t clamp_history;
    float phi_depth;
    float phi_normal; // radians; the shader compares against cos(phi_normal)
    float max_hist_len;
    float clamp_growth;
    uint32_t debug_mode;
};

struct SvgfAtrousPushConstants {
    uint32_t width;
    uint32_t height;
    uint32_t stride;      // 2^iteration
    float phi_color;
    float phi_normal;     // exponent
    float phi_depth;
    uint32_t iteration;
    uint32_t debug_mode;
};

struct SvgfModulatePushConstants {
    float exposure;
    uint32_t tonemap;
    uint32_t use_srgb;
    float gamma;
    uint32_t modulate;
    uint32_t debug_view;
    // The output is a flat structured buffer, so it cannot report its own
    // dimensions - the shader needs these to index it (see svgf_modulate.slang).
    uint32_t width;
    uint32_t height;
    float variance_scale;
    uint32_t use_raw_illum;
};
#pragma pack(pop)

constexpr uint32_t GROUP_SIZE = 8; // must match [numthreads(8, 8, 1)] in every kernel

// Inverts a row-major 4x4. The temporal pass reprojects the G-Buffer's world
// positions through the previous frame's camera, which needs world_to_camera -
// the inverse of the camera_to_world the caller supplies.
std::array<float, 16> invert_transform(const float m[16]) {
    std::array<float, 16> out{};

    Matrix4f M;
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c)
            M(r, c) = m[r * 4 + c];

    // Gauss-Jordan elimination in DOUBLE precision.
    //
    // Two reasons not to use TMatrix::inverse() here:
    //
    // 1. Its `valid` flag is not trustworthy - it reports false for perfectly
    //    well-conditioned projection matrices (it appears to flag the exact
    //    zeros intrinsic to a perspective matrix's sparsity). Acting on that
    //    flag would silently fall back to the input matrix, which for the
    //    camera_to_sample matrix means reprojecting through it BACKWARDS:
    //    history samples land ~1/near_clip pixels away, every validity test
    //    fails, and the denoiser quietly degrades to a single-sample estimate.
    //
    // 2. Even when correct it is not accurate ENOUGH. A projection matrix with
    //    near=1e-4 and far=1e4 has a condition number around 1e8, so a float
    //    inversion leaves ~0.1 relative error - visible as reprojection landing
    //    a fraction of a pixel off, which is enough to bias the validity tests.
    //    Doubles bring that down to ~1e-8.
    double a[4][4], inv[4][4];
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
            a[r][c]   = static_cast<double>(M(r, c));
            inv[r][c] = (r == c) ? 1.0 : 0.0;
        }
    }

    for (int col = 0; col < 4; ++col) {
        // Partial pivoting.
        int pivot = col;
        for (int r = col + 1; r < 4; ++r)
            if (std::fabs(a[r][col]) > std::fabs(a[pivot][col])) pivot = r;

        if (std::fabs(a[pivot][col]) < 1e-18) {
            std::cerr << "[session] warning: singular camera matrix; temporal reprojection disabled for "
                         "this frame" << std::endl;
            std::copy_n(m, 16, out.begin());
            return out;
        }
        if (pivot != col) {
            for (int c = 0; c < 4; ++c) {
                std::swap(a[col][c], a[pivot][c]);
                std::swap(inv[col][c], inv[pivot][c]);
            }
        }

        double d = a[col][col];
        for (int c = 0; c < 4; ++c) {
            a[col][c]   /= d;
            inv[col][c] /= d;
        }
        for (int r = 0; r < 4; ++r) {
            if (r == col) continue;
            double f = a[r][col];
            if (f == 0.0) continue;
            for (int c = 0; c < 4; ++c) {
                a[r][c]   -= f * a[col][c];
                inv[r][c] -= f * inv[col][c];
            }
        }
    }

    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c)
            out[r * 4 + c] = static_cast<float>(inv[r][c]);
    return out;
}

} // namespace

InteractiveSession::InteractiveSession(GPUScene gs, uint32_t width, uint32_t height)
    : m_gs(std::move(gs)), m_width(width), m_height(height), m_ctx(false) {

    if (m_width == 0 || m_height == 0)
        throw std::runtime_error("InteractiveSession: render resolution must be non-zero");

    // Fail loudly on a scene whose Accel never built a BVH. Without this the
    // session would happily trace rays against an empty acceleration structure
    // and display a black image with no indication of why - the same guard
    // render_dispatch() in gpu_renderer.cpp applies to the offline paths.
    if (m_gs.bvh_nodes.empty())
        throw std::runtime_error("InteractiveSession: scene's Accel implementation did not export a flat BVH "
                                 "(only <accelerate type=\"bvh\"> does); did the caller forget "
                                 "Scene::construct()?");

    // The shaders bounds-check thread ids against camera.width/height, so the
    // uploaded camera must describe the render resolution, not whatever the
    // XML happened to specify (resolution scaling is implemented here).
    m_gs.camera.width  = m_width;
    m_gs.camera.height = m_height;

    m_packed        = pack_gpu_scene(m_gs);
    m_scene_buffers = upload_scene_buffers(m_ctx, m_gs, m_packed);

    // The shaders allocate their traversal stack statically; deeper than that
    // and traversal overwrites memory instead of failing. Say so loudly rather
    // than shipping subtly missing geometry.
    if (const int depth = compute_bvh_depth(m_gs); depth > 32)
        throw std::runtime_error("InteractiveSession: BVH depth " + std::to_string(depth) +
                                 " exceeds the shaders' traversal stack (32); raise M_BVH_STACK_DEPTH in "
                                 "shaders/common/scene_common.slang");

    const VkDeviceSize pixels = static_cast<VkDeviceSize>(m_width) * m_height;
    m_accum_buf   = create_empty_buffer(m_ctx, pixels * sizeof(float) * 4, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    m_display_buf = create_empty_buffer(m_ctx, pixels * sizeof(uint32_t), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);

    // ---- Denoiser surfaces (see include/gpu/gpu_image.h for why images) ----
    m_gbuffer_albedo       = create_storage_image(m_ctx, m_width, m_height, VK_FORMAT_R16G16B16A16_SFLOAT);
    m_gbuffer_normal_depth = create_storage_image(m_ctx, m_width, m_height, VK_FORMAT_R32G32B32A32_SFLOAT);
    m_gbuffer_pos_id       = create_storage_image(m_ctx, m_width, m_height, VK_FORMAT_R32G32B32A32_SFLOAT);
    m_illum_cur            = create_storage_image(m_ctx, m_width, m_height, VK_FORMAT_R16G16B16A16_SFLOAT);
    m_moments_cur          = create_storage_image(m_ctx, m_width, m_height, VK_FORMAT_R32G32_SFLOAT);
    for (int i = 0; i < 2; ++i) {
        m_illum_hist[i]   = create_storage_image(m_ctx, m_width, m_height, VK_FORMAT_R16G16B16A16_SFLOAT);
        m_moments_hist[i] = create_storage_image(m_ctx, m_width, m_height, VK_FORMAT_R32G32_SFLOAT);
        m_hist_len[i]     = create_storage_image(m_ctx, m_width, m_height, VK_FORMAT_R32_SFLOAT);
        m_atrous_illum[i] = create_storage_image(m_ctx, m_width, m_height, VK_FORMAT_R16G16B16A16_SFLOAT);
    }

    // Previous-frame camera for the temporal pass: world_to_camera (16 floats)
    // then camera_to_sample (16 floats), plus width/height.
    m_prev_camera_buf = create_empty_buffer(m_ctx, sizeof(float) * 34, VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT);

    // ---- Timestamp queries ----
    // Entirely optional, and deliberately DISABLED on MoltenVK.
    //
    // MoltenVK does report timestampPeriod and accepts the queries, but it does
    // not resolve them at compute-dispatch granularity: the values come out
    // describing a whole command buffer rather than the interval between two
    // writes. Measured here, that produced nonsense - the four a-trous
    // iterations reported 2.3 microseconds (impossible for ~92M texture taps),
    // and one scene reported 115 ms of path tracing while actually running at
    // 54 FPS. Wrong numbers are worse than none, so it is turned off instead.
    //
    // On a native Vulkan implementation these are accurate and the breakdown
    // will populate; the code path is unchanged, only gated.
    {
        VkPhysicalDeviceProperties props{};
        vkGetPhysicalDeviceProperties(m_ctx.physical_device, &props);

        const bool moltenvk = (props.vendorID == 0x106b); // Apple; MoltenVK reports the GPU's ID
        m_timestamp_period_ns = moltenvk ? 0.0f : props.limits.timestampPeriod;

        if (m_timestamp_period_ns > 0.0f) {
            VkQueryPoolCreateInfo q{};
            q.sType      = VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO;
            q.queryType  = VK_QUERY_TYPE_TIMESTAMP;
            q.queryCount = kTimestampCount;
            if (vkCreateQueryPool(m_ctx.device, &q, nullptr, &m_query_pool) != VK_SUCCESS) {
                m_query_pool = VK_NULL_HANDLE;
                m_timestamp_period_ns = 0.0f;
            }
        }
        m_timings.available = (m_query_pool != VK_NULL_HANDLE);
    }

    // One pool serves every pass's descriptor sets.
    VkDescriptorPoolSize pool_sizes[3] = {};
    pool_sizes[0].type                 = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    pool_sizes[0].descriptorCount      = 64;
    pool_sizes[1].type                 = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
    pool_sizes[1].descriptorCount      = 8;
    pool_sizes[2].type                 = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    pool_sizes[2].descriptorCount      = 96;

    VkDescriptorPoolCreateInfo pool_info{};
    pool_info.sType         = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
    pool_info.maxSets       = 16;
    pool_info.poolSizeCount = 3;
    pool_info.pPoolSizes    = pool_sizes;
    VK_CHECK(vkCreateDescriptorPool(m_ctx.device, &pool_info, nullptr, &m_descriptor_pool));

    // Bindings 1..15 are the shared flattened scene (see
    // common/scene_common.slang); every pass that needs the scene uses this
    // exact list, in order.
    const std::pair<VkDescriptorType, VkBuffer> scene[] = {
        { VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, m_scene_buffers.camera_buf.buffer },
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.positions_buf.buffer },
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.normals_buf.buffer },
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.uvs_buf.buffer },
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.indices_buf.buffer },
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.bvh_nodes_buf.buffer },
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.bvh_prims_buf.buffer },
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.tri_mat_buf.buffer },
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.tri_light_buf.buffer },
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.materials_buf.buffer },
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.textures_buf.buffer },
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.texture_pixels_buf.buffer },
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.material_lut_buf.buffer },
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.lights_buf.buffer },
        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.light_cdf_buf.buffer },
    };

    // ---- Path trace pass ----
    // binding 0 = accumulation buffer, 1..15 = scene, 16..20 = G-Buffer +
    // demodulated illumination + moments.
    {
        std::vector<VkDescriptorSetLayoutBinding> bindings;
        std::vector<DescriptorSlot> slots;

        bindings.push_back({ 0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr });
        slots.push_back({ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_accum_buf.buffer, VK_NULL_HANDLE });
        for (uint32_t i = 0; i < 15; ++i) {
            bindings.push_back({ i + 1, scene[i].first, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr });
            slots.push_back({ scene[i].first, scene[i].second, VK_NULL_HANDLE });
        }

        // Bindings 16..20, in rt_path_trace.slang's declared order.
        const VkImageView gbuffer_views[] = {
            m_gbuffer_albedo.view,       // 16
            m_gbuffer_normal_depth.view, // 17
            m_gbuffer_pos_id.view,       // 18
            m_illum_cur.view,            // 19
            m_moments_cur.view,          // 20
        };
        for (uint32_t i = 0; i < 5; ++i) {
            bindings.push_back({ i + 16, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr });
            slots.push_back({ VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, VK_NULL_HANDLE, gbuffer_views[i] });
        }
        // binding 32: packed triangles, declared by scene_common.
        bindings.push_back({ 32, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr });
        slots.push_back({ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_scene_buffers.packed_triangles_buf.buffer,
                          VK_NULL_HANDLE });

        create_pipeline(m_path_trace, "rt_path_trace.spv", sizeof(RtPathTracePushConstants), bindings, { slots });
    }

    // ---- Spatial (a-trous) pass ----
    //
    // Six descriptor sets: three per history pair. Which pair holds this frame's
    // temporal output alternates every frame, and the variance that guides the
    // filter comes from that same pair, so the whole iteration chain has to be
    // bound twice - once for each possible source.
    //
    // For a given source pair `h`:
    //   h*3+0:  out scratch[0] <- in history[h]     (first iteration)
    //   h*3+1:  out scratch[1] <- in scratch[0]     (odd iterations)
    //   h*3+2:  out scratch[0] <- in scratch[1]     (even iterations after the first)
    {
        std::vector<VkDescriptorSetLayoutBinding> bindings;
        const VkDescriptorType img = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
        bindings.push_back({ 0, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // illum_out
        bindings.push_back({ 1, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // illum_in
        bindings.push_back({ 2, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // moments (variance)
        bindings.push_back({ 3, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // normal_depth
        bindings.push_back({ 4, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // pos_id

        std::vector<std::vector<DescriptorSlot>> slot_sets(6);
        for (int h = 0; h < 2; ++h) {
            const VkImageView in_first = m_illum_hist[h].view;      // temporal output
            const VkImageView moments  = m_moments_hist[h].view;    // its variance
            const struct { VkImageView out_view, in_view; } steps[3] = {
                { m_atrous_illum[0].view, in_first },          // first
                { m_atrous_illum[1].view, m_atrous_illum[0].view }, // odd
                { m_atrous_illum[0].view, m_atrous_illum[1].view }, // even
            };
            for (int s = 0; s < 3; ++s) {
                slot_sets[h * 3 + s] = {
                    { img, VK_NULL_HANDLE, steps[s].out_view },
                    { img, VK_NULL_HANDLE, steps[s].in_view },
                    { img, VK_NULL_HANDLE, moments },
                    { img, VK_NULL_HANDLE, m_gbuffer_normal_depth.view },
                    { img, VK_NULL_HANDLE, m_gbuffer_pos_id.view },
                };
            }
        }
        create_pipeline(m_atrous, "svgf_atrous.spv", sizeof(SvgfAtrousPushConstants), bindings, slot_sets);
    }

    // ---- Temporal accumulation pass ----
    // Two descriptor sets: set 0 reads history pair 0 / writes pair 1, set 1 does
    // the reverse. The host picks one by frame parity.
    {
        std::vector<VkDescriptorSetLayoutBinding> bindings;
        const VkDescriptorType img = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;

        bindings.push_back({ 0, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // illum_out
        bindings.push_back({ 1, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // moments_out
        bindings.push_back({ 2, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // hist_len_out
        bindings.push_back({ 3, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // illum_in
        bindings.push_back({ 4, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // moments_in
        bindings.push_back({ 5, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // hist_len_in
        bindings.push_back({ 6, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // illum_cur
        bindings.push_back({ 7, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // moments_cur
        bindings.push_back({ 8, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // pos_id
        bindings.push_back({ 9, img, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr }); // normal_depth
        bindings.push_back({ 10, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT,
                             nullptr }); // prev camera

        std::vector<std::vector<DescriptorSlot>> slot_sets(2);
        for (int dst = 0; dst < 2; ++dst) {
            const int src = 1 - dst;
            slot_sets[dst] = {
                { img, VK_NULL_HANDLE, m_illum_hist[dst].view },   // 0 out
                { img, VK_NULL_HANDLE, m_moments_hist[dst].view }, // 1
                { img, VK_NULL_HANDLE, m_hist_len[dst].view },     // 2
                { img, VK_NULL_HANDLE, m_illum_hist[src].view },   // 3 in
                { img, VK_NULL_HANDLE, m_moments_hist[src].view }, // 4
                { img, VK_NULL_HANDLE, m_hist_len[src].view },     // 5
                { img, VK_NULL_HANDLE, m_illum_cur.view },         // 6 cur
                { img, VK_NULL_HANDLE, m_moments_cur.view },       // 7
                { img, VK_NULL_HANDLE, m_gbuffer_pos_id.view },    // 8
                { img, VK_NULL_HANDLE, m_gbuffer_normal_depth.view }, // 9
                { VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, m_prev_camera_buf.buffer, VK_NULL_HANDLE }, // 10
            };
        }

        create_pipeline(m_temporal, "svgf_temporal.spv", sizeof(SvgfTemporalPushConstants), bindings, slot_sets);
    }

    // ---- Display pass ----
    // Four sets: it must read whatever the last pass produced, and that can be
    // either history pair (spatial pass off) or either scratch image (spatial
    // pass on) - four combinations, all alternating as the frame parity flips.
    //   0: history[0]    1: history[1]    2: scratch[0]    3: scratch[1]
    {
        std::vector<VkDescriptorSetLayoutBinding> bindings;
        bindings.push_back({ 0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr });
        for (uint32_t i = 1; i <= 7; ++i)
            bindings.push_back({ i, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr });

        std::vector<std::vector<DescriptorSlot>> slot_sets(4);
        // Sets 0/1 pair the illumination with history[0]/history[1] (and that
        // pair's moments + history length); sets 2/3 for the scratch images.
        const VkImageView illum_src[4] = { m_illum_hist[0].view, m_illum_hist[1].view, m_atrous_illum[0].view,
                                           m_atrous_illum[1].view };
        const VkImageView mom_src[4]   = { m_moments_hist[0].view, m_moments_hist[1].view, m_moments_hist[0].view,
                                           m_moments_hist[1].view };
        const VkImageView len_src[4]   = { m_hist_len[0].view, m_hist_len[1].view, m_hist_len[0].view,
                                           m_hist_len[1].view };
        for (int k = 0; k < 4; ++k) {
            slot_sets[k] = {
                { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, m_display_buf.buffer, VK_NULL_HANDLE },
                { VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, VK_NULL_HANDLE, illum_src[k] },           // 1 illum
                { VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, VK_NULL_HANDLE, m_gbuffer_albedo.view },  // 2
                { VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, VK_NULL_HANDLE, m_gbuffer_normal_depth.view }, // 3
                { VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, VK_NULL_HANDLE, mom_src[k] },             // 4
                { VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, VK_NULL_HANDLE, len_src[k] },             // 5
                { VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, VK_NULL_HANDLE, m_gbuffer_pos_id.view },  // 6
                { VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, VK_NULL_HANDLE, m_illum_cur.view },       // 7 raw
            };
        }

        create_pipeline(m_modulate, "svgf_modulate.spv", sizeof(SvgfModulatePushConstants), bindings, slot_sets);
    }

    std::cout << "[session] Interactive session ready on '" << m_ctx.device_name << "' (" << m_width << "x"
              << m_height << ")\n";
}

InteractiveSession::~InteractiveSession() {
    // Anything that touches the GPU must be idle first: the destructor may run
    // while a submitted frame is still in flight.
    vkDeviceWaitIdle(m_ctx.device);

    destroy_pipeline(m_modulate);
    destroy_pipeline(m_atrous);
    destroy_pipeline(m_temporal);
    destroy_pipeline(m_path_trace);
    if (m_query_pool)
        vkDestroyQueryPool(m_ctx.device, m_query_pool, nullptr);
    if (m_descriptor_pool)
        vkDestroyDescriptorPool(m_ctx.device, m_descriptor_pool, nullptr);

    destroy_buffer(m_ctx, m_display_buf);
    destroy_buffer(m_ctx, m_accum_buf);
    destroy_buffer(m_ctx, m_prev_camera_buf);

    for (int i = 0; i < 2; ++i) {
        destroy_image(m_ctx, m_hist_len[i]);
        destroy_image(m_ctx, m_moments_hist[i]);
        destroy_image(m_ctx, m_illum_hist[i]);
        destroy_image(m_ctx, m_atrous_illum[i]);
    }
    destroy_image(m_ctx, m_moments_cur);
    destroy_image(m_ctx, m_illum_cur);
    destroy_image(m_ctx, m_gbuffer_pos_id);
    destroy_image(m_ctx, m_gbuffer_normal_depth);
    destroy_image(m_ctx, m_gbuffer_albedo);

    m_scene_buffers.destroy(m_ctx);
}

void InteractiveSession::create_pipeline(Pipeline &out, const std::string &spv_name,
                                         size_t push_constants_size,
                                         const std::vector<VkDescriptorSetLayoutBinding> &bindings,
                                         const std::vector<std::vector<DescriptorSlot>> &slot_sets) {
    VkDescriptorSetLayoutCreateInfo layout_info{};
    layout_info.sType        = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
    layout_info.bindingCount = static_cast<uint32_t>(bindings.size());
    layout_info.pBindings    = bindings.data();
    VK_CHECK(vkCreateDescriptorSetLayout(m_ctx.device, &layout_info, nullptr, &out.set_layout));

    VkPushConstantRange push_range{};
    push_range.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    push_range.offset     = 0;
    push_range.size       = static_cast<uint32_t>(push_constants_size);

    VkPipelineLayoutCreateInfo pipeline_layout_info{};
    pipeline_layout_info.sType                  = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
    pipeline_layout_info.setLayoutCount         = 1;
    pipeline_layout_info.pSetLayouts            = &out.set_layout;
    pipeline_layout_info.pushConstantRangeCount = 1;
    pipeline_layout_info.pPushConstantRanges    = &push_range;
    VK_CHECK(vkCreatePipelineLayout(m_ctx.device, &pipeline_layout_info, nullptr, &out.pipeline_layout));

    out.module = m_ctx.load_shader_module(std::string(M_GPU_SHADER_BINARY_DIR) + "/" + spv_name);

    VkPipelineShaderStageCreateInfo stage_info{};
    stage_info.sType  = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stage_info.stage  = VK_SHADER_STAGE_COMPUTE_BIT;
    stage_info.module = out.module;
    stage_info.pName  = "main";

    VkComputePipelineCreateInfo pipeline_info{};
    pipeline_info.sType  = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
    pipeline_info.stage  = stage_info;
    pipeline_info.layout = out.pipeline_layout;
    VK_CHECK(vkCreateComputePipelines(m_ctx.device, VK_NULL_HANDLE, 1, &pipeline_info, nullptr, &out.pipeline));

    // One descriptor set per entry in `slot_sets`, all sharing this layout.
    out.sets.assign(slot_sets.size(), VK_NULL_HANDLE);
    std::vector<VkDescriptorSetLayout> layouts(slot_sets.size(), out.set_layout);

    VkDescriptorSetAllocateInfo set_alloc_info{};
    set_alloc_info.sType              = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
    set_alloc_info.descriptorPool     = m_descriptor_pool;
    set_alloc_info.descriptorSetCount = static_cast<uint32_t>(layouts.size());
    set_alloc_info.pSetLayouts        = layouts.data();
    VK_CHECK(vkAllocateDescriptorSets(m_ctx.device, &set_alloc_info, out.sets.data()));

    // The three arrays are sized up front and never resized: the write structs
    // below store RAW POINTERS into buf_infos/img_infos, so any reallocation
    // after those pointers are taken would leave them dangling.
    size_t total_slots = 0;
    for (const auto &s : slot_sets) total_slots += s.size();
    std::vector<VkDescriptorBufferInfo> buf_infos(total_slots);
    std::vector<VkDescriptorImageInfo> img_infos(total_slots);
    std::vector<VkWriteDescriptorSet> writes(total_slots);

    size_t k = 0;
    for (size_t set_i = 0; set_i < slot_sets.size(); ++set_i) {
        for (size_t i = 0; i < slot_sets[set_i].size(); ++i, ++k) {
            const DescriptorSlot &slot = slot_sets[set_i][i];
            const bool is_image = (slot.type == VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);

            if (is_image) {
                img_infos[k]             = {};
                img_infos[k].imageView   = slot.image_view;
                // Images are kept permanently in GENERAL (see create_storage_image),
                // so no layout transition is needed between passes.
                img_infos[k].imageLayout = VK_IMAGE_LAYOUT_GENERAL;
            } else {
                buf_infos[k]        = {};
                buf_infos[k].buffer = slot.buffer;
                buf_infos[k].offset = 0;
                buf_infos[k].range  = VK_WHOLE_SIZE;
            }

            writes[k]                 = {};
            writes[k].sType           = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            writes[k].dstSet          = out.sets[set_i];
            // Take the binding number from the layout declaration, NOT from the
            // slot index. They coincide for the denoiser passes (which use a
            // contiguous 0..N), but the path-trace kernel has a gap: it declares
            // the scene on 1..15, its own G-Buffer on 16..20, and then
            // scene_common's packed triangles on 32. Indexing would bind that
            // last one to slot 21 and leave 32 unset - which reads as an
            // all-zero vertex array and renders a completely black image, with
            // nothing anywhere reporting a missing descriptor.
            writes[k].dstBinding      = bindings[i].binding;
            writes[k].descriptorCount = 1;
            writes[k].descriptorType  = slot.type;
            writes[k].pBufferInfo     = is_image ? nullptr : &buf_infos[k];
            writes[k].pImageInfo      = is_image ? &img_infos[k] : nullptr;
        }
    }
    vkUpdateDescriptorSets(m_ctx.device, static_cast<uint32_t>(writes.size()), writes.data(), 0, nullptr);
}

void InteractiveSession::destroy_pipeline(Pipeline &p) {
    if (p.pipeline)
        vkDestroyPipeline(m_ctx.device, p.pipeline, nullptr);
    if (p.pipeline_layout)
        vkDestroyPipelineLayout(m_ctx.device, p.pipeline_layout, nullptr);
    if (p.set_layout)
        vkDestroyDescriptorSetLayout(m_ctx.device, p.set_layout, nullptr);
    if (p.module)
        vkDestroyShaderModule(m_ctx.device, p.module, nullptr);
    p = Pipeline{};
}

// IEEE 754 binary16 -> binary32. The illumination history is stored as
// RGBA16F (half the bandwidth of RGBA32F, and plenty of range for radiance),
// so reading it back to the host needs an explicit widening.
static float half_to_float(uint16_t h) {
    const uint32_t sign = static_cast<uint32_t>(h & 0x8000u) << 16;
    const uint32_t exp  = (h & 0x7C00u) >> 10;
    const uint32_t mant = h & 0x03FFu;

    uint32_t bits;
    if (exp == 0u) {
        if (mant == 0u) {
            bits = sign; // +/-0
        } else {
            // Subnormal half: renormalise into a normal float.
            uint32_t e = 127 - 15 + 1;
            uint32_t m = mant;
            while (!(m & 0x0400u)) {
                m <<= 1;
                --e;
            }
            m &= 0x03FFu;
            bits = sign | (e << 23) | (m << 13);
        }
    } else if (exp == 31u) {
        bits = sign | 0x7F800000u | (mant << 13); // inf / nan
    } else {
        bits = sign | ((exp + 112u) << 23) | (mant << 13);
    }

    float f;
    std::memcpy(&f, &bits, sizeof(f));
    return f;
}

void InteractiveSession::read_pass_timings() {
    m_timings = PassTimings{};
    m_timings.available = (m_query_pool != VK_NULL_HANDLE);
    if (!m_query_pool) return;

    uint64_t ts[kTimestampCount] = {};
    const VkResult r = vkGetQueryPoolResults(m_ctx.device, m_query_pool, 0, kTimestampCount,
                                             sizeof(ts), ts, sizeof(uint64_t),
                                             VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT);
    if (r != VK_SUCCESS) {
        // Not worth failing a frame over: the timings are diagnostic.
        m_timings.available = false;
        return;
    }

    // Timestamp values are ticks; convert with the device's period.
    auto ms = [&](uint32_t a, uint32_t b) {
        const uint64_t delta = (ts[b] > ts[a]) ? (ts[b] - ts[a]) : 0;
        return static_cast<float>(static_cast<double>(delta) * m_timestamp_period_ns / 1e6);
    };

    m_timings.path_trace_ms = ms(0, 1);
    m_timings.temporal_ms   = ms(1, 2);
    m_timings.atrous_ms     = ms(2, 3);
    m_timings.modulate_ms   = ms(3, 4);
    m_timings.total_ms      = ms(0, 4);
}

std::vector<float> InteractiveSession::read_filtered_radiance() {
    const size_t pixels = static_cast<size_t>(m_width) * m_height;
    std::vector<float> out(pixels * 3, 0.0f);

    if (!m_history_valid || m_display_source == nullptr) {
        // Nothing has been written yet, so a readback would return uninitialised
        // memory. Callers use this to compare against an offline render, where
        // silently returning zeros would look like a perfect black image - hence
        // the hard failure.
        throw std::runtime_error("InteractiveSession::read_filtered_radiance: no history yet "
                                 "(the temporal pass has not run since the last reset)");
    }

    // RGBA16F => 8 bytes per pixel.
    GpuBuffer staging = create_empty_buffer(m_ctx, static_cast<VkDeviceSize>(pixels) * 8,
                                            VK_BUFFER_USAGE_TRANSFER_DST_BIT);
    // Whichever image the display pass is showing: the spatial pass's last
    // output, or the temporal output when the spatial pass is disabled.
    const GpuImage &src = *m_display_source;

    m_ctx.submit_and_wait([&](VkCommandBuffer cmd) {
        VkBufferImageCopy region{};
        region.bufferOffset      = 0;
        region.bufferRowLength   = 0;   // tightly packed
        region.bufferImageHeight = 0;
        region.imageSubresource  = { VK_IMAGE_ASPECT_COLOR_BIT, 0, 0, 1 };
        region.imageOffset       = { 0, 0, 0 };
        region.imageExtent       = { m_width, m_height, 1 };
        // The image lives permanently in GENERAL (see create_storage_image), so
        // it is already a valid transfer source - no layout transition needed.
        vkCmdCopyImageToBuffer(cmd, src.image, VK_IMAGE_LAYOUT_GENERAL, staging.buffer, 1, &region);
    });

    // The buffer is host-visible but not guaranteed coherent, so make the
    // device's writes visible before the CPU reads them.
    vmaInvalidateAllocation(m_ctx.allocator, staging.allocation, 0, VK_WHOLE_SIZE);

    const uint16_t *data = static_cast<const uint16_t *>(staging.mapped);
    for (size_t i = 0; i < pixels; ++i) {
        out[i * 3 + 0] = half_to_float(data[i * 4 + 0]);
        out[i * 3 + 1] = half_to_float(data[i * 4 + 1]);
        out[i * 3 + 2] = half_to_float(data[i * 4 + 2]);
    }

    destroy_buffer(m_ctx, staging);
    return out;
}

std::vector<float> InteractiveSession::read_linear_radiance() const {
    const size_t pixels = static_cast<size_t>(m_width) * m_height;
    std::vector<float> out(pixels * 3, 0.0f);
    const float *src = static_cast<const float *>(m_accum_buf.mapped);
    for (size_t i = 0; i < pixels; ++i) {
        float count = std::max(src[i * 4 + 3], 1.0f);
        out[i * 3 + 0] = src[i * 4 + 0] / count;
        out[i * 3 + 1] = src[i * 4 + 1] / count;
        out[i * 3 + 2] = src[i * 4 + 2] / count;
    }
    return out;
}

const uint8_t *InteractiveSession::render_frame(const RenderParams &render, const DenoiseParams &denoise,
                                                const DisplayParams &display) {
    if (!render.camera_to_world || !render.sample_to_camera)
        throw std::runtime_error("InteractiveSession::render_frame: camera matrices are required");

    // Two independent notions of "start over", and conflating them was a real
    // bug: `reset` is the caller's intent (the camera moved, so the running
    // average is stale), while `temporal_reset` additionally covers the case
    // where the ping-pong history has never been written. Folding the latter
    // into the former meant that with the temporal pass DISABLED history was
    // never valid, so every frame reset and the accumulation buffer could never
    // get past one sample.
    const bool reset          = render.reset_history;
    const bool temporal_reset = reset || !m_history_valid;

    // The temporal pass needs the PREVIOUS frame's camera to reproject world
    // positions, so it is snapshotted BEFORE the live camera is overwritten.
    if (m_history_valid && !temporal_reset) {
        float *prev = static_cast<float *>(m_prev_camera_buf.mapped);
        std::memcpy(prev + 0, m_prev_w2c.data(), 16 * sizeof(float));
        std::memcpy(prev + 16, m_prev_c2s.data(), 16 * sizeof(float));
        uint32_t *dims = reinterpret_cast<uint32_t *>(prev + 32);
        dims[0] = m_width;
        dims[1] = m_height;
    }

    // Upload this frame's camera. camera_buf is host-visible and persistently
    // mapped, so this is a straight memcpy - no staging buffer, and no need to
    // rebuild any descriptor. Only the two matrices are rewritten;
    // width/height/near/far (the 4 scalars that follow them in GPUCamera) were
    // baked in at construction and must survive frame updates.
    float *cam = static_cast<float *>(m_scene_buffers.camera_buf.mapped);
    std::memcpy(cam + 0, render.camera_to_world, 16 * sizeof(float));
    std::memcpy(cam + 16, render.sample_to_camera, 16 * sizeof(float));

    m_accumulated_spp = reset ? 0 : m_accumulated_spp;

    RtPathTracePushConstants rtp{};
    rtp.frame_index          = m_frame_index;
    rtp.max_depth            = render.max_depth;
    rtp.rr_depth             = render.rr_depth;
    rtp.num_lights           = static_cast<uint32_t>(m_gs.lights.size());
    rtp.environment_light_id = m_gs.environment_light_id;
    rtp.accumulate           = (reset || m_accumulated_spp == 0) ? 0u : 1u;
    rtp.radiance_clamp       = render.radiance_clamp;
    rtp.spp_per_frame        = render.spp_per_frame < 1u ? 1u : render.spp_per_frame;
    rtp.demodulate           = render.demodulate ? 1u : 0u;

    // The temporal pass READS the pair that currently holds valid history and
    // WRITES the other one; the pair it wrote then becomes the valid one.
    //
    // Getting this backwards is silent and looks exactly like "reprojection is
    // broken": the pass reads a pair that has never been written, every validity
    // test fails against garbage, and the filter falls back to a single-sample
    // estimate every frame - so the image never converges and no error is
    // reported anywhere.
    const uint32_t src = m_history_index;      // holds this frame's input history
    const uint32_t dst = 1u - m_history_index; // receives this frame's output

    // A fixed blend floor is what stops the filter from responding to a moving
    // view, but it also caps how clean a STATIC view can ever get: with
    // alpha_color = 0.05 the EMA freezes after ~20 frames and the image stops
    // improving no matter how long you wait.
    //
    // So the floor decays the longer the view stays still, handing control back
    // to the 1/(n+1) term - which is exactly a running mean, so the image keeps
    // converging towards the offline renderer's quality. This is what makes
    // "stop moving and the picture refines itself" work.
    const float floor_scale = std::exp(-static_cast<float>(m_still_frames) / 60.0f);

    SvgfTemporalPushConstants tmp{};
    tmp.width         = m_width;
    tmp.height        = m_height;
    tmp.reset         = temporal_reset ? 1u : 0u;
    tmp.alpha_color   = denoise.alpha_color * floor_scale;
    tmp.alpha_moments = denoise.alpha_moments * floor_scale;
    tmp.clamp_history = denoise.clamp_history ? 1u : 0u;
    tmp.phi_depth     = denoise.phi_depth;
    tmp.phi_normal    = denoise.phi_normal;
    tmp.max_hist_len  = denoise.max_hist_len;
    tmp.clamp_growth  = denoise.clamp_growth;
    tmp.debug_mode    = denoise.debug_mode ? 1u : 0u;

    SvgfModulatePushConstants mod{};
    mod.exposure        = display.exposure;
    mod.tonemap         = display.tonemap;
    mod.use_srgb        = display.use_srgb ? 1u : 0u;
    mod.gamma           = display.gamma;
    mod.modulate        = render.demodulate ? 1u : 0u;
    mod.debug_view      = display.debug_view;
    mod.width           = m_width;
    mod.height          = m_height;
    mod.variance_scale  = display.variance_scale;
    // With the temporal pass off nothing writes the history pair, so the
    // display pass must read this frame's raw sample instead (binding 7).
    mod.use_raw_illum   = denoise.enabled ? 0u : 1u;

    const uint32_t groups_x = (m_width + GROUP_SIZE - 1) / GROUP_SIZE;
    const uint32_t groups_y = (m_height + GROUP_SIZE - 1) / GROUP_SIZE;

    m_ctx.submit_and_wait([&](VkCommandBuffer cmd) {
        // Timestamp 0: before anything runs.
        if (m_query_pool) {
            vkCmdResetQueryPool(cmd, m_query_pool, 0, kTimestampCount);
            vkCmdWriteTimestamp(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, m_query_pool, 0);
        }

        // ---- Pass 1: path trace + G-Buffer ----
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, m_path_trace.pipeline);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, m_path_trace.pipeline_layout, 0, 1,
                                &m_path_trace.sets[0], 0, nullptr);
        vkCmdPushConstants(cmd, m_path_trace.pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(rtp), &rtp);
        vkCmdDispatch(cmd, groups_x, groups_y, 1);

        // Consecutive dispatches in one command buffer may be reordered or
        // overlapped by the hardware, so every producer/consumer pair needs an
        // explicit barrier even though there is no submission boundary between
        // them: an execution barrier (dispatch N finishes before N+1 starts)
        // plus a memory barrier (N's writes are visible to N+1's reads).
        auto barrier_between_passes = [&]() {
            VkMemoryBarrier barrier{};
            barrier.sType         = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
            barrier.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
            barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
            vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                                 0, 1, &barrier, 0, nullptr, 0, nullptr);
        };

        if (m_query_pool)
            vkCmdWriteTimestamp(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, m_query_pool, 1);

        // ---- Pass 2: temporal accumulation ----
        if (denoise.enabled) {
            barrier_between_passes();

            // slot_sets[dst] reads pair `src` and writes pair `dst`.
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, m_temporal.pipeline);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, m_temporal.pipeline_layout, 0, 1,
                                    &m_temporal.sets[dst], 0, nullptr);
            vkCmdPushConstants(cmd, m_temporal.pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(tmp), &tmp);
            vkCmdDispatch(cmd, groups_x, groups_y, 1);
        }

        if (m_query_pool)
            vkCmdWriteTimestamp(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, m_query_pool, 2);

        // ---- Pass 3: spatial (a-trous) filtering ----
        // Consecutive iterations are read-after-write on the scratch pair, so
        // each one needs its own barrier; unlike the passes above they cannot be
        // batched behind a single one.
        const uint32_t iterations = denoise.enabled ? std::min(denoise.atrous_iterations, 6u) : 0u;
        uint32_t filtered_index = 0; // which scratch image holds the final result
        for (uint32_t it = 0; it < iterations; ++it) {
            barrier_between_passes();

            // set = (source pair)*3 + {0 first, 1 odd, 2 even}
            const uint32_t kind = (it == 0) ? 0u : ((it & 1u) ? 1u : 2u);
            const uint32_t set_index = dst * 3u + kind;

            // Back the spatial filter off as the accumulation takes over.
            //
            // The filter exists to remove noise, and early on it is doing nearly
            // all the work (it cuts error several-fold at 1spp). But once the
            // temporal accumulator has converged, there is no noise left to
            // remove and the filter can only ADD blur - measured on box, running
            // it at full strength took the converged error from 8.7% to 29.5%.
            //
            // Shrinking phi_color makes the luminance term stricter, so a long
            // history keeps its detail while a fresh view still gets cleaned up.
            const float spatial_fade = std::exp(-static_cast<float>(m_still_frames) / 40.0f);

            SvgfAtrousPushConstants at{};
            at.width      = m_width;
            at.height     = m_height;
            at.stride     = 1u << it; // doubling stride = the "a-trous" (with holes) wavelet
            at.phi_color  = denoise.phi_color * spatial_fade;
            at.phi_normal = denoise.phi_normal_a;
            at.phi_depth  = denoise.phi_depth_a;
            at.iteration  = it;

            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, m_atrous.pipeline);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, m_atrous.pipeline_layout, 0, 1,
                                    &m_atrous.sets[set_index], 0, nullptr);
            vkCmdPushConstants(cmd, m_atrous.pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(at), &at);
            vkCmdDispatch(cmd, groups_x, groups_y, 1);

            // steps: first -> scratch[0], odd -> scratch[1], even -> scratch[0]
            filtered_index = (kind == 1u) ? 1u : 0u;
        }

        if (m_query_pool)
            vkCmdWriteTimestamp(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, m_query_pool, 3);

        // ---- Pass 4: remodulate + tonemap + pack to RGBA8 ----
        barrier_between_passes();

        // Which descriptor set to read from:
        //   spatial pass off -> the history pair the temporal pass wrote (dst),
        //     or, if denoising is entirely off, the pair holding current history.
        //   spatial pass on  -> the scratch image holding the final iteration.
        uint32_t display_src;
        if (iterations > 0) {
            display_src = 2u + filtered_index;
            m_display_source = &m_atrous_illum[filtered_index];
        } else {
            display_src = denoise.enabled ? dst : src;
            m_display_source = &m_illum_hist[denoise.enabled ? dst : src];
        }

        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, m_modulate.pipeline);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, m_modulate.pipeline_layout, 0, 1,
                                &m_modulate.sets[display_src], 0, nullptr);
        vkCmdPushConstants(cmd, m_modulate.pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(mod), &mod);
        vkCmdDispatch(cmd, groups_x, groups_y, 1);

        if (m_query_pool)
            vkCmdWriteTimestamp(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, m_query_pool, 4);
    });

    read_pass_timings();

    // Remember this frame's camera as "previous" for the next frame's
    // reprojection. The temporal pass maps a G-Buffer world position back to a
    // previous-frame pixel, so it needs:
    //   world_to_camera  - the inverse of the camera_to_world we were handed
    //   camera_to_sample - the inverse of sample_to_camera, which is named for
    //                      the direction it goes (sample -> camera). Using it
    //                      as-is reprojects with the matrix backwards, which
    //                      silently fetches history from the wrong pixels.
    m_prev_w2c = invert_transform(render.camera_to_world);
    m_prev_c2s = invert_transform(render.sample_to_camera);

    if (denoise.enabled) {
        m_history_index = dst;
        m_history_valid = true;
    } else {
        // Without the temporal pass nothing is written, so the pair holding
        // current history stays where it is and `dst` remains scratch.
        m_history_index = src;
        m_history_valid = false;
    }

    m_still_frames = reset ? 0 : m_still_frames + 1;

    m_frame_index++;
    m_accumulated_spp += rtp.spp_per_frame;

    return static_cast<const uint8_t *>(m_display_buf.mapped);
}

} // namespace tiny_renderer::gpu

// Stage 2 - host side: the wavefront path tracer.
//
// Loads a scene through the same parser/Scene machinery as gpu-path-trace,
// uploads Scene::build_gpu_scene()'s flattened geometry/materials/textures/
// lights/BVH (bindings 0-15, byte-identical to gpu-path-trace's layout - see
// shaders/common/wavefront_common.slang), plus the additional per-pixel
// persistent-path-state / work-queue buffers (bindings 16-31) that the
// megakernel didn't need. Orchestrates the 6-kernel-per-bounce sequence
// described in wavefront_common.slang's file header, using
// vkCmdDispatchIndirect so that later bounces (with fewer surviving paths)
// dispatch proportionally fewer workgroups instead of always re-dispatching
// width*height threads. Renders `spp` samples, each spp's whole bounce loop
// recorded into ONE command buffer (with explicit barriers between kernels)
// and submitted once, then reads back accum_buffer and writes a tonemapped
// PNG - matching gpu-path-trace's output conventions exactly, so the two
// executables' outputs are directly comparable.

#include <chrono>
#include <components/scene.h>
#include <core/gpu_scene.h>
#include <cstring>
#include <filesystem/resolver.h>
#include <gpu/gpu_buffer.h>
#include <gpu/vk_context.h>
#include <iostream>
#include <parse/parser.h>
#include <vector>

#include <stb_image_write.h>

using namespace tiny_renderer;
using namespace tiny_renderer::gpu;

namespace {

constexpr uint32_t GROUP_SIZE = 256; // must match wavefront_common.slang's GROUP_SIZE

// ---- Packed structs, byte-for-byte matching wavefront_common.slang's
// GPU*GPU structs (identical to gpu-path-trace's path_trace.slang layout) ----
#pragma pack(push, 1)
struct GPUTextureGPU {
    uint32_t type;
    uint32_t channels;
    uint32_t width;
    uint32_t height;
    uint32_t pixel_offset;
    float color0_r, color0_g, color0_b;
    float color1_r, color1_g, color1_b;
    float scale_u, scale_v;
};

struct GPUMaterialGPU {
    uint32_t type;
    uint32_t flags;
    uint32_t tex_reflectance;
    uint32_t tex_specular_reflectance;
    uint32_t tex_specular_transmittance;
    uint32_t tex_eta;
    uint32_t tex_k;
    uint32_t tex_alpha;
    uint32_t tex_opacity;
    uint32_t tex_bump;
    uint32_t nested_bsdf;
    uint32_t front_bsdf;
    uint32_t back_bsdf;
    float eta;
    float inv_eta;
    float alpha;
    float scale;
    uint32_t has_reflection;
    uint32_t has_transmission;
    uint32_t nonlinear;
    float fdr_int;
    float fdr_ext;
    float inv_eta_2;
    float specular_sampling_weight;
    uint32_t transmittance_lut_offset;
    uint32_t transmittance_lut_count;
    float internal_reflectance;
};

struct GPULightGPU {
    uint32_t type;
    uint32_t mesh_id;
    uint32_t radiance_tex;
    uint32_t cdf_offset;
    uint32_t cdf_count;
    float inv_total_area;
    uint32_t base_triangle;
    uint32_t spatial_varying;
    float to_world[16];
    float to_world_inv[16];
    float bsc_x, bsc_y, bsc_z;
    float bounding_sphere_radius;
};

struct PushConstants {
    uint32_t spp_index;
    uint32_t max_depth;
    uint32_t rr_depth;
    uint32_t num_lights;
    uint32_t environment_light_id;
    uint32_t depth;
    uint32_t read_is_a;
    uint32_t total_pixels;
};
#pragma pack(pop)

// One pipeline + its (single, fixed) shader module, all sharing the same
// pipeline layout / descriptor set (see main()).
struct Kernel {
    VkPipeline pipeline = VK_NULL_HANDLE;
    VkShaderModule module = VK_NULL_HANDLE;
};

Kernel create_kernel(VulkanContext &ctx, VkPipelineLayout layout, const std::string &spv_path) {
    Kernel k;
    k.module = ctx.load_shader_module(spv_path);

    VkPipelineShaderStageCreateInfo stage{};
    stage.sType  = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stage.stage  = VK_SHADER_STAGE_COMPUTE_BIT;
    stage.module = k.module;
    stage.pName  = "main";

    VkComputePipelineCreateInfo info{};
    info.sType  = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
    info.stage  = stage;
    info.layout = layout;
    VK_CHECK(vkCreateComputePipelines(ctx.device, VK_NULL_HANDLE, 1, &info, nullptr, &k.pipeline));
    return k;
}

// A full compute-to-compute (+indirect-command-read) memory barrier. Used
// between every dispatch in the per-bounce sequence: correctness (avoiding
// a hard-to-debug race between e.g. wavefront_shade's atomic counter writes
// and wavefront_prepare_shadow_indirect's read of them) is prioritized over
// the small amount of pipelining this gives up.
void full_barrier(VkCommandBuffer cmd) {
    VkMemoryBarrier barrier{};
    barrier.sType         = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
    barrier.srcAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;
    barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_INDIRECT_COMMAND_READ_BIT;
    vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT | VK_PIPELINE_STAGE_DRAW_INDIRECT_BIT, 0, 1, &barrier, 0,
                        nullptr, 0, nullptr);
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Syntax: " << argv[0] << " <scene.xml> [output.png] [spp]" << std::endl;
        return -1;
    }

    std::string scene_path  = argv[1];
    std::string output_path = argc >= 3 ? argv[2] : "gpu_wavefront.png";
    int spp                 = argc >= 4 ? std::atoi(argv[3]) : 64;

    try {
        filesystem::path path(scene_path);
        get_file_resolver()->prepend(path.parent_path());

        std::cout << "[gpu-wavefront] Loading " << scene_path << " ..." << std::endl;
        auto root = load_from_xml(scene_path);
        if (root->get_class_type() != Object::EScene) {
            std::cerr << "[gpu-wavefront] Root object is not a <scene>" << std::endl;
            return -1;
        }
        auto scene = std::dynamic_pointer_cast<Scene>(root);
        scene->construct();

        std::cout << "[gpu-wavefront] Building GPUScene ..." << std::endl;
        GPUScene gs = scene->build_gpu_scene();

        if (gs.bvh_nodes.empty()) {
            std::cerr << "[gpu-wavefront] Scene's Accel implementation did not export a flat BVH. Aborting."
                      << std::endl;
            return -1;
        }

        uint32_t width  = gs.camera.width;
        uint32_t height = gs.camera.height;
        uint32_t total_pixels = width * height;
        std::cout << "[gpu-wavefront] vertices=" << gs.vertex_positions.size()
                  << " triangles=" << gs.indices.size() / 3 << " materials=" << gs.materials.size()
                  << " textures=" << gs.textures.size() << " lights=" << gs.lights.size() << " camera=" << width
                  << "x" << height << " spp=" << spp << std::endl;

        // ---- Pack geometry (float4-padded positions/normals) ----
        std::vector<float> packed_positions(gs.vertex_positions.size() * 4);
        std::vector<float> packed_normals(gs.vertex_normals.size() * 4);
        for (size_t i = 0; i < gs.vertex_positions.size(); ++i) {
            const auto &p              = gs.vertex_positions[i];
            packed_positions[i * 4 + 0] = p.x();
            packed_positions[i * 4 + 1] = p.y();
            packed_positions[i * 4 + 2] = p.z();
            packed_positions[i * 4 + 3] = 0.0f;
            const auto &n             = gs.vertex_normals[i];
            packed_normals[i * 4 + 0]  = n.x();
            packed_normals[i * 4 + 1]  = n.y();
            packed_normals[i * 4 + 2]  = n.z();
            packed_normals[i * 4 + 3]  = 0.0f;
        }
        std::vector<float> packed_uvs(gs.vertex_uvs.size() * 2);
        for (size_t i = 0; i < gs.vertex_uvs.size(); ++i) {
            packed_uvs[i * 2 + 0] = gs.vertex_uvs[i].x();
            packed_uvs[i * 2 + 1] = gs.vertex_uvs[i].y();
        }

        // ---- Pack textures + flatten pixel data ----
        std::vector<GPUTextureGPU> gpu_textures(gs.textures.size());
        std::vector<float> texture_pixels;
        for (size_t i = 0; i < gs.textures.size(); ++i) {
            const GPUTexture &t = gs.textures[i];
            GPUTextureGPU &g    = gpu_textures[i];
            g.type              = static_cast<uint32_t>(t.type);
            g.channels          = static_cast<uint32_t>(t.channels);
            g.width              = static_cast<uint32_t>(t.width);
            g.height             = static_cast<uint32_t>(t.height);
            g.color0_r = t.color0(0); g.color0_g = t.color0(1); g.color0_b = t.color0(2);
            g.color1_r = t.color1(0); g.color1_g = t.color1(1); g.color1_b = t.color1(2);
            g.scale_u = t.scale_u; g.scale_v = t.scale_v;
            if (!t.pixels.empty()) {
                g.pixel_offset = static_cast<uint32_t>(texture_pixels.size());
                texture_pixels.insert(texture_pixels.end(), t.pixels.begin(), t.pixels.end());
            } else {
                g.pixel_offset = 0;
            }
        }

        // ---- Pack materials + flatten RoughPlastic's transmittance LUTs ----
        std::vector<GPUMaterialGPU> gpu_materials(gs.materials.size());
        std::vector<float> material_lut;
        for (size_t i = 0; i < gs.materials.size(); ++i) {
            const GPUMaterial &m = gs.materials[i];
            GPUMaterialGPU &g    = gpu_materials[i];
            g.type                        = static_cast<uint32_t>(m.type);
            g.flags                       = m.flags;
            g.tex_reflectance             = m.tex_reflectance;
            g.tex_specular_reflectance    = m.tex_specular_reflectance;
            g.tex_specular_transmittance  = m.tex_specular_transmittance;
            g.tex_eta                     = m.tex_eta;
            g.tex_k                       = m.tex_k;
            g.tex_alpha                   = m.tex_alpha;
            g.tex_opacity                 = m.tex_opacity;
            g.tex_bump                    = m.tex_bump;
            g.nested_bsdf                 = m.nested_bsdf;
            g.front_bsdf                  = m.front_bsdf;
            g.back_bsdf                   = m.back_bsdf;
            g.eta                         = m.eta;
            g.inv_eta                     = m.inv_eta;
            g.alpha                       = m.alpha;
            g.scale                       = m.scale;
            g.has_reflection               = m.has_reflection ? 1u : 0u;
            g.has_transmission             = m.has_transmission ? 1u : 0u;
            g.nonlinear                    = m.nonlinear ? 1u : 0u;
            g.fdr_int                      = m.fdr_int;
            g.fdr_ext                      = m.fdr_ext;
            g.inv_eta_2                   = m.inv_eta_2;
            g.specular_sampling_weight    = m.specular_sampling_weight;
            g.internal_reflectance         = m.internal_reflectance;
            if (!m.external_transmittance.empty()) {
                g.transmittance_lut_offset = static_cast<uint32_t>(material_lut.size());
                g.transmittance_lut_count  = static_cast<uint32_t>(m.external_transmittance.size());
                material_lut.insert(material_lut.end(), m.external_transmittance.begin(),
                                    m.external_transmittance.end());
            } else {
                g.transmittance_lut_offset = 0;
                g.transmittance_lut_count  = 0;
            }
        }
        if (material_lut.empty()) material_lut.push_back(0.0f);
        if (texture_pixels.empty()) texture_pixels.push_back(0.0f);

        // ---- Pack lights ----
        std::vector<GPULightGPU> gpu_lights(gs.lights.size());
        for (size_t i = 0; i < gs.lights.size(); ++i) {
            const GPULight &l = gs.lights[i];
            GPULightGPU &g    = gpu_lights[i];
            g.type            = static_cast<uint32_t>(l.type);
            g.mesh_id         = l.mesh_id;
            g.radiance_tex    = l.radiance_tex;
            g.cdf_offset      = l.cdf_offset;
            g.cdf_count       = l.cdf_count;
            g.inv_total_area  = l.inv_total_area;
            g.base_triangle   = l.base_triangle;
            g.spatial_varying = l.spatial_varying ? 1u : 0u;
            for (int r = 0; r < 4; ++r) {
                for (int c = 0; c < 4; ++c) {
                    g.to_world[r * 4 + c]     = l.to_world.get_transform()(r, c);
                    g.to_world_inv[r * 4 + c] = l.to_world.get_inverse()(r, c);
                }
            }
            g.bsc_x = l.bounding_sphere_center.x();
            g.bsc_y = l.bounding_sphere_center.y();
            g.bsc_z = l.bounding_sphere_center.z();
            g.bounding_sphere_radius = l.bounding_sphere_radius;
        }
        if (gpu_lights.empty()) gpu_lights.emplace_back();
        if (gs.light_triangle_cdf.empty()) const_cast<GPUScene &>(gs).light_triangle_cdf.push_back(0.0f);

        // ============================================================
        VulkanContext ctx;

        // ---- Scene buffers (bindings 0-15), identical to gpu-path-trace ----
        GpuBuffer positions_buf = create_buffer_with_data(ctx, packed_positions.data(), packed_positions.size() * sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer normals_buf   = create_buffer_with_data(ctx, packed_normals.data(), packed_normals.size() * sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer uvs_buf       = create_buffer_with_data(ctx, packed_uvs.data(), packed_uvs.size() * sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer indices_buf   = create_buffer_with_data(ctx, gs.indices.data(), gs.indices.size() * sizeof(uint32_t), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer bvh_nodes_buf = create_buffer_with_data(ctx, gs.bvh_nodes.data(), gs.bvh_nodes.size() * sizeof(GPUBVHNode), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer bvh_prims_buf = create_buffer_with_data(ctx, gs.bvh_primitives.data(), gs.bvh_primitives.size() * sizeof(uint32_t), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer tri_mat_buf   = create_buffer_with_data(ctx, gs.triangle_material_id.data(), gs.triangle_material_id.size() * sizeof(uint32_t), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer tri_light_buf = create_buffer_with_data(ctx, gs.triangle_light_id.data(), gs.triangle_light_id.size() * sizeof(uint32_t), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer materials_buf = create_buffer_with_data(ctx, gpu_materials.data(), gpu_materials.size() * sizeof(GPUMaterialGPU), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer textures_buf  = create_buffer_with_data(ctx, gpu_textures.data(), gpu_textures.size() * sizeof(GPUTextureGPU), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer texture_pixels_buf = create_buffer_with_data(ctx, texture_pixels.data(), texture_pixels.size() * sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer material_lut_buf   = create_buffer_with_data(ctx, material_lut.data(), material_lut.size() * sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer lights_buf    = create_buffer_with_data(ctx, gpu_lights.data(), gpu_lights.size() * sizeof(GPULightGPU), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer light_cdf_buf = create_buffer_with_data(ctx, gs.light_triangle_cdf.data(), gs.light_triangle_cdf.size() * sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer camera_buf    = create_buffer_with_data(ctx, &gs.camera, sizeof(GPUCamera), VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT);

        VkDeviceSize accum_size = static_cast<VkDeviceSize>(total_pixels) * sizeof(float) * 4;
        // VK_BUFFER_USAGE_TRANSFER_DST_BIT is required because we vkCmdFillBuffer
        // this buffer once (to clear it) before the spp loop below.
        GpuBuffer accum_buf = create_empty_buffer(ctx, accum_size,
                                                   VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);

        // ---- Wavefront-specific buffers (bindings 16-31) ----
        VkBufferUsageFlags storage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT;
        GpuBuffer path_ray_o_buf      = create_empty_buffer(ctx, total_pixels * sizeof(float) * 4, storage);
        GpuBuffer path_ray_d_buf      = create_empty_buffer(ctx, total_pixels * sizeof(float) * 4, storage);
        GpuBuffer path_throughput_buf = create_empty_buffer(ctx, total_pixels * sizeof(float) * 4, storage);
        GpuBuffer path_rng_buf        = create_empty_buffer(ctx, total_pixels * sizeof(uint32_t), storage);
        GpuBuffer path_prev_bsdf_buf  = create_empty_buffer(ctx, total_pixels * sizeof(float) * 2, storage);
        GpuBuffer path_prev_p_buf     = create_empty_buffer(ctx, total_pixels * sizeof(float) * 4, storage);
        GpuBuffer hit_tri_buf         = create_empty_buffer(ctx, total_pixels * sizeof(uint32_t), storage);
        GpuBuffer hit_tuv_buf         = create_empty_buffer(ctx, total_pixels * sizeof(float) * 4, storage);
        GpuBuffer active_a_buf        = create_empty_buffer(ctx, total_pixels * sizeof(uint32_t), storage);
        GpuBuffer active_b_buf        = create_empty_buffer(ctx, total_pixels * sizeof(uint32_t), storage);
        GpuBuffer counters_buf        = create_empty_buffer(ctx, 3 * sizeof(uint32_t), storage);
        GpuBuffer shadow_o_buf        = create_empty_buffer(ctx, total_pixels * sizeof(float) * 4, storage);
        GpuBuffer shadow_d_buf        = create_empty_buffer(ctx, total_pixels * sizeof(float) * 4, storage);
        GpuBuffer shadow_pixel_buf    = create_empty_buffer(ctx, total_pixels * sizeof(uint32_t), storage);
        GpuBuffer shadow_contrib_buf  = create_empty_buffer(ctx, total_pixels * sizeof(float) * 4, storage);
        GpuBuffer indirect_args_buf = create_empty_buffer(
            ctx, 6 * sizeof(uint32_t), storage | VK_BUFFER_USAGE_INDIRECT_BUFFER_BIT);

        // ---- Descriptor set layout: 32 bindings shared by all 6 kernels ----
        struct B {
            uint32_t binding;
            VkDescriptorType type;
            VkBuffer buffer;
            VkDeviceSize size;
        };
        std::vector<B> bindings = {
            { 0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, accum_buf.buffer, accum_buf.size },
            { 1, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, camera_buf.buffer, camera_buf.size },
            { 2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, positions_buf.buffer, positions_buf.size },
            { 3, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, normals_buf.buffer, normals_buf.size },
            { 4, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, uvs_buf.buffer, uvs_buf.size },
            { 5, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, indices_buf.buffer, indices_buf.size },
            { 6, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, bvh_nodes_buf.buffer, bvh_nodes_buf.size },
            { 7, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, bvh_prims_buf.buffer, bvh_prims_buf.size },
            { 8, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, tri_mat_buf.buffer, tri_mat_buf.size },
            { 9, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, tri_light_buf.buffer, tri_light_buf.size },
            { 10, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, materials_buf.buffer, materials_buf.size },
            { 11, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, textures_buf.buffer, textures_buf.size },
            { 12, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, texture_pixels_buf.buffer, texture_pixels_buf.size },
            { 13, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, material_lut_buf.buffer, material_lut_buf.size },
            { 14, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, lights_buf.buffer, lights_buf.size },
            { 15, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, light_cdf_buf.buffer, light_cdf_buf.size },
            { 16, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, path_ray_o_buf.buffer, path_ray_o_buf.size },
            { 17, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, path_ray_d_buf.buffer, path_ray_d_buf.size },
            { 18, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, path_throughput_buf.buffer, path_throughput_buf.size },
            { 19, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, path_rng_buf.buffer, path_rng_buf.size },
            { 20, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, path_prev_bsdf_buf.buffer, path_prev_bsdf_buf.size },
            { 21, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, path_prev_p_buf.buffer, path_prev_p_buf.size },
            { 22, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, hit_tri_buf.buffer, hit_tri_buf.size },
            { 23, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, hit_tuv_buf.buffer, hit_tuv_buf.size },
            { 24, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, active_a_buf.buffer, active_a_buf.size },
            { 25, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, active_b_buf.buffer, active_b_buf.size },
            { 26, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, counters_buf.buffer, counters_buf.size },
            { 27, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, shadow_o_buf.buffer, shadow_o_buf.size },
            { 28, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, shadow_d_buf.buffer, shadow_d_buf.size },
            { 29, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, shadow_pixel_buf.buffer, shadow_pixel_buf.size },
            { 30, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, shadow_contrib_buf.buffer, shadow_contrib_buf.size },
            { 31, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, indirect_args_buf.buffer, indirect_args_buf.size },
        };

        std::vector<VkDescriptorSetLayoutBinding> layout_bindings(bindings.size());
        for (size_t i = 0; i < bindings.size(); ++i) {
            layout_bindings[i]                  = {};
            layout_bindings[i].binding          = bindings[i].binding;
            layout_bindings[i].descriptorType   = bindings[i].type;
            layout_bindings[i].descriptorCount  = 1;
            layout_bindings[i].stageFlags        = VK_SHADER_STAGE_COMPUTE_BIT;
        }

        VkDescriptorSetLayoutCreateInfo layout_info{};
        layout_info.sType        = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
        layout_info.bindingCount = static_cast<uint32_t>(layout_bindings.size());
        layout_info.pBindings    = layout_bindings.data();

        VkDescriptorSetLayout descriptor_set_layout;
        VK_CHECK(vkCreateDescriptorSetLayout(ctx.device, &layout_info, nullptr, &descriptor_set_layout));

        VkPushConstantRange push_range{};
        push_range.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        push_range.offset     = 0;
        push_range.size       = sizeof(PushConstants);

        VkPipelineLayoutCreateInfo pipeline_layout_info{};
        pipeline_layout_info.sType                  = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        pipeline_layout_info.setLayoutCount         = 1;
        pipeline_layout_info.pSetLayouts            = &descriptor_set_layout;
        pipeline_layout_info.pushConstantRangeCount = 1;
        pipeline_layout_info.pPushConstantRanges    = &push_range;

        VkPipelineLayout pipeline_layout;
        VK_CHECK(vkCreatePipelineLayout(ctx.device, &pipeline_layout_info, nullptr, &pipeline_layout));

        std::string shader_dir = std::string(M_GPU_SHADER_BINARY_DIR) + "/";
        Kernel k_raygen  = create_kernel(ctx, pipeline_layout, shader_dir + "wavefront_raygen.spv");
        Kernel k_extend  = create_kernel(ctx, pipeline_layout, shader_dir + "wavefront_extend.spv");
        Kernel k_shade   = create_kernel(ctx, pipeline_layout, shader_dir + "wavefront_shade.spv");
        Kernel k_prep_sh = create_kernel(ctx, pipeline_layout, shader_dir + "wavefront_prepare_shadow_indirect.spv");
        Kernel k_shadow  = create_kernel(ctx, pipeline_layout, shader_dir + "wavefront_shadow.spv");
        Kernel k_advance = create_kernel(ctx, pipeline_layout, shader_dir + "wavefront_advance.spv");

        VkDescriptorPoolSize pool_sizes[2] = {};
        pool_sizes[0].type            = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        pool_sizes[0].descriptorCount = 31;
        pool_sizes[1].type            = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
        pool_sizes[1].descriptorCount = 1;

        VkDescriptorPoolCreateInfo pool_info{};
        pool_info.sType         = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        pool_info.maxSets       = 1;
        pool_info.poolSizeCount = 2;
        pool_info.pPoolSizes    = pool_sizes;

        VkDescriptorPool descriptor_pool;
        VK_CHECK(vkCreateDescriptorPool(ctx.device, &pool_info, nullptr, &descriptor_pool));

        VkDescriptorSetAllocateInfo set_alloc_info{};
        set_alloc_info.sType              = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
        set_alloc_info.descriptorPool     = descriptor_pool;
        set_alloc_info.descriptorSetCount = 1;
        set_alloc_info.pSetLayouts        = &descriptor_set_layout;

        VkDescriptorSet descriptor_set;
        VK_CHECK(vkAllocateDescriptorSets(ctx.device, &set_alloc_info, &descriptor_set));

        std::vector<VkDescriptorBufferInfo> buf_infos(bindings.size());
        std::vector<VkWriteDescriptorSet> writes(bindings.size());
        for (size_t i = 0; i < bindings.size(); ++i) {
            buf_infos[i]        = {};
            buf_infos[i].buffer = bindings[i].buffer;
            buf_infos[i].offset = 0;
            buf_infos[i].range  = bindings[i].size;

            writes[i]                 = {};
            writes[i].sType           = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            writes[i].dstSet          = descriptor_set;
            writes[i].dstBinding      = bindings[i].binding;
            writes[i].descriptorCount = 1;
            writes[i].descriptorType  = bindings[i].type;
            writes[i].pBufferInfo     = &buf_infos[i];
        }
        vkUpdateDescriptorSets(ctx.device, static_cast<uint32_t>(writes.size()), writes.data(), 0, nullptr);

        uint32_t max_depth = static_cast<uint32_t>(scene->get_integrator()->get_max_depth());
        uint32_t rr_depth  = static_cast<uint32_t>(scene->get_integrator()->get_rr_depth());
        uint32_t num_lights = static_cast<uint32_t>(gs.lights.size());
        uint32_t env_id      = gs.environment_light_id;

        // Zero accum_buffer ONCE before the spp loop (all subsequent updates
        // are +=, matching gpu-path-trace's incremental-accumulation model).
        ctx.submit_and_wait([&](VkCommandBuffer cmd) {
            vkCmdFillBuffer(cmd, accum_buf.buffer, 0, VK_WHOLE_SIZE, 0);
            VkMemoryBarrier barrier{};
            barrier.sType         = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
            barrier.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
            barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;
            vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1,
                                &barrier, 0, nullptr, 0, nullptr);
        });

        std::cout << "[gpu-wavefront] Rendering " << spp << " spp (max_depth=" << max_depth
                  << ", rr_depth=" << rr_depth << ") on " << ctx.device_name << " ..." << std::endl;

        auto t_start = std::chrono::steady_clock::now();

        auto bind_and_push = [&](VkCommandBuffer cmd, VkPipeline pipeline, const PushConstants &pc) {
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_layout, 0, 1, &descriptor_set, 0,
                                    nullptr);
            vkCmdPushConstants(cmd, pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(PushConstants), &pc);
        };

        for (int s = 0; s < spp; ++s) {
            ctx.submit_and_wait([&](VkCommandBuffer cmd) {
                PushConstants pc{};
                pc.spp_index            = static_cast<uint32_t>(s);
                pc.max_depth            = max_depth;
                pc.rr_depth             = rr_depth;
                pc.num_lights           = num_lights;
                pc.environment_light_id = env_id;
                pc.depth                = 0;
                pc.read_is_a            = 1; // raygen always seeds active_indices_a
                pc.total_pixels         = total_pixels;

                bind_and_push(cmd, k_raygen.pipeline, pc);
                vkCmdDispatch(cmd, (total_pixels + GROUP_SIZE - 1) / GROUP_SIZE, 1, 1);
                full_barrier(cmd);

                bool read_is_a = true;
                for (uint32_t depth = 0; depth < max_depth; ++depth) {
                    pc.depth     = depth;
                    pc.read_is_a = read_is_a ? 1u : 0u;

                    bind_and_push(cmd, k_extend.pipeline, pc);
                    vkCmdDispatchIndirect(cmd, indirect_args_buf.buffer, 0);
                    full_barrier(cmd);

                    bind_and_push(cmd, k_shade.pipeline, pc);
                    vkCmdDispatchIndirect(cmd, indirect_args_buf.buffer, 0);
                    full_barrier(cmd);

                    bind_and_push(cmd, k_prep_sh.pipeline, pc);
                    vkCmdDispatch(cmd, 1, 1, 1);
                    full_barrier(cmd);

                    bind_and_push(cmd, k_shadow.pipeline, pc);
                    vkCmdDispatchIndirect(cmd, indirect_args_buf.buffer, 3 * sizeof(uint32_t));
                    full_barrier(cmd);

                    bind_and_push(cmd, k_advance.pipeline, pc);
                    vkCmdDispatch(cmd, 1, 1, 1);
                    full_barrier(cmd);

                    read_is_a = !read_is_a;
                }
            });

            if ((s + 1) % 16 == 0 || s + 1 == spp) {
                std::cout << "[gpu-wavefront]   " << (s + 1) << "/" << spp << " spp done" << std::endl;
            }
        }

        auto t_end     = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(t_end - t_start).count();
        std::cout << "[gpu-wavefront] Done in " << elapsed << "s" << std::endl;

        // ---- Readback: average, then EXACTLY Bitmap::save_png()'s tonemap
        // (i.e. NO Reinhard - direct linear radiance -> sRGB transfer
        // function -> clamp to [0,255]), identical to gpu-path-trace's
        // readback, so this output is byte-for-byte comparable to both the
        // megakernel's and the CPU renderer's PNG for the same scene/spp
        // (see include/core/spectrum.h's TRGBSpectrum::to_srgb()). ----
        const auto *pixels = static_cast<const float *>(accum_buf.mapped);
        std::vector<unsigned char> ldr(static_cast<size_t>(width) * height * 3);
        for (uint32_t i = 0; i < total_pixels; ++i) {
            float count = std::max(pixels[i * 4 + 3], 1.0f);
            for (int c = 0; c < 3; ++c) {
                float v = pixels[i * 4 + c] / count;
                float srgb;
                if (v <= 0.0031308f) {
                    srgb = 12.92f * v;
                } else {
                    srgb = 1.055f * std::pow(v, 1.0f / 2.4f) - 0.055f;
                }
                ldr[i * 3 + c] = static_cast<unsigned char>(std::min(std::max(srgb, 0.0f), 1.0f) * 255.0f + 0.5f);
            }
        }
        stbi_write_png(output_path.c_str(), static_cast<int>(width), static_cast<int>(height), 3, ldr.data(),
                       static_cast<int>(width) * 3);
        std::cout << "[gpu-wavefront] Wrote " << output_path << std::endl;

        // ---- Cleanup ----
        auto destroy_kernel = [&](Kernel &k) {
            vkDestroyPipeline(ctx.device, k.pipeline, nullptr);
            vkDestroyShaderModule(ctx.device, k.module, nullptr);
        };
        destroy_kernel(k_raygen);
        destroy_kernel(k_extend);
        destroy_kernel(k_shade);
        destroy_kernel(k_prep_sh);
        destroy_kernel(k_shadow);
        destroy_kernel(k_advance);
        vkDestroyDescriptorPool(ctx.device, descriptor_pool, nullptr);
        vkDestroyPipelineLayout(ctx.device, pipeline_layout, nullptr);
        vkDestroyDescriptorSetLayout(ctx.device, descriptor_set_layout, nullptr);

        for (GpuBuffer *b : { &positions_buf, &normals_buf, &uvs_buf, &indices_buf, &bvh_nodes_buf, &bvh_prims_buf,
                              &tri_mat_buf, &tri_light_buf, &materials_buf, &textures_buf, &texture_pixels_buf,
                              &material_lut_buf, &lights_buf, &light_cdf_buf, &camera_buf, &accum_buf,
                              &path_ray_o_buf, &path_ray_d_buf, &path_throughput_buf, &path_rng_buf,
                              &path_prev_bsdf_buf, &path_prev_p_buf, &hit_tri_buf, &hit_tuv_buf, &active_a_buf,
                              &active_b_buf, &counters_buf, &shadow_o_buf, &shadow_d_buf, &shadow_pixel_buf,
                              &shadow_contrib_buf, &indirect_args_buf }) {
            destroy_buffer(ctx, *b);
        }

    } catch (const std::exception &e) {
        std::cerr << "[gpu-wavefront] FAILED: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

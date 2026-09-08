// GPU backend entry point for the MAIN tiny-renderer executable - see
// include/gpu/gpu_renderer.h for the rationale. This is host-side glue only:
// it packs Scene::build_gpu_scene()'s flattened GPUScene into the exact same
// GPU buffer layouts as the standalone gpu-path-trace / gpu-wavefront tools
// (see those files' headers for the full design notes on why each layout
// looks the way it does - bindings 0-15 are shared scene data, bindings
// 16-31 are wavefront-only persistent path state / work queues), runs the
// same dispatch sequence, and returns the result as a linear-radiance
// Bitmap instead of writing an already-tonemapped PNG.
//
// Deliberately NOT sharing translation units with path_trace_main.cpp /
// wavefront_main.cpp: those are standalone executables that recompile all
// of tiny-renderer's core sources into their own object files (see
// src/gpu/CMakeLists.txt's M_CORE_SOURCES glob) so they can run without the
// main tiny-renderer executable. This file is linked directly into
// tiny-renderer instead, where those core symbols already exist - pulling
// in that same glob here would produce duplicate-symbol link errors.

#include <gpu/gpu_renderer.h>

#include <chrono>
#include <core/gpu_scene.h>
#include <cstring>
#include <gpu/gpu_buffer.h>
#include <gpu/gpu_scene_upload.h>
#include <gpu/vk_context.h>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace tiny_renderer::gpu;

namespace tiny_renderer::gpu {
namespace {

// The packed scene structs (GPUTextureGPU/GPUMaterialGPU/GPULightGPU),
// PackedScene/pack_gpu_scene() and SceneBufferSet/upload_scene_buffers() used
// to live here. They now live in gpu_scene_upload.{h,cpp} so the real-time
// interactive backend (gpu/gpu_session.h) can share them instead of carrying a
// second copy that would silently drift from the shader-side layout. Only the
// push-constant layouts below are still backend-specific and stay here.

#pragma pack(push, 1)
struct MegakernelPushConstants {
    uint32_t spp_index;
    uint32_t max_depth;
    uint32_t rr_depth;
    uint32_t num_lights;
    uint32_t environment_light_id;
    uint32_t accumulate;
};

struct WavefrontPushConstants {
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

constexpr uint32_t WAVEFRONT_GROUP_SIZE = 256; // must match wavefront_common.slang's GROUP_SIZE

// Builds a Bitmap from accum_buffer's raw float4 (rgb + sample-count)
// per-pixel data, averaging by the count - LINEAR radiance, no tonemap
// (the caller applies Bitmap::save_png()'s sRGB transfer function).
std::shared_ptr<Bitmap> accum_buffer_to_bitmap(const float *pixels, uint32_t width, uint32_t height) {
    auto bitmap = std::make_shared<Bitmap>(static_cast<int>(height), static_cast<int>(width), 3);
    for (uint32_t y = 0; y < height; ++y) {
        for (uint32_t x = 0; x < width; ++x) {
            uint32_t i    = y * width + x;
            float count   = std::max(pixels[i * 4 + 3], 1.0f);
            (*bitmap)(static_cast<int>(y), static_cast<int>(x), 0) = pixels[i * 4 + 0] / count;
            (*bitmap)(static_cast<int>(y), static_cast<int>(x), 1) = pixels[i * 4 + 1] / count;
            (*bitmap)(static_cast<int>(y), static_cast<int>(x), 2) = pixels[i * 4 + 2] / count;
        }
    }
    return bitmap;
}

// Turns a final accumulation buffer plus a snapshot of that same buffer taken
// at the halfway point of the spp loop into the denoiser input set.
//
// The accumulation buffer holds a running (sum_rgb, count) per pixel, so the
// first half's mean is just the snapshot averaged by its own count, and the
// second half's mean is (final_sum - snapshot_sum) / (final_count -
// snapshot_count). The two halves are built from disjoint sets of samples and
// are therefore statistically independent - which is exactly the property the
// variance estimate 1/4 (A - B)^2 needs (see FrameBufferSet::variance).
//
// `half` may be empty, in which case only `color` is produced.
FrameBufferSet accum_buffers_to_framebuffers(const float *final_pixels, const std::vector<float> &half,
                                             uint32_t width, uint32_t height) {
    FrameBufferSet fbs;
    fbs.color = accum_buffer_to_bitmap(final_pixels, width, height);

    if (half.size() != static_cast<size_t>(width) * height * 4)
        return fbs;

    const int h = static_cast<int>(height), w = static_cast<int>(width);
    fbs.color_a  = std::make_shared<Bitmap>(h, w, 3);
    fbs.color_b  = std::make_shared<Bitmap>(h, w, 3);
    fbs.variance = std::make_shared<Bitmap>(h, w, 3);

    std::vector<float> raw(static_cast<size_t>(h) * w * 3, 0.0f);

    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            size_t i = static_cast<size_t>(y) * w + x;

            float count_a = half[i * 4 + 3];
            float count_t = final_pixels[i * 4 + 3];
            float count_b = count_t - count_a;

            for (int c = 0; c < 3; ++c) {
                float sum_a = half[i * 4 + c];
                float sum_b = final_pixels[i * 4 + c] - sum_a;

                float a = count_a > 0.0f ? sum_a / count_a : 0.0f;
                float b = count_b > 0.0f ? sum_b / count_b : 0.0f;

                (*fbs.color_a)(y, x, c) = a;
                (*fbs.color_b)(y, x, c) = b;

                // Var(mean of all samples) ~= 1/4 (A - B)^2, valid because A
                // and B average disjoint, equally-sized sample sets.
                float d                                            = a - b;
                raw[(static_cast<size_t>(y) * w + x) * 3 + c] = 0.25f * d * d;
            }
        }
    }

    // A raw half-buffer difference has one degree of freedom and is far too
    // noisy to drive a filter's edge-stopping weights directly; smooth it over
    // a 3x3 box first, matching what ImageBlock::to_framebuffers() does on the
    // CPU side so both paths hand the denoisers comparable variance buffers.
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            for (int c = 0; c < 3; ++c) {
                float sum = 0.0f;
                int n     = 0;
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        int yy = y + dy, xx = x + dx;
                        if (yy < 0 || xx < 0 || yy >= h || xx >= w)
                            continue;
                        sum += raw[(static_cast<size_t>(yy) * w + xx) * 3 + c];
                        ++n;
                    }
                }
                (*fbs.variance)(y, x, c) = n > 0 ? sum / static_cast<float>(n) : 0.0f;
            }
        }
    }

    return fbs;
}

// ----------------------------------------------------------------------
// Megakernel backend (mirrors path_trace_main.cpp, minus argv/PNG output)
// ----------------------------------------------------------------------
FrameBufferSet render_megakernel(const std::shared_ptr<Scene> &scene, const GPUScene &gs, const PackedScene &packed,
                                 int spp, const ProgressCallback &progress, bool want_aov) {
    uint32_t width  = gs.camera.width;
    uint32_t height = gs.camera.height;

    VulkanContext ctx;
    SceneBufferSet sb = upload_scene_buffers(ctx, gs, packed);

    VkDeviceSize accum_size = static_cast<VkDeviceSize>(width) * height * sizeof(float) * 4;
    GpuBuffer accum_buf     = create_empty_buffer(ctx, accum_size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);

    struct B {
        uint32_t binding;
        VkDescriptorType type;
        VkBuffer buffer;
        VkDeviceSize size;
    };
    std::vector<B> bindings = {
        { 0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, accum_buf.buffer, accum_buf.size },
        { 1, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, sb.camera_buf.buffer, sb.camera_buf.size },
        { 2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.positions_buf.buffer, sb.positions_buf.size },
        { 3, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.normals_buf.buffer, sb.normals_buf.size },
        { 4, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.uvs_buf.buffer, sb.uvs_buf.size },
        { 5, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.indices_buf.buffer, sb.indices_buf.size },
        { 6, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.bvh_nodes_buf.buffer, sb.bvh_nodes_buf.size },
        { 7, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.bvh_prims_buf.buffer, sb.bvh_prims_buf.size },
        { 8, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.tri_mat_buf.buffer, sb.tri_mat_buf.size },
        { 9, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.tri_light_buf.buffer, sb.tri_light_buf.size },
        { 10, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.materials_buf.buffer, sb.materials_buf.size },
        { 11, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.textures_buf.buffer, sb.textures_buf.size },
        { 12, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.texture_pixels_buf.buffer, sb.texture_pixels_buf.size },
        { 13, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.material_lut_buf.buffer, sb.material_lut_buf.size },
        { 14, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.lights_buf.buffer, sb.lights_buf.size },
        { 15, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.light_cdf_buf.buffer, sb.light_cdf_buf.size },
        // 16..31 are taken by the wavefront backend's per-path state; 32 is the
        // first slot every kernel that includes scene_common agrees on.
        { 32, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.packed_triangles_buf.buffer, sb.packed_triangles_buf.size },
    };

    std::vector<VkDescriptorSetLayoutBinding> layout_bindings(bindings.size());
    for (size_t i = 0; i < bindings.size(); ++i) {
        layout_bindings[i]                 = {};
        layout_bindings[i].binding         = bindings[i].binding;
        layout_bindings[i].descriptorType  = bindings[i].type;
        layout_bindings[i].descriptorCount = 1;
        layout_bindings[i].stageFlags       = VK_SHADER_STAGE_COMPUTE_BIT;
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
    push_range.size       = sizeof(MegakernelPushConstants);

    VkPipelineLayoutCreateInfo pipeline_layout_info{};
    pipeline_layout_info.sType                  = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
    pipeline_layout_info.setLayoutCount         = 1;
    pipeline_layout_info.pSetLayouts            = &descriptor_set_layout;
    pipeline_layout_info.pushConstantRangeCount = 1;
    pipeline_layout_info.pPushConstantRanges    = &push_range;

    VkPipelineLayout pipeline_layout;
    VK_CHECK(vkCreatePipelineLayout(ctx.device, &pipeline_layout_info, nullptr, &pipeline_layout));

    VkShaderModule shader_module = ctx.load_shader_module(std::string(M_GPU_SHADER_BINARY_DIR) + "/path_trace.spv");

    VkPipelineShaderStageCreateInfo stage_info{};
    stage_info.sType  = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stage_info.stage  = VK_SHADER_STAGE_COMPUTE_BIT;
    stage_info.module = shader_module;
    stage_info.pName  = "main";

    VkComputePipelineCreateInfo pipeline_info{};
    pipeline_info.sType  = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
    pipeline_info.stage  = stage_info;
    pipeline_info.layout = pipeline_layout;

    VkPipeline pipeline;
    VK_CHECK(vkCreateComputePipelines(ctx.device, VK_NULL_HANDLE, 1, &pipeline_info, nullptr, &pipeline));

    VkDescriptorPoolSize pool_sizes[2] = {};
    pool_sizes[0].type            = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    pool_sizes[0].descriptorCount = 16; // 15 scene slots + packed_triangles
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

    constexpr uint32_t group_size = 16;
    uint32_t groups_x             = (width + group_size - 1) / group_size;
    uint32_t groups_y             = (height + group_size - 1) / group_size;

    uint32_t max_depth = static_cast<uint32_t>(scene->get_integrator()->get_max_depth());
    uint32_t rr_depth  = static_cast<uint32_t>(scene->get_integrator()->get_rr_depth());

    std::cout << "[gpu] Dispatching " << spp << " spp passes (megakernel, " << groups_x << "x" << groups_y
              << " groups each) on " << ctx.device_name << " ..." << std::endl;

    // Snapshot of the accumulation buffer taken right after the first half of
    // the spp passes, from which the two independent half-buffers are derived
    // (see accum_buffers_to_framebuffers). Empty when AOVs were not requested.
    std::vector<float> half_snapshot;
    const int half_spp = spp / 2;

    for (int s = 0; s < spp; ++s) {
        MegakernelPushConstants pc{};
        pc.spp_index            = static_cast<uint32_t>(s);
        pc.max_depth            = max_depth;
        pc.rr_depth             = rr_depth;
        pc.num_lights           = static_cast<uint32_t>(gs.lights.size());
        pc.environment_light_id = gs.environment_light_id;
        pc.accumulate           = s == 0 ? 0u : 1u;

        ctx.submit_and_wait([&](VkCommandBuffer cmd) {
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_layout, 0, 1, &descriptor_set, 0,
                                    nullptr);
            vkCmdPushConstants(cmd, pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(MegakernelPushConstants), &pc);
            vkCmdDispatch(cmd, groups_x, groups_y, 1);
        });

        // submit_and_wait() above already blocked until this dispatch
        // completed, so the persistently-mapped accumulation buffer holds
        // exactly the first (s + 1) samples' running sums right now.
        if (want_aov && half_spp > 0 && s + 1 == half_spp) {
            const float *p = static_cast<const float *>(accum_buf.mapped);
            half_snapshot.assign(p, p + static_cast<size_t>(width) * height * 4);
        }

        if (progress) {
            // accum_buf.mapped is a persistently-mapped, host-visible pointer
            // (see include/gpu/gpu_buffer.h) - submit_and_wait() above already
            // blocked until the GPU finished this dispatch, so it is safe to
            // read it here for a live progress snapshot.
            auto snapshot = accum_buffer_to_bitmap(static_cast<const float *>(accum_buf.mapped), width, height);
            progress(s + 1, spp, snapshot->get_data()->get_data(), width, height);
        }
    }

    auto fbs = accum_buffers_to_framebuffers(static_cast<const float *>(accum_buf.mapped), half_snapshot, width, height);

    vkDestroyDescriptorPool(ctx.device, descriptor_pool, nullptr);
    vkDestroyPipeline(ctx.device, pipeline, nullptr);
    vkDestroyPipelineLayout(ctx.device, pipeline_layout, nullptr);
    vkDestroyShaderModule(ctx.device, shader_module, nullptr);
    vkDestroyDescriptorSetLayout(ctx.device, descriptor_set_layout, nullptr);
    sb.destroy(ctx);
    destroy_buffer(ctx, accum_buf);

    return fbs;
}

// ----------------------------------------------------------------------
// Wavefront backend (mirrors wavefront_main.cpp, minus argv/PNG output)
// ----------------------------------------------------------------------
struct Kernel {
    VkPipeline pipeline    = VK_NULL_HANDLE;
    VkShaderModule module  = VK_NULL_HANDLE;
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

void full_barrier(VkCommandBuffer cmd) {
    VkMemoryBarrier barrier{};
    barrier.sType         = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
    barrier.srcAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;
    barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_INDIRECT_COMMAND_READ_BIT;
    vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT | VK_PIPELINE_STAGE_DRAW_INDIRECT_BIT, 0, 1, &barrier, 0,
                        nullptr, 0, nullptr);
}

FrameBufferSet render_wavefront(const std::shared_ptr<Scene> &scene, const GPUScene &gs, const PackedScene &packed,
                                int spp, const ProgressCallback &progress, bool want_aov) {
    uint32_t width        = gs.camera.width;
    uint32_t height       = gs.camera.height;
    uint32_t total_pixels = width * height;

    VulkanContext ctx;
    SceneBufferSet sb = upload_scene_buffers(ctx, gs, packed);

    VkDeviceSize accum_size = static_cast<VkDeviceSize>(total_pixels) * sizeof(float) * 4;
    GpuBuffer accum_buf     = create_empty_buffer(ctx, accum_size,
                                               VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT);

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
    GpuBuffer indirect_args_buf   = create_empty_buffer(ctx, 6 * sizeof(uint32_t),
                                                       storage | VK_BUFFER_USAGE_INDIRECT_BUFFER_BIT);

    struct B {
        uint32_t binding;
        VkDescriptorType type;
        VkBuffer buffer;
        VkDeviceSize size;
    };
    std::vector<B> bindings = {
        { 0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, accum_buf.buffer, accum_buf.size },
        { 1, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, sb.camera_buf.buffer, sb.camera_buf.size },
        { 2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.positions_buf.buffer, sb.positions_buf.size },
        { 3, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.normals_buf.buffer, sb.normals_buf.size },
        { 4, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.uvs_buf.buffer, sb.uvs_buf.size },
        { 5, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.indices_buf.buffer, sb.indices_buf.size },
        { 6, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.bvh_nodes_buf.buffer, sb.bvh_nodes_buf.size },
        { 7, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.bvh_prims_buf.buffer, sb.bvh_prims_buf.size },
        { 8, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.tri_mat_buf.buffer, sb.tri_mat_buf.size },
        { 9, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.tri_light_buf.buffer, sb.tri_light_buf.size },
        { 10, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.materials_buf.buffer, sb.materials_buf.size },
        { 11, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.textures_buf.buffer, sb.textures_buf.size },
        { 12, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.texture_pixels_buf.buffer, sb.texture_pixels_buf.size },
        { 13, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.material_lut_buf.buffer, sb.material_lut_buf.size },
        { 14, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.lights_buf.buffer, sb.lights_buf.size },
        { 15, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.light_cdf_buf.buffer, sb.light_cdf_buf.size },
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
        { 32, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, sb.packed_triangles_buf.buffer, sb.packed_triangles_buf.size },
    };

    std::vector<VkDescriptorSetLayoutBinding> layout_bindings(bindings.size());
    for (size_t i = 0; i < bindings.size(); ++i) {
        layout_bindings[i]                 = {};
        layout_bindings[i].binding         = bindings[i].binding;
        layout_bindings[i].descriptorType  = bindings[i].type;
        layout_bindings[i].descriptorCount = 1;
        layout_bindings[i].stageFlags       = VK_SHADER_STAGE_COMPUTE_BIT;
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
    push_range.size       = sizeof(WavefrontPushConstants);

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
    pool_sizes[0].descriptorCount = 32; // 15 scene + 16 path-state + packed_triangles
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

    uint32_t max_depth  = static_cast<uint32_t>(scene->get_integrator()->get_max_depth());
    uint32_t rr_depth   = static_cast<uint32_t>(scene->get_integrator()->get_rr_depth());
    uint32_t num_lights = static_cast<uint32_t>(gs.lights.size());
    uint32_t env_id     = gs.environment_light_id;

    ctx.submit_and_wait([&](VkCommandBuffer cmd) {
        vkCmdFillBuffer(cmd, accum_buf.buffer, 0, VK_WHOLE_SIZE, 0);
        VkMemoryBarrier barrier{};
        barrier.sType         = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
        barrier.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
        barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;
        vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1,
                            &barrier, 0, nullptr, 0, nullptr);
    });

    std::cout << "[gpu] Dispatching " << spp << " spp passes (wavefront, max_depth=" << max_depth
              << ", rr_depth=" << rr_depth << ") on " << ctx.device_name << " ..." << std::endl;

    auto bind_and_push = [&](VkCommandBuffer cmd, VkPipeline pipeline, const WavefrontPushConstants &pc) {
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_layout, 0, 1, &descriptor_set, 0,
                                nullptr);
        vkCmdPushConstants(cmd, pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(WavefrontPushConstants), &pc);
    };

    // See render_megakernel(): a single mid-loop readback of the running
    // accumulation buffer is all that is needed to reconstruct two independent
    // half-sample estimates afterwards.
    std::vector<float> half_snapshot;
    const int half_spp = spp / 2;

    for (int s = 0; s < spp; ++s) {
        ctx.submit_and_wait([&](VkCommandBuffer cmd) {
            WavefrontPushConstants pc{};
            pc.spp_index            = static_cast<uint32_t>(s);
            pc.max_depth            = max_depth;
            pc.rr_depth             = rr_depth;
            pc.num_lights           = num_lights;
            pc.environment_light_id = env_id;
            pc.depth                = 0;
            pc.read_is_a            = 1;
            pc.total_pixels         = total_pixels;

            bind_and_push(cmd, k_raygen.pipeline, pc);
            vkCmdDispatch(cmd, (total_pixels + WAVEFRONT_GROUP_SIZE - 1) / WAVEFRONT_GROUP_SIZE, 1, 1);
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

        if (want_aov && half_spp > 0 && s + 1 == half_spp) {
            const float *p = static_cast<const float *>(accum_buf.mapped);
            half_snapshot.assign(p, p + static_cast<size_t>(width) * height * 4);
        }

        if (progress) {
            // Same reasoning as the megakernel backend above: accum_buf's
            // mapped pointer is safe to read here since submit_and_wait()
            // already blocked until this spp's dispatches all completed.
            auto snapshot = accum_buffer_to_bitmap(static_cast<const float *>(accum_buf.mapped), width, height);
            progress(s + 1, spp, snapshot->get_data()->get_data(), width, height);
        }
    }

    auto fbs = accum_buffers_to_framebuffers(static_cast<const float *>(accum_buf.mapped), half_snapshot, width, height);

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

    sb.destroy(ctx);
    for (GpuBuffer *b : { &accum_buf, &path_ray_o_buf, &path_ray_d_buf, &path_throughput_buf, &path_rng_buf,
                          &path_prev_bsdf_buf, &path_prev_p_buf, &hit_tri_buf, &hit_tuv_buf, &active_a_buf,
                          &active_b_buf, &counters_buf, &shadow_o_buf, &shadow_d_buf, &shadow_pixel_buf,
                          &shadow_contrib_buf, &indirect_args_buf }) {
        destroy_buffer(ctx, *b);
    }

    return fbs;
}

} // namespace
} // namespace tiny_renderer::gpu

namespace tiny_renderer::gpu {

namespace {

FrameBufferSet render_dispatch(const std::shared_ptr<Scene> &scene, GPUBackend backend, int spp,
                               const ProgressCallback &progress, bool want_aov) {
    if (spp <= 0) {
        spp = static_cast<int>(scene->get_sampler()->get_sample_count());
    }

    GPUScene gs = scene->build_gpu_scene();
    if (gs.bvh_nodes.empty()) {
        throw std::runtime_error(
            "[gpu] Scene's Accel implementation did not export a flat BVH (only <accelerate type=\"bvh\"> does); "
            "cannot render on GPU.");
    }

    PackedScene packed = pack_gpu_scene(gs);

    switch (backend) {
    case GPUBackend::Megakernel:
        return render_megakernel(scene, gs, packed, spp, progress, want_aov);
    case GPUBackend::Wavefront:
        return render_wavefront(scene, gs, packed, spp, progress, want_aov);
    }
    throw std::runtime_error("[gpu] Unknown GPUBackend");
}

} // namespace

std::shared_ptr<Bitmap> render_gpu(const std::shared_ptr<Scene> &scene, GPUBackend backend, int spp,
                                   const ProgressCallback &progress) {
    return render_dispatch(scene, backend, spp, progress, false).color;
}

FrameBufferSet render_gpu_with_aov(const std::shared_ptr<Scene> &scene, GPUBackend backend, int spp,
                                   const ProgressCallback &progress) {
    return render_dispatch(scene, backend, spp, progress, true);
}

} // namespace tiny_renderer::gpu

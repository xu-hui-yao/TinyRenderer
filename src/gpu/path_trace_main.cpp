// Stage 1, checkpoint 3 - host side.
//
// Loads a real scene, builds Scene::build_gpu_scene()'s GPUScene, packs the
// materials/textures/lights into GPU-uploadable structs (see path_trace.slang
// for the matching layouts), and runs the full megakernel path tracer,
// accumulating one dispatch per sample-per-pixel and averaging on readback.
// Output is written as a tone-mapped PNG for visual/quantitative comparison
// against the CPU renderer's output for the same scene.

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

// ---- Packed structs, byte-for-byte matching path_trace.slang's GPU*GPU structs ----

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
    uint32_t accumulate;
};
#pragma pack(pop)

// GPUMaterialType/GPULightType enum values must match path_trace.slang's MAT_*/LIGHT_* constants
// (both are defined from the same GPUMaterialType/GPULightType enum order in gpu_scene.h).

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Syntax: " << argv[0] << " <scene.xml> [output.png] [spp]" << std::endl;
        return -1;
    }

    std::string scene_path  = argv[1];
    std::string output_path = argc >= 3 ? argv[2] : "gpu_path_trace.png";
    int spp                 = argc >= 4 ? std::atoi(argv[3]) : 64;

    try {
        filesystem::path path(scene_path);
        get_file_resolver()->prepend(path.parent_path());

        std::cout << "[gpu-path-trace] Loading " << scene_path << " ..." << std::endl;
        auto root = load_from_xml(scene_path);
        if (root->get_class_type() != Object::EScene) {
            std::cerr << "[gpu-path-trace] Root object is not a <scene>" << std::endl;
            return -1;
        }
        auto scene = std::dynamic_pointer_cast<Scene>(root);
        scene->construct();

        std::cout << "[gpu-path-trace] Building GPUScene ..." << std::endl;
        GPUScene gs = scene->build_gpu_scene();

        if (gs.bvh_nodes.empty()) {
            std::cerr << "[gpu-path-trace] Scene's Accel implementation did not export a flat BVH. Aborting."
                      << std::endl;
            return -1;
        }

        uint32_t width  = gs.camera.width;
        uint32_t height = gs.camera.height;
        std::cout << "[gpu-path-trace] vertices=" << gs.vertex_positions.size()
                  << " triangles=" << gs.indices.size() / 3 << " materials=" << gs.materials.size()
                  << " textures=" << gs.textures.size() << " lights=" << gs.lights.size() << " camera=" << width
                  << "x" << height << " spp=" << spp << std::endl;

        // ---- Pack geometry (float4-padded positions/normals, as in checkpoint 2) ----
        std::vector<float> packed_positions(gs.vertex_positions.size() * 4);
        std::vector<float> packed_normals(gs.vertex_normals.size() * 4);
        for (size_t i = 0; i < gs.vertex_positions.size(); ++i) {
            const auto &p             = gs.vertex_positions[i];
            packed_positions[i * 4 + 0] = p.x();
            packed_positions[i * 4 + 1] = p.y();
            packed_positions[i * 4 + 2] = p.z();
            packed_positions[i * 4 + 3] = 0.0f;
            const auto &n             = gs.vertex_normals[i];
            packed_normals[i * 4 + 0]   = n.x();
            packed_normals[i * 4 + 1]   = n.y();
            packed_normals[i * 4 + 2]   = n.z();
            packed_normals[i * 4 + 3]   = 0.0f;
        }
        std::vector<float> packed_uvs(gs.vertex_uvs.size() * 2);
        for (size_t i = 0; i < gs.vertex_uvs.size(); ++i) {
            packed_uvs[i * 2 + 0] = gs.vertex_uvs[i].x();
            packed_uvs[i * 2 + 1] = gs.vertex_uvs[i].y();
        }

        // ---- Pack textures + flatten their pixel data into one shared buffer ----
        std::vector<GPUTextureGPU> gpu_textures(gs.textures.size());
        std::vector<float> texture_pixels;
        for (size_t i = 0; i < gs.textures.size(); ++i) {
            const GPUTexture &t = gs.textures[i];
            GPUTextureGPU &g    = gpu_textures[i];
            g.type              = static_cast<uint32_t>(t.type);
            g.channels           = static_cast<uint32_t>(t.channels);
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
                // RoughPlastic always has a non-empty LUT; other material types simply
                // never read transmittance_lut_* at all, so a zero-size placeholder is fine.
                g.transmittance_lut_offset = 0;
                g.transmittance_lut_count  = 0;
            }
        }
        if (material_lut.empty()) {
            material_lut.push_back(0.0f); // avoid a zero-sized buffer
        }
        if (texture_pixels.empty()) {
            texture_pixels.push_back(0.0f);
        }

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
        if (gpu_lights.empty()) {
            gpu_lights.emplace_back(); // avoid a zero-sized buffer (num_lights push constant will be 0)
        }
        if (gs.light_triangle_cdf.empty()) {
            const_cast<GPUScene &>(gs).light_triangle_cdf.push_back(0.0f);
        }

        // ============================================================
        VulkanContext ctx;

        GpuBuffer positions_buf = create_buffer_with_data(ctx, packed_positions.data(),
                                                           packed_positions.size() * sizeof(float),
                                                           VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer normals_buf = create_buffer_with_data(ctx, packed_normals.data(), packed_normals.size() * sizeof(float),
                                                        VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer uvs_buf = create_buffer_with_data(ctx, packed_uvs.data(), packed_uvs.size() * sizeof(float),
                                                    VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer indices_buf = create_buffer_with_data(ctx, gs.indices.data(), gs.indices.size() * sizeof(uint32_t),
                                                        VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer bvh_nodes_buf = create_buffer_with_data(ctx, gs.bvh_nodes.data(), gs.bvh_nodes.size() * sizeof(GPUBVHNode),
                                                          VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer bvh_prims_buf = create_buffer_with_data(ctx, gs.bvh_primitives.data(),
                                                          gs.bvh_primitives.size() * sizeof(uint32_t),
                                                          VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer tri_mat_buf = create_buffer_with_data(ctx, gs.triangle_material_id.data(),
                                                        gs.triangle_material_id.size() * sizeof(uint32_t),
                                                        VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer tri_light_buf = create_buffer_with_data(ctx, gs.triangle_light_id.data(),
                                                          gs.triangle_light_id.size() * sizeof(uint32_t),
                                                          VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer materials_buf = create_buffer_with_data(ctx, gpu_materials.data(),
                                                          gpu_materials.size() * sizeof(GPUMaterialGPU),
                                                          VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer textures_buf = create_buffer_with_data(ctx, gpu_textures.data(),
                                                         gpu_textures.size() * sizeof(GPUTextureGPU),
                                                         VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer texture_pixels_buf = create_buffer_with_data(ctx, texture_pixels.data(),
                                                               texture_pixels.size() * sizeof(float),
                                                               VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer material_lut_buf = create_buffer_with_data(ctx, material_lut.data(), material_lut.size() * sizeof(float),
                                                             VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer lights_buf = create_buffer_with_data(ctx, gpu_lights.data(), gpu_lights.size() * sizeof(GPULightGPU),
                                                       VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer light_cdf_buf = create_buffer_with_data(ctx, gs.light_triangle_cdf.data(),
                                                          gs.light_triangle_cdf.size() * sizeof(float),
                                                          VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer camera_buf = create_buffer_with_data(ctx, &gs.camera, sizeof(GPUCamera), VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT);

        VkDeviceSize accum_size = static_cast<VkDeviceSize>(width) * height * sizeof(float) * 4;
        GpuBuffer accum_buf     = create_empty_buffer(ctx, accum_size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);

        // ---- Descriptor set layout: 16 bindings, matching path_trace.slang exactly ----
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
        push_range.size       = sizeof(PushConstants);

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
        pool_sizes[0].descriptorCount = 15;
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

            writes[i]                = {};
            writes[i].sType          = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            writes[i].dstSet         = descriptor_set;
            writes[i].dstBinding     = bindings[i].binding;
            writes[i].descriptorCount = 1;
            writes[i].descriptorType = bindings[i].type;
            writes[i].pBufferInfo    = &buf_infos[i];
        }
        vkUpdateDescriptorSets(ctx.device, static_cast<uint32_t>(writes.size()), writes.data(), 0, nullptr);

        // ---- Dispatch: one pass per spp, accumulating into accum_buf ----
        constexpr uint32_t group_size = 16;
        uint32_t groups_x             = (width + group_size - 1) / group_size;
        uint32_t groups_y             = (height + group_size - 1) / group_size;

        // Mirror the scene's actual <integrator> settings exactly (rather than
        // hardcoding values) so max_depth/rr_depth match the CPU renderer for
        // a fair comparison.
        uint32_t max_depth = static_cast<uint32_t>(scene->get_integrator()->get_max_depth());
        uint32_t rr_depth  = static_cast<uint32_t>(scene->get_integrator()->get_rr_depth());

        std::cout << "[gpu-path-trace] Dispatching " << spp << " spp passes (" << groups_x << "x" << groups_y
                  << " groups each) on " << ctx.device_name << " ..." << std::endl;

        auto t_start = std::chrono::steady_clock::now();
        for (int s = 0; s < spp; ++s) {
            PushConstants pc{};
            pc.spp_index             = static_cast<uint32_t>(s);
            pc.max_depth             = max_depth;
            pc.rr_depth              = rr_depth;
            pc.num_lights            = static_cast<uint32_t>(gs.lights.size());
            pc.environment_light_id  = gs.environment_light_id;
            pc.accumulate            = s == 0 ? 0u : 1u;

            ctx.submit_and_wait([&](VkCommandBuffer cmd) {
                vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
                vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_layout, 0, 1, &descriptor_set, 0,
                                        nullptr);
                vkCmdPushConstants(cmd, pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(PushConstants), &pc);
                vkCmdDispatch(cmd, groups_x, groups_y, 1);
            });

            if ((s + 1) % 16 == 0 || s + 1 == spp) {
                std::cout << "[gpu-path-trace]   " << (s + 1) << "/" << spp << " spp done" << std::endl;
            }
        }
        auto t_end = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(t_end - t_start).count();
        std::cout << "[gpu-path-trace] Done in " << elapsed << "s" << std::endl;

        // ---- Readback: average, then EXACTLY Bitmap::save_png()'s tonemap
        // (i.e. NO Reinhard - direct linear radiance -> sRGB transfer
        // function -> clamp to [0,255]), so this output is byte-for-byte
        // comparable to the CPU renderer's PNG for the same scene/spp
        // (see include/core/spectrum.h's TRGBSpectrum::to_srgb()). ----
        const auto *pixels = static_cast<const float *>(accum_buf.mapped);
        std::vector<unsigned char> ldr(static_cast<size_t>(width) * height * 3);
        for (uint32_t i = 0; i < width * height; ++i) {
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
        std::cout << "[gpu-path-trace] Wrote " << output_path << std::endl;

        // ---- Cleanup ----
        vkDestroyDescriptorPool(ctx.device, descriptor_pool, nullptr);
        vkDestroyPipeline(ctx.device, pipeline, nullptr);
        vkDestroyPipelineLayout(ctx.device, pipeline_layout, nullptr);
        vkDestroyShaderModule(ctx.device, shader_module, nullptr);
        vkDestroyDescriptorSetLayout(ctx.device, descriptor_set_layout, nullptr);
        destroy_buffer(ctx, positions_buf);
        destroy_buffer(ctx, normals_buf);
        destroy_buffer(ctx, uvs_buf);
        destroy_buffer(ctx, indices_buf);
        destroy_buffer(ctx, bvh_nodes_buf);
        destroy_buffer(ctx, bvh_prims_buf);
        destroy_buffer(ctx, tri_mat_buf);
        destroy_buffer(ctx, tri_light_buf);
        destroy_buffer(ctx, materials_buf);
        destroy_buffer(ctx, textures_buf);
        destroy_buffer(ctx, texture_pixels_buf);
        destroy_buffer(ctx, material_lut_buf);
        destroy_buffer(ctx, lights_buf);
        destroy_buffer(ctx, light_cdf_buf);
        destroy_buffer(ctx, camera_buf);
        destroy_buffer(ctx, accum_buf);

    } catch (const std::exception &e) {
        std::cerr << "[gpu-path-trace] FAILED: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

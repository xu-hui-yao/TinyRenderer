// Stage 1, checkpoint 2 - host side.
//
// Loads a real scene (.xml) through the SAME parser/Scene machinery used by
// tiny-renderer, calls Scene::build_gpu_scene() (Stage 0's flattening export
// layer) to obtain a GPUScene, uploads its geometry + flat BVH to the GPU,
// dispatches raytrace_debug.slang (primary rays + BVH traversal + geometric
// normal visualization), reads the result back, and writes a PNG. This is a
// visual/structural correctness check for the geometry+BVH pipeline, not a
// pixel-exact comparison against the CPU renderer (no materials/lights/
// shading are involved yet - that is the next checkpoint).

#include <components/scene.h>
#include <core/gpu_scene.h>
#include <cstring>
#include <filesystem/resolver.h>
#include <gpu/gpu_buffer.h>
#include <gpu/vk_context.h>
#include <iostream>
#include <parse/parser.h>
#include <vector>

// NOTE: stb_image_write's implementation is already compiled once into
// src/textures/bitmap.cpp (which this target also builds, see
// src/gpu/CMakeLists.txt's M_CORE_SOURCES glob), so it must NOT be defined
// again here - doing so caused duplicate-symbol link errors.
#include <stb_image_write.h>

using namespace tiny_renderer;
using namespace tiny_renderer::gpu;

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Syntax: " << argv[0] << " <scene.xml> [output.png]" << std::endl;
        return -1;
    }

    std::string scene_path  = argv[1];
    std::string output_path = argc >= 3 ? argv[2] : "gpu_raytrace_debug.png";

    try {
        filesystem::path path(scene_path);
        get_file_resolver()->prepend(path.parent_path());

        std::cout << "[gpu-raytrace-debug] Loading " << scene_path << " ..." << std::endl;
        auto root = load_from_xml(scene_path);
        if (root->get_class_type() != Object::EScene) {
            std::cerr << "[gpu-raytrace-debug] Root object is not a <scene>" << std::endl;
            return -1;
        }
        auto scene = std::dynamic_pointer_cast<Scene>(root);
        scene->construct();

        std::cout << "[gpu-raytrace-debug] Building GPUScene ..." << std::endl;
        GPUScene gpu_scene = scene->build_gpu_scene();

        std::cout << "[gpu-raytrace-debug] vertices=" << gpu_scene.vertex_positions.size()
                  << " triangles=" << gpu_scene.indices.size() / 3 << " bvh_nodes=" << gpu_scene.bvh_nodes.size()
                  << " camera=" << gpu_scene.camera.width << "x" << gpu_scene.camera.height << std::endl;

        if (gpu_scene.bvh_nodes.empty()) {
            std::cerr << "[gpu-raytrace-debug] Scene's Accel implementation did not export a flat BVH "
                         "(only <accelerate type=\"bvh\"/> currently supports this). Aborting."
                      << std::endl;
            return -1;
        }

        // ---- Pack vertex positions as float4 (w unused): sidesteps any
        // ambiguity around how a shading language's StructuredBuffer<float3>
        // pads array elements; see the discussion in gpu_scene.h/the shader. ----
        std::vector<float> packed_positions(gpu_scene.vertex_positions.size() * 4);
        for (size_t i = 0; i < gpu_scene.vertex_positions.size(); ++i) {
            const auto &p            = gpu_scene.vertex_positions[i];
            packed_positions[i * 4 + 0] = p.x();
            packed_positions[i * 4 + 1] = p.y();
            packed_positions[i * 4 + 2] = p.z();
            packed_positions[i * 4 + 3] = 0.0f;
        }

        uint32_t width  = gpu_scene.camera.width;
        uint32_t height = gpu_scene.camera.height;

        VulkanContext ctx;

        GpuBuffer positions_buf =
            create_buffer_with_data(ctx, packed_positions.data(), packed_positions.size() * sizeof(float),
                                    VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer indices_buf =
            create_buffer_with_data(ctx, gpu_scene.indices.data(), gpu_scene.indices.size() * sizeof(uint32_t),
                                    VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer bvh_nodes_buf =
            create_buffer_with_data(ctx, gpu_scene.bvh_nodes.data(), gpu_scene.bvh_nodes.size() * sizeof(GPUBVHNode),
                                    VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        GpuBuffer bvh_prims_buf = create_buffer_with_data(ctx, gpu_scene.bvh_primitives.data(),
                                                          gpu_scene.bvh_primitives.size() * sizeof(uint32_t),
                                                          VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
        // See GPUCamera's comment in core/gpu_scene.h: its raw byte layout is
        // relied upon to exactly match the shader's GPUCameraGPU constant
        // buffer struct, so this is a direct memcpy with no repacking.
        GpuBuffer camera_buf = create_buffer_with_data(ctx, &gpu_scene.camera, sizeof(GPUCamera),
                                                       VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT);

        VkDeviceSize output_size_bytes = static_cast<VkDeviceSize>(width) * height * sizeof(float) * 4;
        GpuBuffer output_buf           = create_empty_buffer(ctx, output_size_bytes, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);

        // ---- Descriptor set layout: 4 storage buffers (positions, indices,
        // bvh_nodes, bvh_primitives) + 1 storage buffer (output) + 1 uniform
        // buffer (camera). Bindings must match raytrace_debug.slang exactly. ----
        struct BindingSpec {
            uint32_t binding;
            VkDescriptorType type;
        };
        const BindingSpec binding_specs[] = {
            { 0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER }, // output_image
            { 1, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER }, // camera
            { 2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER }, // positions
            { 3, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER }, // indices
            { 4, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER }, // bvh_nodes
            { 5, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER }, // bvh_primitives
        };
        constexpr uint32_t binding_count = sizeof(binding_specs) / sizeof(binding_specs[0]);

        std::vector<VkDescriptorSetLayoutBinding> layout_bindings(binding_count);
        for (uint32_t i = 0; i < binding_count; ++i) {
            layout_bindings[i]                 = {};
            layout_bindings[i].binding         = binding_specs[i].binding;
            layout_bindings[i].descriptorType  = binding_specs[i].type;
            layout_bindings[i].descriptorCount = 1;
            layout_bindings[i].stageFlags       = VK_SHADER_STAGE_COMPUTE_BIT;
        }

        VkDescriptorSetLayoutCreateInfo layout_info{};
        layout_info.sType        = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
        layout_info.bindingCount = binding_count;
        layout_info.pBindings    = layout_bindings.data();

        VkDescriptorSetLayout descriptor_set_layout;
        VK_CHECK(vkCreateDescriptorSetLayout(ctx.device, &layout_info, nullptr, &descriptor_set_layout));

        VkPipelineLayoutCreateInfo pipeline_layout_info{};
        pipeline_layout_info.sType          = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        pipeline_layout_info.setLayoutCount = 1;
        pipeline_layout_info.pSetLayouts    = &descriptor_set_layout;

        VkPipelineLayout pipeline_layout;
        VK_CHECK(vkCreatePipelineLayout(ctx.device, &pipeline_layout_info, nullptr, &pipeline_layout));

        VkShaderModule shader_module =
            ctx.load_shader_module(std::string(M_GPU_SHADER_BINARY_DIR) + "/raytrace_debug.spv");

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
        pool_sizes[0].descriptorCount = 5; // output + positions + indices + bvh_nodes + bvh_primitives
        pool_sizes[1].type            = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
        pool_sizes[1].descriptorCount = 1; // camera

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

        struct BufferBindingInfo {
            uint32_t binding;
            VkDescriptorType type;
            VkBuffer buffer;
            VkDeviceSize size;
        };
        BufferBindingInfo buffer_bindings[] = {
            { 0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, output_buf.buffer, output_buf.size },
            { 1, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, camera_buf.buffer, camera_buf.size },
            { 2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, positions_buf.buffer, positions_buf.size },
            { 3, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, indices_buf.buffer, indices_buf.size },
            { 4, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, bvh_nodes_buf.buffer, bvh_nodes_buf.size },
            { 5, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, bvh_prims_buf.buffer, bvh_prims_buf.size },
        };
        constexpr uint32_t num_buffer_bindings = sizeof(buffer_bindings) / sizeof(buffer_bindings[0]);

        std::vector<VkDescriptorBufferInfo> descriptor_buffer_infos(num_buffer_bindings);
        std::vector<VkWriteDescriptorSet> writes(num_buffer_bindings);
        for (uint32_t i = 0; i < num_buffer_bindings; ++i) {
            descriptor_buffer_infos[i]        = {};
            descriptor_buffer_infos[i].buffer = buffer_bindings[i].buffer;
            descriptor_buffer_infos[i].offset = 0;
            descriptor_buffer_infos[i].range  = buffer_bindings[i].size;

            writes[i]                = {};
            writes[i].sType          = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            writes[i].dstSet         = descriptor_set;
            writes[i].dstBinding     = buffer_bindings[i].binding;
            writes[i].descriptorCount = 1;
            writes[i].descriptorType = buffer_bindings[i].type;
            writes[i].pBufferInfo    = &descriptor_buffer_infos[i];
        }
        vkUpdateDescriptorSets(ctx.device, num_buffer_bindings, writes.data(), 0, nullptr);

        constexpr uint32_t group_size = 16;
        uint32_t groups_x             = (width + group_size - 1) / group_size;
        uint32_t groups_y             = (height + group_size - 1) / group_size;

        std::cout << "[gpu-raytrace-debug] Dispatching " << groups_x << "x" << groups_y << " groups on "
                  << ctx.device_name << " ..." << std::endl;

        ctx.submit_and_wait([&](VkCommandBuffer cmd) {
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_layout, 0, 1, &descriptor_set, 0,
                                    nullptr);
            vkCmdDispatch(cmd, groups_x, groups_y, 1);
        });

        // ---- Readback + tone-mapped (simple clamp) PNG write ----
        const auto *pixels = static_cast<const float *>(output_buf.mapped);
        std::vector<unsigned char> ldr(static_cast<size_t>(width) * height * 3);
        for (uint32_t i = 0; i < width * height; ++i) {
            for (int c = 0; c < 3; ++c) {
                float v        = pixels[i * 4 + c];
                ldr[i * 3 + c] = static_cast<unsigned char>(std::min(std::max(v, 0.0f), 1.0f) * 255.0f + 0.5f);
            }
        }
        stbi_write_png(output_path.c_str(), static_cast<int>(width), static_cast<int>(height), 3, ldr.data(),
                       static_cast<int>(width) * 3);
        std::cout << "[gpu-raytrace-debug] Wrote " << output_path << std::endl;

        // ---- Cleanup ----
        vkDestroyDescriptorPool(ctx.device, descriptor_pool, nullptr);
        vkDestroyPipeline(ctx.device, pipeline, nullptr);
        vkDestroyPipelineLayout(ctx.device, pipeline_layout, nullptr);
        vkDestroyShaderModule(ctx.device, shader_module, nullptr);
        vkDestroyDescriptorSetLayout(ctx.device, descriptor_set_layout, nullptr);
        destroy_buffer(ctx, positions_buf);
        destroy_buffer(ctx, indices_buf);
        destroy_buffer(ctx, bvh_nodes_buf);
        destroy_buffer(ctx, bvh_prims_buf);
        destroy_buffer(ctx, camera_buf);
        destroy_buffer(ctx, output_buf);

    } catch (const std::exception &e) {
        std::cerr << "[gpu-raytrace-debug] FAILED: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

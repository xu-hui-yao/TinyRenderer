// GPU compute smoke test (Stage 1, checkpoint 1) - host side.
//
// Validates the full Vulkan/MoltenVK/Slang toolchain end-to-end on this
// machine: instance+device creation (with the MoltenVK portability flags),
// shader module loading from a Slang-compiled SPIR-V binary, a compute
// pipeline with a storage buffer + push constants, dispatch, GPU->CPU
// readback via VMA, and finally a PNG write so the result can be inspected
// visually. This does not depend on Scene/GPUScene yet - that wiring is the
// next checkpoint once this baseline is confirmed working.

#include <cstdint>
#include <cstdlib>
#include <gpu/vk_context.h>
#include <iostream>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

using namespace tiny_renderer::gpu;

namespace {

struct PushConstants {
    uint32_t width;
    uint32_t height;
};

} // namespace

int main() {
    const uint32_t width  = 512;
    const uint32_t height = 512;

    try {
        VulkanContext ctx;

        // ---- Output buffer: one float4 per pixel, host-visible so we can
        // map and read it back directly (fine for a correctness smoke test;
        // the real renderer will use a device-local buffer + explicit
        // staging copy for performance). ----
        VkDeviceSize buffer_size = static_cast<VkDeviceSize>(width) * height * sizeof(float) * 4;

        VkBufferCreateInfo buffer_info{};
        buffer_info.sType       = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
        buffer_info.size        = buffer_size;
        buffer_info.usage       = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT;
        buffer_info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

        VmaAllocationCreateInfo alloc_info{};
        alloc_info.usage = VMA_MEMORY_USAGE_AUTO;
        alloc_info.flags = VMA_ALLOCATION_CREATE_HOST_ACCESS_RANDOM_BIT | VMA_ALLOCATION_CREATE_MAPPED_BIT;

        VkBuffer output_buffer;
        VmaAllocation output_allocation;
        VmaAllocationInfo output_alloc_info;
        VK_CHECK(vmaCreateBuffer(ctx.allocator, &buffer_info, &alloc_info, &output_buffer, &output_allocation,
                                 &output_alloc_info));

        // ---- Descriptor set layout: binding 0 = storage buffer ----
        VkDescriptorSetLayoutBinding binding{};
        binding.binding         = 0;
        binding.descriptorType  = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        binding.descriptorCount = 1;
        binding.stageFlags      = VK_SHADER_STAGE_COMPUTE_BIT;

        VkDescriptorSetLayoutCreateInfo layout_info{};
        layout_info.sType        = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
        layout_info.bindingCount = 1;
        layout_info.pBindings    = &binding;

        VkDescriptorSetLayout descriptor_set_layout;
        VK_CHECK(vkCreateDescriptorSetLayout(ctx.device, &layout_info, nullptr, &descriptor_set_layout));

        // ---- Pipeline layout: descriptor set + push constants ----
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

        // ---- Compute pipeline ----
        VkShaderModule shader_module = ctx.load_shader_module(std::string(M_GPU_SHADER_BINARY_DIR) + "/smoke_test.spv");

        VkPipelineShaderStageCreateInfo stage_info{};
        stage_info.sType  = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        stage_info.stage  = VK_SHADER_STAGE_COMPUTE_BIT;
        stage_info.module = shader_module;
        stage_info.pName  = "main";

        VkComputePipelineCreateInfo compute_pipeline_info{};
        compute_pipeline_info.sType  = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
        compute_pipeline_info.stage  = stage_info;
        compute_pipeline_info.layout = pipeline_layout;

        VkPipeline pipeline;
        VK_CHECK(vkCreateComputePipelines(ctx.device, VK_NULL_HANDLE, 1, &compute_pipeline_info, nullptr, &pipeline));

        // ---- Descriptor pool + set ----
        VkDescriptorPoolSize pool_size{};
        pool_size.type            = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        pool_size.descriptorCount = 1;

        VkDescriptorPoolCreateInfo pool_info{};
        pool_info.sType         = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        pool_info.maxSets       = 1;
        pool_info.poolSizeCount = 1;
        pool_info.pPoolSizes    = &pool_size;

        VkDescriptorPool descriptor_pool;
        VK_CHECK(vkCreateDescriptorPool(ctx.device, &pool_info, nullptr, &descriptor_pool));

        VkDescriptorSetAllocateInfo set_alloc_info{};
        set_alloc_info.sType              = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
        set_alloc_info.descriptorPool     = descriptor_pool;
        set_alloc_info.descriptorSetCount = 1;
        set_alloc_info.pSetLayouts        = &descriptor_set_layout;

        VkDescriptorSet descriptor_set;
        VK_CHECK(vkAllocateDescriptorSets(ctx.device, &set_alloc_info, &descriptor_set));

        VkDescriptorBufferInfo buffer_descriptor{};
        buffer_descriptor.buffer = output_buffer;
        buffer_descriptor.offset = 0;
        buffer_descriptor.range  = buffer_size;

        VkWriteDescriptorSet write{};
        write.sType           = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        write.dstSet          = descriptor_set;
        write.dstBinding      = 0;
        write.descriptorCount = 1;
        write.descriptorType  = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        write.pBufferInfo     = &buffer_descriptor;
        vkUpdateDescriptorSets(ctx.device, 1, &write, 0, nullptr);

        // ---- Dispatch ----
        PushConstants pc{ width, height };
        constexpr uint32_t group_size = 16;
        uint32_t groups_x             = (width + group_size - 1) / group_size;
        uint32_t groups_y             = (height + group_size - 1) / group_size;

        std::cout << "[gpu-smoketest] Dispatching " << groups_x << "x" << groups_y << " groups (" << width << "x"
                  << height << " pixels) on " << ctx.device_name << " ..." << std::endl;

        ctx.submit_and_wait([&](VkCommandBuffer cmd) {
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_layout, 0, 1, &descriptor_set, 0,
                                    nullptr);
            vkCmdPushConstants(cmd, pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(PushConstants), &pc);
            vkCmdDispatch(cmd, groups_x, groups_y, 1);
        });

        // ---- Readback + PNG write ----
        const auto *pixels = reinterpret_cast<const float *>(output_alloc_info.pMappedData);
        std::vector<unsigned char> ldr(static_cast<size_t>(width) * height * 3);
        for (uint32_t i = 0; i < width * height; ++i) {
            for (int c = 0; c < 3; ++c) {
                float v      = pixels[i * 4 + c];
                ldr[i * 3 + c] = static_cast<unsigned char>(std::min(std::max(v, 0.0f), 1.0f) * 255.0f + 0.5f);
            }
        }

        std::string out_path = "gpu_smoketest.png";
        stbi_write_png(out_path.c_str(), static_cast<int>(width), static_cast<int>(height), 3, ldr.data(),
                       static_cast<int>(width) * 3);
        std::cout << "[gpu-smoketest] Wrote " << out_path << ". If it shows a colorful gradient with a ring,"
                  << " the Vulkan/MoltenVK/Slang toolchain works end-to-end on this machine." << std::endl;

        // ---- Cleanup ----
        vkDestroyDescriptorPool(ctx.device, descriptor_pool, nullptr);
        vkDestroyPipeline(ctx.device, pipeline, nullptr);
        vkDestroyPipelineLayout(ctx.device, pipeline_layout, nullptr);
        vkDestroyShaderModule(ctx.device, shader_module, nullptr);
        vkDestroyDescriptorSetLayout(ctx.device, descriptor_set_layout, nullptr);
        vmaDestroyBuffer(ctx.allocator, output_buffer, output_allocation);

    } catch (const std::exception &e) {
        std::cerr << "[gpu-smoketest] FAILED: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

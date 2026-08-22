// Standalone diagnostic: verifies whether Slang's InterlockedAdd compiles
// to a real atomic on this Vulkan/MoltenVK target (see the wavefront
// energy-loss investigation in the conversation history). Dispatches N
// threads that each InterlockedAdd(counter, 1); the readback value MUST
// equal N exactly if atomics are correctly implemented.
#include <gpu/vk_context.h>
#include <iostream>

using namespace tiny_renderer::gpu;

int main() {
    constexpr uint32_t N = 1'000'000;

    try {
        VulkanContext ctx;

        VkBufferCreateInfo buffer_info{};
        buffer_info.sType       = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
        buffer_info.size        = sizeof(uint32_t);
        buffer_info.usage       = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT;
        buffer_info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

        VmaAllocationCreateInfo alloc_info{};
        alloc_info.usage = VMA_MEMORY_USAGE_AUTO;
        alloc_info.flags = VMA_ALLOCATION_CREATE_HOST_ACCESS_RANDOM_BIT | VMA_ALLOCATION_CREATE_MAPPED_BIT;

        VkBuffer buf;
        VmaAllocation alloc;
        VmaAllocationInfo out_info;
        VK_CHECK(vmaCreateBuffer(ctx.allocator, &buffer_info, &alloc_info, &buf, &alloc, &out_info));
        *static_cast<uint32_t *>(out_info.pMappedData) = 0;

        VkDescriptorSetLayoutBinding binding{};
        binding.binding         = 0;
        binding.descriptorType  = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        binding.descriptorCount = 1;
        binding.stageFlags      = VK_SHADER_STAGE_COMPUTE_BIT;

        VkDescriptorSetLayoutCreateInfo layout_info{};
        layout_info.sType        = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
        layout_info.bindingCount = 1;
        layout_info.pBindings    = &binding;
        VkDescriptorSetLayout set_layout;
        VK_CHECK(vkCreateDescriptorSetLayout(ctx.device, &layout_info, nullptr, &set_layout));

        VkPipelineLayoutCreateInfo pl_info{};
        pl_info.sType          = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        pl_info.setLayoutCount = 1;
        pl_info.pSetLayouts    = &set_layout;
        VkPipelineLayout pipeline_layout;
        VK_CHECK(vkCreatePipelineLayout(ctx.device, &pl_info, nullptr, &pipeline_layout));

        VkShaderModule module = ctx.load_shader_module(std::string(M_GPU_SHADER_BINARY_DIR) + "/atomic_test.spv");
        VkPipelineShaderStageCreateInfo stage{};
        stage.sType  = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        stage.stage  = VK_SHADER_STAGE_COMPUTE_BIT;
        stage.module = module;
        stage.pName  = "main";
        VkComputePipelineCreateInfo cp_info{};
        cp_info.sType  = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
        cp_info.stage  = stage;
        cp_info.layout = pipeline_layout;
        VkPipeline pipeline;
        VK_CHECK(vkCreateComputePipelines(ctx.device, VK_NULL_HANDLE, 1, &cp_info, nullptr, &pipeline));

        VkDescriptorPoolSize pool_size{};
        pool_size.type            = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        pool_size.descriptorCount = 1;
        VkDescriptorPoolCreateInfo pool_info{};
        pool_info.sType         = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        pool_info.maxSets       = 1;
        pool_info.poolSizeCount = 1;
        pool_info.pPoolSizes    = &pool_size;
        VkDescriptorPool pool;
        VK_CHECK(vkCreateDescriptorPool(ctx.device, &pool_info, nullptr, &pool));

        VkDescriptorSetAllocateInfo set_alloc{};
        set_alloc.sType              = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
        set_alloc.descriptorPool     = pool;
        set_alloc.descriptorSetCount = 1;
        set_alloc.pSetLayouts        = &set_layout;
        VkDescriptorSet set;
        VK_CHECK(vkAllocateDescriptorSets(ctx.device, &set_alloc, &set));

        VkDescriptorBufferInfo buf_info{};
        buf_info.buffer = buf;
        buf_info.offset = 0;
        buf_info.range  = sizeof(uint32_t);
        VkWriteDescriptorSet write{};
        write.sType           = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        write.dstSet          = set;
        write.dstBinding      = 0;
        write.descriptorCount = 1;
        write.descriptorType  = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        write.pBufferInfo     = &buf_info;
        vkUpdateDescriptorSets(ctx.device, 1, &write, 0, nullptr);

        uint32_t groups = (N + 255) / 256;
        uint32_t expected = groups * 256; // shader has no bounds check, every dispatched thread increments once
        ctx.submit_and_wait([&](VkCommandBuffer cmd) {
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_layout, 0, 1, &set, 0, nullptr);
            vkCmdDispatch(cmd, groups, 1, 1);
        });

        uint32_t result = *static_cast<const uint32_t *>(out_info.pMappedData);
        std::cout << "Dispatched " << groups << " groups x 256 threads = " << expected
                  << " expected increments; actual = " << result << std::endl;
        std::cout << (result == expected
                          ? "PASS: atomics are correct.\n"
                          : "FAIL: lost updates detected -- InterlockedAdd is NOT atomic on this target!\n");

    } catch (const std::exception &e) {
        std::cerr << "FAILED: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

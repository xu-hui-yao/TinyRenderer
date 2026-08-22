#include <cstring>
#include <gpu/gpu_buffer.h>

namespace tiny_renderer::gpu {

namespace {

GpuBuffer create_buffer(VulkanContext &ctx, VkDeviceSize size, VkBufferUsageFlags usage,
                        VmaAllocationCreateFlags vma_flags) {
    GpuBuffer buf;
    buf.size = size;

    VkBufferCreateInfo buffer_info{};
    buffer_info.sType       = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
    buffer_info.size        = size;
    buffer_info.usage       = usage;
    buffer_info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

    VmaAllocationCreateInfo alloc_info{};
    alloc_info.usage = VMA_MEMORY_USAGE_AUTO;
    alloc_info.flags = vma_flags;

    VmaAllocationInfo out_info;
    VK_CHECK(vmaCreateBuffer(ctx.allocator, &buffer_info, &alloc_info, &buf.buffer, &buf.allocation, &out_info));
    buf.mapped = out_info.pMappedData;
    return buf;
}

} // namespace

GpuBuffer create_buffer_with_data(VulkanContext &ctx, const void *data, VkDeviceSize size_bytes,
                                  VkBufferUsageFlags usage) {
    VkDeviceSize alloc_size = size_bytes == 0 ? 4 : size_bytes;
    GpuBuffer buf = create_buffer(ctx, alloc_size, usage,
                                  VMA_ALLOCATION_CREATE_HOST_ACCESS_SEQUENTIAL_WRITE_BIT |
                                      VMA_ALLOCATION_CREATE_MAPPED_BIT);
    if (data != nullptr && size_bytes > 0 && buf.mapped != nullptr) {
        std::memcpy(buf.mapped, data, size_bytes);
    }
    return buf;
}

GpuBuffer create_empty_buffer(VulkanContext &ctx, VkDeviceSize size_bytes, VkBufferUsageFlags usage) {
    VkDeviceSize alloc_size = size_bytes == 0 ? 4 : size_bytes;
    return create_buffer(ctx, alloc_size, usage,
                         VMA_ALLOCATION_CREATE_HOST_ACCESS_RANDOM_BIT | VMA_ALLOCATION_CREATE_MAPPED_BIT);
}

void destroy_buffer(VulkanContext &ctx, GpuBuffer &buffer) {
    if (buffer.buffer != VK_NULL_HANDLE) {
        vmaDestroyBuffer(ctx.allocator, buffer.buffer, buffer.allocation);
        buffer.buffer     = VK_NULL_HANDLE;
        buffer.allocation = VK_NULL_HANDLE;
        buffer.mapped     = nullptr;
    }
}

} // namespace tiny_renderer::gpu

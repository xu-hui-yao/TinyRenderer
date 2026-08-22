#pragma once

#include <gpu/vk_context.h>

namespace tiny_renderer::gpu {

// A single VMA-backed Vulkan buffer plus its allocation, bundled for
// convenient bookkeeping in the Stage 1 GPU prototype code. `mapped` is a
// persistently-mapped host pointer (valid because every buffer created here
// uses VMA_ALLOCATION_CREATE_MAPPED_BIT), so callers can read/write it
// directly without extra map/unmap calls. Ownership: the caller is
// responsible for calling destroy_buffer() when done.
struct GpuBuffer {
    VkBuffer buffer          = VK_NULL_HANDLE;
    VmaAllocation allocation = VK_NULL_HANDLE;
    VkDeviceSize size         = 0;
    void *mapped              = nullptr;
};

// Creates a buffer with the given usage flags and immediately fills it with
// `data` (size `size_bytes`) via a host-visible, mapped VMA allocation. If
// `size_bytes` is 0 (e.g. an optional/empty scene array), a minimal 4-byte
// placeholder is allocated instead, since Vulkan disallows zero-sized
// buffers - this keeps descriptor binding code uniform regardless of
// whether a given array happens to be empty for a particular scene.
//
// This is meant for one-time scene upload during setup, not the hot render
// loop (which will eventually want device-local memory plus an explicit
// staging-buffer copy for performance - out of scope for this
// correctness-focused checkpoint).
GpuBuffer create_buffer_with_data(VulkanContext &ctx, const void *data, VkDeviceSize size_bytes,
                                  VkBufferUsageFlags usage);

// Creates an empty, host-visible+mapped buffer of the given size (e.g. for
// a compute shader's output image, later read back via GpuBuffer::mapped).
GpuBuffer create_empty_buffer(VulkanContext &ctx, VkDeviceSize size_bytes, VkBufferUsageFlags usage);

void destroy_buffer(VulkanContext &ctx, GpuBuffer &buffer);

} // namespace tiny_renderer::gpu

#pragma once

// ============================================================================
// VMA-backed 2D images + image views, for the denoiser's screen-sized buffers.
//
// Why images and not the GpuBuffer storage buffers used elsewhere:
//
//   MoltenVK (the Vulkan-on-Metal translation layer this renderer runs through
//   on macOS, and the only Vulkan implementation available there) exposes Metal's
//   per-stage argument-table limits. Buffers are the scarce resource there -
//   roughly 31 usable slots per stage - while textures get 128. The scene alone
//   already occupies bindings 1..15 as buffers, and the denoiser needs a dozen
//   more screen-sized surfaces (history ping-pong, moments, G-Buffer). Binding
//   those as buffers would risk hitting the ceiling, and the failure mode is an
//   opaque MoltenVK device-lost rather than a clean validation error.
//
//   Images are also the better fit on the merits: the a-trous filter (P3) is a
//   2D gather loop, and texture samplers hit cache far better than strided
//   buffer reads for that access pattern.
//
// All images here are STORAGE (read+write from compute), TRANSFER_SRC (so a
// debug view can be read back to the host) and live in VK_IMAGE_LAYOUT_GENERAL
// from creation, matching what the shaders' RWTexture2D bindings expect.
// ============================================================================

#include <gpu/vk_context.h>

namespace tiny_renderer::gpu {

struct GpuImage {
    VkImage image               = VK_NULL_HANDLE;
    VkImageView view            = VK_NULL_HANDLE;
    VmaAllocation allocation    = VK_NULL_HANDLE;
    VkFormat format             = VK_FORMAT_UNDEFINED;
    uint32_t width              = 0;
    uint32_t height             = 0;
};

// Creates a device-local 2D storage image of `width`x`height` in `format`.
// Device-local (not host-visible): these buffers are written and read only by
// the GPU in the hot loop, and making them host-visible would push them onto
// the PCIe/BAR path on discrete GPUs or force a less optimal tiling on
// integrated ones. Anything the host needs goes through an explicit staging
// copy (see read_image_rgba8()).
GpuImage create_storage_image(VulkanContext &ctx, uint32_t width, uint32_t height, VkFormat format);

void destroy_image(VulkanContext &ctx, GpuImage &image);

} // namespace tiny_renderer::gpu

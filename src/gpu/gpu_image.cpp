#include <gpu/gpu_image.h>

#include <stdexcept>

namespace tiny_renderer::gpu {

GpuImage create_storage_image(VulkanContext &ctx, uint32_t width, uint32_t height, VkFormat format) {
    if (width == 0 || height == 0)
        throw std::runtime_error("create_storage_image: dimensions must be non-zero");

    GpuImage img;
    img.format = format;
    img.width  = width;
    img.height = height;

    VkImageCreateInfo image_info{};
    image_info.sType         = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
    image_info.imageType     = VK_IMAGE_TYPE_2D;
    image_info.extent        = { width, height, 1 };
    image_info.mipLevels     = 1;
    image_info.arrayLayers   = 1;
    image_info.format        = format;
    image_info.tiling        = VK_IMAGE_TILING_OPTIMAL;
    image_info.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    // STORAGE   - read+write from the compute passes.
    // TRANSFER_SRC/DST - lets a debug view be copied out to a host-visible
    //                    staging buffer, and lets a buffer be blitted in.
    image_info.usage       = VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT |
                       VK_IMAGE_USAGE_TRANSFER_DST_BIT;
    image_info.samples     = VK_SAMPLE_COUNT_1_BIT;
    image_info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

    VmaAllocationCreateInfo alloc_info{};
    alloc_info.usage = VMA_MEMORY_USAGE_AUTO; // device-local for OPTIMAL tiling

    VK_CHECK(vmaCreateImage(ctx.allocator, &image_info, &alloc_info, &img.image, &img.allocation, nullptr));

    VkImageViewCreateInfo view_info{};
    view_info.sType                       = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
    view_info.image                       = img.image;
    view_info.viewType                    = VK_IMAGE_VIEW_TYPE_2D;
    view_info.format                      = format;
    view_info.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    view_info.subresourceRange.baseMipLevel = 0;
    view_info.subresourceRange.levelCount   = 1;
    view_info.subresourceRange.baseArrayLayer = 0;
    view_info.subresourceRange.layerCount     = 1;
    VK_CHECK(vkCreateImageView(ctx.device, &view_info, nullptr, &img.view));

    // Images are born UNDEFINED, but storage-image access requires GENERAL.
    // Doing the transition here means every pass can assume GENERAL and no
    // per-frame barrier is ever needed for these images.
    ctx.submit_and_wait([&](VkCommandBuffer cmd) {
        VkImageMemoryBarrier barrier{};
        barrier.sType               = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;
        barrier.oldLayout           = VK_IMAGE_LAYOUT_UNDEFINED;
        barrier.newLayout           = VK_IMAGE_LAYOUT_GENERAL;
        barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.image               = img.image;
        barrier.subresourceRange    = { VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1 };
        barrier.srcAccessMask       = 0;
        barrier.dstAccessMask       = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;

        vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0,
                             0, nullptr, 0, nullptr, 1, &barrier);
    });

    return img;
}

void destroy_image(VulkanContext &ctx, GpuImage &image) {
    if (image.view)
        vkDestroyImageView(ctx.device, image.view, nullptr);
    if (image.image)
        vmaDestroyImage(ctx.allocator, image.image, image.allocation);
    image = GpuImage{};
}

} // namespace tiny_renderer::gpu

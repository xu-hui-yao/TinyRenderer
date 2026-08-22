#pragma once

// ============================================================================
// Minimal Vulkan compute context (Stage 1, checkpoint 1 of the CPU -> GPU
// path tracer port).
//
// Design notes / rationale (see the conversation history for the full
// API survey):
//   - Vulkan compute is used instead of a native per-platform API so that a
//     single implementation runs on Windows, Linux, AND macOS (via MoltenVK).
//   - MoltenVK does NOT implement VK_KHR_acceleration_structure / ray_query,
//     so this renderer's BVH traversal is (and will remain) a hand-written
//     compute shader loop over the flattened BVH from Scene::build_gpu_scene(),
//     not a hardware ray-tracing pipeline. This keeps the same shader code
//     path working uniformly on every platform.
//   - No vk-bootstrap dependency is used here: the instance/device selection
//     logic needed for a single dedicated compute queue is small enough to
//     write directly, and doing so keeps full, explicit control over the
//     MoltenVK portability enumeration flags (see create_instance()).
//   - Memory allocation goes through AMD's VulkanMemoryAllocator (VMA)
//     rather than hand-rolled vkAllocateMemory bookkeeping.
// ============================================================================

#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

// VK_KHR_portability_subset (needed to create a logical device on top of
// MoltenVK, see vk_context.cpp) lives in the "beta" extension header and
// requires this guard to be defined before including vulkan.h.
#define VK_ENABLE_BETA_EXTENSIONS
#include <vulkan/vulkan.h>

// VMA emits a large number of -Wnullability-completeness warnings under
// recent Apple Clang; they are internal to the (vendored) header and not
// actionable from this project, so they are silenced locally.
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnullability-completeness"
#endif
#include <vk_mem_alloc.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#define VK_CHECK(expr)                                                                                               \
    do {                                                                                                             \
        VkResult _vk_result = (expr);                                                                                \
        if (_vk_result != VK_SUCCESS) {                                                                             \
            throw std::runtime_error(std::string("Vulkan call failed (") + #expr + "): VkResult=" +                 \
                                     std::to_string(static_cast<int>(_vk_result)) + " at " + __FILE__ + ":" +       \
                                     std::to_string(__LINE__));                                                      \
        }                                                                                                            \
    } while (0)

namespace tiny_renderer::gpu {

// Owns an Vulkan instance + a single physical/logical device pair exposing
// one dedicated compute queue, plus a VMA allocator bound to that device.
// This is intentionally minimal: no swapchain, no graphics queue, no
// window - this project's renderer is offline / headless.
class VulkanContext {
public:
    // If `enable_validation` is true and validation layers are present on
    // the system, VK_LAYER_KHRONOS_validation is enabled and a debug
    // messenger prints validation messages to stderr.
    explicit VulkanContext(bool enable_validation = true);
    ~VulkanContext();

    VulkanContext(const VulkanContext &)            = delete;
    VulkanContext &operator=(const VulkanContext &) = delete;

    // Allocates a one-shot command buffer, invokes `record` to fill it in,
    // submits it to the compute queue, and blocks until it has completed.
    // Intended for setup/readback code, NOT the hot render loop (Stage 2's
    // wavefront renderer will need persistent command buffers instead).
    template <typename Fn> void submit_and_wait(Fn &&record) {
        VkCommandBufferAllocateInfo alloc_info{};
        alloc_info.sType              = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
        alloc_info.commandPool        = command_pool;
        alloc_info.level              = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
        alloc_info.commandBufferCount = 1;

        VkCommandBuffer cmd;
        VK_CHECK(vkAllocateCommandBuffers(device, &alloc_info, &cmd));

        VkCommandBufferBeginInfo begin_info{};
        begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
        begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
        VK_CHECK(vkBeginCommandBuffer(cmd, &begin_info));

        record(cmd);

        VK_CHECK(vkEndCommandBuffer(cmd));

        VkSubmitInfo submit_info{};
        submit_info.sType              = VK_STRUCTURE_TYPE_SUBMIT_INFO;
        submit_info.commandBufferCount = 1;
        submit_info.pCommandBuffers    = &cmd;
        VK_CHECK(vkQueueSubmit(compute_queue, 1, &submit_info, fence));

        VK_CHECK(vkWaitForFences(device, 1, &fence, VK_TRUE, UINT64_MAX));
        VK_CHECK(vkResetFences(device, 1, &fence));
        vkFreeCommandBuffers(device, command_pool, 1, &cmd);
    }

    // Loads a SPIR-V binary from disk and creates a shader module.
    [[nodiscard]] VkShaderModule load_shader_module(const std::string &spv_path) const;

    VkInstance instance                 = VK_NULL_HANDLE;
    VkPhysicalDevice physical_device    = VK_NULL_HANDLE;
    VkDevice device                     = VK_NULL_HANDLE;
    VkQueue compute_queue                = VK_NULL_HANDLE;
    uint32_t compute_queue_family        = UINT32_MAX;
    VkCommandPool command_pool          = VK_NULL_HANDLE;
    VkFence fence                       = VK_NULL_HANDLE;
    VmaAllocator allocator              = VK_NULL_HANDLE;
    VkDebugUtilsMessengerEXT debug_messenger = VK_NULL_HANDLE;
    std::string device_name;

private:
    void create_instance(bool enable_validation);
    void select_physical_device();
    void create_logical_device();
    void create_allocator();
};

} // namespace tiny_renderer::gpu

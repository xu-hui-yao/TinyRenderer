#define VMA_IMPLEMENTATION
#include <gpu/vk_context.h>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>

namespace tiny_renderer::gpu {

namespace {

VKAPI_ATTR VkBool32 VKAPI_CALL debug_callback(VkDebugUtilsMessageSeverityFlagBitsEXT severity,
                                              VkDebugUtilsMessageTypeFlagsEXT /*type*/,
                                              const VkDebugUtilsMessengerCallbackDataEXT *callback_data,
                                              void * /*user_data*/) {
    if (severity >= VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT) {
        std::cerr << "[Vulkan] " << callback_data->pMessage << std::endl;
    }
    return VK_FALSE;
}

bool layer_available(const char *name) {
    uint32_t count = 0;
    vkEnumerateInstanceLayerProperties(&count, nullptr);
    std::vector<VkLayerProperties> layers(count);
    vkEnumerateInstanceLayerProperties(&count, layers.data());
    return std::any_of(layers.begin(), layers.end(), [&](const VkLayerProperties &l) { return std::strcmp(l.layerName, name) == 0; });
}

bool instance_extension_available(const char *name) {
    uint32_t count = 0;
    vkEnumerateInstanceExtensionProperties(nullptr, &count, nullptr);
    std::vector<VkExtensionProperties> exts(count);
    vkEnumerateInstanceExtensionProperties(nullptr, &count, exts.data());
    return std::any_of(exts.begin(), exts.end(), [&](const VkExtensionProperties &e) { return std::strcmp(e.extensionName, name) == 0; });
}

bool device_extension_available(VkPhysicalDevice pd, const char *name) {
    uint32_t count = 0;
    vkEnumerateDeviceExtensionProperties(pd, nullptr, &count, nullptr);
    std::vector<VkExtensionProperties> exts(count);
    vkEnumerateDeviceExtensionProperties(pd, nullptr, &count, exts.data());
    return std::any_of(exts.begin(), exts.end(), [&](const VkExtensionProperties &e) { return std::strcmp(e.extensionName, name) == 0; });
}

} // namespace

VulkanContext::VulkanContext(bool enable_validation) {
    create_instance(enable_validation);
    select_physical_device();
    create_logical_device();
    create_allocator();

    VkCommandPoolCreateInfo pool_info{};
    pool_info.sType            = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    pool_info.flags            = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    pool_info.queueFamilyIndex = compute_queue_family;
    VK_CHECK(vkCreateCommandPool(device, &pool_info, nullptr, &command_pool));

    VkFenceCreateInfo fence_info{};
    fence_info.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
    VK_CHECK(vkCreateFence(device, &fence_info, nullptr, &fence));

    std::cout << "[VulkanContext] Ready. Device: " << device_name << std::endl;
}

VulkanContext::~VulkanContext() {
    if (device != VK_NULL_HANDLE) {
        vkDeviceWaitIdle(device);
    }
    if (fence != VK_NULL_HANDLE) {
        vkDestroyFence(device, fence, nullptr);
    }
    if (command_pool != VK_NULL_HANDLE) {
        vkDestroyCommandPool(device, command_pool, nullptr);
    }
    if (allocator != VK_NULL_HANDLE) {
        vmaDestroyAllocator(allocator);
    }
    if (device != VK_NULL_HANDLE) {
        vkDestroyDevice(device, nullptr);
    }
    if (debug_messenger != VK_NULL_HANDLE) {
        auto destroy_fn =
            reinterpret_cast<PFN_vkDestroyDebugUtilsMessengerEXT>(vkGetInstanceProcAddr(instance, "vkDestroyDebugUtilsMessengerEXT"));
        if (destroy_fn) {
            destroy_fn(instance, debug_messenger, nullptr);
        }
    }
    if (instance != VK_NULL_HANDLE) {
        vkDestroyInstance(instance, nullptr);
    }
}

void VulkanContext::create_instance(bool enable_validation) {
    VkApplicationInfo app_info{};
    app_info.sType              = VK_STRUCTURE_TYPE_APPLICATION_INFO;
    app_info.pApplicationName   = "TinyRenderer GPU Backend";
    app_info.applicationVersion = VK_MAKE_VERSION(1, 0, 0);
    app_info.pEngineName        = "TinyRenderer";
    app_info.apiVersion         = VK_API_VERSION_1_2;

    std::vector<const char *> instance_extensions;
    // Required to enumerate MoltenVK (a non-conformant "portability")
    // physical device on macOS; a no-op on Windows/Linux where the
    // extension simply is not advertised.
    bool has_portability_enum = instance_extension_available(VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME);
    if (has_portability_enum) {
        instance_extensions.push_back(VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME);
    }

    std::vector<const char *> layers;
    bool validation_available = enable_validation && layer_available("VK_LAYER_KHRONOS_validation");
    if (validation_available) {
        layers.push_back("VK_LAYER_KHRONOS_validation");
    }
    bool has_debug_utils = instance_extension_available(VK_EXT_DEBUG_UTILS_EXTENSION_NAME);
    if (validation_available && has_debug_utils) {
        instance_extensions.push_back(VK_EXT_DEBUG_UTILS_EXTENSION_NAME);
    }

    VkInstanceCreateInfo create_info{};
    create_info.sType                   = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    create_info.pApplicationInfo        = &app_info;
    create_info.enabledExtensionCount   = static_cast<uint32_t>(instance_extensions.size());
    create_info.ppEnabledExtensionNames = instance_extensions.data();
    create_info.enabledLayerCount       = static_cast<uint32_t>(layers.size());
    create_info.ppEnabledLayerNames     = layers.data();
    if (has_portability_enum) {
        create_info.flags |= VK_INSTANCE_CREATE_ENUMERATE_PORTABILITY_BIT_KHR;
    }

    VK_CHECK(vkCreateInstance(&create_info, nullptr, &instance));

    if (validation_available && has_debug_utils) {
        auto create_fn =
            reinterpret_cast<PFN_vkCreateDebugUtilsMessengerEXT>(vkGetInstanceProcAddr(instance, "vkCreateDebugUtilsMessengerEXT"));
        if (create_fn) {
            VkDebugUtilsMessengerCreateInfoEXT dbg_info{};
            dbg_info.sType = VK_STRUCTURE_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT;
            dbg_info.messageSeverity = VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT |
                                       VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT;
            dbg_info.messageType = VK_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT |
                                   VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT;
            dbg_info.pfnUserCallback = debug_callback;
            create_fn(instance, &dbg_info, nullptr, &debug_messenger);
        }
    } else if (enable_validation) {
        std::cerr << "[VulkanContext] Validation layer requested but not available on this system; continuing "
                     "without it."
                  << std::endl;
    }
}

void VulkanContext::select_physical_device() {
    uint32_t count = 0;
    VK_CHECK(vkEnumeratePhysicalDevices(instance, &count, nullptr));
    if (count == 0) {
        throw std::runtime_error("No Vulkan-capable physical devices found");
    }
    std::vector<VkPhysicalDevice> devices(count);
    VK_CHECK(vkEnumeratePhysicalDevices(instance, &count, devices.data()));

    // Prefer a discrete GPU, then integrated, then anything else. Every
    // candidate is required to expose a queue family with VK_QUEUE_COMPUTE_BIT.
    auto find_compute_family = [](VkPhysicalDevice pd) -> int {
        uint32_t qf_count = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(pd, &qf_count, nullptr);
        std::vector<VkQueueFamilyProperties> qfs(qf_count);
        vkGetPhysicalDeviceQueueFamilyProperties(pd, &qf_count, qfs.data());
        for (uint32_t i = 0; i < qf_count; ++i) {
            if (qfs[i].queueFlags & VK_QUEUE_COMPUTE_BIT) {
                return static_cast<int>(i);
            }
        }
        return -1;
    };

    VkPhysicalDevice best           = VK_NULL_HANDLE;
    int best_family                 = -1;
    int best_score                  = -1;
    for (auto pd : devices) {
        int family = find_compute_family(pd);
        if (family < 0) {
            continue;
        }
        VkPhysicalDeviceProperties props;
        vkGetPhysicalDeviceProperties(pd, &props);
        int score = 0;
        if (props.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU) {
            score = 3;
        } else if (props.deviceType == VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU) {
            score = 2;
        } else {
            score = 1;
        }
        if (score > best_score) {
            best_score  = score;
            best        = pd;
            best_family = family;
        }
    }

    if (best == VK_NULL_HANDLE) {
        throw std::runtime_error("No physical device with a compute-capable queue family found");
    }

    physical_device       = best;
    compute_queue_family  = static_cast<uint32_t>(best_family);

    VkPhysicalDeviceProperties props;
    vkGetPhysicalDeviceProperties(physical_device, &props);
    device_name = props.deviceName;
}

void VulkanContext::create_logical_device() {
    float queue_priority = 1.0f;
    VkDeviceQueueCreateInfo queue_info{};
    queue_info.sType            = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
    queue_info.queueFamilyIndex = compute_queue_family;
    queue_info.queueCount       = 1;
    queue_info.pQueuePriorities = &queue_priority;

    std::vector<const char *> device_extensions;
    // Required on macOS: MoltenVK only exposes a "portability subset" of
    // Vulkan, and the spec mandates enabling this extension when present.
    if (device_extension_available(physical_device, VK_KHR_PORTABILITY_SUBSET_EXTENSION_NAME)) {
        device_extensions.push_back(VK_KHR_PORTABILITY_SUBSET_EXTENSION_NAME);
    }

    VkDeviceCreateInfo device_info{};
    device_info.sType                   = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
    device_info.queueCreateInfoCount    = 1;
    device_info.pQueueCreateInfos       = &queue_info;
    device_info.enabledExtensionCount   = static_cast<uint32_t>(device_extensions.size());
    device_info.ppEnabledExtensionNames = device_extensions.data();

    VK_CHECK(vkCreateDevice(physical_device, &device_info, nullptr, &device));
    vkGetDeviceQueue(device, compute_queue_family, 0, &compute_queue);
}

void VulkanContext::create_allocator() {
    VmaAllocatorCreateInfo alloc_info{};
    alloc_info.physicalDevice = physical_device;
    alloc_info.device         = device;
    alloc_info.instance       = instance;
    alloc_info.vulkanApiVersion = VK_API_VERSION_1_2;
    VK_CHECK(vmaCreateAllocator(&alloc_info, &allocator));
}

VkShaderModule VulkanContext::load_shader_module(const std::string &spv_path) const {
    std::ifstream file(spv_path, std::ios::binary | std::ios::ate);
    if (!file) {
        throw std::runtime_error("Failed to open SPIR-V shader: " + spv_path);
    }
    auto size = static_cast<size_t>(file.tellg());
    std::vector<char> buffer(size);
    file.seekg(0);
    file.read(buffer.data(), static_cast<std::streamsize>(size));

    VkShaderModuleCreateInfo create_info{};
    create_info.sType    = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
    create_info.codeSize = buffer.size();
    create_info.pCode    = reinterpret_cast<const uint32_t *>(buffer.data());

    VkShaderModule module;
    VK_CHECK(vkCreateShaderModule(device, &create_info, nullptr, &module));
    return module;
}

} // namespace tiny_renderer::gpu

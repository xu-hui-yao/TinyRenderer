# ----------------------------------------------------------------------------
# GPU backend build rules (Stage 1 checkpoint 1: Vulkan/MoltenVK/Slang smoke test)
# ----------------------------------------------------------------------------

find_package(Vulkan REQUIRED)

include(FetchContent)
FetchContent_Declare(
        vma
        GIT_REPOSITORY https://github.com/GPUOpen-LibrariesAndSDKs/VulkanMemoryAllocator.git
        GIT_TAG v3.1.0
        GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(vma)

# Locate slangc (Slang -> SPIR-V compiler). It ships as part of the Vulkan
# SDK installation used by this project (same install as glslc/vulkaninfo).
find_program(SLANGC_EXECUTABLE slangc REQUIRED)
message(STATUS "Using slangc: ${SLANGC_EXECUTABLE}")

# ---- Shader compilation: shaders/*.slang -> ${CMAKE_CURRENT_BINARY_DIR}/shaders/*.spv ----
set(M_GPU_SHADER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/shaders")
set(M_GPU_SHADER_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/shaders")
file(MAKE_DIRECTORY ${M_GPU_SHADER_BINARY_DIR})

file(GLOB M_GPU_SHADER_SOURCES "${M_GPU_SHADER_DIR}/*.slang")

# Headers #include'd by the kernels above (they live in a subdirectory so
# they stay out of the non-recursive glob - they have no [shader] entry
# point and must not be compiled standalone). They still have to be listed
# as DEPENDS, otherwise editing e.g. common/wavefront_common.slang silently
# leaves every wavefront_*.spv stale.
file(GLOB_RECURSE M_GPU_SHADER_HEADERS "${M_GPU_SHADER_DIR}/common/*.slang")

set(M_GPU_SHADER_OUTPUTS "")
foreach (shader_src ${M_GPU_SHADER_SOURCES})
    get_filename_component(shader_name ${shader_src} NAME_WE)
    set(shader_out "${M_GPU_SHADER_BINARY_DIR}/${shader_name}.spv")
    add_custom_command(
            OUTPUT ${shader_out}
            COMMAND ${SLANGC_EXECUTABLE} ${shader_src}
            -profile sm_6_0
            -target spirv
            -entry main
            -o ${shader_out}
            DEPENDS ${shader_src} ${M_GPU_SHADER_HEADERS}
            COMMENT "Compiling Slang shader ${shader_name}.slang -> SPIR-V"
            VERBATIM
    )
    list(APPEND M_GPU_SHADER_OUTPUTS ${shader_out})
endforeach ()

add_custom_target(gpu_shaders DEPENDS ${M_GPU_SHADER_OUTPUTS})

# ---- gpu-smoketest executable ----
add_executable(gpu-smoketest
        smoke_test_main.cpp
        vk_context.cpp
)
add_dependencies(gpu-smoketest gpu_shaders)

target_include_directories(gpu-smoketest PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/stb
)
target_link_libraries(gpu-smoketest PRIVATE Vulkan::Vulkan GPUOpen::VulkanMemoryAllocator)

target_compile_definitions(gpu-smoketest PRIVATE
        M_GPU_SHADER_BINARY_DIR="${M_GPU_SHADER_BINARY_DIR}"
        $<$<CONFIG:Debug>:M_DEBUG>
)

set_target_properties(gpu-smoketest PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
)

# ---- gpu-raytrace-debug executable ----
#
# Stage 1 checkpoint 2: loads a real scene through the same parser/Scene code
# tiny-renderer uses, uploads Scene::build_gpu_scene()'s flattened geometry +
# BVH, and runs a BVH-traversal + normal-visualization compute shader.
#
# This needs tiny-renderer's core sources (everything except main/main.cpp),
# so they are compiled a second time into this target via a glob, rather
# than restructuring tiny-renderer's existing per-directory CMakeLists.txt
# files into a shared static library - lower risk (zero changes to the
# working tiny-renderer target) at the cost of a modest extra compile.
file(GLOB_RECURSE M_CORE_SOURCES
        ${CMAKE_SOURCE_DIR}/src/core/*.cpp
        ${CMAKE_SOURCE_DIR}/src/parse/*.cpp
        ${CMAKE_SOURCE_DIR}/src/components/*.cpp
        ${CMAKE_SOURCE_DIR}/src/shapes/*.cpp
        ${CMAKE_SOURCE_DIR}/src/bsdfs/*.cpp
        ${CMAKE_SOURCE_DIR}/src/sensors/*.cpp
        ${CMAKE_SOURCE_DIR}/src/sampler/*.cpp
        ${CMAKE_SOURCE_DIR}/src/rfilters/*.cpp
        ${CMAKE_SOURCE_DIR}/src/integrators/*.cpp
        ${CMAKE_SOURCE_DIR}/src/render/*.cpp
        ${CMAKE_SOURCE_DIR}/src/emitters/*.cpp
        ${CMAKE_SOURCE_DIR}/src/textures/*.cpp
        ${CMAKE_SOURCE_DIR}/src/accelerate/*.cpp
)

add_executable(gpu-raytrace-debug
        raytrace_debug_main.cpp
        vk_context.cpp
        gpu_buffer.cpp
        ${M_CORE_SOURCES}
)
add_dependencies(gpu-raytrace-debug gpu_shaders)

target_include_directories(gpu-raytrace-debug PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/stb
        ${CMAKE_SOURCE_DIR}/ext/pugixml
        ${CMAKE_SOURCE_DIR}/ext/tinyexr
)
target_link_libraries(gpu-raytrace-debug PRIVATE Vulkan::Vulkan GPUOpen::VulkanMemoryAllocator pugixml)

find_package(Threads REQUIRED)
target_link_libraries(gpu-raytrace-debug PRIVATE Threads::Threads)

target_compile_definitions(gpu-raytrace-debug PRIVATE
        M_GPU_SHADER_BINARY_DIR="${M_GPU_SHADER_BINARY_DIR}"
        M_PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        $<$<CONFIG:Debug>:M_DEBUG>
)

set_target_properties(gpu-raytrace-debug PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
)

# ---- gpu-path-trace executable ----
#
# Stage 1 checkpoint 3: the full megakernel path tracer (materials, textures,
# lights, MIS, multi-bounce) - see shaders/path_trace.slang. Reuses the same
# core-sources glob as gpu-raytrace-debug.
add_executable(gpu-path-trace
        path_trace_main.cpp
        vk_context.cpp
        gpu_buffer.cpp
        ${M_CORE_SOURCES}
)
add_dependencies(gpu-path-trace gpu_shaders)

target_include_directories(gpu-path-trace PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/stb
        ${CMAKE_SOURCE_DIR}/ext/pugixml
        ${CMAKE_SOURCE_DIR}/ext/tinyexr
)
target_link_libraries(gpu-path-trace PRIVATE Vulkan::Vulkan GPUOpen::VulkanMemoryAllocator pugixml Threads::Threads)

target_compile_definitions(gpu-path-trace PRIVATE
        M_GPU_SHADER_BINARY_DIR="${M_GPU_SHADER_BINARY_DIR}"
        M_PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        $<$<CONFIG:Debug>:M_DEBUG>
)

set_target_properties(gpu-path-trace PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
)

# ---- gpu-wavefront executable ----
#
# Stage 2: the wavefront path tracer - the same rendering algorithm as
# gpu-path-trace's megakernel, restructured into 5 separate kernels
# (raygen/extend/shade/shadow + 2 tiny indirect-dispatch-arg housekeeping
# kernels) connected by persistent per-pixel state + compacted work queues.
# See shaders/common/wavefront_common.slang for the shared code and
# shaders/wavefront_*.slang for each stage. Reuses the same core-sources
# glob and scene-buffer layout (bindings 0-15) as gpu-path-trace, so the
# same GPUScene upload code works for both.
add_executable(gpu-wavefront
        wavefront_main.cpp
        vk_context.cpp
        gpu_buffer.cpp
        ${M_CORE_SOURCES}
)
add_dependencies(gpu-wavefront gpu_shaders)

target_include_directories(gpu-wavefront PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ext/stb
        ${CMAKE_SOURCE_DIR}/ext/pugixml
        ${CMAKE_SOURCE_DIR}/ext/tinyexr
)
target_link_libraries(gpu-wavefront PRIVATE Vulkan::Vulkan GPUOpen::VulkanMemoryAllocator pugixml Threads::Threads)

target_compile_definitions(gpu-wavefront PRIVATE
        M_GPU_SHADER_BINARY_DIR="${M_GPU_SHADER_BINARY_DIR}"
        M_PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        $<$<CONFIG:Debug>:M_DEBUG>
)

set_target_properties(gpu-wavefront PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
)

# ---- atomic-test executable (diagnostic, not part of the renderer) ----
add_executable(atomic-test atomic_test_main.cpp vk_context.cpp)
add_dependencies(atomic-test gpu_shaders)
target_include_directories(atomic-test PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(atomic-test PRIVATE Vulkan::Vulkan GPUOpen::VulkanMemoryAllocator)
target_compile_definitions(atomic-test PRIVATE M_GPU_SHADER_BINARY_DIR="${M_GPU_SHADER_BINARY_DIR}")
set_target_properties(atomic-test PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

# ----------------------------------------------------------------------------
# gpu-renderer-lib: the GPU backend as a LIBRARY, linked directly into the
# main tiny-renderer executable (see include/gpu/gpu_renderer.h).
#
# Deliberately does NOT compile the M_CORE_SOURCES glob used by the
# gpu-raytrace-debug/gpu-path-trace/gpu-wavefront executables above: those
# are standalone binaries with no other source of Scene/BSDF/etc. symbols,
# so they need their own copy. tiny-renderer already compiles all of those
# sources itself (via add_subdirectory(src) in the top-level CMakeLists.txt)
# - linking gpu-renderer-lib's object files alongside them only ADDS the
# genuinely new GPU-only symbols (VulkanContext, GpuBuffer, render_gpu()),
# so no duplicate-symbol link errors result.
# ----------------------------------------------------------------------------
add_library(gpu-renderer-lib STATIC
        gpu_renderer.cpp
        vk_context.cpp
        gpu_buffer.cpp
)
add_dependencies(gpu-renderer-lib gpu_shaders)

target_include_directories(gpu-renderer-lib PUBLIC
        ${CMAKE_SOURCE_DIR}/include
)
target_link_libraries(gpu-renderer-lib PUBLIC Vulkan::Vulkan GPUOpen::VulkanMemoryAllocator)

target_compile_definitions(gpu-renderer-lib PUBLIC
        M_GPU_SHADER_BINARY_DIR="${M_GPU_SHADER_BINARY_DIR}"
        M_PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        M_ENABLE_GPU_BACKEND
)

set_target_properties(gpu-renderer-lib PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
)


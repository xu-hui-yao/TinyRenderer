cmake_minimum_required(VERSION 3.19)

project(tiny-renderer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_CUDA_STANDARD 17)
#set(CMAKE_CUDA_STANDARD_REQUIRED ON)

#if (MSVC)
#    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /fsanitize=address /Zi /MDd")
#    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /INCREMENTAL:NO /DEBUG")
#    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /fsanitize=address")
#endif ()
if(MSVC)
    add_compile_definitions(NOMINMAX)
endif()

#set(CUDA_TOOLKIT_ROOT_DIR "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.6")
#set(CUDA_INCLUDE_DIRS "${CUDA_TOOLKIT_ROOT_DIR}/include")
#set(CUDA_LIBRARIES "${CUDA_TOOLKIT_ROOT_DIR}/lib/x64/cudart.lib")
#set(CUBLAS_LIBRARIES "${CUDA_TOOLKIT_ROOT_DIR}/lib/x64/cublas.lib")

add_subdirectory(src)
add_subdirectory(ext)

set(THREADS_PREFER_PTHREAD_FLAG OFF) # use system thread library
find_package(Threads REQUIRED)
target_link_libraries(tiny-renderer PRIVATE Threads::Threads)

target_include_directories(tiny-renderer PRIVATE ${CMAKE_SOURCE_DIR}/include)

#target_include_directories(tiny-renderer PRIVATE ${CUDA_INCLUDE_DIRS})
#target_link_libraries(tiny-renderer PRIVATE ${CUDA_LIBRARIES} ${CUBLAS_LIBRARIES})

#set_target_properties(tiny-renderer PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Define M_DEBUG only for Debug builds
target_compile_definitions(tiny-renderer PRIVATE
        $<$<CONFIG:Debug>:M_DEBUG>
        M_PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)

# OpenMP
#find_package(OpenMP REQUIRED)
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(tiny-renderer PRIVATE OpenMP::OpenMP_CXX)
#endif()

# ----------------------------------------------------------------------------
# GPU backend (Stage 1 of the CPU -> GPU path tracer port).
#
# Uses Vulkan compute (portable to Windows/Linux/macOS; on macOS this runs
# through MoltenVK, which does not implement hardware ray tracing extensions,
# so BVH traversal is done in-shader rather than via VK_KHR_ray_query). Shaders
# are written in Slang and cross-compiled to SPIR-V, so the same shader source
# can later be retargeted to other backends (e.g. native Metal) without a
# host-side rewrite. See include/gpu/vk_context.h for the rationale.
#
# This is opt-in (M_ENABLE_GPU). It builds a handful of standalone debug/
# comparison executables (gpu-smoketest, gpu-raytrace-debug, gpu-path-trace,
# gpu-wavefront, atomic-test - see src/gpu/CMakeLists.txt), AND a
# `gpu-renderer-lib` static library that is linked directly into the main
# `tiny-renderer` executable below, powering its `--gpu[=megakernel|wavefront]`
# CLI flag (see src/main/main.cpp and include/gpu/gpu_renderer.h).
# ----------------------------------------------------------------------------
option(M_ENABLE_GPU "Build the experimental Vulkan GPU backend" ON)
if (M_ENABLE_GPU)
    add_subdirectory(src/gpu)

    # Link the GPU backend library (see src/gpu/CMakeLists.txt's
    # gpu-renderer-lib target) directly into the main tiny-renderer
    # executable, enabling its --gpu[=megakernel|wavefront] CLI flag
    # (see src/main/main.cpp). gpu-renderer-lib PUBLIC-propagates the
    # M_ENABLE_GPU_BACKEND compile definition and its include directories,
    # so no further target_compile_definitions/target_include_directories
    # call is needed here.
    target_link_libraries(tiny-renderer PRIVATE gpu-renderer-lib)
endif ()

# ----------------------------------------------------------------------------
# Optional live preview window (GLFW + Dear ImGui) powering --progress's
# on-screen display (see include/render/progress.h and src/main/main.cpp).
# Purely a display convenience, completely decoupled from both the CPU and
# GPU (headless Vulkan) render paths: if this is OFF, or glfw3/OpenGL are
# unavailable on the build machine, --progress still works and simply falls
# back to console-only progress reporting.
# ----------------------------------------------------------------------------
option(M_ENABLE_PREVIEW_GUI "Build the ImGui live preview window used by --progress" ON)
if (M_ENABLE_PREVIEW_GUI)
    find_package(glfw3 QUIET)
    find_package(OpenGL QUIET)
    if (glfw3_FOUND AND OpenGL_FOUND)
        add_subdirectory(src/gui)
        target_link_libraries(tiny-renderer PRIVATE preview-gui-lib)
    else ()
        message(STATUS "M_ENABLE_PREVIEW_GUI is ON but glfw3/OpenGL were not found - "
                       "building without the live preview window (--progress will still "
                       "print console progress).")
    endif ()
endif ()


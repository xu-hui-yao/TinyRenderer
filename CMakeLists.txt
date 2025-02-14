cmake_minimum_required(VERSION 3.19)

project(tiny-renderer LANGUAGES CUDA CXX)

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

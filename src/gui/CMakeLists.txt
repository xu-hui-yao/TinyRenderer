# ----------------------------------------------------------------------------
# preview-gui-lib: the optional live-preview window used by --progress (see
# include/render/progress.h / include/gui/preview_window.h). Only added from
# the top-level CMakeLists.txt when M_ENABLE_PREVIEW_GUI is ON AND both glfw3
# and OpenGL were found - see there for the fallback (console-only progress)
# when they are not.
#
# Dear ImGui itself is fetched the same way src/gpu/CMakeLists.txt already
# fetches VulkanMemoryAllocator (FetchContent, pinned tag) rather than being
# vendored under ext/ - it has no build system of its own to speak of, just
# a handful of source files compiled directly into this library below.
# ----------------------------------------------------------------------------

include(FetchContent)
FetchContent_Declare(
        imgui
        GIT_REPOSITORY https://github.com/ocornut/imgui.git
        GIT_TAG v1.92.9b
        GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(imgui)

add_library(preview-gui-lib STATIC
        preview_window.cpp
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)

target_include_directories(preview-gui-lib PUBLIC
        ${CMAKE_SOURCE_DIR}/include
        ${imgui_SOURCE_DIR}
        ${imgui_SOURCE_DIR}/backends
)

target_link_libraries(preview-gui-lib PUBLIC glfw OpenGL::GL)

target_compile_definitions(preview-gui-lib PUBLIC M_ENABLE_PREVIEW_GUI)

set_target_properties(preview-gui-lib PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
)

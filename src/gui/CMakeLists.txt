# ----------------------------------------------------------------------------
# preview-gui-lib: the optional live-preview window used by --progress (see
# include/render/progress.h / include/gui/preview_window.h). Only added from
# the top-level CMakeLists.txt when M_ENABLE_PREVIEW_GUI is ON AND both glfw3
# and OpenGL were found - see there for the fallback (console-only progress)
# when they are not.
#
# Dear ImGui is vendored as source under ext/imgui (it has no build system of
# its own to speak of, just a handful of source files compiled directly into
# this library below). Previously this was fetched via FetchContent, but on
# some platforms (e.g. Windows) that could not be relied on to work out of
# the box, so it is now built directly from the local checkout instead.
# ----------------------------------------------------------------------------

set(imgui_SOURCE_DIR ${CMAKE_SOURCE_DIR}/ext/imgui)

add_library(preview-gui-lib STATIC
        preview_window.cpp
        interactive_window.cpp
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

# core/array.h -> core/common.h -> filesystem/resolver.h needs
# M_PROJECT_SOURCE_DIR unconditionally (it seeds the file resolver's search
# path). The top-level CMakeLists.txt only defines it PRIVATE for the
# tiny-renderer target, so this library needs its own copy now that
# interactive_window.cpp pulls in core math types (via camera_controller.h).
target_compile_definitions(preview-gui-lib PUBLIC
        M_ENABLE_PREVIEW_GUI
        M_PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)

set_target_properties(preview-gui-lib PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
)

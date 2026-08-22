#pragma once

// Optional, decoupled live-preview window for the --progress CLI flag (see
// include/render/progress.h). Uses GLFW + OpenGL3 + Dear ImGui to display a
// texture that is periodically refreshed from the render's in-progress
// linear-radiance buffer, alongside a progress bar / status line.
//
// Deliberately independent from the Vulkan compute backend
// (include/gpu/vk_context.h): the GPU path tracer stays fully headless as
// before. This window just reads the ALREADY-COMPUTED accumulation buffer
// (its host-visible, persistently-mapped pointer - see
// include/gpu/gpu_buffer.h) or the CPU ImageBlock's snapshot from the host
// side, once per progress tick, and blits it through a completely separate
// GL context - zero coupling to how the pixels got produced.
//
// Only compiled/linked when M_ENABLE_PREVIEW_GUI is defined (GLFW/OpenGL
// found at CMake configure time - see the top-level CMakeLists.txt and
// src/gui/CMakeLists.txt). Callers that may run without it (i.e.
// include/render/progress.h) guard all use of this header behind the same
// macro.

#include <string>

namespace tiny_renderer::gui {

class PreviewWindow {
public:
    // `width`/`height` size the preview image (and initial window client
    // area); `title` is the OS window title. Throws std::runtime_error if
    // GLFW/OpenGL/ImGui setup fails (e.g. no display available).
    PreviewWindow(int width, int height, const std::string &title);
    ~PreviewWindow();

    PreviewWindow(const PreviewWindow &)            = delete;
    PreviewWindow &operator=(const PreviewWindow &) = delete;

    // Call once per progress tick (e.g. once per finished CPU block, or
    // once per completed GPU spp pass): pumps GLFW/OS events, re-uploads
    // `rgb` (linear radiance, row-major height*width*3 floats; may be
    // nullptr to skip the image refresh and just redraw the progress bar)
    // as the preview texture, and redraws the ImGui frame (a progress bar
    // at `progress` in [0,1] plus the `status` text line). Returns false
    // once the user has closed the window - the caller should stop calling
    // update() after that (rendering itself is NOT aborted).
    bool update(const float *rgb, int width, int height, float progress, const std::string &status);

    [[nodiscard]] bool is_open() const;

private:
    struct Impl;
    Impl *m_impl;
};

} // namespace tiny_renderer::gui

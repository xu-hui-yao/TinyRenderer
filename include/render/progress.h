#pragma once

// Progress reporting for the --progress CLI flag (see src/main/main.cpp):
// prints a live console progress bar (percentage / elapsed / ETA) for both
// the CPU and GPU render paths, and - when the optional GLFW+ImGui preview
// window is available (M_ENABLE_PREVIEW_GUI, see include/gui/preview_window.h
// and the top-level CMakeLists.txt) - also pops up a window showing the
// in-progress linear-radiance image, refreshed every time update() is
// called (once per finished CPU block, or once per completed GPU spp pass -
// see src/main/main.cpp and src/gpu/gpu_renderer.cpp for the call sites).
//
// Deliberately NOT a virtual/abstract interface: there is only ever one
// reporting policy (console + optional GUI), so a single concrete class is
// used directly by both render paths instead of introducing polymorphism.
//
// If M_ENABLE_PREVIEW_GUI is not defined (GLFW/OpenGL unavailable at CMake
// configure time), this class still works and behaves exactly the same
// minus the popup window - --progress then only prints the console bar.

#include <chrono>
#include <core/common.h>
#include <memory>
#include <string>

#ifdef M_ENABLE_PREVIEW_GUI
namespace tiny_renderer::gui {
class PreviewWindow; // fwd decl; full definition included where needed (see progress.cpp)
}
#endif

M_NAMESPACE_BEGIN

class ProgressReporter {
public:
    // `enabled` mirrors --progress: if false, every method below is a
    // no-op, so call sites can construct this unconditionally instead of
    // branching on the flag themselves. `width`/`height` size the optional
    // preview window's initial image; `label` is used both as the window
    // title and the "[label] " prefix on the console progress line (e.g.
    // "CPU", "GPU megakernel", "GPU wavefront").
    ProgressReporter(bool enabled, int width, int height, std::string label);
    ~ProgressReporter();

    ProgressReporter(const ProgressReporter &)            = delete;
    ProgressReporter &operator=(const ProgressReporter &) = delete;

    // Reports that `done` out of `total` units of work (CPU blocks, or GPU
    // spp passes) have completed. `rgb`, if non-null, is a row-major
    // rgb_height*rgb_width*3 linear-radiance snapshot of the image-so-far
    // and refreshes the preview window's texture; pass nullptr to just
    // update the progress bar/percentage text without touching the image.
    void update(int done, int total, const float *rgb = nullptr, int rgb_width = 0, int rgb_height = 0);

    // Prints a final newline on the console so subsequent output starts on
    // its own line. Safe (no-op) to call when disabled.
    void finish();

    // False once the user has closed the optional preview window (only
    // ever true if a window was actually created - i.e. GUI enabled AND
    // available). Rendering itself is never auto-cancelled by this; the
    // window simply stops refreshing once the caller stops calling
    // update() in response to this returning false.
    [[nodiscard]] bool window_open() const;

private:
    bool m_enabled;
    std::string m_label;
    std::chrono::steady_clock::time_point m_start;
    int m_last_percent = -1;
    // Only a data member when the GUI backend is actually being built for
    // THIS target (see the header doc comment above and progress.cpp) - not
    // just guarding the class' behavior, but its very layout, so that
    // translation units which never define M_ENABLE_PREVIEW_GUI (e.g. the
    // standalone src/gpu/*_main.cpp debug tools, which glob this very file
    // via src/gpu/CMakeLists.txt's M_CORE_SOURCES but never link
    // preview-gui-lib) never need gui::PreviewWindow's full definition.
#ifdef M_ENABLE_PREVIEW_GUI
    std::unique_ptr<gui::PreviewWindow> m_window;
#endif
};

M_NAMESPACE_END

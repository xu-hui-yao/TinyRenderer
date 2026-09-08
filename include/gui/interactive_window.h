#pragma once

// ============================================================================
// InteractiveWindow - GLFW + OpenGL3 + Dear ImGui window driving a real-time
// render loop (the --interactive CLI flag).
//
// Complements (and deliberately does NOT replace) PreviewWindow, which stays
// as the offline renderer's passive progress display. The difference is that
// this window OWNS the loop:
//
//   - it captures mouse/keyboard input each frame and feeds a CameraController
//   - it blits the GPU's RGBA8 output (already tonemapped on the GPU, so the
//     host does zero per-pixel work - see gpu_session.h)
//   - it hosts the settings panel and reports when a change invalidates the
//     temporal accumulation
//
// Like PreviewWindow it is fully decoupled from how pixels are produced: it is
// handed a pointer to RGBA8 data and knows nothing about Vulkan.
//
// Only compiled/linked when M_ENABLE_PREVIEW_GUI is defined (see
// src/gui/CMakeLists.txt) - interactive rendering additionally needs the GPU
// backend, which the top-level CMakeLists.txt checks.
// ============================================================================

#include <cstdint>
#include <gpu/gpu_timings.h>
#include <string>

namespace tiny_renderer::gui {

class CameraController;

// Mutable settings owned by the window and read by the render loop each frame.
struct InteractiveSettings {
    // ---- Render ----
    int spp_per_frame   = 1;    // 1..8
    int max_depth       = 3;    // 1..8
    int rr_depth        = 2;    // 1..8
    float radiance_clamp = 10.0f; // 0 disables
    bool demodulate     = true; // divide out albedo before filtering

    // ---- Denoise, temporal (SVGF) ----
    bool denoise       = true;
    bool clamp_history = true;  // rectify history into the current 3x3 range
    float clamp_growth = 0.25f; // how fast rectification relaxes with history length
    float alpha_color   = 0.05f; // minimum blend weight for new samples
    float alpha_moments = 0.05f;
    float phi_depth     = 1.0f;  // depth-rejection tolerance
    float phi_normal    = 0.5f;  // normal-rejection threshold, radians

    // ---- Denoise, spatial (a-trous) ----
    int atrous_iterations = 4;   // 0 disables the spatial pass entirely
    float phi_color    = 4.0f;   // luminance edge-stopping scale
    float phi_normal_a = 128.0f; // normal edge-stopping EXPONENT
    float phi_depth_a  = 1.0f;   // depth edge-stopping scale

    // ---- Display ----
    float exposure  = 1.0f;
    int tonemap     = 2;    // 0 = clamp, 1 = Reinhard, 2 = ACES
    bool use_srgb   = true; // matches Bitmap::save_png()
    float gamma     = 2.2f; // only when use_srgb == false
    int debug_view  = 0;    // 0 = final; see svgf_modulate.slang's DEBUG_*

    // ---- Camera ----
    int camera_mode = 0;    // 0 = Orbit, 1 = Fly
    float fov       = 30.0f;
    float move_speed = 1.0f;
};

class InteractiveWindow {
public:
    // `width`/`height` are the RENDER resolution; the window is clamped to fit
    // the display and the image is scaled to fit the window.
    InteractiveWindow(int width, int height, const std::string &title,
                      const InteractiveSettings &initial = {});
    ~InteractiveWindow();

    InteractiveWindow(const InteractiveWindow &)            = delete;
    InteractiveWindow &operator=(const InteractiveWindow &) = delete;

    // ---- Per-frame sequence ----

    // Pumps OS events, starts a new ImGui frame, and updates the FPS counter.
    // Returns false once the user closed the window.
    bool begin_frame();

    // Applies held mouse buttons / keys to `camera`. `dt` is seconds since the
    // previous frame. No-ops while the cursor is over an ImGui widget, so
    // dragging a slider doesn't also spin the camera.
    void update_camera(CameraController &camera, float dt);

    // Uploads and draws the RGBA8 render output as a full-window background.
    void draw_image(const uint8_t *rgba, int width, int height);

    // Draws the floating settings/stats panel. `accumulated_spp` is shown in
    // the stats section, as is the per-pass GPU timing breakdown when available.
    void draw_panel(uint32_t accumulated_spp, const gpu::PassTimings &timings);

    void end_frame();

    [[nodiscard]] bool is_open() const;

    InteractiveSettings &settings() { return m_settings; }

    // True if a render setting changed since the last call, meaning the
    // accumulation is stale and history must be reset. Display-only settings
    // (exposure/tonemap/gamma) deliberately do NOT set this.
    [[nodiscard]] bool consume_render_reset();

    [[nodiscard]] bool consume_view_reset();

    [[nodiscard]] float fps() const { return m_fps; }
    [[nodiscard]] float frame_time_ms() const { return m_frame_ms; }

    // Seconds elapsed since the previous begin_frame(). This is the value to
    // feed update_camera(), NOT frame_time_ms(): the latter is a display
    // readout smoothed over ~0.5s and stays exactly 0 until the first window
    // elapses, which made Fly-mode WASD dead on startup.
    [[nodiscard]] float delta_time() const;

private:
    struct Impl;
    Impl *m_impl;
    InteractiveSettings m_settings;
    bool m_render_reset = false;
    bool m_view_reset   = false;
    float m_fps         = 0.0f;
    float m_frame_ms    = 0.0f;
};

} // namespace tiny_renderer::gui

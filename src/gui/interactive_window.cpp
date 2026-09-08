#include <gui/interactive_window.h>

#include <gui/camera_controller.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>

#if defined(__APPLE__)
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl3.h>
#else
#include <GL/gl.h>
#endif

#include <GLFW/glfw3.h>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

namespace tiny_renderer::gui {

struct InteractiveWindow::Impl {
    GLFWwindow *window     = nullptr;
    GLuint texture         = 0;
    int tex_width          = 0;
    int tex_height         = 0;
    bool imgui_created     = false;
    bool has_last_cursor   = false;
    double last_cursor_x   = 0.0;
    double last_cursor_y   = 0.0;
    double last_time       = 0.0;
    double delta_time      = 0.0; // seconds since the previous begin_frame()
    double fps_accum       = 0.0;
    int fps_frames         = 0;
    // Wheel input arrives via a GLFW callback (there is no poll-based API for
    // it), so the callback accumulates here and update_camera() drains it.
    double scroll_delta    = 0.0;

    // Declared inside Impl (rather than as a free function) so it can name the
    // private nested type without needing befriending.
    static void scroll_callback(GLFWwindow *window, double xoffset, double yoffset) {
        // Chain to ImGui's own handler first. ImGui_ImplGlfw_InitForOpenGL was
        // asked to install callbacks, and registering ours afterwards REPLACES
        // its scroll callback - without forwarding, scrolling inside the
        // Controls panel does nothing.
        ImGui_ImplGlfw_ScrollCallback(window, xoffset, yoffset);

        if (auto *impl = static_cast<Impl *>(glfwGetWindowUserPointer(window)))
            impl->scroll_delta += yoffset;
    }
};

namespace {

// Exponential moving average over ~0.5s so the readout is stable but responsive.
constexpr double FPS_SMOOTH_WINDOW = 0.5;

int clamp_window_dim(int value, int lo, int hi) { return std::clamp(value, lo, hi); }

} // namespace

InteractiveWindow::InteractiveWindow(int width, int height, const std::string &title,
                                     const InteractiveSettings &initial)
    : m_impl(new Impl()), m_settings(initial) {

    if (!glfwInit())
        throw std::runtime_error("InteractiveWindow: glfwInit() failed (no display available?)");

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);

    // Fit the window to the display: a 4K render should still open a usable
    // window. The image is scaled to fit in draw_image() regardless.
    int win_w = width, win_h = height;
    if (GLFWmonitor *monitor = glfwGetPrimaryMonitor()) {
        const GLFWvidmode *mode = glfwGetVideoMode(monitor);
        if (mode) {
            float scale = std::min(1.0f, std::min(static_cast<float>(mode->width) * 0.9f / width,
                                                  static_cast<float>(mode->height) * 0.9f / height));
            win_w = clamp_window_dim(static_cast<int>(width * scale), 320, mode->width);
            win_h = clamp_window_dim(static_cast<int>(height * scale), 240, mode->height);
        }
    }

    m_impl->window = glfwCreateWindow(win_w, win_h, title.c_str(), nullptr, nullptr);
    if (!m_impl->window) {
        glfwTerminate();
        throw std::runtime_error("InteractiveWindow: glfwCreateWindow() failed");
    }

    glfwMakeContextCurrent(m_impl->window);
    // vsync ON: the render loop is GPU-bound and there is no benefit to
    // spinning faster than the display; it also keeps input latency low.
    glfwSwapInterval(1);

    // Deliberately NOT GLFW_STICKY_MOUSE_BUTTONS: the camera polls the held
    // state with glfwGetMouseButton() every frame, and sticky mode latches a
    // release so the next poll still reports PRESS - which drags the view one
    // extra frame after the button is let go.
    glfwSetInputMode(m_impl->window, GLFW_STICKY_MOUSE_BUTTONS, GLFW_FALSE);

    glfwSetWindowUserPointer(m_impl->window, m_impl);
    glfwSetScrollCallback(m_impl->window, &Impl::scroll_callback);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();

    if (!ImGui_ImplGlfw_InitForOpenGL(m_impl->window, true) || !ImGui_ImplOpenGL3_Init("#version 150")) {
        ImGui::DestroyContext();
        glfwDestroyWindow(m_impl->window);
        glfwTerminate();
        throw std::runtime_error("InteractiveWindow: Dear ImGui GLFW/OpenGL3 backend init failed");
    }
    m_impl->imgui_created = true;

    // Don't let ImGui persist layout to disk - interacting with the panel
    // shouldn't leave state behind between runs.
    ImGui::GetIO().IniFilename = nullptr;

    glGenTextures(1, &m_impl->texture);
    glBindTexture(GL_TEXTURE_2D, m_impl->texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    m_impl->last_time = glfwGetTime();
}

InteractiveWindow::~InteractiveWindow() {
    if (m_impl) {
        if (m_impl->texture)
            glDeleteTextures(1, &m_impl->texture);
        if (m_impl->imgui_created) {
            ImGui_ImplOpenGL3_Shutdown();
            ImGui_ImplGlfw_Shutdown();
            ImGui::DestroyContext();
        }
        if (m_impl->window)
            glfwDestroyWindow(m_impl->window);
        glfwTerminate();
        delete m_impl;
    }
}

bool InteractiveWindow::is_open() const { return m_impl->window && !glfwWindowShouldClose(m_impl->window); }

float InteractiveWindow::delta_time() const { return static_cast<float>(m_impl->delta_time); }

bool InteractiveWindow::begin_frame() {
    if (!m_impl->window || glfwWindowShouldClose(m_impl->window))
        return false;

    glfwMakeContextCurrent(m_impl->window);
    glfwPollEvents();

    double now = glfwGetTime();
    double dt  = now - m_impl->last_time;
    m_impl->last_time = now;

    // Guard against absurd dt (e.g. the first frame, or after a stall in a
    // debugger) which would make both the camera motion and the stats jump.
    dt = std::min(dt, 0.25);
    m_impl->delta_time = dt;

    m_impl->fps_accum += dt;
    m_impl->fps_frames++;
    if (m_impl->fps_accum >= FPS_SMOOTH_WINDOW) {
        m_fps = static_cast<float>(m_impl->fps_frames / m_impl->fps_accum);
        m_frame_ms = static_cast<float>(m_impl->fps_accum * 1000.0 / m_impl->fps_frames);
        m_impl->fps_accum = 0.0;
        m_impl->fps_frames = 0;
    }

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    return true;
}

void InteractiveWindow::update_camera(CameraController &camera, float dt) {
    const ImGuiIO &io = ImGui::GetIO();

    // Never drive the camera while the user is interacting with the panel, or
    // dragging a slider would spin the view.
    if (io.WantCaptureMouse) {
        m_impl->has_last_cursor = false;
        // Drop, rather than keep, wheel input received over the panel. Leaving
        // it queued means the accumulated notches all fire at once the moment
        // the cursor leaves the panel, snapping the zoom.
        m_impl->scroll_delta = 0.0;
        return;
    }

    double x = 0.0, y = 0.0;
    glfwGetCursorPos(m_impl->window, &x, &y);

    if (m_impl->has_last_cursor) {
        float dx = static_cast<float>(x - m_impl->last_cursor_x);
        float dy = static_cast<float>(y - m_impl->last_cursor_y);

        if (glfwGetMouseButton(m_impl->window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
            camera.orbit(dx, dy);
        }
        if (glfwGetMouseButton(m_impl->window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS ||
            glfwGetMouseButton(m_impl->window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
            camera.pan(dx, dy);
        }
    }
    if (m_impl->scroll_delta != 0.0) {
        camera.dolly(static_cast<float>(m_impl->scroll_delta));
        m_impl->scroll_delta = 0.0;
    }

    m_impl->last_cursor_x = x;
    m_impl->last_cursor_y = y;
    m_impl->has_last_cursor = true;

    // Fly-mode translation. Shift to sprint.
    // Skipped while ImGui owns the keyboard, so typing into a widget does not
    // also fly the camera.
    if (camera.mode() == CameraController::Mode::Fly && !io.WantCaptureKeyboard) {
        float sprint = (glfwGetKey(m_impl->window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                        glfwGetKey(m_impl->window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS)
                           ? 4.0f
                           : 1.0f;
        auto held = [this](int key) { return glfwGetKey(m_impl->window, key) == GLFW_PRESS; };

        float right   = (held(GLFW_KEY_D) ? 1.0f : 0.0f) - (held(GLFW_KEY_A) ? 1.0f : 0.0f);
        float up      = (held(GLFW_KEY_E) ? 1.0f : 0.0f) - (held(GLFW_KEY_Q) ? 1.0f : 0.0f);
        float forward = (held(GLFW_KEY_W) ? 1.0f : 0.0f) - (held(GLFW_KEY_S) ? 1.0f : 0.0f);
        if (right != 0.0f || up != 0.0f || forward != 0.0f)
            camera.move_local(right, up, forward, dt * sprint);
    }
}

void InteractiveWindow::draw_image(const uint8_t *rgba, int width, int height) {
    if (!rgba || width <= 0 || height <= 0)
        return;

    glBindTexture(GL_TEXTURE_2D, m_impl->texture);
    if (width != m_impl->tex_width || height != m_impl->tex_height) {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgba);
        m_impl->tex_width  = width;
        m_impl->tex_height = height;
    } else {
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, rgba);
    }

    // Full-window background quad, with the panel drawn on top afterwards.
    const ImGuiIO &io = ImGui::GetIO();
    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::SetNextWindowSize(io.DisplaySize);
    ImGui::Begin("viewport", nullptr,
                 ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                     ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoScrollbar |
                     ImGuiWindowFlags_NoBringToFrontOnFocus |
                     // NoInputs is not cosmetic - it is what makes the camera
                     // controllable at all. This window covers the ENTIRE
                     // client area, so without it ImGui reports
                     // io.WantCaptureMouse == true for every cursor position,
                     // and update_camera() (which correctly refuses to drive
                     // the view while the user is on a widget) then bails on
                     // every single frame: drags, wheel and WASD all do
                     // nothing. With it, only the Controls panel captures.
                     ImGuiWindowFlags_NoInputs);

    if (m_impl->tex_width > 0 && m_impl->tex_height > 0) {
        float avail_x = io.DisplaySize.x, avail_y = io.DisplaySize.y;
        float aspect  = static_cast<float>(m_impl->tex_width) / static_cast<float>(m_impl->tex_height);
        float w = avail_x, h = w / aspect;
        if (h > avail_y) {
            h = avail_y;
            w = h * aspect;
        }
        ImGui::SetCursorPos(ImVec2((avail_x - w) * 0.5f, (avail_y - h) * 0.5f));
        // ImGui images are UV-origin top-left; our buffer's first row is also
        // the top row, so no flip is needed.
        ImGui::Image(static_cast<ImTextureID>(static_cast<intptr_t>(m_impl->texture)), ImVec2(w, h));
    }
    ImGui::End();
}

void InteractiveWindow::draw_panel(uint32_t accumulated_spp, const gpu::PassTimings &timings) {
    ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(330, 470), ImGuiCond_FirstUseEver);
    ImGui::Begin("Controls");

    // ---- Stats (top: the thing you watch while tuning) ----
    if (ImGui::CollapsingHeader("Stats", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Text("FPS: %.1f (%.2f ms)", static_cast<double>(m_fps), static_cast<double>(m_frame_ms));
        ImGui::Text("Accumulated spp: %u", accumulated_spp);

        // Per-pass GPU breakdown. This is the only way to tell "the path tracer
        // is the bottleneck" from "the denoiser is", which matters because the
        // fixes are completely different (BVH/triangle work vs filter params).
        if (timings.available) {
            ImGui::Separator();
            const float total = timings.total_ms;
            auto row = [&](const char *name, float v) {
                const float pct = (total > 0.0f) ? (v / total * 100.0f) : 0.0f;
                ImGui::Text("%-11s %6.2f ms  %4.1f%%", name, static_cast<double>(v), static_cast<double>(pct));
            };
            row("path trace", timings.path_trace_ms);
            row("temporal", timings.temporal_ms);
            row("a-trous", timings.atrous_ms);
            row("display", timings.modulate_ms);
            ImGui::Text("total        %6.2f ms", static_cast<double>(total));
        } else {
            ImGui::TextDisabled("(GPU timings unavailable)");
        }
        ImGui::Separator();
    }

    // ---- Camera ----
    if (ImGui::CollapsingHeader("Camera", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::Combo("Mode", &m_settings.camera_mode, "Orbit\0Fly\0"))
            m_view_reset = true;
        if (ImGui::SliderFloat("FOV", &m_settings.fov, 5.0f, 120.0f, "%.1f"))
            m_render_reset = true;
        if (ImGui::SliderFloat("Move speed", &m_settings.move_speed, 0.05f, 8.0f, "%.2f")) {
            // motion-only, no reset needed
        }
        if (ImGui::Button("Reset view"))
            m_view_reset = true;
        ImGui::TextDisabled("LMB drag = orbit, MMB/RMB = pan, wheel = zoom");
        ImGui::TextDisabled("Fly mode: WASD/QE, Shift = sprint");
    }

    // ---- Render ----
    if (ImGui::CollapsingHeader("Render", ImGuiTreeNodeFlags_DefaultOpen)) {
        if (ImGui::SliderInt("spp / frame", &m_settings.spp_per_frame, 1, 8))
            m_render_reset = true;
        if (ImGui::SliderInt("Max depth", &m_settings.max_depth, 1, 8))
            m_render_reset = true;
        if (ImGui::SliderInt("RR depth", &m_settings.rr_depth, 1, 8))
            m_render_reset = true;
        if (ImGui::SliderFloat("Radiance clamp", &m_settings.radiance_clamp, 0.0f, 50.0f, "%.1f"))
            m_render_reset = true;
        ImGui::TextDisabled("clamp = 0 disables");
        // Demodulation changes what the stored illumination MEANS, so existing
        // history is no longer comparable and has to be discarded.
        if (ImGui::Checkbox("Demodulate albedo", &m_settings.demodulate))
            m_render_reset = true;
        if (ImGui::Button("Reset accumulation"))
            m_render_reset = true;
    }

    // ---- Denoise (SVGF) ----
    // None of these reset the accumulation either: they only change how the
    // NEXT frame blends in, so the EMA absorbs them over a frame or two.
    if (ImGui::CollapsingHeader("Denoise", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Checkbox("Enable temporal", &m_settings.denoise);
        ImGui::Checkbox("Demodulate albedo", &m_settings.demodulate);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Divide out albedo before filtering, multiply back at display.\n"
                              "Off = filter raw radiance (blurs texture detail).");
        ImGui::Checkbox("Clamp history", &m_settings.clamp_history);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Rectify history into the current frame's 3x3 range.\nRemoves ghosting and residual fireflies.");
        ImGui::SliderFloat("Alpha color", &m_settings.alpha_color, 0.01f, 1.0f, "%.3f");
        ImGui::SliderFloat("Alpha moments", &m_settings.alpha_moments, 0.01f, 1.0f, "%.3f");
        ImGui::SliderFloat("Phi depth", &m_settings.phi_depth, 0.1f, 20.0f, "%.2f");
        ImGui::SliderFloat("Phi normal", &m_settings.phi_normal, 0.05f, 1.5f, "%.3f");
        ImGui::SliderFloat("Clamp growth", &m_settings.clamp_growth, 0.0f, 2.0f, "%.3f");
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("How fast history rectification relaxes as history grows.\n"
                              "0 = fixed radius forever (stops convergence);\n"
                              "higher = lets a static view keep refining.");

        ImGui::SeparatorText("Spatial (a-trous)");
        ImGui::SliderInt("Iterations", &m_settings.atrous_iterations, 0, 6);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Each doubles the gather stride. 0 disables the spatial pass,\n"
                              "which is the quickest way to see what the temporal pass alone gives.");
        ImGui::SliderFloat("Phi color (spatial)", &m_settings.phi_color, 0.1f, 32.0f, "%.2f");
        ImGui::SliderFloat("Phi normal (spatial)", &m_settings.phi_normal_a, 1.0f, 512.0f, "%.1f");
        ImGui::SliderFloat("Phi depth (spatial)", &m_settings.phi_depth_a, 0.01f, 8.0f, "%.3f");
    }

    // ---- Display ----
    // None of these invalidate the accumulation: they are applied by the
    // modulate pass over already-accumulated radiance, so changing exposure is
    // free and instant even at 1000 spp.
    if (ImGui::CollapsingHeader("Display")) {
        ImGui::SliderFloat("Exposure", &m_settings.exposure, -4.0f, 4.0f, "%.2f");
        ImGui::Combo("Tonemap", &m_settings.tonemap, "Clamp\0Reinhard\0ACES\0");
        ImGui::Checkbox("sRGB transfer", &m_settings.use_srgb);
        if (!m_settings.use_srgb)
            ImGui::SliderFloat("Gamma", &m_settings.gamma, 1.0f, 3.0f, "%.2f");

        // The debug views are how the filter actually gets tuned - SVGF's
        // parameters are strongly coupled, and "wrong settings" looks identical
        // to "wrong G-Buffer" until you look at the intermediates.
        ImGui::Combo("Debug view", &m_settings.debug_view,
                     "Final\0Illum (raw)\0Illum (filtered)\0Albedo\0Normal\0Depth\0Variance\0Hist len\0Position\0");
    }

    ImGui::End();
}

void InteractiveWindow::end_frame() {
    ImGui::Render();

    int fb_w = 0, fb_h = 0;
    glfwGetFramebufferSize(m_impl->window, &fb_w, &fb_h);
    glViewport(0, 0, fb_w, fb_h);
    glClearColor(0.05f, 0.05f, 0.06f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    glfwSwapBuffers(m_impl->window);
}

bool InteractiveWindow::consume_render_reset() {
    bool r = m_render_reset;
    m_render_reset = false;
    return r;
}

bool InteractiveWindow::consume_view_reset() {
    bool r = m_view_reset;
    m_view_reset = false;
    return r;
}

} // namespace tiny_renderer::gui

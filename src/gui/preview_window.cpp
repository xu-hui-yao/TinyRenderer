#include <gui/preview_window.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#if defined(__APPLE__)
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl3.h>
#else
// The Windows SDK GL/gl.h needs windows.h first (WINGDIAPI/APIENTRY); glfw3.h
// does not pull it in for us.
#if defined(_WIN32)
#include <windows.h>
#endif
#include <GL/gl.h>

// The Windows SDK only ships OpenGL 1.1 headers, so this 1.2 constant is
// missing there.
#ifndef GL_CLAMP_TO_EDGE
#define GL_CLAMP_TO_EDGE 0x812F
#endif
#endif

#include <GLFW/glfw3.h>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

namespace tiny_renderer::gui {

struct PreviewWindow::Impl {
    GLFWwindow *window   = nullptr;
    GLuint texture        = 0;
    int tex_width         = 0;
    int tex_height        = 0;
    bool imgui_created    = false;
    std::vector<unsigned char> ldr_scratch; // reused across update() calls to avoid reallocating every tick
};

namespace {
// Converts a linear-radiance RGB buffer to display-ready RGBA8 (simple
// clamp + 1/2.2 gamma - purely a PREVIEW tonemap, independent of
// Bitmap::save_png()'s ToneMapMode used for the actual saved output).
void convert_to_ldr(const float *rgb, int width, int height, std::vector<unsigned char> &out) {
    out.resize(static_cast<size_t>(width) * height * 4);
    for (int i = 0; i < width * height; ++i) {
        for (int c = 0; c < 3; ++c) {
            float v                 = rgb[i * 3 + c];
            v                        = std::clamp(v, 0.0f, 1.0f);
            v                        = std::pow(v, 1.0f / 2.2f);
            out[i * 4 + c]           = static_cast<unsigned char>(v * 255.0f + 0.5f);
        }
        out[i * 4 + 3] = 255;
    }
}
} // namespace

PreviewWindow::PreviewWindow(int width, int height, const std::string &title) : m_impl(new Impl()) {
    if (!glfwInit()) {
        delete m_impl;
        throw std::runtime_error("PreviewWindow: glfwInit() failed (no display available?)");
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);

    // Clamp the initial window size so very large renders (e.g. 4K) still
    // open a window that fits on a normal display; the image itself is
    // scaled to fit inside the window in update() regardless.
    int win_w = std::clamp(width, 320, 1600);
    int win_h = std::clamp(height, 240, 1000) + 70; // +70 for the status/progress-bar text above the image

    m_impl->window = glfwCreateWindow(win_w, win_h, title.c_str(), nullptr, nullptr);
    if (!m_impl->window) {
        glfwTerminate();
        delete m_impl;
        throw std::runtime_error("PreviewWindow: glfwCreateWindow() failed");
    }

    glfwMakeContextCurrent(m_impl->window);
    glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();

    if (!ImGui_ImplGlfw_InitForOpenGL(m_impl->window, true) || !ImGui_ImplOpenGL3_Init("#version 150")) {
        ImGui::DestroyContext();
        glfwDestroyWindow(m_impl->window);
        glfwTerminate();
        delete m_impl;
        throw std::runtime_error("PreviewWindow: Dear ImGui GLFW/OpenGL3 backend init failed");
    }
    m_impl->imgui_created = true;

    glGenTextures(1, &m_impl->texture);
    glBindTexture(GL_TEXTURE_2D, m_impl->texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}

PreviewWindow::~PreviewWindow() {
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

bool PreviewWindow::is_open() const { return m_impl->window && !glfwWindowShouldClose(m_impl->window); }

bool PreviewWindow::update(const float *rgb, int width, int height, float progress, const std::string &status) {
    if (glfwWindowShouldClose(m_impl->window))
        return false;

    glfwMakeContextCurrent(m_impl->window);
    glfwPollEvents();

    if (rgb != nullptr && width > 0 && height > 0) {
        convert_to_ldr(rgb, width, height, m_impl->ldr_scratch);
        glBindTexture(GL_TEXTURE_2D, m_impl->texture);
        if (width != m_impl->tex_width || height != m_impl->tex_height) {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE,
                        m_impl->ldr_scratch.data());
            m_impl->tex_width  = width;
            m_impl->tex_height = height;
        } else {
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE,
                            m_impl->ldr_scratch.data());
        }
    }

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    const ImGuiIO &io = ImGui::GetIO();
    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::SetNextWindowSize(io.DisplaySize);
    ImGui::Begin("preview", nullptr,
                 ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                     ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings);

    ImGui::TextUnformatted(status.c_str());
    ImGui::ProgressBar(progress, ImVec2(-1, 0));

    ImVec2 avail = ImGui::GetContentRegionAvail();
    if (m_impl->tex_width > 0 && m_impl->tex_height > 0 && avail.x > 1 && avail.y > 1) {
        float aspect = static_cast<float>(m_impl->tex_width) / static_cast<float>(m_impl->tex_height);
        float w = avail.x, h = w / aspect;
        if (h > avail.y) {
            h = avail.y;
            w = h * aspect;
        }
        ImGui::SetCursorPosX((avail.x - w) * 0.5f + ImGui::GetCursorPosX());
        ImGui::Image(static_cast<ImTextureID>(static_cast<intptr_t>(m_impl->texture)), ImVec2(w, h));
    }

    ImGui::End();
    ImGui::Render();

    int fb_w, fb_h;
    glfwGetFramebufferSize(m_impl->window, &fb_w, &fb_h);
    glViewport(0, 0, fb_w, fb_h);
    glClearColor(0.08f, 0.08f, 0.09f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    glfwSwapBuffers(m_impl->window);

    return !glfwWindowShouldClose(m_impl->window);
}

} // namespace tiny_renderer::gui

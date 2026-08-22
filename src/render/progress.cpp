#include <render/progress.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>

#ifdef M_ENABLE_PREVIEW_GUI
#include <gui/preview_window.h>
#endif

M_NAMESPACE_BEGIN

namespace {
std::string format_duration(double seconds) {
    if (!std::isfinite(seconds) || seconds < 0)
        seconds = 0;
    int total = static_cast<int>(seconds + 0.5);
    int h     = total / 3600;
    int m     = (total % 3600) / 60;
    int s     = total % 60;
    char buf[32];
    if (h > 0)
        std::snprintf(buf, sizeof(buf), "%dh%02dm%02ds", h, m, s);
    else if (m > 0)
        std::snprintf(buf, sizeof(buf), "%dm%02ds", m, s);
    else
        std::snprintf(buf, sizeof(buf), "%ds", s);
    return std::string(buf);
}
} // namespace

ProgressReporter::ProgressReporter(bool enabled, int width, int height, std::string label)
    : m_enabled(enabled), m_label(std::move(label)), m_start(std::chrono::steady_clock::now()) {
    if (!m_enabled)
        return;

#ifdef M_ENABLE_PREVIEW_GUI
    try {
        m_window = std::make_unique<gui::PreviewWindow>(width, height, "TinyRenderer - " + m_label);
    } catch (const std::exception &e) {
        std::cerr << "[progress] Live preview window unavailable (" << e.what()
                  << "); continuing with console progress only." << std::endl;
        m_window.reset();
    }
#else
    (void)width;
    (void)height;
#endif
}

ProgressReporter::~ProgressReporter() = default;

void ProgressReporter::update(int done, int total, const float *rgb, int rgb_width, int rgb_height) {
    if (!m_enabled)
        return;

    total       = std::max(total, 1);
    done        = std::clamp(done, 0, total);
    float frac  = static_cast<float>(done) / static_cast<float>(total);
    int percent = static_cast<int>(frac * 100.0f + 0.5f);

    double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - m_start).count();
    double eta     = (frac > 0.0f) ? elapsed * (1.0 / frac - 1.0) : 0.0;

    std::string status = "[" + m_label + "] " + std::to_string(percent) + "% (" + std::to_string(done) + "/" +
                         std::to_string(total) + ") elapsed " + format_duration(elapsed) + ", ETA " +
                         format_duration(eta);

    if (percent != m_last_percent || done == total) {
        m_last_percent = percent;
        std::cout << "\r" << status << "        " << std::flush;
    }

#ifdef M_ENABLE_PREVIEW_GUI
    if (m_window) {
        if (!m_window->update(rgb, rgb_width, rgb_height, frac, status)) {
            // User closed the window: stop trying to refresh it further
            // (rendering itself keeps going - see the header doc comment).
            m_window.reset();
        }
    }
#else
    (void)rgb;
    (void)rgb_width;
    (void)rgb_height;
#endif
}

void ProgressReporter::finish() {
    if (!m_enabled)
        return;
    std::cout << std::endl;
}

bool ProgressReporter::window_open() const {
#ifdef M_ENABLE_PREVIEW_GUI
    return static_cast<bool>(m_window);
#else
    return false;
#endif
}

M_NAMESPACE_END

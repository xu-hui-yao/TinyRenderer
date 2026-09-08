#pragma once

// Per-pass GPU timings, shared between the renderer (which measures them) and
// the GUI (which displays them).
//
// Its own header so the GUI can depend on it without pulling in Vulkan or VMA:
// preview-gui-lib does not have those include paths, and gpu_session.h needs
// both.

namespace tiny_renderer::gpu {

struct PassTimings {
    float path_trace_ms = 0.0f;
    float temporal_ms   = 0.0f;
    float atrous_ms     = 0.0f;
    float modulate_ms   = 0.0f;
    float total_ms      = 0.0f;
    // False when the device exposes no timestamp queries; all fields stay 0.
    bool available = false;
};

} // namespace tiny_renderer::gpu

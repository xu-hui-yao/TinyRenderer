#pragma once

#include <algorithm>
#include <components/bitmap.h>
#include <render/framebuffer.h>
#include <thread>
#include <vector>

M_NAMESPACE_BEGIN

/**
 * \brief Small helpers shared by the denoiser implementations in
 * src/denoisers/.
 *
 * Everything here operates on plain float buffers in row-major
 * `(y * cols + x) * channels + c` layout - the same layout \ref TTensor (and
 * therefore \ref Bitmap) uses - so a denoiser can read straight out of a
 * Bitmap's data pointer without any copying.
 */
namespace denoise {

// Floor for the albedo used as a demodulation divisor. Pure-black surfaces
// (albedo 0) carry no signal to recover anyway, and dividing by them would
// manufacture infinities; 1e-3 keeps the round trip exact to well below 8-bit
// quantization while staying numerically harmless.
constexpr float k_albedo_floor = 1e-3f;

// ---------------------------------------------------------------------------
// Parallelism
// ---------------------------------------------------------------------------

/**
 * \brief Run `body(y_begin, y_end)` over a partition of `[0, rows)` across
 * `thread_count` threads.
 *
 * Every denoiser here is a *gather* operation: each output pixel reads a
 * neighborhood of the (immutable) input and writes only itself. Splitting the
 * output by horizontal bands is therefore embarrassingly parallel with no
 * locking and no false sharing beyond band boundaries.
 */
template <typename Body> void parallel_bands(int rows, int thread_count, const Body &body) {
    thread_count = std::max(1, thread_count);
    if (thread_count == 1 || rows <= 1) {
        body(0, rows);
        return;
    }

    thread_count = std::min(thread_count, rows);

    std::vector<std::thread> threads;
    threads.reserve(thread_count);
    int base = rows / thread_count, extra = rows % thread_count;
    int y = 0;
    for (int t = 0; t < thread_count; ++t) {
        int count = base + (t < extra ? 1 : 0);
        threads.emplace_back([&body, y, count] { body(y, y + count); });
        y += count;
    }
    for (auto &th : threads)
        th.join();
}

// ---------------------------------------------------------------------------
// Buffer access
// ---------------------------------------------------------------------------

inline int clamp_i(int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); }

// Read `channels`-channel buffer at a position clamped into the image, which
// extends the image by edge replication. Cheaper and visually better behaved
// at the border than renormalizing a truncated filter footprint.
inline float fetch(const std::vector<float> &buf, int rows, int cols, int channels, int y, int x, int c) {
    y = clamp_i(y, 0, rows - 1);
    x = clamp_i(x, 0, cols - 1);
    return buf[(static_cast<size_t>(y) * cols + x) * channels + c];
}

// Copy a Bitmap into a flat float vector. `channels` must match the bitmap.
inline std::vector<float> to_vector(const Bitmap &bitmap, int channels) {
    const int rows = bitmap.get_rows(), cols = bitmap.get_cols();
    std::vector<float> out(static_cast<size_t>(rows) * cols * channels);
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < cols; ++x)
            for (int c = 0; c < channels; ++c)
                out[(static_cast<size_t>(y) * cols + x) * channels + c] = bitmap(y, x, c);
    return out;
}

inline std::shared_ptr<Bitmap> to_bitmap(const std::vector<float> &buf, int rows, int cols) {
    auto out = std::make_shared<Bitmap>(rows, cols, 3);
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < cols; ++x)
            for (int c = 0; c < 3; ++c)
                (*out)(y, x, c) = buf[(static_cast<size_t>(y) * cols + x) * 3 + c];
    return out;
}

// ---------------------------------------------------------------------------
// Separable box blur (used to prefilter variance and error maps)
// ---------------------------------------------------------------------------

inline void box_blur(const std::vector<float> &src, std::vector<float> &dst, int rows, int cols, int channels,
                     int radius) {
    if (radius <= 0) {
        dst = src;
        return;
    }

    dst.assign(src.size(), 0.0f);
    std::vector<float> tmp(src.size(), 0.0f);
    const float inv = 1.0f / static_cast<float>(2 * radius + 1);

    // Horizontal pass
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < cols; ++x)
            for (int c = 0; c < channels; ++c) {
                float sum = 0.0f;
                for (int d = -radius; d <= radius; ++d)
                    sum += src[(static_cast<size_t>(y) * cols + clamp_i(x + d, 0, cols - 1)) * channels + c];
                tmp[(static_cast<size_t>(y) * cols + x) * channels + c] = sum * inv;
            }

    // Vertical pass
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < cols; ++x)
            for (int c = 0; c < channels; ++c) {
                float sum = 0.0f;
                for (int d = -radius; d <= radius; ++d)
                    sum += tmp[(static_cast<size_t>(clamp_i(y + d, 0, rows - 1)) * cols + x) * channels + c];
                dst[(static_cast<size_t>(y) * cols + x) * channels + c] = sum * inv;
            }
}

// ---------------------------------------------------------------------------
// Variance
// ---------------------------------------------------------------------------

/**
 * \brief Per-pixel, per-channel variance of `color` (the full-sample mean).
 *
 * Prefers the precomputed `input.variance` buffer. If that is missing but the
 * half-buffers are present, recompute 1/4 (A - B)^2 with a 3x3 box prefilter
 * (see FrameBufferSet::variance for the derivation). If neither is available -
 * e.g. a GPU render, or an integrator without li_aov() support - fall back to
 * a local spatial variance estimate over a 3x3 window. The latter conflates
 * genuine image detail with noise and is therefore markedly worse, but it lets
 * the feature-guided filters still run instead of failing outright.
 */
inline std::vector<float> estimate_variance(const DenoiseInput &input, int rows, int cols) {
    std::vector<float> var(static_cast<size_t>(rows) * cols * 3, 0.0f);

    if (input.variance) {
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < cols; ++x)
                for (int c = 0; c < 3; ++c)
                    var[(static_cast<size_t>(y) * cols + x) * 3 + c] = (*input.variance)(y, x, c);
        return var;
    }

    if (input.has_half_buffers()) {
        std::vector<float> raw(static_cast<size_t>(rows) * cols * 3);
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < cols; ++x)
                for (int c = 0; c < 3; ++c) {
                    float d = (*input.color_a)(y, x, c) - (*input.color_b)(y, x, c);
                    raw[(static_cast<size_t>(y) * cols + x) * 3 + c] = 0.25f * d * d;
                }
        box_blur(raw, var, rows, cols, 3, 1);
        return var;
    }

    // Last resort: local sample variance of the color buffer itself.
    const auto &color = *input.color;
    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < cols; ++x) {
            for (int c = 0; c < 3; ++c) {
                float sum = 0.0f, sum_sq = 0.0f;
                int n = 0;
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        int yy = clamp_i(y + dy, 0, rows - 1), xx = clamp_i(x + dx, 0, cols - 1);
                        float v = color(yy, xx, c);
                        sum += v;
                        sum_sq += v * v;
                        ++n;
                    }
                }
                float mean = sum / static_cast<float>(n);
                float v    = sum_sq / static_cast<float>(n) - mean * mean;
                var[(static_cast<size_t>(y) * cols + x) * 3 + c] = v > 0.0f ? v : 0.0f;
            }
        }
    }
    return var;
}

// ---------------------------------------------------------------------------
// Albedo demodulation
// ---------------------------------------------------------------------------

/**
 * \brief In-place `color /= albedo`, turning radiance into (approximate)
 * irradiance.
 *
 * The point is that irradiance is smooth and low-frequency while the albedo
 * carries all the high-frequency texture detail. Filtering the former and
 * multiplying the latter back in afterwards (\ref modulate) therefore removes
 * noise without ever touching texture - which a filter operating directly on
 * radiance cannot do, since it has no way to tell a texture edge apart from a
 * noise spike.
 *
 * `variance` must be scaled consistently: Var(c/a) = Var(c)/a^2.
 */
inline void demodulate(std::vector<float> &color, std::vector<float> &variance, const Bitmap &albedo, int rows,
                       int cols) {
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < cols; ++x)
            for (int c = 0; c < 3; ++c) {
                float a         = std::max(albedo(y, x, c), k_albedo_floor);
                size_t i        = (static_cast<size_t>(y) * cols + x) * 3 + c;
                color[i]    /= a;
                variance[i] /= a * a;
            }
}

inline void modulate(std::vector<float> &color, const Bitmap &albedo, int rows, int cols) {
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < cols; ++x)
            for (int c = 0; c < 3; ++c) {
                float a  = std::max(albedo(y, x, c), k_albedo_floor);
                size_t i = (static_cast<size_t>(y) * cols + x) * 3 + c;
                color[i] *= a;
            }
}

} // namespace denoise

M_NAMESPACE_END

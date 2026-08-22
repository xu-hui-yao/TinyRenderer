#include <algorithm>
#include <cmath>
#include <components/denoiser.h>
#include <render/denoise_utils.h>
#include <thread>

M_NAMESPACE_BEGIN

/**
 * \brief Edge-avoiding a-trous wavelet transform, i.e. a joint (cross)
 * bilateral filter applied iteratively with an exponentially growing stride.
 *
 * This is the spatial half of SVGF (Schied et al. 2017), which in turn builds
 * on Dammertz et al. 2010's edge-avoiding a-trous wavelet transform, and it is
 * the de-facto standard for real-time Monte Carlo denoising.
 *
 * ## Why a-trous instead of one large joint bilateral kernel
 *
 * A plain joint bilateral filter wide enough to remove low-frequency noise
 * (say 65x65) costs ~4225 taps per pixel. The a-trous decomposition gets the
 * same reach by applying a small 5x5 kernel five times, doubling the tap
 * SPACING each iteration (1, 2, 4, 8, 16) while keeping the tap COUNT at 25:
 *
 *     iteration 0:  x x x x x        (stride 1,  footprint  5)
 *     iteration 1:  x _ x _ x ...    (stride 2,  footprint  9)
 *     iteration 4:                   (stride 16, footprint 65)
 *
 * That is 125 taps instead of 4225 - roughly 34x cheaper - and because each
 * iteration re-evaluates the edge-stopping weights on the PROGRESSIVELY
 * SMOOTHED image, the weights get less noisy (and hence more reliable) as the
 * footprint grows, which a single-pass filter cannot do.
 *
 * ## The edge-stopping weights
 *
 * For center pixel p and tap q, the weight is the product of a fixed B3-spline
 * kernel h and four data-dependent terms:
 *
 *     w = h(q-p)
 *       * exp( -|l_p - l_q| / (sigma_c * sqrt(Var_p) + eps) )     <- color
 *       * max(0, n_p . n_q)^sigma_n                                <- normal
 *       * exp( -|d_p - d_q| / sigma_d )                            <- depth
 *       * exp( -|a_p - a_q|^2 / sigma_a )                          <- albedo
 *
 * The color term is the crucial one, and the reason per-pixel variance had to
 * be estimated in the first place: normalizing the LUMINANCE difference by the
 * pixel's own standard deviation is what lets the SAME sigma_c behave
 * correctly in a converged dark region and in a noisy bright one - it makes
 * sigma_c a threshold measured in standard deviations, so "is this difference
 * bigger than the noise?" is asked in units of the noise itself. Note that
 * both numerator and denominator must be FIRST order in radiance for that to
 * hold; a squared difference over a standard deviation would silently make the
 * effective threshold proportional to brightness, over-blurring dark areas and
 * under-blurring bright ones.
 *
 * A single scalar (luminance) weight is applied to the whole RGB triple rather
 * than three independent per-channel weights, which would let the channels
 * drift apart and introduce color fringes at edges.
 *
 * The variance itself is filtered alongside the color, weighted by w^2 (the
 * variance of a weighted mean sum(w*c)/sum(w) is sum(w^2*Var)/(sum(w))^2), so
 * that later iterations know how much noise is actually left.
 *
 * When feature buffers are absent (GPU renders, integrators without li_aov()),
 * the geometric terms are simply dropped and this degrades into a pure
 * variance-guided color bilateral filter - noticeably blurrier at geometric
 * edges, but still a large improvement over the raw render.
 *
 * XML properties:
 *   - `iterations`  (int,   default 5)     number of a-trous passes
 *   - `sigma_c`     (float, default 4.0)   color bandwidth (in std-devs)
 *   - `sigma_n`     (float, default 128.0) normal exponent (larger = sharper)
 *   - `sigma_d`     (float, default 1.0)   depth bandwidth, relative units
 *   - `sigma_a`     (float, default 0.05)  albedo bandwidth
 *   - `demodulate`  (bool,  default true)  filter radiance/albedo, not radiance
 *   - `threads`     (int,   default 0)     0 = hardware concurrency
 */
class AtrousDenoiser : public Denoiser {
public:
    explicit AtrousDenoiser(const PropertyList &properties) {
        m_iterations = properties.get_integer("iterations", 5);
        m_sigma_c    = properties.get_float("sigma_c", 4.0f);
        m_sigma_n    = properties.get_float("sigma_n", 128.0f);
        m_sigma_d    = properties.get_float("sigma_d", 1.0f);
        m_sigma_a    = properties.get_float("sigma_a", 0.05f);
        m_demodulate = properties.get_boolean("demodulate", true);
        m_threads    = properties.get_integer("threads", 0);
        m_name       = properties.get_string("name", "atrous");

        if (m_iterations < 1)
            throw std::runtime_error("AtrousDenoiser: iterations must be >= 1");
        if (m_sigma_c <= 0.0f || m_sigma_d <= 0.0f || m_sigma_a <= 0.0f)
            throw std::runtime_error("AtrousDenoiser: sigma_c / sigma_d / sigma_a must be > 0");
    }

    void construct() override {
        Denoiser::construct();
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << " (atrous)" << std::endl;
#endif
    }

    [[nodiscard]] std::shared_ptr<Bitmap> denoise(const DenoiseInput &input) const override {
        auto prefiltered    = run_prefilter(input);
        const Bitmap &color = prefiltered ? *prefiltered : *input.color;

        const int rows = color.get_rows(), cols = color.get_cols();
        const int threads = m_threads > 0 ? m_threads : static_cast<int>(std::thread::hardware_concurrency());

        std::vector<float> c = denoise::to_vector(color, 3);
        std::vector<float> v = denoise::estimate_variance(input, rows, cols);

        // Demodulate: filter irradiance (radiance / albedo) so that texture
        // detail, which lives entirely in the albedo, is bypassed by the
        // filter instead of being smoothed along with the noise.
        const bool demodulate = m_demodulate && input.albedo != nullptr;
        if (demodulate)
            denoise::demodulate(c, v, *input.albedo, rows, cols);

        // Feature buffers are optional; without them the geometric terms drop
        // out of the weight and this becomes a plain variance-guided bilateral.
        const bool use_features = input.normal != nullptr && input.depth != nullptr;
        const bool use_albedo_w = input.albedo != nullptr;

        // Depth differences must be compared against the scene's own scale, or
        // sigma_d would mean something completely different for a scene
        // measured in millimeters and one measured in kilometers. Normalize by
        // the mean depth of the visible geometry.
        float depth_scale = 1.0f;
        if (use_features) {
            double sum = 0.0;
            int n      = 0;
            for (int y = 0; y < rows; ++y)
                for (int x = 0; x < cols; ++x) {
                    float d = (*input.depth)(y, x, 0);
                    if (d > 0.0f && d < M_MAX_FLOAT) {
                        sum += d;
                        ++n;
                    }
                }
            if (n > 0 && sum > 0.0)
                depth_scale = static_cast<float>(sum / n);
        }

        // B3-spline (1/16, 1/4, 3/8, 1/4, 1/16) - the standard a-trous kernel.
        static constexpr float h[5] = { 1.f / 16.f, 1.f / 4.f, 3.f / 8.f, 1.f / 4.f, 1.f / 16.f };

        std::vector<float> c_next(c.size()), v_next(v.size());

        for (int it = 0; it < m_iterations; ++it) {
            const int stride = 1 << it;

            denoise::parallel_bands(rows, threads, [&](int y_begin, int y_end) {
                for (int y = y_begin; y < y_end; ++y) {
                    for (int x = 0; x < cols; ++x) {
                        const size_t ip = (static_cast<size_t>(y) * cols + x) * 3;

                        // Standard deviation of the center pixel's LUMINANCE.
                        // Var(l) = sum (k_ch^2 * Var_ch) for l = sum k_ch*c_ch,
                        // the channels being independent estimates here.
                        constexpr float kr = 0.212671f, kg = 0.715160f, kb = 0.072169f;
                        float var_p = kr * kr * v[ip] + kg * kg * v[ip + 1] + kb * kb * v[ip + 2];
                        float sd_p  = std::sqrt(std::max(var_p, 0.0f));
                        float lum_p = kr * c[ip] + kg * c[ip + 1] + kb * c[ip + 2];

                        Normal3f n_p(0.f);
                        float d_p = 0.f;
                        if (use_features) {
                            n_p = Normal3f((*input.normal)(y, x, 0), (*input.normal)(y, x, 1),
                                           (*input.normal)(y, x, 2));
                            d_p = (*input.depth)(y, x, 0);
                        }

                        float sum_w = 0.0f;
                        float acc_c[3] = { 0.f, 0.f, 0.f };
                        float acc_v[3] = { 0.f, 0.f, 0.f };

                        for (int ky = 0; ky < 5; ++ky) {
                            for (int kx = 0; kx < 5; ++kx) {
                                int yy = y + (ky - 2) * stride;
                                int xx = x + (kx - 2) * stride;
                                if (yy < 0 || xx < 0 || yy >= rows || xx >= cols)
                                    continue; // truncate + renormalize at borders

                                const size_t iq = (static_cast<size_t>(yy) * cols + xx) * 3;

                                float w = h[ky] * h[kx];

                                // ------------------ color term ------------------
                                // Absolute luminance difference over the noise
                                // level: both are FIRST order in radiance, so
                                // sigma_c is a pure "number of standard
                                // deviations" threshold, independent of how
                                // bright this part of the image happens to be.
                                float lum_q = kr * c[iq] + kg * c[iq + 1] + kb * c[iq + 2];
                                // eps keeps a fully converged pixel (sd == 0)
                                // from producing a zero denominator, in which
                                // case any color difference at all is treated
                                // as a real edge - which is exactly right.
                                w *= std::exp(-std::abs(lum_p - lum_q) / (m_sigma_c * sd_p + 1e-6f));

                                // ---------------- geometric terms ----------------
                                if (use_features) {
                                    Normal3f n_q((*input.normal)(yy, xx, 0), (*input.normal)(yy, xx, 1),
                                                 (*input.normal)(yy, xx, 2));
                                    float nd = n_p.dot(n_q);
                                    // Pixels with no recorded normal (escaped
                                    // rays, silhouette averages) have n == 0,
                                    // making nd == 0 and hence w == 0: they
                                    // simply do not participate, rather than
                                    // contaminating the result.
                                    w *= std::pow(std::max(nd, 0.0f), m_sigma_n);

                                    float d_q = (*input.depth)(yy, xx, 0);
                                    w *= std::exp(-std::abs(d_p - d_q) / (m_sigma_d * depth_scale + 1e-8f));
                                }

                                if (use_albedo_w) {
                                    float da = 0.0f;
                                    for (int ch = 0; ch < 3; ++ch) {
                                        float diff = (*input.albedo)(y, x, ch) - (*input.albedo)(yy, xx, ch);
                                        da += diff * diff;
                                    }
                                    w *= std::exp(-da / m_sigma_a);
                                }

                                if (w <= 0.0f)
                                    continue;

                                sum_w += w;
                                for (int ch = 0; ch < 3; ++ch) {
                                    acc_c[ch] += w * c[iq + ch];
                                    // Var(sum w*c / sum w) = sum w^2 Var(c) / (sum w)^2
                                    acc_v[ch] += w * w * v[iq + ch];
                                }
                            }
                        }

                        if (sum_w > 0.0f) {
                            float inv = 1.0f / sum_w;
                            for (int ch = 0; ch < 3; ++ch) {
                                c_next[ip + ch] = acc_c[ch] * inv;
                                v_next[ip + ch] = acc_v[ch] * inv * inv;
                            }
                        } else {
                            // Every neighbor was rejected - keep the pixel as
                            // is (the center tap alone would give the same).
                            for (int ch = 0; ch < 3; ++ch) {
                                c_next[ip + ch] = c[ip + ch];
                                v_next[ip + ch] = v[ip + ch];
                            }
                        }
                    }
                }
            });

            c.swap(c_next);
            v.swap(v_next);
        }

        if (demodulate)
            denoise::modulate(c, *input.albedo, rows, cols);

        return denoise::to_bitmap(c, rows, cols);
    }

    [[nodiscard]] std::string to_string() const override {
        return "AtrousDenoiser[\n  iterations = " + std::to_string(m_iterations) +
               ",\n  sigma_c = " + std::to_string(m_sigma_c) + ",\n  sigma_n = " + std::to_string(m_sigma_n) +
               ",\n  sigma_d = " + std::to_string(m_sigma_d) + ",\n  sigma_a = " + std::to_string(m_sigma_a) +
               ",\n  demodulate = " + (m_demodulate ? "true" : "false") +
               ",\n  prefilter = " + (m_prefilter ? indent(m_prefilter->to_string(), 2) : "null") + "\n]";
    }

private:
    int m_iterations;
    float m_sigma_c;
    float m_sigma_n;
    float m_sigma_d;
    float m_sigma_a;
    bool m_demodulate;
    int m_threads;
};

REGISTER_CLASS(AtrousDenoiser, "atrous")

M_NAMESPACE_END

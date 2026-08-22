#include <algorithm>
#include <cmath>
#include <components/denoiser.h>
#include <render/denoise_utils.h>
#include <thread>

M_NAMESPACE_BEGIN

/**
 * \brief Feature-guided Non-Local Means with half-buffer cross-validation.
 *
 * This is the offline / high-quality tier, modelled on Rousselle et al.'s
 * "Robust Denoising using Feature and Colour Information" (RDFC, 2013) - the
 * approach that dominated offline MC denoising before neural networks, and
 * still the strongest option that needs no trained weights.
 *
 * ## Why NLM rather than a bilateral filter
 *
 * A bilateral filter decides whether two pixels may be averaged by comparing
 * two single values, c_p and c_q. At low sample counts those values ARE mostly
 * noise, so the decision is mostly noise too. NLM instead compares whole
 * PATCHES around p and q: averaging (2f+1)^2 squared differences suppresses the
 * noise in the *distance metric itself* by that same factor, so the filter can
 * tell "these two neighborhoods are the same surface under the same lighting"
 * apart from "these two pixels happened to land on similar noise" far more
 * reliably. That is what lets NLM use a much larger search window (21x21 here)
 * without smearing structure.
 *
 * ## The variance-cancelling distance
 *
 * A naive patch distance is biased: even for two pixels with the SAME true
 * mean, E[(u_p - u_q)^2] = Var_p + Var_q > 0, so noisy regions look
 * systematically dissimilar and get under-filtered exactly where filtering is
 * needed most. Following Rousselle et al., the known per-pixel variances are
 * subtracted out and the residual is normalized by them:
 *
 *     d^2(p,q) = [ (u_p - u_q)^2 - alpha*(Var_p + min(Var_p, Var_q)) ]
 *                / ( eps + k^2 * (Var_p + Var_q) )
 *
 * The `min(Var_p, Var_q)` (rather than `Var_q`) is deliberate: it prevents an
 * outlier pixel with a hugely overestimated variance from cancelling away an
 * arbitrarily large true difference and thereby dragging its neighborhood
 * towards itself.
 *
 * ## Cross-validation
 *
 * `k` controls the bandwidth: small k preserves detail but leaves noise, large
 * k does the opposite, and the best value differs per pixel (flat walls want
 * large k, textured/detailed areas want small k). Rather than exposing this as
 * a parameter to be hand-tuned per scene, the filter is run at TWO bandwidths
 * and the choice is made per pixel by cross-validation on the half-buffers:
 *
 *   - filter half-buffer A, compare against the *unfiltered* half-buffer B
 *   - filter half-buffer B, compare against the *unfiltered* half-buffer A
 *
 * Because A and B are statistically independent, the resulting squared error is
 * an unbiased estimate of (filtering bias)^2 + (residual variance) up to a
 * constant offset - i.e. of the actual MSE, without any reference image. The
 * error maps are smoothed (a per-pixel decision from one noisy estimate would
 * itself be noise) and whichever bandwidth wins locally is used for the final
 * full-buffer output, blended smoothly to avoid visible seams between regions
 * that picked different bandwidths.
 *
 * If half-buffers are unavailable, cross-validation is skipped and the filter
 * runs once at `k`.
 *
 * ## Cost
 *
 * O(rows * cols * (2*r+1)^2 * (2*f+1)^2) per filter run, and up to five runs
 * with cross-validation. This is deliberately the "quality, not speed" option -
 * expect seconds, not milliseconds. Use `atrous` for interactive work.
 *
 * XML properties:
 *   - `radius`        (int,   default 10)   search window half-width
 *   - `patch_radius`  (int,   default 3)    patch half-width
 *   - `k`             (float, default 0.45) bandwidth (detail-preserving)
 *   - `k_smooth`      (float, default 1.0)  second bandwidth (noise-removing)
 *   - `alpha`         (float, default 1.0)  variance cancellation strength
 *   - `cross_validate`(bool,  default true) per-pixel bandwidth selection
 *   - `demodulate`    (bool,  default true) filter radiance/albedo
 *   - `sigma_n`       (float, default 0.8)  normal feature bandwidth
 *   - `sigma_d`       (float, default 0.6)  depth feature bandwidth
 *   - `threads`       (int,   default 0)    0 = hardware concurrency
 */
class NLMDenoiser : public Denoiser {
public:
    explicit NLMDenoiser(const PropertyList &properties) {
        m_radius         = properties.get_integer("radius", 10);
        m_patch_radius   = properties.get_integer("patch_radius", 3);
        m_k              = properties.get_float("k", 0.45f);
        m_k_smooth       = properties.get_float("k_smooth", 1.0f);
        m_alpha          = properties.get_float("alpha", 1.0f);
        m_cross_validate = properties.get_boolean("cross_validate", true);
        m_demodulate     = properties.get_boolean("demodulate", true);
        m_sigma_n        = properties.get_float("sigma_n", 0.8f);
        m_sigma_d        = properties.get_float("sigma_d", 0.6f);
        m_threads        = properties.get_integer("threads", 0);
        m_name           = properties.get_string("name", "nlm");

        if (m_radius < 1 || m_patch_radius < 0)
            throw std::runtime_error("NLMDenoiser: radius must be >= 1 and patch_radius >= 0");
        if (m_k <= 0.0f || m_k_smooth <= 0.0f)
            throw std::runtime_error("NLMDenoiser: k and k_smooth must be > 0");
    }

    void construct() override {
        Denoiser::construct();
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << " (nlm)" << std::endl;
#endif
    }

    [[nodiscard]] std::shared_ptr<Bitmap> denoise(const DenoiseInput &input) const override {
        auto prefiltered    = run_prefilter(input);
        const Bitmap &color = prefiltered ? *prefiltered : *input.color;

        const int rows = color.get_rows(), cols = color.get_cols();
        const int threads = m_threads > 0 ? m_threads : static_cast<int>(std::thread::hardware_concurrency());

        std::vector<float> c = denoise::to_vector(color, 3);
        std::vector<float> v = denoise::estimate_variance(input, rows, cols);

        const bool demodulate = m_demodulate && input.albedo != nullptr;
        if (demodulate)
            denoise::demodulate(c, v, *input.albedo, rows, cols);

        Guides guides = build_guides(input, rows, cols);

        std::vector<float> out;

        if (m_cross_validate && input.has_half_buffers()) {
            // Half-buffers, demodulated consistently with the full buffer so
            // that the cross-validation error is measured in the same space the
            // filter actually operates in.
            std::vector<float> a = denoise::to_vector(*input.color_a, 3);
            std::vector<float> b = denoise::to_vector(*input.color_b, 3);
            // Each half holds ~half the samples, so its variance is ~2x the
            // full buffer's.
            std::vector<float> v_half(v.size());
            for (size_t i = 0; i < v.size(); ++i)
                v_half[i] = v[i] * 2.0f;

            if (demodulate) {
                // demodulate() scales color and variance together, so give each
                // call its own variance copy and keep one of the (identical)
                // results.
                std::vector<float> v_a = v_half, v_b = v_half;
                denoise::demodulate(a, v_a, *input.albedo, rows, cols);
                denoise::demodulate(b, v_b, *input.albedo, rows, cols);
                v_half = v_a;
            }

            // Per-pixel MSE proxy for each candidate bandwidth. Filtering A and
            // scoring it against the INDEPENDENT B (and vice versa) makes the
            // squared difference an unbiased MSE estimate up to a constant.
            std::vector<float> err_sharp = cross_validation_error(a, b, v_half, guides, rows, cols, threads, m_k);
            std::vector<float> err_smooth =
                cross_validation_error(a, b, v_half, guides, rows, cols, threads, m_k_smooth);

            // A per-pixel argmin over two single-sample error estimates would
            // itself be noisy and produce a speckled patchwork of bandwidths.
            // Smoothing the error maps first makes the selection spatially
            // coherent.
            std::vector<float> es, em;
            denoise::box_blur(err_sharp, es, rows, cols, 1, 5);
            denoise::box_blur(err_smooth, em, rows, cols, 1, 5);

            std::vector<float> filt_sharp  = filter(c, c, v, guides, rows, cols, threads, m_k);
            std::vector<float> filt_smooth = filter(c, c, v, guides, rows, cols, threads, m_k_smooth);

            out.assign(c.size(), 0.0f);
            for (int y = 0; y < rows; ++y) {
                for (int x = 0; x < cols; ++x) {
                    size_t p  = static_cast<size_t>(y) * cols + x;
                    float e_s = es[p], e_m = em[p];
                    // Soft (rather than hard) selection: interpolate by the
                    // relative errors, so neighboring pixels that disagree
                    // about the best bandwidth blend instead of forming a
                    // visible boundary between two differently-filtered
                    // regions.
                    float total = e_s + e_m;
                    float t     = total > 0.0f ? e_s / total : 0.5f; // weight of the SMOOTH result
                    for (int ch = 0; ch < 3; ++ch)
                        out[p * 3 + ch] = (1.0f - t) * filt_sharp[p * 3 + ch] + t * filt_smooth[p * 3 + ch];
                }
            }
        } else {
            out = filter(c, c, v, guides, rows, cols, threads, m_k);
        }

        if (demodulate)
            denoise::modulate(out, *input.albedo, rows, cols);

        return denoise::to_bitmap(out, rows, cols);
    }

    [[nodiscard]] std::string to_string() const override {
        return "NLMDenoiser[\n  radius = " + std::to_string(m_radius) +
               ",\n  patch_radius = " + std::to_string(m_patch_radius) + ",\n  k = " + std::to_string(m_k) +
               ",\n  k_smooth = " + std::to_string(m_k_smooth) + ",\n  alpha = " + std::to_string(m_alpha) +
               ",\n  cross_validate = " + (m_cross_validate ? "true" : "false") +
               ",\n  demodulate = " + (m_demodulate ? "true" : "false") +
               ",\n  prefilter = " + (m_prefilter ? indent(m_prefilter->to_string(), 2) : "null") + "\n]";
    }

private:
    // Normalized, optional geometric guides. Empty vectors mean "not available",
    // in which case the corresponding term drops out of the weight.
    struct Guides {
        std::vector<float> normal; // 3 channels, unit length (or zero)
        std::vector<float> depth;  // 1 channel, divided by the mean scene depth
        bool has_normal = false;
        bool has_depth  = false;
    };

    static Guides build_guides(const DenoiseInput &input, int rows, int cols) {
        Guides g;
        if (input.normal) {
            g.normal     = denoise::to_vector(*input.normal, 3);
            g.has_normal = true;
        }
        if (input.depth) {
            g.depth = denoise::to_vector(*input.depth, 1);
            // Same scale normalization as the a-trous filter: sigma_d must mean
            // the same thing regardless of the units the scene was modelled in.
            double sum = 0.0;
            int n      = 0;
            for (float d : g.depth)
                if (d > 0.0f && d < M_MAX_FLOAT) {
                    sum += d;
                    ++n;
                }
            float scale = (n > 0 && sum > 0.0) ? static_cast<float>(sum / n) : 1.0f;
            for (float &d : g.depth)
                d /= scale;
            g.has_depth = true;
        }
        (void) rows;
        (void) cols;
        return g;
    }

    /**
     * \brief NLM-filter `src` using patch distances computed on `guide_color`.
     *
     * Splitting these two lets cross-validation filter half-buffer A while
     * driving the weights from A itself, and lets the final pass filter the
     * full buffer - all through one code path.
     */
    [[nodiscard]] std::vector<float> filter(const std::vector<float> &src, const std::vector<float> &guide_color,
                                            const std::vector<float> &var, const Guides &guides, int rows, int cols,
                                            int threads, float k) const {
        std::vector<float> dst(src.size(), 0.0f);
        const float k2 = k * k;

        denoise::parallel_bands(rows, threads, [&](int y_begin, int y_end) {
            for (int y = y_begin; y < y_end; ++y) {
                for (int x = 0; x < cols; ++x) {
                    const size_t p = static_cast<size_t>(y) * cols + x;

                    float sum_w    = 0.0f;
                    float acc[3]   = { 0.f, 0.f, 0.f };

                    for (int dy = -m_radius; dy <= m_radius; ++dy) {
                        for (int dx = -m_radius; dx <= m_radius; ++dx) {
                            int qy = y + dy, qx = x + dx;
                            if (qy < 0 || qx < 0 || qy >= rows || qx >= cols)
                                continue;
                            const size_t q = static_cast<size_t>(qy) * cols + qx;

                            // ---------------- patch distance ----------------
                            float d_sum = 0.0f;
                            int d_count = 0;
                            for (int py = -m_patch_radius; py <= m_patch_radius; ++py) {
                                for (int px = -m_patch_radius; px <= m_patch_radius; ++px) {
                                    int ay = denoise::clamp_i(y + py, 0, rows - 1);
                                    int ax = denoise::clamp_i(x + px, 0, cols - 1);
                                    int by = denoise::clamp_i(qy + py, 0, rows - 1);
                                    int bx = denoise::clamp_i(qx + px, 0, cols - 1);

                                    size_t ia = static_cast<size_t>(ay) * cols + ax;
                                    size_t ib = static_cast<size_t>(by) * cols + bx;

                                    for (int ch = 0; ch < 3; ++ch) {
                                        float diff = guide_color[ia * 3 + ch] - guide_color[ib * 3 + ch];
                                        float va   = var[ia * 3 + ch];
                                        float vb   = var[ib * 3 + ch];

                                        // Rousselle et al. 2013 eq. 2: subtract
                                        // the variance that the difference is
                                        // EXPECTED to contain even for two
                                        // identical means, and normalize by the
                                        // remaining noise level. min(va, vb)
                                        // guards against an outlier's inflated
                                        // variance cancelling a real edge.
                                        float num = diff * diff - m_alpha * (va + std::min(va, vb));
                                        float den = 1e-10f + k2 * (va + vb);
                                        d_sum += num / den;
                                        ++d_count;
                                    }
                                }
                            }
                            float d2 = d_count > 0 ? d_sum / static_cast<float>(d_count) : 0.0f;

                            float w = std::exp(-std::max(d2, 0.0f));

                            // ---------------- geometric guides ----------------
                            if (guides.has_normal) {
                                Normal3f n_p(guides.normal[p * 3], guides.normal[p * 3 + 1],
                                             guides.normal[p * 3 + 2]);
                                Normal3f n_q(guides.normal[q * 3], guides.normal[q * 3 + 1],
                                             guides.normal[q * 3 + 2]);
                                float nd = 1.0f - std::max(n_p.dot(n_q), 0.0f);
                                w *= std::exp(-nd * nd / (m_sigma_n * m_sigma_n));
                            }
                            if (guides.has_depth) {
                                float dd = guides.depth[p] - guides.depth[q];
                                w *= std::exp(-dd * dd / (m_sigma_d * m_sigma_d));
                            }

                            if (w <= 0.0f)
                                continue;

                            sum_w += w;
                            for (int ch = 0; ch < 3; ++ch)
                                acc[ch] += w * src[q * 3 + ch];
                        }
                    }

                    if (sum_w > 0.0f) {
                        float inv = 1.0f / sum_w;
                        for (int ch = 0; ch < 3; ++ch)
                            dst[p * 3 + ch] = acc[ch] * inv;
                    } else {
                        for (int ch = 0; ch < 3; ++ch)
                            dst[p * 3 + ch] = src[p * 3 + ch];
                    }
                }
            }
        });

        return dst;
    }

    /**
     * \brief Per-pixel, reference-free MSE proxy for bandwidth `k`.
     *
     * Filters A and measures it against the untouched, statistically
     * independent B, and symmetrically. Independence is what makes this work:
     * the noise in the "reference" is uncorrelated with the noise the filter
     * saw, so it contributes only a constant offset that is identical for every
     * candidate bandwidth and therefore cancels out of the comparison.
     */
    [[nodiscard]] std::vector<float> cross_validation_error(const std::vector<float> &a, const std::vector<float> &b,
                                                            const std::vector<float> &var_half, const Guides &guides,
                                                            int rows, int cols, int threads, float k) const {
        std::vector<float> fa = filter(a, a, var_half, guides, rows, cols, threads, k);
        std::vector<float> fb = filter(b, b, var_half, guides, rows, cols, threads, k);

        std::vector<float> err(static_cast<size_t>(rows) * cols, 0.0f);
        for (size_t p = 0; p < err.size(); ++p) {
            float e = 0.0f;
            for (int ch = 0; ch < 3; ++ch) {
                float d1 = fa[p * 3 + ch] - b[p * 3 + ch];
                float d2 = fb[p * 3 + ch] - a[p * 3 + ch];
                e += d1 * d1 + d2 * d2;
            }
            err[p] = e;
        }
        return err;
    }

    int m_radius;
    int m_patch_radius;
    float m_k;
    float m_k_smooth;
    float m_alpha;
    bool m_cross_validate;
    bool m_demodulate;
    float m_sigma_n;
    float m_sigma_d;
    int m_threads;
};

REGISTER_CLASS(NLMDenoiser, "nlm")

M_NAMESPACE_END

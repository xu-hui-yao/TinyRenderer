#include <algorithm>
#include <components/denoiser.h>
#include <render/denoise_utils.h>
#include <thread>

M_NAMESPACE_BEGIN

/**
 * \brief Neighborhood clamping / firefly rejection ("outlier removal").
 *
 * Strictly speaking this is not a denoiser but a *heavy-tail suppressor*, and
 * it is best understood as a variance-reduction step rather than a filter.
 *
 * Path tracing occasionally produces samples whose contribution is orders of
 * magnitude larger than the pixel's true mean - a caustic path found by pure
 * luck, a near-degenerate BSDF sample divided by a tiny pdf, a tiny bright
 * emitter hit head-on. Those isolated white pixels ("fireflies") converge
 * eventually, but at practical sample counts they dominate the error, and they
 * are actively harmful to any subsequent smoothing filter: a single spike gets
 * smeared into a bright blob across the entire filter footprint, turning one
 * bad pixel into dozens.
 *
 * The fix is a robust local statistic: compare each pixel's luminance against
 * the *second largest* luminance in its (2r+1)^2 neighborhood - which by
 * construction ignores the center pixel's own spike as well as one additional
 * outlier - and rescale the pixel down if it exceeds `threshold` times that
 * value. Rescaling is done on the RGB triple as a whole (a single luminance
 * ratio), so hue is preserved exactly.
 *
 * This introduces BIAS: energy is genuinely removed, and a legitimately tiny,
 * legitimately very bright highlight will be dimmed. That is why `threshold`
 * defaults to a fairly permissive 2.0 and why this filter is most useful as a
 * *pre-filter* nested inside one of the real denoisers rather than as an
 * output filter of its own.
 *
 * Needs no auxiliary buffers whatsoever, so it also works on GPU renders.
 *
 * XML properties:
 *   - `radius`    (int,   default 1)   half-width of the neighborhood
 *   - `threshold` (float, default 2.0) allowed multiple of the robust maximum
 *   - `threads`   (int,   default 0)   0 = use hardware concurrency
 */
class OutlierDenoiser : public Denoiser {
public:
    explicit OutlierDenoiser(const PropertyList &properties) {
        m_radius    = properties.get_integer("radius", 1);
        m_threshold = properties.get_float("threshold", 2.0f);
        m_threads   = properties.get_integer("threads", 0);
        m_name      = properties.get_string("name", "outlier");

        if (m_radius < 1)
            throw std::runtime_error("OutlierDenoiser: radius must be >= 1");
        if (m_threshold <= 1.0f)
            throw std::runtime_error("OutlierDenoiser: threshold must be > 1 (a threshold of 1 or less would clamp "
                                     "essentially every pixel and destroy the image)");
    }

    void construct() override {
        Denoiser::construct();
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << " (outlier)" << std::endl;
#endif
    }

    [[nodiscard]] std::shared_ptr<Bitmap> denoise(const DenoiseInput &input) const override {
        // A nested pre-filter is unusual here but harmless, and keeps the
        // composition rules uniform across all denoisers.
        auto prefiltered = run_prefilter(input);
        const Bitmap &color = prefiltered ? *prefiltered : *input.color;

        const int rows = color.get_rows(), cols = color.get_cols();
        std::vector<float> src = denoise::to_vector(color, 3);
        std::vector<float> dst = src;

        int threads = m_threads > 0 ? m_threads : static_cast<int>(std::thread::hardware_concurrency());

        denoise::parallel_bands(rows, threads, [&](int y_begin, int y_end) {
            // Luminances of the neighborhood, reused across pixels.
            std::vector<float> lum;
            lum.reserve(static_cast<size_t>(2 * m_radius + 1) * (2 * m_radius + 1));

            for (int y = y_begin; y < y_end; ++y) {
                for (int x = 0; x < cols; ++x) {
                    size_t i = (static_cast<size_t>(y) * cols + x) * 3;
                    Color3f c({ src[i], src[i + 1], src[i + 2] });
                    float center = c.luminance();
                    if (center <= 0.0f)
                        continue;

                    lum.clear();
                    for (int dy = -m_radius; dy <= m_radius; ++dy) {
                        for (int dx = -m_radius; dx <= m_radius; ++dx) {
                            if (dy == 0 && dx == 0)
                                continue; // the center is what we are judging
                            int yy = denoise::clamp_i(y + dy, 0, rows - 1);
                            int xx = denoise::clamp_i(x + dx, 0, cols - 1);
                            size_t j = (static_cast<size_t>(yy) * cols + xx) * 3;
                            lum.push_back(Color3f({ src[j], src[j + 1], src[j + 2] }).luminance());
                        }
                    }
                    if (lum.size() < 2)
                        continue;

                    // Second largest neighbor luminance: a robust upper bound
                    // on "how bright it is reasonable to be here". Using the
                    // plain maximum would let a PAIR of adjacent fireflies
                    // shield each other from detection.
                    std::nth_element(lum.begin(), lum.begin() + 1, lum.end(), std::greater<float>());
                    float robust_max = lum[1];

                    float limit = m_threshold * robust_max;
                    if (center > limit) {
                        // Uniform luminance rescale - preserves chromaticity,
                        // so a clamped firefly desaturates towards nothing
                        // rather than shifting hue.
                        float scale  = limit / center;
                        dst[i]     = src[i] * scale;
                        dst[i + 1] = src[i + 1] * scale;
                        dst[i + 2] = src[i + 2] * scale;
                    }
                }
            }
        });

        return denoise::to_bitmap(dst, rows, cols);
    }

    [[nodiscard]] std::string to_string() const override {
        return "OutlierDenoiser[\n  radius = " + std::to_string(m_radius) +
               ",\n  threshold = " + std::to_string(m_threshold) +
               ",\n  prefilter = " + (m_prefilter ? indent(m_prefilter->to_string(), 2) : "null") + "\n]";
    }

private:
    int m_radius;
    float m_threshold;
    int m_threads;
};

REGISTER_CLASS(OutlierDenoiser, "outlier")

M_NAMESPACE_END

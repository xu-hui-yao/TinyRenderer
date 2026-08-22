#pragma once

#include <components/object.h>
#include <render/framebuffer.h>

M_NAMESPACE_BEGIN

/**
 * \brief Abstract single-frame Monte Carlo denoiser.
 *
 * A denoiser consumes the linear-HDR output of a finished render together with
 * whatever auxiliary buffers it came with (see \ref DenoiseInput) and returns
 * a filtered image in the same linear-HDR space. It runs strictly *after*
 * light transport and *before* tonemapping / sRGB encoding, which means the
 * exact same implementation serves both the CPU and the GPU render paths.
 *
 * Denoisers are configured from the scene XML as a child of <scene>:
 *
 * \code
 * <denoiser type="atrous">
 *     <integer name="iterations" value="5"/>
 *     <float name="sigma_c" value="4.0"/>
 *     <!-- optional pre-filter, applied to `color` before the main filter -->
 *     <denoiser type="outlier"/>
 * </denoiser>
 * \endcode
 *
 * A denoiser may accept a single nested denoiser as a child, which acts as a
 * pre-filter (see \ref m_prefilter): its output replaces `color` before the
 * outer filter runs. The canonical use is firefly / outlier rejection, whose
 * heavy-tailed spikes would otherwise leak across a whole filter footprint.
 */
class Denoiser : public Object {
public:
    ~Denoiser() override = default;

    /**
     * \brief Filter `input.color` and return the result.
     *
     * The returned bitmap always has 3 channels and the same resolution as
     * `input.color`, and stays in LINEAR radiance space.
     *
     * Implementations must degrade gracefully when the optional members of
     * `input` are absent (see \ref DenoiseInput): a render produced without
     * AOVs - e.g. by the GPU backend, or by an integrator that does not
     * override \ref Integrator::li_aov() - carries neither feature buffers nor
     * half-buffers, and denoising must still produce a sane image rather than
     * fail.
     */
    [[nodiscard]] virtual std::shared_ptr<Bitmap> denoise(const DenoiseInput &input) const = 0;

    void add_child(const std::shared_ptr<Object> &child) override {
        if (child->get_class_type() != EDenoiser)
            throw std::runtime_error("Denoiser::add_child(<" + class_type_name(child->get_class_type()) +
                                     ">) is not supported!");
        if (m_prefilter)
            throw std::runtime_error("Denoiser: tried to register more than one nested pre-filter!");
        m_prefilter = std::dynamic_pointer_cast<Denoiser>(child);
    }

    void construct() override {
        if (m_prefilter)
            m_prefilter->construct();
    }

    [[nodiscard]] EClassType get_class_type() const override { return EDenoiser; }

    [[nodiscard]] std::string to_string() const override = 0;

    [[nodiscard]] const std::shared_ptr<Denoiser> &get_prefilter() const { return m_prefilter; }

protected:
    /**
     * \brief Run the nested pre-filter (if any) on `input`, returning the
     * pre-filtered color buffer, or null if there is no pre-filter.
     *
     * Call this at the top of denoise() and substitute the result for
     * `input.color`. Note that only `color` is replaced - the variance and
     * feature buffers describe the original render and stay valid.
     */
    [[nodiscard]] std::shared_ptr<Bitmap> run_prefilter(const DenoiseInput &input) const {
        return m_prefilter ? m_prefilter->denoise(input) : nullptr;
    }

    std::shared_ptr<Denoiser> m_prefilter = nullptr;
};

M_NAMESPACE_END

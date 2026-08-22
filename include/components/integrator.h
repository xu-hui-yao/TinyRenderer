#pragma once

#include <components/object.h>
#include <render/aov.h>

M_NAMESPACE_BEGIN
/**
 * \brief Abstract integrator (i.e. a rendering technique)
 *
 * In Nori, the different rendering techniques are collectively referred to as
 * integrators, since they perform integration over a high-dimensional
 * space. Each integrator represents a specific approach for solving
 * the light transport equation---usually favored in certain scenarios, but
 * at the same time affected by its own set of intrinsic limitations.
 */
class Integrator : public Object {
public:
    /// Release all memory
    ~Integrator() override = default;

    void construct() override = 0;

    void add_child(const std::shared_ptr<Object> &child) override = 0;

    /// Perform an (optional) preprocess step
    virtual void preprocess(const std::shared_ptr<Scene> &scene) {}

    /**
     * \brief Sample the incident radiance along a ray
     *
     * \param scene
     *    A pointer to the underlying scene
     * \param sampler
     *    A pointer to a sample generator
     * \param ray
     *    The ray in question
     * \param valid
     *	  The valid
     * \return
     *    A (usually) unbiased estimate of the radiance in this direction
     */
    [[nodiscard]] virtual Color3f li(const std::shared_ptr<Scene> &scene, std::shared_ptr<Sampler> sampler,
                                     const Ray3f &ray, bool &valid) const = 0;

    /**
     * \brief Sample the incident radiance along a ray AND report the auxiliary
     * G-buffer features of the primary hit (see \ref AOVSample).
     *
     * This is what the denoising pipeline calls instead of \ref li(): the
     * feature buffers it fills are noise-free surface attributes that a
     * denoiser needs in order to tell "these two neighboring pixels look
     * different because of noise" apart from "...because there is an actual
     * edge here".
     *
     * The default implementation simply forwards to \ref li() and leaves
     * `aov` at its neutral defaults (albedo 1, no geometric feature), so
     * integrators that have no meaningful notion of a primary hit keep
     * working unchanged - they just cannot benefit from feature-guided
     * filtering, and the denoisers fall back to purely color-based weights.
     *
     * \param aov
     *    Output parameter, populated with the primary hit's albedo / shading
     *    normal / depth.
     */
    [[nodiscard]] virtual Color3f li_aov(const std::shared_ptr<Scene> &scene, const std::shared_ptr<Sampler> &sampler,
                                         const Ray3f &ray, bool &valid, AOVSample &aov) const {
        (void) aov;
        return li(scene, sampler, ray, valid);
    }

    /**
     * \brief Return the type of object (i.e. Mesh/BSDF/etc.)
     * provided by this instance
     * */
    [[nodiscard]] EClassType get_class_type() const override { return EIntegrator; }

    // Exposes the integrator's max_depth/rr_depth so the GPU port's host
    // code (path_trace_main.cpp) can mirror the exact same values instead of
    // hardcoding its own, for an apples-to-apples CPU/GPU comparison.
    // Default implementation returns Path's own defaults (see path.cpp);
    // integrators without a notion of depth may leave these unimplemented.
    [[nodiscard]] virtual int get_max_depth() const { return 10; }
    [[nodiscard]] virtual int get_rr_depth() const { return 5; }
};

M_NAMESPACE_END

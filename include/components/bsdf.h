#pragma once

#include <components/object.h>

M_NAMESPACE_BEGIN
class BSDF : public Object {
public:
    /**
     * \brief Importance sample the BSDF model
     *
     * The function returns a sample data structure along with the importance
     * weight, which is the value of the BSDF divided by the probability
     * density, and multiplied by the cosine foreshortening factor (if needed
     * --- it is omitted for degenerate BSDFs like smooth mirrors/dielectrics).
     *
     * If the supplied context data structures selects subset of components in
     * a multi-lobe BRDF model, the sampling is restricted to this subset.
     * Depending on the provided transport type, either the BSDF or its joint
     * version is sampled.
     *
     * When sampling a continuous/non-delta component, this method also
     * multiplies by the cosine foreshortening factor with respect to the
     * sampled direction.
     *
     * \param si
     *     A surface interaction data structure describing the underlying
     *     surface position. The incident direction is obtained from
     *     the field <tt>si.wi</tt>.
     *
     * \param sample1
     *     A uniformly distributed sample on \f$[0,1]\f$. It is used
     *     to select the BSDF lobe in multi-lobe models.
     *
     * \param sample2
     *     A uniformly distributed sample on \f$[0,1]^2\f$. It is
     *     used to generate the sampled direction.
     *
     * \param active
     *
     * \return A pair (bs, value) consisting of
     *
     *     bs:    Sampling record, indicating the sampled direction, PDF values
     *            and other information. The contents are undefined if sampling
     *            failed.
     *
     *     value: The BSDF value divided by the probability multiplied by the
     *            cosine foreshortening factor when a non-delta component is
     *            sampled. A zero spectrum indicates that sampling failed.
     */
    [[nodiscard]] virtual std::pair<BSDFSample3f, Color3f> sample(const SurfaceIntersection3f &si, float sample1,
                                                                  const Point2f &sample2, bool active) const = 0;

    /**
     * \brief Evaluate the BSDF f(wi, wo) or its joint version f^{*}(wi, wo)
     * and multiply by the cosine foreshortening term.
     *
     * Based on the information in the supplied query context \c ctx, this
     * method will either evaluate the entire BSDF or query individual
     * components (e.g. the diffuse lobe). Only smooth (i.e. non Dirac-delta)
     * components are supported: calling ``eval()`` on a perfectly specular
     * material will return zero.
     *
     * Note that the incident direction does not need to be explicitly
     * specified. It is obtained from the field <tt>si.wi</tt>.
     *
     * \param si
     *     A surface interaction data structure describing the underlying
     *     surface position. The incident direction is obtained from
     *     the field <tt>si.wi</tt>.
     *
     * \param wo
     *     The outgoing direction
     *
     * \param active
     */
    [[nodiscard]] virtual Color3f eval(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const = 0;

    /**
     * \brief Compute the probability per unit solid angle of sampling a
     * given direction
     *
     * This method provides access to the probability density that would result
     * when supplying the same BSDF context and surface interaction data
     * structures to the \ref sample() method. It correctly handles changes in
     * probability when only a subset of the components is chosen for sampling
     * (this can be done using the \ref BSDFContext::component and \ref
     * BSDFContext::type_mask fields).
     *
     * Note that the incident direction does not need to be explicitly
     * specified. It is obtained from the field <tt>si.wi</tt>.
     *
     * \param si
     *     A surface interaction data structure describing the underlying
     *     surface position. The incident direction is obtained from
     *     the field <tt>si.wi</tt>.
     *
     * \param wo
     *     The outgoing direction
     *
     * \param active
     */
    [[nodiscard]] virtual float pdf(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const = 0;

    void add_child(const std::shared_ptr<Object> &child) override = 0;

    void construct() override = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/BSDF/etc.)
     * provided by this instance
     */
    [[nodiscard]] EClassType get_class_type() const override { return EBSDF; }

    // Return a human-readable summary of this instance
    [[nodiscard]] std::string to_string() const override = 0;
};

M_NAMESPACE_END

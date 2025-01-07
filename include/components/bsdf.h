#pragma once

#include <components/object.h>

M_NAMESPACE_BEGIN

enum BSDFFlags : uint32_t {
    // BSDF lobe types
    // No flags set (default value)
    EEmpty = 0x00000,

    // 'null' scattering event, i.e. particles do not undergo deflection
    ENull = 0x00001,

    // Ideally diffuse reflection
    EDiffuseReflection = 0x00002,

    // Ideally diffuse transmission
    EDiffuseTransmission = 0x00004,

    // Glossy reflection
    EGlossyReflection = 0x00008,

    // Glossy transmission
    EGlossyTransmission = 0x00010,

    // Reflection into a discrete set of directions
    EDeltaReflection = 0x00020,

    // Transmission into a discrete set of directions
    EDeltaTransmission = 0x00040,

    // Compound lobe attributes
    // Any reflection component (scattering into discrete, 1D, or 2D set of directions)
    EReflection = EDiffuseReflection | EDeltaReflection | EGlossyReflection,

    // Any transmission component (scattering into discrete, 1D, or 2D set of directions)
    ETransmission = EDiffuseTransmission | EDeltaTransmission | EGlossyTransmission | ENull,

    // Diffuse scattering into a 2D set of directions
    EDiffuse = EDiffuseReflection | EDiffuseTransmission,

    // Non-diffuse scattering into a 2D set of directions
    EGlossy = EGlossyReflection | EGlossyTransmission,

    // Scattering into a 2D set of directions
    ESmooth = EDiffuse | EGlossy,

    // Scattering into a discrete set of directions
    EDelta = ENull | EDeltaReflection | EDeltaTransmission,

    // Any kind of scattering
    EAll = EDiffuse | EGlossy | EDelta
};

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

    [[nodiscard]] EClassType get_class_type() const override;

    bool has_flag(BSDFFlags bsdf_flags) const;

    [[nodiscard]] BSDFFlags get_flag() const;

    [[nodiscard]] std::string to_string() const override = 0;

protected:
    BSDFFlags m_flags = EEmpty;
};

M_NAMESPACE_END

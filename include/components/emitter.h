#pragma once

#include <components/object.h>

M_NAMESPACE_BEGIN
class Emitter : public Object {
public:
    void add_child(const std::shared_ptr<Object> &child) override = 0;

    void construct() override = 0;
    // =============================================================
    //! @{ \name Direction sampling interface
    // =============================================================

    /**
     * \brief Given a reference point in the scene, sample a direction from the
     * reference point towards the endpoint (ideally proportional to the
     * emission/sensitivity profile)
     *
     * The default implementation throws an exception.
     *
     * return A \ref DirectionSample instance describing the generated sample
     * along with a spectral importance weight.
     */
    [[nodiscard]] virtual std::pair<DirectionSample3f, Color3f>
    sample_direction(const Intersection3f &ref, const Point2f &sample, bool &active) const = 0;

    /**
     * \brief Evaluate the probability density of the \a direct sampling
     * method implemented by the \ref sample_direction() method.
     *
     * The returned probability will always be zero when the
     * emission/sensitivity profile contains a Dirac delta term (e.g. point or
     * directional emitters/sensors).
     */
    [[nodiscard]] virtual float pdf_direction(const Intersection3f &ref, const DirectionSample3f &ds,
                                              bool &active) const = 0;

    // =============================================================
    //! @{ \name Position sampling interface
    // =============================================================

    /**
     * \brief Importance sample the spatial component of the
     * emission or importance profile of the endpoint.
     *
     * return A \ref PositionSample instance describing the generated sample
     * along with an importance weight.
     */
    [[nodiscard]] virtual std::pair<PositionSample3f, float> sample_position(const Point2f &sample,
                                                                             bool &active) const = 0;

    /**
     * \brief Evaluate the probability density of the position sampling
     * method implemented by \ref sample_position().
     *
     * In simple cases, this will be the reciprocal of the endpoint's
     * surface area.
     */
    [[nodiscard]] virtual float pdf_position(const PositionSample3f &ps, bool &active) const = 0;

    // =============================================================
    //! @{ \name Other query functions
    // =============================================================

    /**
     * \brief Given a ray-surface intersection, return the emitted
     * radiance or importance traveling along the reverse direction
     *
     * This function is e.g. used when an area light source has been hit by a
     * ray in a path tracing-style integrator, and it subsequently needs to be
     * queried for the emitted radiance along the negative ray direction.
     *
     * return The emitted radiance or importance
     */
    [[nodiscard]] virtual Color3f eval(const SurfaceIntersection3f &si, bool &active) const = 0;

    [[nodiscard]] EClassType get_class_type() const override { return EEmitter; }

    [[nodiscard]] std::string to_string() const override;

    [[nodiscard]] virtual bool is_spatial_varying() const = 0;
};

M_NAMESPACE_END

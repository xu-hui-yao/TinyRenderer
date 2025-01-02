#pragma once

#include <components/object.h>

// Reconstruction filters will be tabulated at this resolution
#define M_FILTER_RESOLUTION 32

M_NAMESPACE_BEGIN
/**
 * \brief Generic radially symmetric image reconstruction filter
 *
 * When adding radiance-valued samples to the rendered image, Nori
 * first convolve them with a so-called image reconstruction filter.
 *
 * To learn more about reconstruction filters and sampling theory
 * in general, take a look at the excellent chapter 7 of PBRT,
 * which is freely available at:
 *
 * http://graphics.stanford.edu/~mmp/chapters/pbrt_chapter7.pdf
 */
class ReconstructionFilter : public Object {
public:
    /// Return the filter radius in fractional pixels
    [[nodiscard]] float get_radius() const { return m_radius; }

    /// Evaluate the filter function
    [[nodiscard]] virtual float eval(float x) const = 0;

    void construct() override = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/Camera/etc.)
     * provided by this instance
     * */
    [[nodiscard]] EClassType get_class_type() const override { return EReconstructionFilter; }

protected:
    float m_radius{};
};

M_NAMESPACE_END

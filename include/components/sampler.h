#pragma once

#include <components/object.h>
#include <memory>
#include <render/block.h>

M_NAMESPACE_BEGIN
/**
 * \brief Abstract sample generator
 *
 * A sample generator is responsible for generating the random number stream
 * that will be passed an \ref Integrator implementation as it computes the
 * radiance incident along a specified ray.
 *
 * The most simple conceivable sample generator is just a wrapper around the
 * Mersenne-Twister random number generator and is implemented in
 * <tt>independent.cpp</tt> (it is named this way because it generates
 * statistically independent random numbers).
 *
 * Fancier samplers might use stratification or low-discrepancy sequences
 * (e.g. Halton, Hammersley, or Sobol point sets) for improved convergence.
 * Another use of this class is in producing intentionally correlated
 * random numbers, e.g. as part of a Metropolis-Hastings integration scheme.
 *
 * The general interface between a sampler and a rendering algorithm is as
 * follows: Before beginning to render a pixel, the rendering algorithm calls
 * \ref generate(). The first pixel sample can now be computed, after which
 * \ref advance() needs to be invoked. This repeats until all pixel samples have
 * been exhausted.  While computing a pixel sample, the rendering
 * algorithm requests (pseudo-) random numbers using the \ref next1D() and
 * \ref next2D() functions.
 *
 * Conceptually, the right way of thinking of this goes as follows:
 * For each sample in a pixel, a sample generator produces a (hypothetical)
 * point in an infinite dimensional random number hypercube. A rendering
 * algorithm can then request subsequent 1D or 2D components of this point
 * using the \ref next1D() and \ref next2D() functions. Fancy implementations
 * of this class make certain guarantees about the stratification of the
 * first n components with respect to the other points that are sampled
 * within a pixel.
 */
class Sampler : public Object {
public:
    // Release all memory
    ~Sampler() override = default;

    void construct() override = 0;

    void add_child(const std::shared_ptr<Object> &child) override = 0;

    // Create an exact clone of the current instance
    [[nodiscard]] virtual std::shared_ptr<Sampler> clone() const = 0;

    /**
     * \brief Prepare to render a new image block
     *
     * This function is called when the sampler begins rendering
     * a new image block. This can be used to deterministically
     * initialize the sampler so that repeated program runs
     * always create the same image.
     */
    virtual void prepare(const ImageBlock &block) = 0;

    /**
     * \brief Prepare to generate new samples
     *
     * This function is called initially and every time the
     * integrator starts rendering a new pixel.
     */
    virtual void generate() = 0;

    // Advance to the next sample
    virtual void advance() = 0;

    // Retrieve the next component value from the current sample
    virtual float next1d() = 0;

    // Retrieve the next two component values from the current sample
    virtual Point2f next2d() = 0;

    // Return the number of configured pixel samples
    [[nodiscard]] virtual size_t get_sample_count() const { return m_sample_count; }

    /**
     * \brief Return the type of object (i.e. Mesh/Sampler/etc.)
     * provided by this instance
     * */
    [[nodiscard]] EClassType get_class_type() const override { return ESampler; }

protected:
    unsigned long long m_sample_count{};
};

M_NAMESPACE_END

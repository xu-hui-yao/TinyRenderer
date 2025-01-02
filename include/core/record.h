#pragma once

#include <components/emitter.h>
#include <core/common.h>

M_NAMESPACE_BEGIN
/**
 * \brief Generic sampling record for positions
 *
 * This sampling record is used to implement techniques that draw a position
 * from a point, line, surface, or volume domain in 3D and furthermore provide
 * auxiliary information about the sample.
 *
 * Apart from returning the position and (optionally) the surface normal, the
 * responsible sampling method must annotate the record with the associated
 * probability density and delta.
 */
template <typename Scalar_> class TPositionSample {
public:
    typedef Scalar_ Scalar;
    typedef TArray<Scalar, 3, ArrayType::Vector> VectorType;
    typedef TArray<Scalar, 3, ArrayType::Point> PointType;
    typedef TArray<Scalar, 2, ArrayType::Point> PointType2;
    typedef TArray<Scalar, 3, ArrayType::Normal> NormalType;
    typedef TFrame<Scalar> FrameType;

    PointType p;        // Position of the surface intersection
    NormalType n;       // Sampled surface normal (if applicable)
    PointType2 uv;      // Optional: 2D sample position associated with the record
    Scalar t;           // Un-occluded distance along the ray
    Scalar pdf;         // Probability density at the sample
    bool delta = false; // Set if the sample was drawn from a degenerate (Dirac delta) distribution

    TPositionSample() = default;

    explicit TPositionSample(const TDirectionSample<Scalar> other) {
        p     = other.p;
        n     = other.n;
        uv    = other.uv;
        t     = other.t;
        pdf   = other.pdf;
        delta = other.delta;
    };

    explicit TPositionSample(const TSurfaceIntersection<Scalar> &si)
        : p(si.p), n(si.n), uv(si.uv), t(si.t), pdf(0), delta(false) {}

    explicit TPositionSample(const PointType &_p, const NormalType &_n, const PointType2 &_uv, Scalar _t, Scalar _pdf,
                             bool _delta)
        : p(_p), n(_n), uv(_uv), t(_t), pdf(_pdf), delta(_delta) {}

    // Return a human-readable string summary of this frame
    [[nodiscard]] std::string to_string() const {
        std::ostringstream oss;
        oss << "Position Sample[\n"
               "  p = "
            << p.to_string()
            << ",\n"
               "  n = "
            << n.to_string()
            << ",\n"
               "  uv = "
            << uv.to_string()
            << "\n"
               "  t = "
            << t
            << ",\n"
               "  pdf = "
            << pdf
            << "\n"
               "  delta = "
            << delta
            << "\n"
               "]\n";
        return oss.str();
    }
};

template <typename Scalar_> class TDirectionSample : public TPositionSample<Scalar_> {
public:
    typedef Scalar_ Scalar;
    typedef TArray<Scalar, 3, ArrayType::Vector> VectorType;
    typedef TArray<Scalar, 3, ArrayType::Point> PointType;
    typedef TArray<Scalar, 2, ArrayType::Point> PointType2;
    typedef TArray<Scalar, 3, ArrayType::Normal> NormalType;
    typedef TFrame<Scalar> FrameType;
    typedef TPositionSample<Scalar> Base;

    VectorType d; // Unit direction from the reference point to the target shape
    Scalar dist;  // Distance from the reference point to the target shape
    std::shared_ptr<Emitter> emitter = nullptr;

    TDirectionSample() = default;

    /**
     * \brief Create a direct sampling record, which can be used to \a query
     * the density of a surface position with respect to a given reference
     * position.
     *
     * \param si
     *     Surface interaction
     *
     * \param ref
     *     Reference position
     */
    explicit TDirectionSample(const TSurfaceIntersection<Scalar> &si, const TIntersection<Scalar> &ref) : Base(si) {
        VectorType rel = si.p - ref.p;
        dist           = rel.magnitude();
        d              = si.is_valid() ? rel / dist : -si.wi;
        emitter        = si.mesh->get_emitter();
    }

    explicit TDirectionSample(const Point3f &p, const Normal3f &n, const Point2f &uv, const Scalar &t,
                              const Scalar &pdf, const bool &delta, const Vector3f &d, const Scalar &dist,
                              const std::shared_ptr<Emitter> &emitter)
        : Base(p, n, uv, t, pdf, delta), d(d), dist(dist), emitter(emitter) {}

    explicit TDirectionSample(const Base &base) : Base(base) {}

    // Return a human-readable string summary of this frame
    [[nodiscard]] std::string to_string() const {
        std::ostringstream oss;
        oss << "Position Sample[\n"
               "  p = "
            << Base::p.to_string()
            << ",\n"
               "  n = "
            << Base::n.to_string()
            << ",\n"
               "  uv = "
            << Base::uv.to_string()
            << "\n"
               "  t = "
            << Base::t
            << ",\n"
               "  pdf = "
            << Base::pdf
            << "\n"
               "  delta = "
            << Base::delta
            << "\n"
               "  d = "
            << d.to_string()
            << "\n"
               "  dist = "
            << dist
            << "\n"
               "  emitter = "
            << indent(emitter ? emitter->to_string() : "null")
            << "\n"
               "]\n";
        return oss.str();
    }
};

template <typename Scalar_> class TBSDFSample {
public:
    typedef Scalar_ Scalar;
    typedef TArray<Scalar, 3, ArrayType::Vector> VectorType;
    typedef TArray<Scalar, 3, ArrayType::Point> PointType;
    typedef TArray<Scalar, 2, ArrayType::Point> PointType2;
    typedef TArray<Scalar, 3, ArrayType::Normal> NormalType;
    typedef TFrame<Scalar> FrameType;

    VectorType wo; // Normalized outgoing direction in local coordinates
    Scalar pdf;    // Probability density at the sample
    Scalar eta;    // Relative index of refraction in the sampled direction
    bool delta = false;

    explicit TBSDFSample() = default;

    explicit TBSDFSample(const Vector3f &wo) : wo(wo), pdf(0.f), eta(1.f) {}
};

M_NAMESPACE_END

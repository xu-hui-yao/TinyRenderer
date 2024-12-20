#pragma once

#include <core/common.h>
#include <core/fresnel.h>
#include <core/quad.h>
#include <core/warp.h>

M_NAMESPACE_BEGIN
template <typename Scalar_, DeviceType Device_> class TMicrofacetDistribution {
public:
    typedef Scalar_ Scalar;
    static constexpr DeviceType Device = Device_;
    typedef TArray<Scalar, 3, ArrayType::Vector, Device> Vector;
    typedef TArray<Scalar, 2, ArrayType::Vector, Device> Vector2;
    typedef TArray<Scalar, 2, ArrayType::Point, Device> Point2;
    typedef TArray<Scalar, 3, ArrayType::Normal, Device> Normal;
    typedef TFrame<Scalar, Device> Frame;

    M_HOST_DEVICE explicit TMicrofacetDistribution(
        const TMicrofacetDistribution<Scalar, Device> &other) {
        m_alpha = other.m_alpha;
        configure();
    }

    M_HOST_DEVICE explicit TMicrofacetDistribution(Scalar alpha)
        : m_alpha(alpha) {
        configure();
    }

    M_HOST_DEVICE Scalar alpha() const { return m_alpha; }

    M_HOST_DEVICE void set_alpha(Scalar alpha) {
        m_alpha = M_MAX(alpha, Scalar(1e-4));
    }

    M_HOST_DEVICE Scalar eval(const Normal &m) const {
        bool active        = true;
        Scalar alpha_uv    = m_alpha * m_alpha;
        Scalar cos_theta   = Frame::cos_theta(m, active);
        Scalar cos_theta_2 = cos_theta * cos_theta;
        Scalar result;

        Scalar temp1 = m.x() / m_alpha * m.x() / m_alpha;
        Scalar temp2 = m.y() / m_alpha * m.y() / m_alpha;
        Scalar temp3 = temp1 + temp2 + m.z() * m.z();
        result       = Scalar(1) / (M_PI * alpha_uv * temp3 * temp3);

        return result * cos_theta > Scalar(0) ? result : Scalar(0);
    }

    M_HOST_DEVICE Scalar pdf(const Vector &wi, const Normal &m) const {
        bool active   = true;
        Scalar result = eval(m);
        result *=
            smith_g1(wi, m) * abs(wi.dot(m)) / Frame::cos_theta(wi, active);
        return result;
    }

    M_HOST_DEVICE std::pair<Normal, Scalar>
    sample(const Vector &wi, const Point2 &sample, bool &valid) const {
        // Visible normal sampling.
        Scalar sin_phi, cos_phi, cos_theta;

        // Step 1: stretch wi
        Vector wi_p =
            Vector({ m_alpha * wi.x(), m_alpha * wi.y(), wi.z() }).norm(valid);

        std::tie(sin_phi, cos_phi) = Frame::sincos_phi(wi_p);
        cos_theta                  = Frame::cos_theta(wi_p, valid);

        // Step 2: simulate P22_{wi}(slope.x, slope.y, 1, 1)
        Vector2 slope = sample_visible_11(cos_theta, sample);

        // Step 3: rotate & un stretch
        Vector2 slope1(
            { (cos_phi * slope.x() - sin_phi * slope.y()) * m_alpha,
              (sin_phi * slope.x() + cos_phi * slope.y()) * m_alpha });

        // Step 4: compute normal & PDF
        Normal m = Normal({ -slope1.x(), -slope1.y(), Scalar(1) }).norm(valid);
        Scalar pdf = eval(m) * smith_g1(wi, m) * abs(wi.dot(m)) /
                     Frame::cos_theta(wi, valid);

        return { m, pdf };
    }

    // Smith's separable shadowing-masking approximation
    M_HOST_DEVICE Scalar g(const Vector &wi, const Vector &wo,
                           const Normal &m) const {
        return smith_g1(wi, m) * smith_g1(wo, m);
    }

    /**
     * \brief Smith's shadowing-masking function for a single direction
     *
     * \param v
     *     An arbitrary direction
     * \param m
     *     The microfacet normal
     */
    M_HOST_DEVICE Scalar smith_g1(const Vector &v, const Normal &m) const {
        Scalar xy_alpha_2 = m_alpha * v.x() * m_alpha * v.x() +
                            m_alpha * v.y() * m_alpha * v.y();
        Scalar tan_theta_alpha_2 = xy_alpha_2 / v.z() * v.z();
        Scalar result;

        result = Scalar(2) / (Scalar(1) + sqrt(Scalar(1) + tan_theta_alpha_2));

        // Perpendicular incidence -- no shadowing/masking
        if (xy_alpha_2 == Scalar(0)) {
            result = Scalar(1);
        }

        bool active = true;
        // Ensure consistent orientation (can't see the back of the microfacet
        // from the front and vice versa)
        if (v.dot(m) * Frame::cos_theta(v, active) <= Scalar(0)) {
            result = Scalar(0);
        }

        return result;
    }

    M_HOST_DEVICE void scale_alpha(Scalar value) { m_alpha *= value; }

    M_HOST_DEVICE Scalar eval_reflectance(const Vector &wi, Scalar eta) {
        int res = eta > 1 ? 32 : 128;

        auto [nodes, weights] = gauss_legendre<float>(res);

        Scalar result = Scalar(0);

        for (int i = 0; i < res; ++i) {
            for (int j = 0; j < res; ++j) {
                // Map to [0, 1] range
                Scalar u = Scalar(0.5) * (nodes[i] + Scalar(1));
                Scalar v = Scalar(0.5) * (nodes[j] + Scalar(1));
                Point2 sample2({ u, v });
                Scalar weight = weights[i] * weights[j] * Scalar(0.25);

                bool valid   = true;
                auto m       = std::get<0>(sample(wi, sample2, valid));
                Vector wo    = reflect(wi, m);
                Scalar f     = std::get<0>(fresnel(wi.dot(m), eta));
                Scalar smith = smith_g1(wo, m) * f;

                if (wo.z() <= 0.0f || wi.z() <= 0.0f) {
                    smith = 0.0f;
                }
                result += smith * weight;
            }
        }
        return result;
    }

    M_HOST_DEVICE Scalar eval_transmittance(const Vector &wi, Scalar eta) {
        int res = eta > 1 ? 32 : 128;

        auto [nodes, weights] = gauss_legendre<float>(res);

        Scalar result = Scalar(0);

        for (int i = 0; i < res; ++i) {
            for (int j = 0; j < res; ++j) {
                // Map to [0, 1] range
                Scalar u = Scalar(0.5) * (nodes[i] + Scalar(1));
                Scalar v = Scalar(0.5) * (nodes[j] + Scalar(1));
                Point2 sample2({ u, v });
                Scalar weight = weights[i] * weights[j] * Scalar(0.25);

                bool valid = true;
                auto m     = std::get<0>(sample(wi, sample2, valid));
                auto [f, cos_theta_t, eta_it, eta_ti] = fresnel(wi.dot(m), eta);

                Vector wo    = refract(wi, m, cos_theta_t, eta_ti);
                Scalar smith = smith_g1(wo, m) * (1.0f - f);

                if (wo.z() * wi.z() >= Scalar(0)) {
                    smith = 0.0f;
                }

                result += weight * smith;
            }
        }
        return result;
    }

private:
    M_HOST_DEVICE void configure() { m_alpha = M_MAX(m_alpha, Scalar(1e-4)); }

    // \brief Visible normal sampling code for the alpha=1 case
    M_HOST_DEVICE Vector2 sample_visible_11(Scalar cos_theta_i,
                                            const Point2 &sample) const {
        // Choose a projection direction and re-scale the sample
        Point2 p = square_to_uniform_disk_concentric(sample);

        Scalar s = Scalar(0.5) * (Scalar(1) + cos_theta_i);
        p.y()    = lerp(safe_sqrt(Scalar(1) - p.x() * p.x()), p.y(), s);

        // Project onto chosen side of the hemisphere
        Scalar x = p.x(), y = p.y();
        Scalar z = safe_sqrt(Scalar(1) - p.square_magnitude());

        // Convert to slope
        Scalar sin_theta_i = safe_sqrt(Scalar(1) - cos_theta_i * cos_theta_i);
        Scalar norm        = Scalar(1) / (sin_theta_i * y + cos_theta_i * z);
        return Vector2({ cos_theta_i * y - sin_theta_i * z, x }) * norm;
    }

    Scalar m_alpha;
};

M_NAMESPACE_END

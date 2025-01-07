#include <components/bsdf.h>
#include <components/texture.h>
#include <core/frame.h>
#include <core/fresnel.h>
#include <core/intersection.h>
#include <core/microfacet.h>
#include <core/record.h>

M_NAMESPACE_BEGIN
class RoughDielectric : public BSDF {
public:
    explicit RoughDielectric(const PropertyList &properties)
        : m_alpha(nullptr), m_specular_reflectance(nullptr), m_specular_transmittance(nullptr) {
        m_name           = properties.get_string("name", "rough dielectric");
        float ext_ior    = properties.get_float("ext_ior", 1.0);
        float int_ior    = properties.get_float("int_ior", 1.33);
        has_reflection   = properties.get_boolean("has_reflection", true);
        has_transmission = properties.get_boolean("has_transmission", true);
        m_eta            = int_ior / ext_ior;
        m_inv_eta        = ext_ior / int_ior;
    }

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
        if (!m_alpha) {
            throw std::runtime_error("Rough conductor: alpha does not exist");
        }

        m_flags =
            static_cast<BSDFFlags>(static_cast<uint32_t>(EGlossyReflection) | static_cast<uint32_t>(EGlossyTransmission));
    }

    void add_child(const std::shared_ptr<Object> &child) override {
        switch (child->get_class_type()) {
            case ETexture:
                if (!this->m_alpha) {
                    m_alpha = std::dynamic_pointer_cast<Texture>(child);
                } else if (!this->m_specular_reflectance) {
                    m_specular_reflectance = std::dynamic_pointer_cast<Texture>(child);
                } else if (!this->m_specular_transmittance) {
                    m_specular_transmittance = std::dynamic_pointer_cast<Texture>(child);
                } else {
                    throw std::runtime_error(
                        "Rough dielectric only supports alpha and specular_reflectance and specular_transmittance");
                }
                break;

            default:
                throw std::runtime_error("BSDF::add_child(<" + class_type_name(child->get_class_type()) +
                                         ">) is not supported!");
        }
    }

    [[nodiscard]] std::pair<BSDFSample3f, Color3f> sample(const SurfaceIntersection3f &si, float sample1,
                                                          const Point2f &sample2, bool active) const override {
        float cos_theta_i = Frame3f::cos_theta(si.wi, active);
        BSDFSample3f bs(Vector3f({ 0, 0, 0 }));
        active &= cos_theta_i != 0.f;

        MicrofacetDistribution1f distribution(m_alpha->eval(si, active)(0));
        // Trick by Walter et al.: slightly scale the roughness values to reduce importance sampling weights. Not needed
        // for the Heitz and D'Eon sampling technique.
        MicrofacetDistribution1f sample_distribution(distribution);
        sample_distribution.scale_alpha(1.2f - 0.2f * sqrt(abs(cos_theta_i)));

        // Sample the microfacet normal
        Normal3f m;
        std::tie(m, bs.pdf) = sample_distribution.sample(mulsign(si.wi, cos_theta_i), sample2, active);
        active &= bs.pdf != 0.0f;

        auto [f, cos_theta_t, eta_it, eta_ti] = fresnel(si.wi.dot(m), m_eta);

        // Select the lobe to be sampled
        float weight;
        bool selected_r = false, selected_t = false;
        if (has_reflection && has_transmission) {
            selected_r = sample1 <= f && active;
            weight     = 1.0f;
            bs.pdf *= selected_r ? f : 1.0f - f;
        } else if (has_reflection || has_transmission) {
            selected_r = has_reflection && active;
            weight     = has_reflection ? f : 1.0f - f;
        } else {
            return { bs, Color3f(0.0f) };
        }

        selected_t = !selected_r && active;

        bs.eta        = selected_r ? 1.0f : eta_it;
        float dwh_dwo = 0.0f;

        auto result = Color3f(weight);
        if (selected_r) {
            bs.wo   = reflect(si.wi, m);
            dwh_dwo = 1.0f / (4.0f * bs.wo.dot(m));
            if (m_specular_reflectance) {
                result *= m_specular_reflectance->eval(si, active);
            }
        }

        if (selected_t) {
            bs.wo        = refract(si.wi, m, cos_theta_t, eta_ti);
            float temp   = si.wi.dot(m) + bs.eta * bs.wo.dot(m);
            dwh_dwo      = bs.eta * bs.eta * bs.wo.dot(m) / (temp * temp);
            float factor = eta_ti;
            result *= factor * factor;
            if (m_specular_transmittance) {
                result *= m_specular_transmittance->eval(si, active);
            }
        }

        result *= distribution.smith_g1(bs.wo, m);

        bs.pdf *= abs(dwh_dwo);
        bs.delta = false;

        return { bs, result };
    }

    [[nodiscard]] Color3f eval(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        float cos_theta_i = Frame3f::cos_theta(si.wi, active), cos_theta_o = Frame3f::cos_theta(wo, active);

        active &= cos_theta_i != 0.0f;

        bool reflect = cos_theta_i * cos_theta_o > 0.0f;

        // Determine the relative index of refraction
        float eta     = cos_theta_i > 0.f ? m_eta : m_inv_eta;
        float inv_eta = cos_theta_i > 0.f ? m_inv_eta : m_eta;

        // Compute the half-vector
        Normal3f m = Normal3f(si.wi + wo * (reflect ? 1.0f : eta)).norm(active);

        // Ensure that the half-vector points into the same hemisphere as the macrosurface normal
        m = mulsign(m, Frame3f::cos_theta(m, active));

        MicrofacetDistribution1f distribution(m_alpha->eval(si, active)(0));

        // Evaluate the microfacet normal distribution
        float d = distribution.eval(m);

        // Fresnel factor
        float f = std::get<0>(fresnel(si.wi.dot(m), m_eta));

        // Smith's shadow-masking function
        float g = distribution.g(si.wi, wo, m);

        Color3f result(0);

        bool eval_r = has_reflection && reflect && active;
        bool eval_t = has_transmission && !reflect && active;

        if (eval_r) {
            Color3f value(f * d * g / (4.0f * abs(cos_theta_i)));
            result = value;
            if (m_specular_reflectance) {
                result *= m_specular_reflectance->eval(si, active);
            }
        }

        if (eval_t) {
            /* Missing term in the original paper: account for the solid angle
               compression when tracing radiance -- this is necessary for
               bidirectional methods. */
            float scale = inv_eta * inv_eta;
            float temp  = si.wi.dot(m) + eta * wo.dot(m);

            Color3f value(
                abs(scale * (1.0f - f) * d * g * eta * eta * si.wi.dot(m) * wo.dot(m) / (cos_theta_i * temp * temp)));
            result = value;
            if (m_specular_transmittance) {
                result *= m_specular_transmittance->eval(si, active);
            }
        }

        return result;
    }

    [[nodiscard]] float pdf(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        float cos_theta_i = Frame3f::cos_theta(si.wi, active), cos_theta_o = Frame3f::cos_theta(wo, active);

        active &= cos_theta_i != 0.0f;

        bool reflect = cos_theta_i * cos_theta_o > 0.0f;
        active &= (has_reflection && reflect) || (has_transmission && !reflect);

        // Determine the relative index of refraction
        float eta = cos_theta_i > 0.f ? m_eta : m_inv_eta;

        // Compute the half-vector
        Normal3f m = Normal3f(si.wi + wo * (reflect ? 1.0f : eta)).norm(active);

        // Ensure that the half-vector points into the same hemisphere as the macrosurface normal
        m = mulsign(m, Frame3f::cos_theta(m, active));

        /* Filter cases where the micro/macro-surface don't agree on the side.
           This logic is evaluated in smith_g1() called as part of the eval()
           and sample() methods and needs to be replicated in the probability
           density computation as well. */
        active &=
            si.wi.dot(m) * Frame3f::cos_theta(si.wi, active) > 0.f && wo.dot(m) * Frame3f::cos_theta(wo, active) > 0.f;

        // Jacobian of the half-direction mapping
        float temp    = si.wi.dot(m) + eta * wo.dot(m);
        float dwh_dwo = reflect ? 1.0f / (4.0f * wo.dot(m)) : eta * eta * wo.dot(m) / (temp * temp);

        /* Construct the microfacet distribution matching the
           roughness values at the current surface position. */
        MicrofacetDistribution1f sample_distribution(m_alpha->eval(si, active)(0));

        /* Trick by Walter et al.: slightly scale the roughness values to
           reduce importance sampling weights. Not needed for the
           Heitz and D'Eon sampling technique. */
        sample_distribution.scale_alpha(1.2f - 0.2f * sqrt(abs(Frame3f::cos_theta(si.wi, active))));

        // Evaluate the microfacet model sampling density function
        float prob = sample_distribution.pdf(mulsign(si.wi, Frame3f::cos_theta(si.wi, active)), m);

        if (has_transmission && has_reflection) {
            float f = std::get<0>(fresnel(si.wi.dot(m), m_eta));
            prob *= reflect ? f : 1.0f - f;
        }

        return active ? prob * abs(dwh_dwo) : 0.0f;
    }

    // Return a human-readable summary
    [[nodiscard]] std::string to_string() const override {
        return "RoughDielectric[\n"
               "  alpha = " +
               indent(m_alpha->to_string(), 2) +
               "\n"
               "  specular_reflectance = " +
               indent(m_specular_reflectance->to_string(), 2) +
               "\n"
               "  specular_transmittance = " +
               indent(m_specular_transmittance->to_string(), 2) +
               "\n"
               "  eta = " +
               std::to_string(m_eta) +
               "\n"
               "]";
    }

    [[nodiscard]] EClassType get_class_type() const override { return EBSDF; }

private:
    float m_eta, m_inv_eta;
    std::shared_ptr<Texture> m_alpha;
    std::shared_ptr<Texture> m_specular_reflectance;
    std::shared_ptr<Texture> m_specular_transmittance;
    bool has_reflection, has_transmission;
};

REGISTER_CLASS(RoughDielectric, "roughdielectric");

M_NAMESPACE_END

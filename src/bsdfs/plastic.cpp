#include <components/bsdf.h>
#include <components/texture.h>
#include <core/frame.h>
#include <core/fresnel.h>
#include <core/intersection.h>
#include <core/record.h>
#include <core/warp.h>

M_NAMESPACE_BEGIN
class SmoothPlastic : public BSDF {
public:
    explicit SmoothPlastic(const PropertyList &properties)
        : m_diffuse_reflectance(nullptr), m_specular_reflectance(nullptr) {
        m_name        = properties.get_string("name", "smooth plastic");
        float ext_ior = properties.get_float("ext_ior", 1.0f);
        float int_ior = properties.get_float("int_ior", 1.49f); // Default to polypropylene
        m_eta         = int_ior / ext_ior;
        m_nonlinear   = properties.get_boolean("nonlinear", false);

        m_specular_sampling_weight = 0.5f; // Default to equal sampling weight
        m_inv_eta_2                = 1.0f / (m_eta * m_eta);

        m_fdr_int = fresnel_diffuse_reflectance(1.0f / m_eta);
        m_fdr_ext = fresnel_diffuse_reflectance(m_eta);
    }

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
        if (!m_diffuse_reflectance) {
            throw std::runtime_error("Rough plastic requires diffuse reflectance");
        }
    }

    void add_child(const std::shared_ptr<Object> &child) override {
        switch (child->get_class_type()) {
            case ETexture:
                if (!m_diffuse_reflectance) {
                    m_diffuse_reflectance = std::dynamic_pointer_cast<Texture>(child);
                } else if (!m_specular_reflectance) {
                    m_specular_reflectance = std::dynamic_pointer_cast<Texture>(child);
                } else {
                    throw std::runtime_error(
                        "Rough plastic only supports diffuse_reflectance and specular_reflectance");
                }
                break;

            default:
                throw std::runtime_error("BSDF::add_child(<" + class_type_name(child->get_class_type()) +
                                         ">) is not supported!");
        }
    }

    [[nodiscard]] std::pair<BSDFSample3f, Color3f> sample(const SurfaceIntersection3f &si, float sample1,
                                                          const Point2f &sample2, bool active) const override {
        // Compute the cosine of the angle of incidence
        float cos_theta_i = Frame3f::cos_theta(si.wi, active);
        active &= cos_theta_i > 0.f;

        BSDFSample3f bs(Vector3f({0, 0, 0}));
        Color3f result(0.f);

        // Early exit if the sample is invalid
        if (!active) {
            return { bs, result };
        }

        // Fresnel term for the incident angle
        float f_i = std::get<0>(fresnel(cos_theta_i, m_eta));

        float prob_specular = f_i * m_specular_sampling_weight;
        float prob_diffuse  = (1.0f - f_i) * (1.0f - m_specular_sampling_weight);
        prob_specular       = prob_specular / (prob_specular + prob_diffuse);
        prob_diffuse        = 1.0f - prob_specular;

        // Specular Reflection Sampling
        bool sample_specular = active && sample1 < prob_specular;
        bool sample_diffuse  = active && !sample_specular;

        bs.eta   = 1.f; // No refraction
        bs.pdf   = 0.f;
        bs.delta = true;

        // Sample specular component
        if (sample_specular) {
            bs.wo  = reflect(si.wi); // Perfect reflection
            bs.pdf = prob_specular;
            Color3f value(f_i / bs.pdf);
            if (m_specular_reflectance) {
                value *= m_specular_reflectance->eval(si, active);
            }
            result = value;
        }

        // Sample diffuse component
        if (sample_diffuse) {
            bs.wo         = square_to_cosine_hemisphere(sample2); // Diffuse reflection
            bs.pdf        = prob_diffuse * square_to_cosine_hemisphere_pdf(bs.wo);
            float f_o     = std::get<0>(fresnel(Frame3f::cos_theta(bs.wo, active), m_eta));
            Color3f value = m_diffuse_reflectance->eval(si, active);
            value /= Color3f(1.0f) - (m_nonlinear ? value * m_fdr_int : Color3f(m_fdr_int));
            value *= m_inv_eta_2 * (1.0f - f_i) * (1.0f - f_o) / prob_diffuse;
            result = value;
        }

        return { bs, active ? result : Color3f(0.f) };
    }

    [[nodiscard]] Color3f eval(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        float cos_theta_i = Frame3f::cos_theta(si.wi, active), cos_theta_o = Frame3f::cos_theta(wo, active);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        if (!active) {
            return Color3f(0);
        }

        float f_i = std::get<0>(fresnel(cos_theta_i, m_eta));
        float f_o = std::get<0>(fresnel(cos_theta_o, m_eta));

        Color3f diffuse = m_diffuse_reflectance->eval(si, active);
        diffuse /= Color3f(1.0f) - (m_nonlinear ? diffuse * m_fdr_int : Color3f(m_fdr_int));

        diffuse *= square_to_cosine_hemisphere_pdf(wo) * m_inv_eta_2 * (1.0f - f_i) * (1.0f - f_o);

        return active ? diffuse : Color3f(0.f);
    }

    [[nodiscard]] float pdf(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        float cos_theta_i = Frame3f::cos_theta(si.wi, active), cos_theta_o = Frame3f::cos_theta(wo, active);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        if (!active) {
            return 0.0f;
        }

        float f_i           = std::get<0>(fresnel(cos_theta_i, m_eta));
        float prob_specular = f_i * m_specular_sampling_weight;
        float prob_diffuse  = (1.0f - f_i) * (1.0f - m_specular_sampling_weight);
        prob_diffuse        = prob_diffuse / (prob_diffuse + prob_specular);

        float pdf = square_to_cosine_hemisphere_pdf(wo) * prob_diffuse;

        return active ? pdf : 0.0f;
    }

    [[nodiscard]] std::string to_string() const override {
        std::ostringstream oss;
        oss << "SmoothPlastic[" << std::endl
            << "  diffuse_reflectance = " << indent(m_diffuse_reflectance ? m_diffuse_reflectance->to_string() : "null")
            << "," << std::endl
            << "  specular_reflectance = "
            << indent(m_specular_reflectance ? m_specular_reflectance->to_string() : "null") << "," << std::endl
            << "  eta = " << m_eta << "," << std::endl
            << "  nonlinear = " << m_nonlinear << std::endl
            << "]";
        return oss.str();
    }

    [[nodiscard]] EClassType get_class_type() const override { return EBSDF; }

private:
    std::shared_ptr<Texture> m_diffuse_reflectance;
    std::shared_ptr<Texture> m_specular_reflectance;
    float m_eta;
    float m_inv_eta_2;
    float m_specular_sampling_weight;
    bool m_nonlinear;
    float m_fdr_int;
    float m_fdr_ext;
};

REGISTER_CLASS(SmoothPlastic, "plastic");

M_NAMESPACE_END

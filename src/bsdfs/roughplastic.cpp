#include <components/bsdf.h>
#include <components/texture.h>
#include <core/frame.h>
#include <core/fresnel.h>
#include <core/intersection.h>
#include <core/microfacet.h>
#include <core/record.h>

M_NAMESPACE_BEGIN
class RoughPlastic : public BSDF {
public:
    explicit RoughPlastic(const PropertyList &properties)
        : m_diffuse_reflectance(nullptr), m_specular_reflectance(nullptr) {
        m_name                     = properties.get_string("name", "rough plastic");
        float ext_ior              = properties.get_float("ext_ior", 1.0);
        float int_ior              = properties.get_float("int_ior", 1.49);
        m_eta                      = int_ior / ext_ior;
        m_nonlinear                = properties.get_boolean("nonlinear", false);
        m_alpha                    = properties.get_float("alpha", 0.1);
        m_rough_transmittance_res  = properties.get_integer("rough_transmittance_res", 64);
        m_inv_eta_2                = 1.0f / (m_eta * m_eta);
        m_internal_reflectance     = 0.0f;
        m_specular_sampling_weight = 0.0f;
    }

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
        if (!m_diffuse_reflectance) {
            throw std::runtime_error("Rough plastic requires diffuse reflectance");
        }

        float d_mean = m_diffuse_reflectance->mean().luminance();
        float s_mean = 1.0f;

        if (m_specular_reflectance) {
            s_mean = m_specular_reflectance->mean().luminance();
        }

        m_specular_sampling_weight = s_mean / (d_mean + s_mean);

        // Create a Microfacet distribution
        MicrofacetDistribution1f distribution(m_alpha);

        // Precompute rough reflectance (vectorized approach using std::vector)
        std::vector<float> mu(m_rough_transmittance_res);
        std::vector<float> zero(m_rough_transmittance_res, 0);
        std::vector<Vector3f> wi(m_rough_transmittance_res);

        // Initialize mu, representing cosine of the incident angle
        float step = 1.0f / static_cast<float>(m_rough_transmittance_res - 1);
        for (int i = 0; i < m_rough_transmittance_res; ++i) {
            mu[i] = M_MAX(1e-6f, i * step); // Ensure no zero value for cosine
            wi[i] = Vector3f(std::sqrt(1.0f - mu[i] * mu[i]), zero[i], mu[i]);
        }

        // Compute the external transmittance for each direction
        m_external_transmittance.resize(m_rough_transmittance_res);
        for (int i = 0; i < m_rough_transmittance_res; ++i) {
            m_external_transmittance[i] = distribution.eval_transmittance(wi[i], m_eta);
        }

        // Compute internal reflectance
        m_internal_reflectance = 0.0f;
        for (int i = 0; i < m_rough_transmittance_res; ++i) {
            m_internal_reflectance += distribution.eval_reflectance(wi[i], 1.0f / m_eta) * wi[i].z();
        }

        m_internal_reflectance = m_internal_reflectance * 2.0f / static_cast<float>(m_rough_transmittance_res);

        m_flags =
            static_cast<BSDFFlags>(static_cast<uint32_t>(EGlossyReflection) | static_cast<uint32_t>(EDiffuseReflection));
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

        // Interpolate external transmittance for the incoming direction
        float t_i = lerp_gather(m_external_transmittance, cos_theta_i, m_rough_transmittance_res);

        // Determine the probabilities for specular and diffuse reflection
        float prob_specular = (1.f - t_i) * m_specular_sampling_weight;
        float prob_diffuse  = t_i * (1.f - m_specular_sampling_weight);

        // Assuming has both specular and diffuse
        prob_specular /= prob_specular + prob_diffuse;

        // Sample specular or diffuse based on probabilities
        bool sample_specular = active && sample1 < prob_specular;
        bool sample_diffuse  = active && !sample_specular;

        bs.eta = 1.f;

        // Sample the specular component using the microfacet distribution
        if (sample_specular) {
            MicrofacetDistribution1f distribution(m_alpha);
            Normal3f m = std::get<0>(distribution.sample(si.wi, sample2, active));

            // Reflect the incoming direction around the microfacet normal
            bs.wo = reflect(si.wi, m);
        }

        // Sample the diffuse component
        if (sample_diffuse) {
            bs.wo = square_to_cosine_hemisphere(sample2);
        }

        // Compute the PDF for the sampled direction
        bs.pdf = pdf(si, bs.wo, active);
        active &= bs.pdf > 0.f;
        result   = eval(si, bs.wo, active);
        bs.delta = false;

        return { bs, active ? result / bs.pdf : Color3f(0.f) };
    }

    [[nodiscard]] Color3f eval(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        float cos_theta_i = Frame3f::cos_theta(si.wi, active), cos_theta_o = Frame3f::cos_theta(wo, active);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        if (!active) {
            return Color3f(0);
        }

        Color3f value(0.0f);

        // Specular
        MicrofacetDistribution1f distribution(m_alpha);
        Normal3f h = Normal3f(wo + si.wi).norm(active);
        float d    = distribution.eval(h);
        float f    = std::get<0>(fresnel(si.wi.dot(h), m_eta));
        float g    = distribution.g(si.wi, wo, h);
        value      = Color3f(f * d * g / (4.0f * cos_theta_i));
        if (m_specular_reflectance) {
            value *= m_specular_reflectance->eval(si, active);
        }

        // Diffuse
        float t_i       = lerp_gather(m_external_transmittance, cos_theta_i, m_rough_transmittance_res);
        float t_o       = lerp_gather(m_external_transmittance, cos_theta_o, m_rough_transmittance_res);
        Color3f diffuse = m_diffuse_reflectance->eval(si, active);
        diffuse /= Color3f(1.0f) - (m_nonlinear ? diffuse * m_internal_reflectance : Color3f(m_internal_reflectance));
        value += diffuse * (static_cast<float>(M_INV_PI) * m_inv_eta_2 * cos_theta_o * t_i * t_o);

        return active ? value : Color3f(0.f);
    }

    [[nodiscard]] float pdf(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        float cos_theta_i = Frame3f::cos_theta(si.wi, active), cos_theta_o = Frame3f::cos_theta(wo, active);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        if (!active) {
            return 0.0f;
        }

        float t_i = lerp_gather(m_external_transmittance, cos_theta_i, m_rough_transmittance_res);

        float prob_specular = (1.f - t_i) * m_specular_sampling_weight;
        float prob_diffuse  = t_i * (1.f - m_specular_sampling_weight);
        prob_specular       = prob_specular / (prob_specular + prob_diffuse);
        prob_diffuse        = 1.f - prob_specular;

        Normal3f h = Normal3f(wo + si.wi).norm(active);

        MicrofacetDistribution1f distribution(m_alpha);
        float result = distribution.eval(h) * distribution.smith_g1(si.wi, h) / (4.0f * cos_theta_i);
        result *= prob_specular;

        result += prob_diffuse * square_to_cosine_hemisphere_pdf(wo);

        return result;
    }

    [[nodiscard]] std::string to_string() const override {
        std::ostringstream oss;
        oss << "RoughPlastic[" << std::endl
            << "  diffuse_reflectance = " << indent(m_diffuse_reflectance ? m_diffuse_reflectance->to_string() : "null")
            << "," << std::endl
            << "  specular_reflectance = "
            << indent(m_specular_reflectance ? m_specular_reflectance->to_string() : "null") << "," << std::endl
            << "  alpha = " << m_alpha << "," << std::endl
            << "  eta = " << m_eta << "," << std::endl
            << "  nonlinear = " << m_nonlinear << std::endl
            << "]";
        return oss.str();
    }

    [[nodiscard]] EClassType get_class_type() const override { return EBSDF; }

private:
    [[nodiscard]] static float lerp_gather(const std::vector<float> &data, float x, int size) {
        // Scale `x` to fit into the range of the data array (0 to size-1)
        x *= static_cast<float>(size - 1);
        int index = M_MIN(static_cast<int>(x), size - 2);

        // Gather the two closest values
        float v0 = data[index];
        float v1 = data[index + 1];

        // Return the linearly interpolated value
        return lerp(v0, v1, x - static_cast<float>(index));
    }

    std::shared_ptr<Texture> m_diffuse_reflectance;
    std::shared_ptr<Texture> m_specular_reflectance;
    float m_eta;
    float m_inv_eta_2;
    float m_alpha;
    float m_specular_sampling_weight;
    std::vector<float> m_external_transmittance;
    float m_internal_reflectance;
    bool m_nonlinear;
    int m_rough_transmittance_res;
};

REGISTER_CLASS(RoughPlastic, "roughplastic");

M_NAMESPACE_END

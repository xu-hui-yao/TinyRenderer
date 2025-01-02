#include <components/bsdf.h>
#include <components/texture.h>
#include <core/frame.h>
#include <core/fresnel.h>
#include <core/intersection.h>
#include <core/microfacet.h>
#include <core/record.h>

M_NAMESPACE_BEGIN
class ThinDielectric : public BSDF {
public:
    explicit ThinDielectric(const PropertyList &properties) {
        m_name           = properties.get_string("name", "rough dielectric");
        float ext_ior    = properties.get_float("ext_ior", 1.0);
        float int_ior    = properties.get_float("int_ior", 1.33);
        has_reflection   = properties.get_boolean("has_reflection", true);
        has_transmission = properties.get_boolean("has_transmission", true);
        m_eta            = int_ior / ext_ior;
    }

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
    }

    void add_child(const std::shared_ptr<Object> &child) override {
        switch (child->get_class_type()) {
            case ETexture:
                if (!this->m_specular_reflectance) {
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
        float r = std::get<0>(fresnel(abs(Frame3f::cos_theta(si.wi, active)), m_eta));

        // Account for internal reflections: r = r + trt + tr^3t + ..
        r *= 2.f / (1.f + r);
        float t = 1.0f - r;

        // Select the lobe to be sampled
        BSDFSample3f bs(Vector3f({0, 0, 0}));
        bool selected_r;
        float weight = 0.0f;
        if (has_reflection && has_transmission) {
            selected_r = sample1 <= r && active;
            bs.pdf     = selected_r ? r : t;
            weight     = 1.0f;
        } else if (has_reflection || has_transmission) {
            selected_r = has_reflection && active;
            bs.pdf     = 1.0f;
            weight     = has_reflection ? r : t;
        } else {
            return { bs, Color3f(0.0f) };
        }

        Color3f result(weight);
        if (selected_r) {
            if (m_specular_reflectance) {
                result *= m_specular_reflectance->eval(si, active);
            }
        } else {
            if (m_specular_transmittance) {
                result *= m_specular_transmittance->eval(si, active);
            }
        }

        bs.wo    = selected_r ? reflect(si.wi) : -si.wi;
        bs.eta   = 1.0f;
        bs.delta = true;

        return { bs, active ? result : Color3f(0.0f) };
    }

    [[nodiscard]] Color3f eval(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        return Color3f(0.0f);
    }

    [[nodiscard]] float pdf(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        return 0.0f;
    }

    // Return a human-readable summary
    [[nodiscard]] std::string to_string() const override {
        return "ThinDielectric[\n"
               "  eta = " +
               std::to_string(m_eta) +
               "\n"
               "]";
    }

    [[nodiscard]] EClassType get_class_type() const override { return EBSDF; }

private:
    float m_eta;
    bool has_reflection, has_transmission;
    std::shared_ptr<Texture> m_specular_reflectance;
    std::shared_ptr<Texture> m_specular_transmittance;
};

REGISTER_CLASS(ThinDielectric, "thindielectric");

M_NAMESPACE_END

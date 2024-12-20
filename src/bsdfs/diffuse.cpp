#include <components/bsdf.h>
#include <components/texture.h>
#include <core/frame.h>
#include <core/intersection.h>
#include <core/record.h>
#include <core/warp.h>

M_NAMESPACE_BEGIN
class Diffuse : public BSDF {
public:
    explicit Diffuse(const Color3f &color) {
        m_name        = "diffuse";
        m_reflectance = std::make_shared<Constant>(color);
    }

    explicit Diffuse(const PropertyList &properties) : m_reflectance(nullptr) {
        m_name = properties.get_string("name", "diffuse");
    }

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif

        if (!m_reflectance) {
            throw std::runtime_error("Diffuse: No reflection provided");
        }
    }

    void add_child(const std::shared_ptr<Object> &child) override {
        switch (child->get_class_type()) {
            case ETexture:
                if (this->m_reflectance)
                    throw std::runtime_error("Diffuse: tried to register multiple Texture instances!");
                m_reflectance = std::dynamic_pointer_cast<Texture>(child);
                break;

            default:
                throw std::runtime_error("BSDF::add_child(<" + class_type_name(child->get_class_type()) +
                                         ">) is not supported!");
        }
    }

    [[nodiscard]] std::pair<BSDFSample3f, Color3f> sample(const SurfaceIntersection3f &si, float /* sample1 */,
                                                          const Point2f &sample2, bool active) const override {
        float cos_theta_i = Frame3f::cos_theta(si.wi, active);
        BSDFSample3f bs(Vector3f(0));

        active &= cos_theta_i > 0.f;

        if (!active) {
            return { bs, Color3f(0) };
        }

        bs.wo    = square_to_cosine_hemisphere(sample2);
        bs.pdf   = square_to_cosine_hemisphere_pdf(bs.wo);
        bs.eta   = 1.f;
        bs.delta = false;

        Color3f value = m_reflectance->eval(si, active);

        return { bs, active && bs.pdf > 0.f ? value : Color3f(0) };
    }

    [[nodiscard]] Color3f eval(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        float cos_theta_i = Frame3f::cos_theta(si.wi, active), cos_theta_o = Frame3f::cos_theta(wo, active);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        Color3f value = m_reflectance->eval(si, active) * M_INV_PI * cos_theta_o;

        return active ? value : Color3f(0);
    }

    [[nodiscard]] float pdf(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        float cos_theta_i = Frame3f::cos_theta(si.wi, active), cos_theta_o = Frame3f::cos_theta(wo, active);

        float pdf = square_to_cosine_hemisphere_pdf(wo);

        return cos_theta_i > 0.f && cos_theta_o > 0.f ? pdf : 0.f;
    }

    // Return a human-readable summary
    [[nodiscard]] std::string to_string() const override {
        return "Diffuse[\n"
               "  albedo = " +
               indent(m_reflectance->to_string(), 2) +
               "\n"
               "]";
    }

    [[nodiscard]] EClassType get_class_type() const override { return EBSDF; }

private:
    std::shared_ptr<Texture> m_reflectance;
};
REGISTER_CLASS(Diffuse, "diffuse");

M_NAMESPACE_END

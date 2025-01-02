#include <components/bsdf.h>
#include <components/texture.h>
#include <core/frame.h>
#include <core/intersection.h>
#include <core/record.h>

M_NAMESPACE_BEGIN
class BumpMap : public BSDF {
public:
    explicit BumpMap(const PropertyList &properties) {
        m_nested_texture = nullptr;
        m_nested_bsdf    = nullptr;
        m_scale          = properties.get_float("scale", 1.0f);
        m_name           = properties.get_string("name", "bumpmap");
    }

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
        if (!m_nested_texture)
            throw std::runtime_error("BumpMap: Texture is not provided.");
        if (!m_nested_bsdf)
            throw std::runtime_error("BumpMap: Nested BSDF is not provided.");
    }

    void add_child(const std::shared_ptr<Object> &child) override {
        switch (child->get_class_type()) {
            case ETexture:
                if (!m_nested_texture)
                    m_nested_texture = std::dynamic_pointer_cast<Texture>(child);
                else
                    throw std::runtime_error("BumpMap: Texture already provided.");
                break;
            case EBSDF:
                if (!m_nested_bsdf)
                    m_nested_bsdf = std::dynamic_pointer_cast<BSDF>(child);
                else
                    throw std::runtime_error("BumpMap: Nested BSDF already provided.");
                break;
            default:
                throw std::runtime_error("BumpMap::add_child(<" + class_type_name(child->get_class_type()) +
                                         ">) is not supported!");
        }
    }

    [[nodiscard]] std::pair<BSDFSample3f, Color3f> sample(const SurfaceIntersection3f &si, float sample1,
                                                          const Point2f &sample2, bool active) const override {
        SurfaceIntersection3f perturbed_si = si;
        perturbed_si.shading_frame         = frame(si, active);
        perturbed_si.wi                    = perturbed_si.to_local(si.wi);

        auto [bs, value] = m_nested_bsdf->sample(perturbed_si, sample1, sample2, active);
        active &= (value != 0.0f).any();
        Vector3f perturbed_wo = perturbed_si.to_world(bs.wo);
        active &= Frame3f::cos_theta(perturbed_wo, active) * Frame3f::cos_theta(bs.wo, active) > 0;
        bs.wo = perturbed_wo;

        return { bs, active ? value : Color3f(0.0f) };
    }

    [[nodiscard]] Color3f eval(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        SurfaceIntersection3f perturbed_si = si;
        perturbed_si.shading_frame         = frame(si, active);
        perturbed_si.wi                    = perturbed_si.to_local(si.wi);
        Vector3f perturbed_wo              = perturbed_si.to_local(wo);

        active &= Frame3f::cos_theta(perturbed_wo, active) * Frame3f::cos_theta(wo, active) > 0;

        return active ? m_nested_bsdf->eval(perturbed_si, perturbed_wo, active) : Color3f(0.0f);
    }

    [[nodiscard]] float pdf(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        SurfaceIntersection3f perturbed_si = si;
        perturbed_si.shading_frame         = frame(si, active);
        perturbed_si.wi                    = perturbed_si.to_local(si.wi);
        Vector3f perturbed_wo              = perturbed_si.to_local(wo);

        active &= Frame3f::cos_theta(perturbed_wo, active) * Frame3f::cos_theta(wo, active) > 0;

        return m_nested_bsdf->pdf(perturbed_si, perturbed_wo, active);
    }

    [[nodiscard]] std::string to_string() const override {
        return "BumpMap[\n"
               "  nested_texture = " +
               indent(m_nested_texture ? m_nested_texture->to_string() : "null") +
               "\n"
               "  nested_bsdf = " +
               indent(m_nested_bsdf ? m_nested_bsdf->to_string() : "null") +
               "\n"
               "  scale = " +
               std::to_string(m_scale) + "\n]";
    }

    [[nodiscard]] EClassType get_class_type() const override { return EBSDF; }

private:
    [[nodiscard]] Frame3f frame(const SurfaceIntersection3f &si, bool &active) const {
        // Evaluate texture gradient
        Vector2f grad_uv = m_nested_texture->eval_1_grad(si, active) * m_scale;

        // Compute perturbed differential geometry
        Vector3f dp_du = si.dp_du + si.shading_frame.n * (grad_uv.x() - si.shading_frame.n.dot(si.dp_du));
        Vector3f dp_dv = si.dp_dv + si.shading_frame.n * (grad_uv.y() - si.shading_frame.n.dot(si.dp_dv));

        // Convert to small rotation from original shading frame
        Frame3f result;
        result.n = dp_du.cross(dp_dv).norm(active);

        // Gram-schmidt orthogonality to compute local shading frame
        result.s = (si.dp_du - result.n * dp_du.dot(result.n)).norm(active);
        result.t = result.n.cross(result.s);

        return result;
    }

    std::shared_ptr<Texture> m_nested_texture;
    std::shared_ptr<BSDF> m_nested_bsdf;
    float m_scale;
};

REGISTER_CLASS(BumpMap, "bumpmap");

M_NAMESPACE_END

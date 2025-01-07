#include <components/bsdf.h>
#include <components/texture.h>
#include <core/intersection.h>
#include <core/record.h>

M_NAMESPACE_BEGIN
class Mask : public BSDF {
public:
    explicit Mask(const PropertyList &properties) {
        m_opacity     = nullptr;
        m_nested_bsdf = nullptr;
        m_name        = properties.get_string("name", "mask");
    }

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
        if (!m_opacity)
            throw std::runtime_error("Mask: Opacity texture is not provided.");
        if (!m_nested_bsdf)
            throw std::runtime_error("Mask: Nested BSDF is not provided.");

        m_flags = m_nested_bsdf->get_flag();
    }

    void add_child(const std::shared_ptr<Object> &child) override {
        switch (child->get_class_type()) {
            case ETexture:
                if (!m_opacity)
                    m_opacity = std::dynamic_pointer_cast<Texture>(child);
                else
                    throw std::runtime_error("Mask: Opacity texture already provided.");
                break;
            case EBSDF:
                if (!m_nested_bsdf)
                    m_nested_bsdf = std::dynamic_pointer_cast<BSDF>(child);
                else
                    throw std::runtime_error("Mask: Nested BSDF already provided.");
                break;
            default:
                throw std::runtime_error("Mask::add_child(<" + class_type_name(child->get_class_type()) +
                                         ">) is not supported!");
        }
    }

    [[nodiscard]] std::pair<BSDFSample3f, Color3f> sample(const SurfaceIntersection3f &si, float sample1,
                                                          const Point2f &sample2, bool active) const override {
        BSDFSample3f bs;
        Color3f value(1.0f);

        float opacity = m_opacity->eval(si, active)(0);

        bs.wo  = -si.wi;
        bs.eta = 1.0f;
        bs.pdf = 1.0f - opacity;

        if (active && sample1 < opacity) {
            sample1 /= opacity;
            auto [nested_bs, nested_value] = m_nested_bsdf->sample(si, sample1 / opacity, sample2, active);
            bs                             = nested_bs;
            bs.pdf *= opacity;
            value = nested_value;
        }

        return { bs, value };
    }

    [[nodiscard]] Color3f eval(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        float opacity = m_opacity->eval(si, active)(0);
        return m_nested_bsdf->eval(si, wo, active) * opacity;
    }

    [[nodiscard]] float pdf(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        float opacity = m_opacity->eval(si, active)(0);
        return m_nested_bsdf->pdf(si, wo, active) * opacity;
    }

    [[nodiscard]] std::string to_string() const override {
        return "Mask[\n"
               "  opacity=" +
               indent(m_opacity ? m_opacity->to_string() : "null") +
               "\n"
               "  nested_bsdf=" +
               indent(m_nested_bsdf ? m_nested_bsdf->to_string() : "null") + "\n]";
    }

    [[nodiscard]] EClassType get_class_type() const override { return EBSDF; }

private:
    std::shared_ptr<Texture> m_opacity;
    std::shared_ptr<BSDF> m_nested_bsdf;
};

REGISTER_CLASS(Mask, "mask");

M_NAMESPACE_END

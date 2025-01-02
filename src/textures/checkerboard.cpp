#include <components/texture.h>
#include <core/intersection.h>

M_NAMESPACE_BEGIN
class CheckerBoard : public Texture {
public:
    explicit CheckerBoard(const PropertyList &properties) : Texture(true) {
        m_name    = properties.get_string("name", "checkerboard");
        m_color0  = properties.get_color("color0", Color3f(0.4f));
        m_color1  = properties.get_color("color1", Color3f(0.2f));
        m_scale_u = properties.get_float("scale_u", 1.0f);
        m_scale_v = properties.get_float("scale_v", 1.0f);
    }

    void add_child(const std::shared_ptr<Object> &child) override {}

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct" << class_type_name(get_class_type()) << std::endl;
#endif
    }

    Color3f eval(const SurfaceIntersection3f &si, bool &active) override {
        Point2f uv = si.uv;
        uv.x() *= m_scale_u;
        uv.y() *= m_scale_v;
        auto mask = uv - uv.get_floor() > 0.5f;

        bool m0 = mask.x() == mask.y() && active;
        bool m1 = !m0 && active;

        if (m0) {
            return m_color0;
        }

        if (m1) {
            return m_color1;
        }

        return Color3f(0.0f);
    }

    float eval_1(const SurfaceIntersection3f &si, bool &active) override { return eval(si, active).luminance(); }

    Vector2f eval_1_grad(const SurfaceIntersection3f &si, bool &active) override { return Vector2f({ 0.0f, 0.0f }); }

    Color3f mean() override { return (m_color0 + m_color1) / 2.0f; }

    [[nodiscard]] std::string to_string() const override {
        return "CheckerBoard[\n"
               "  color0 = " +
               m_color0.to_string() +
               "\n"
               "  color1 = " +
               m_color1.to_string() +
               "\n"
               "]";
    }

private:
    Color3f m_color0;
    Color3f m_color1;
    float m_scale_u, m_scale_v;
};

REGISTER_CLASS(CheckerBoard, "checkerboard")

M_NAMESPACE_END

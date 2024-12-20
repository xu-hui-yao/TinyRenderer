#include <components/texture.h>

M_NAMESPACE_BEGIN
class Constant : public Texture {
public:
    explicit Constant(const Color3f &color) : Texture(false) {
        m_name     = "constant";
        m_channels = 3;
        m_color    = color;
        m_color_1  = color.luminance();
    }

    explicit Constant(const PropertyList &properties) : Texture(false) {
        m_name     = properties.get_string("name", "constant");
        m_channels = properties.get_integer("channels", 3);
        if (m_channels == 3) {
            m_color = properties.get_color("color", Color3f({ 1.0f, 1.0f, 1.0f }));
        } else if (m_channels == 1) {
            m_color_1 = properties.get_float("color", 1.0f);
        } else {
            throw std::runtime_error("Unsupported channels");
        }
    }

    void add_child(const std::shared_ptr<Object> &child) override {}

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct" << class_type_name(get_class_type()) << std::endl;
#endif
    }

    Color3f eval(const SurfaceIntersection3f &si, bool &active) override {
        if (m_channels == 3) {
            return m_color;
        } else {
            return Color3f(m_color_1);
        }
    }

    float eval_1(const SurfaceIntersection3f &si, bool &active) override {
        if (m_channels == 1) {
            return m_color_1;
        } else {
            return m_color.luminance();
        }
    }

    Color3f mean() override { return m_color; }

    Vector2f eval_1_grad(const SurfaceIntersection3f &si, bool &active) override { return Vector2f({ 0.0f, 0.0f }); }

    [[nodiscard]] std::string to_string() const override {
        return "Constant[\n"
               "  color = " +
               m_color.to_string() +
               "\n"
               "]";
    }

private:
    Color3f m_color;
    float m_color_1;
    int m_channels;
};

REGISTER_CLASS(Constant, "constant")

M_NAMESPACE_END

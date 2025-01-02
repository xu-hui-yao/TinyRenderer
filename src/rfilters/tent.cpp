#include <components/rfilter.h>
#include <parse/property_list.h>

M_NAMESPACE_BEGIN

class TentFilter : public ReconstructionFilter {
public:
    explicit TentFilter(const PropertyList &property_list) {
        m_radius     = property_list.get_float("radius");
        m_inv_radius = 1.0f / m_radius;
        m_name       = property_list.get_string("name", "tent");
    }

    [[nodiscard]] float eval(float x) const override { return M_MAX(0.0f, 1.0f - abs(x * m_inv_radius)); }

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
    }

    [[nodiscard]] std::string to_string() const override {
        return std::string("TentFilter[radius=") + std::to_string(m_radius) + std::string("]");
    }

private:
    float m_inv_radius;
};

REGISTER_CLASS(TentFilter, "tent")

M_NAMESPACE_END
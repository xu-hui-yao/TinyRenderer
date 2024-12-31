#include <components/emitter.h>
#include <components/mesh.h>
#include <components/texture.h>
#include <core/common.h>
#include <core/intersection.h>
#include <core/record.h>
#include <core/spectrum.h>

M_NAMESPACE_BEGIN
/**
 * This plugin implements an area light, i.e. a light source that emits
 * diffuse illumination from the exterior of an arbitrary shape.
 * Since the emission profile of an area light is completely diffuse, it
 * has the same apparent brightness regardless of the observer's viewing
 * direction. Furthermore, since it occupies a nonzero amount of space, an
 * area light generally causes scene objects to cast soft shadows.
 */
class Area : public Emitter {
public:
    explicit Area(const PropertyList &properties) : m_radiance(nullptr) {
        m_name = properties.get_string("name", "area");
    }

    void add_child(const std::shared_ptr<Object> &child) override {
        switch (child->get_class_type()) {
            case ETexture: {
                m_radiance = std::dynamic_pointer_cast<Texture>(child);
            } break;
            default:
                throw std::runtime_error("Area::add_child(<" + class_type_name(child->get_class_type()) +
                                         ">) is not supported!");
        }
    }

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
    }

    [[nodiscard]] std::pair<DirectionSample3f, Color3f>
    sample_direction(const Intersection3f &it, const Point2f &sample, bool &active) const override {
        DirectionSample3f ds = std::dynamic_pointer_cast<Mesh>(parent)->sample_direction(it, sample, active);
        active &= ds.d.dot(ds.n) < 0.f && ds.pdf != 0.f;
        SurfaceIntersection3f si = SurfaceIntersection3f(ds);
        auto spec                = m_radiance->eval(si, active) / ds.pdf;
        ds.emitter               = std::dynamic_pointer_cast<Emitter>(std::make_shared<Area>(*this));

        return { ds, active ? spec : Color3f(0.0f) };
    }

    [[nodiscard]] float pdf_direction(const Intersection3f &it, const DirectionSample3f &ds,
                                      bool &active) const override {
        float dp = ds.d.dot(ds.n);
        active &= dp < 0.0f;
        float value = std::dynamic_pointer_cast<Mesh>(parent)->pdf_direction(it, ds, active);
        return active ? value : 0.0f;
    }

    [[nodiscard]] std::pair<PositionSample3f, float> sample_position(const Point2f &sample,
                                                                     bool &active) const override {
        PositionSample3f ps = std::dynamic_pointer_cast<Mesh>(parent)->sample_position(sample, active);
        float weight        = active && ps.pdf > 0.f ? 1.0f / ps.pdf : 0.0f;
        return { ps, weight };
    }

    [[nodiscard]] float pdf_position(const PositionSample3f &ps, bool &active) const override {
        float pdf = std::dynamic_pointer_cast<Mesh>(parent)->pdf_position(ps, active);
        return active ? pdf : 0.0f;
    }

    [[nodiscard]] Color3f eval(const SurfaceIntersection3f &si, bool &active) const override {
        return Frame3f::cos_theta(si.wi, active) > 0.0f ? m_radiance->eval(si, active) : Color3f(0.0f);
    }

    [[nodiscard]] std::string to_string() const override {
        return std::string("Area[\n  radiance = ") + indent(m_radiance->to_string()) + "\n]";
    }

    [[nodiscard]] bool is_spatial_varying() const override { return m_radiance->is_spatial_varying(); }

private:
    std::shared_ptr<Texture> m_radiance;
};

REGISTER_CLASS(Area, "area")

M_NAMESPACE_END

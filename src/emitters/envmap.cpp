#include <components/bitmap.h>
#include <components/emitter.h>
#include <components/texture.h>
#include <core/common.h>
#include <core/distribution.h>
#include <core/intersection.h>
#include <core/record.h>
#include <core/spectrum.h>
#include <core/warp.h>

M_NAMESPACE_BEGIN
class EnvironmentMap : public Emitter {
public:
    explicit EnvironmentMap(const PropertyList &properties) : m_radiance(nullptr) {
        m_name                   = properties.get_string("name", "envmap");
        m_to_world               = properties.get_transform("to_world", Transform4f());
        m_bounding_sphere_center = properties.get_point("center", Point3f(0, 0, 0));
        m_bounding_sphere_radius = properties.get_float("radius", 1000.0f);
    }

    ~EnvironmentMap() override = default;

    void add_child(const std::shared_ptr<Object> &child) override {
        switch (child->get_class_type()) {
            case ETexture: {
                m_radiance = std::dynamic_pointer_cast<Texture>(child);
            } break;
            default:
                throw std::runtime_error("EnvMap::add_child(<" + class_type_name(child->get_class_type()) +
                                         ">) is not supported!");
        }
    }

    void set_scene(const std::shared_ptr<Scene> &scene) override {
        const auto &bounding_box = scene->get_accel()->get_bounding_box();
        m_bounding_sphere_center = bounding_box.get_center();
        float length             = (bounding_box.get_max() - bounding_box.get_min()).magnitude();
        m_bounding_sphere_radius = length / 2.0f * (1.0f + static_cast<float>(M_EPSILON));
    }

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
        if (is_spatial_varying()) {
            auto bitmap = std::dynamic_pointer_cast<Bitmap>(m_radiance)->get_data();
            int height  = bitmap->get_rows();
            int width   = bitmap->get_cols();
            std::vector<std::vector<float>> data;
            data.resize(height);
            for (auto &row : data) {
                row.resize(width);
            }
            for (int i = 0; i < height; i++) {
                for (int j = 0; j < width; j++) {
                    Color3f color(
                        { bitmap->operator()(i, j, 0), bitmap->operator()(i, j, 1), bitmap->operator()(i, j, 2) });
                    data[i][j] = color.luminance();
                }
            }
            m_distribution = std::make_shared<HierarchicalDistribution2f>(data);
        }
    }
    /**
     *
     *             \phi|
     *              \  |
     *               \ |
     * x              \|
     * <---------------+---------------
     *                 |
     *                 |
     *                 |
     *                 | z
     */
    [[nodiscard]] std::pair<DirectionSample3f, Color3f> sample_direction(const Intersection3f &it,
                                                                         const Point2f &sample, bool active) override {
        if (is_spatial_varying()) {
            auto [uv, pdf] = m_distribution->sample(sample);
            active &= pdf > 0.0f;

            DirectionSample3f ds;
            ds.uv.x() = uv.x(); // col
            ds.uv.y() = uv.y(); // row

            // The half-texel shift is purely a UV <-> direction mapping
            // convention (it is what eval()/pdf_direction() undo with their
            // `uv.x() -= 0.5/(width-1)`), so it belongs ONLY in the
            // uv -> direction conversion below. It must NOT be baked into
            // ds.uv, because `pdf` was evaluated by the distribution at the
            // UNSHIFTED uv: sampling the radiance at uv + delta while
            // dividing by the density at uv makes the estimator
            // radiance/pdf inconsistent. On a smooth envmap that is
            // invisible, but on a real one (where the sun spans only a few
            // texels) it biases the result by >10%.
            int width       = std::dynamic_pointer_cast<Bitmap>(m_radiance)->get_cols();
            float u_shifted = uv.x() + 0.5f / static_cast<float>(width - 1);

            float theta = ds.uv.y() * static_cast<float>(M_PI);   // [0, pi]
            float phi   = u_shifted * static_cast<float>(M_PI) * 2.0f; // [0, 2pi]

            Vector3f d({ sin(phi) * sin(theta), cos(theta), -cos(phi) * sin(theta) });

            float radius = M_MAX(m_bounding_sphere_radius, (it.p - m_bounding_sphere_center).magnitude());
            float dist   = 2.0f * radius;

            float inv_sin_theta = 1.0f / safe_sqrt(M_MAX(d.x() * d.x() + d.z() * d.z(), static_cast<float>(M_EPSILON)));
            d                   = m_to_world * d;

            ds.p     = it.p + d * dist;
            ds.n     = -d;
            ds.t     = it.t;
            ds.pdf   = active ? pdf * inv_sin_theta * static_cast<float>(1.0f / (2.0f * M_PI * M_PI)) : 0.0f;
            ds.delta = false;
            ds.d     = d;
            ds.dist  = dist;

            SurfaceIntersection3f its;
            its.uv         = ds.uv;
            Color3f weight = m_radiance->eval(its, active);
            ds.emitter     = std::dynamic_pointer_cast<Emitter>(shared_from_this());
            return { ds, active ? weight / ds.pdf : Color3f(0.0f) };
        } else {
            DirectionSample3f ds;
            float radius = M_MAX(m_bounding_sphere_radius, (it.p - m_bounding_sphere_center).magnitude());
            float dist   = 2.0f * radius;

            Vector3f d = square_to_uniform_sphere(sample);
            ds.p       = it.p + d * dist;
            ds.n       = -d;
            ds.t       = it.t;
            ds.pdf     = square_to_uniform_sphere_pdf(d);
            ds.delta   = false;
            ds.d       = d;
            ds.dist    = dist;

            SurfaceIntersection3f its;
            Color3f weight = m_radiance->eval(its, active);
            ds.emitter     = std::dynamic_pointer_cast<Emitter>(shared_from_this());
            return { ds, active ? weight / ds.pdf : Color3f(0.0f) };
        }
    }

    [[nodiscard]] float pdf_direction(const Intersection3f &it, const DirectionSample3f &ds,
                                      bool active) const override {
        if (is_spatial_varying()) {
            Vector3f d = m_to_world.inverse() * ds.d;
            auto phi   = atan2(d.x(), -d.z());
            auto theta = safe_acos(d.y());
            if (phi < 0.0f) {
                phi += 2 * M_PI;
            }
            auto uv   = Point2f({ phi * static_cast<float>(M_INV_TWOPI), theta * static_cast<float>(M_INV_PI) });
            int width = std::dynamic_pointer_cast<Bitmap>(m_radiance)->get_cols();
            uv.x() -= 0.5f / static_cast<float>(width - 1);
            uv -= uv.get_floor();

            float inv_sin_theta = 1.0f / safe_sqrt(M_MAX(d.x() * d.x() + d.z() * d.z(), static_cast<float>(M_EPSILON)));

            // Must be pdf() (the normalized DENSITY), not eval() (the
            // un-normalized cell value): this has to return the same density
            // that sample_direction() divides by, otherwise MIS weights this
            // strategy against a differently-scaled quantity.
            return m_distribution->pdf(uv) * inv_sin_theta * static_cast<float>(1.0f / (2.0f * M_PI * M_PI));
        } else {
            Vector3f d = m_to_world.inverse() * ds.d;
            return square_to_uniform_sphere_pdf(d);
        }
    }

    [[nodiscard]] std::pair<PositionSample3f, float> sample_position(const Point2f &sample,
                                                                     bool active) const override {
        throw std::runtime_error("Not implemented");
    }

    [[nodiscard]] float pdf_position(const PositionSample3f &ps, bool active) const override {
        throw std::runtime_error("Not implemented");
    }

    [[nodiscard]] Color3f eval(const SurfaceIntersection3f &si, bool active) const override {
        if (is_spatial_varying()) {
            Vector3f d = m_to_world.inverse() * -si.wi;
            auto phi   = atan2(d.x(), -d.z());
            auto theta = safe_acos(d.y());
            auto uv = Point2f({ phi * static_cast<float>(M_INV_TWOPI), theta * static_cast<float>(M_INV_PI) });
            auto bitmap = std::dynamic_pointer_cast<Bitmap>(m_radiance)->get_data();
            int width   = bitmap->get_cols();
            uv.x() -= 0.5f / static_cast<float>(width - 1);
            uv -= uv.get_floor();

            SurfaceIntersection3f its;
            its.uv      = uv;
            auto result = m_radiance->eval(its, active);
            return active ? result : Color3f(0.0f);
        } else {
            SurfaceIntersection3f its;
            auto result = m_radiance->eval(its, active);
            return active ? result : Color3f(0.0f);
        }
    }

    [[nodiscard]] std::string to_string() const override {
        return std::string("EnvironmentMap[\n  radiance = ") + indent(m_radiance->to_string()) + "\n]";
    }

    [[nodiscard]] bool is_spatial_varying() const override { return m_radiance->is_spatial_varying(); }

    [[nodiscard]] GPULight to_gpu_light(GPUSceneBuilder &builder, uint32_t /*mesh_id*/) const override {
        GPULight light;
        light.type                   = GPULightType::Envmap;
        light.radiance_tex           = builder.add_texture(m_radiance);
        light.to_world               = m_to_world;
        light.bounding_sphere_center = m_bounding_sphere_center;
        light.bounding_sphere_radius = m_bounding_sphere_radius;
        light.spatial_varying        = is_spatial_varying();
        if (light.spatial_varying && m_distribution) {
            light.distribution_width      = m_distribution->get_cols();
            light.distribution_height     = m_distribution->get_rows();
            light.distribution_luminance  = m_distribution->level0_data();
        }
        return light;
    }

private:
    std::shared_ptr<Texture> m_radiance;
    Transform4f m_to_world;
    std::shared_ptr<HierarchicalDistribution2f> m_distribution;
    Point3f m_bounding_sphere_center;
    float m_bounding_sphere_radius;
};
REGISTER_CLASS(EnvironmentMap, "envmap")

M_NAMESPACE_END

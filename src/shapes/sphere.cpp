#include <array>
#include <components/bsdf.h>
#include <components/mesh.h>
#include <core/transform.h>

M_NAMESPACE_BEGIN
class Sphere : public Mesh {
public:
    explicit Sphere(const PropertyList &properties) {
        bool valid   = true;
        this->m_name = properties.get_string("name", "sphere");

        float m_radius   = properties.get_float("radius", 1.0f);
        Point3f m_center = properties.get_point("center", Point3f(0.0f, 0.0f, 0.0f));
        int resolution   = properties.get_integer("resolution", 72);

        if (resolution < 3)
            throw std::runtime_error("Resolution must be at least 3");

        // Generate vertices, normals, and UVs
        for (int i = 0; i <= resolution; ++i) {
            // Latitude angle
            float theta = static_cast<float>(i) / static_cast<float>(resolution) * static_cast<float>(M_PI);

            for (int j = 0; j <= resolution; ++j) {
                // Longitude angle
                float phi = static_cast<float>(j) / static_cast<float>(resolution) * 2.0f * static_cast<float>(M_PI);

                float x = m_radius * sin(theta) * cos(phi);
                float y = m_radius * cos(theta);
                float z = m_radius * sin(theta) * sin(phi);

                Point3f p  = m_center + Point3f(x, y, z);
                Normal3f n = Normal3f(x, y, z).norm(valid);
                Point2f uv(static_cast<float>(j) / static_cast<float>(resolution),
                           static_cast<float>(i) / static_cast<float>(resolution));

                this->vertices.push_back(p);
                this->normals.push_back(n);
                this->uvs.push_back(uv);
                this->bounding_box.expand_by(this->vertices.back());
            }
        }

        // Generate faces (triangles)
        for (int i = 0; i < resolution; ++i) {
            for (int j = 0; j < resolution; ++j) {
                int p0 = i * (resolution + 1) + j;
                int p1 = p0 + 1;
                int p2 = (i + 1) * (resolution + 1) + j;
                int p3 = p2 + 1;

                this->faces.emplace_back(p0, p1, p2);
                this->faces.emplace_back(p1, p3, p2);
            }
        }
    }

    explicit Sphere(const Point3f &center, float radius, int resolution, bool flip_normal) {
        bool valid   = true;
        this->m_name = "sphere";

        if (resolution < 3)
            throw std::runtime_error("Resolution must be at least 3");

        // Generate vertices, normals, and UVs
        for (int i = 0; i <= resolution; ++i) {
            // Latitude angle
            float theta = static_cast<float>(i) / static_cast<float>(resolution) * static_cast<float>(M_PI);

            for (int j = 0; j <= resolution; ++j) {
                // Longitude angle
                float phi = static_cast<float>(j) / static_cast<float>(resolution) * 2.0f * static_cast<float>(M_PI);

                float x = radius * sin(theta) * cos(phi);
                float y = radius * cos(theta);
                float z = radius * sin(theta) * sin(phi);

                Point3f p  = center + Point3f(x, y, z);
                Normal3f n = Normal3f(x, y, z).norm(valid);
                if (flip_normal) {
                    n = -n;
                }
                Point2f uv(static_cast<float>(j) / static_cast<float>(resolution),
                           static_cast<float>(i) / static_cast<float>(resolution));

                this->vertices.push_back(p);
                this->normals.push_back(n);
                this->uvs.push_back(uv);
                this->bounding_box.expand_by(this->vertices.back());
            }
        }

        // Generate faces (triangles)
        for (int i = 0; i < resolution; ++i) {
            for (int j = 0; j < resolution; ++j) {
                int p0 = i * (resolution + 1) + j;
                int p1 = p0 + 1;
                int p2 = (i + 1) * (resolution + 1) + j;
                int p3 = p2 + 1;

                this->faces.emplace_back(p0, p1, p2);
                this->faces.emplace_back(p1, p3, p2);
            }
        }
    }

    [[nodiscard]] std::string to_string() const override {
        std::ostringstream oss;
        oss << "Sphere[" << std::endl
            << "  vertices = " << vertices.size() << "," << std::endl
            << "  faces = " << faces.size() << "," << std::endl
            << "  normals = " << normals.size() << "," << std::endl
            << std::string("  BSDF = ") + (bsdf ? indent(bsdf->to_string()) : std::string("null")) << std::endl
            << std::string("  Emitter = ") + (emitter ? indent(emitter->to_string()) : std::string("null")) << "\n]";
        return oss.str();
    }
};

REGISTER_CLASS(Sphere, "sphere");

M_NAMESPACE_END

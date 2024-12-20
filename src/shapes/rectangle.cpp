#include <array>
#include <components/mesh.h>
#include <components/bsdf.h>
#include <core/transform.h>

M_NAMESPACE_BEGIN

class Rectangle : public Mesh {
public:
    explicit Rectangle(const PropertyList &properties) {
        this->m_name = properties.get_string("name", "rectangle");

        bool flip_normal     = properties.get_boolean("flip_normal", false);
        Transform4f to_world = properties.get_transform("to_world", Transform4f());

        std::vector<Point3f> vertices = {
            { -1.f, -1.f, 0.f }, { 1.f, -1.f, 0.f }, { 1.f, 1.f, 0.f }, { -1.f, 1.f, 0.f }
        };
        std::vector<Normal3f> normals = { { 0.f, 0.f, 1.f }, { 0.f, 0.f, 1.f }, { 0.f, 0.f, 1.f }, { 0.f, 0.f, 1.f } };
        if (flip_normal) {
            for (auto &n : normals) {
                n = -n;
            }
        }
        std::vector<Point2f> tex_coords = { { 0.f, 0.f }, { 1.f, 0.f }, { 1.f, 1.f }, { 0.f, 1.f } };
        std::vector<Point3i> faces      = { { 0, 1, 2 }, { 0, 2, 3 } };

        this->vertices.resize(vertices.size());
        this->normals.resize(normals.size());
        this->uvs.resize(tex_coords.size());
        this->faces.resize(faces.size());

        bool valid = true;
        for (size_t i = 0; i < vertices.size(); ++i) {
            this->vertices[i] = to_world * vertices[i];
            this->normals[i]  = (to_world * normals[i]).norm(valid);
            this->uvs[i]      = tex_coords[i];
            this->bounding_box.expand_by(this->vertices[i]);
        }

        for (size_t i = 0; i < faces.size(); ++i) {
            this->faces[i] = faces[i];
        }
    }

    [[nodiscard]] std::string to_string() const override {
        std::ostringstream oss;
        oss << "Rectangle[" << std::endl
            << "  vertices = " << vertices.size() << "," << std::endl
            << "  faces = " << faces.size() << "," << std::endl
            << std::string("  BSDF = ") + (bsdf ? indent(bsdf->to_string()) : std::string("null")) << std::endl
            << std::string("  Emitter = ") + (emitter ? indent(emitter->to_string()) : std::string("null")) << "\n]";
        return oss.str();
    }
};

REGISTER_CLASS(Rectangle, "rectangle");

M_NAMESPACE_END

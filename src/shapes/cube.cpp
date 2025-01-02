#include <array>
#include <components/mesh.h>
#include <core/transform.h>

M_NAMESPACE_BEGIN
class Cube : public Mesh {
public:
    explicit Cube(const PropertyList &properties) {
        this->m_name = properties.get_string("name", "cube");

        Transform4f to_world = properties.get_transform("to_world", Transform4f());

        const std::vector<Point3f> vertices = { { 1, -1, -1 }, { 1, -1, 1 },  { -1, -1, 1 },  { -1, -1, -1 },
                                                { 1, 1, -1 },  { -1, 1, -1 }, { -1, 1, 1 },   { 1, 1, 1 },
                                                { 1, -1, -1 }, { 1, 1, -1 },  { 1, 1, 1 },    { 1, -1, 1 },
                                                { 1, -1, 1 },  { 1, 1, 1 },   { -1, 1, 1 },   { -1, -1, 1 },
                                                { -1, -1, 1 }, { -1, 1, 1 },  { -1, 1, -1 },  { -1, -1, -1 },
                                                { 1, 1, -1 },  { 1, -1, -1 }, { -1, -1, -1 }, { -1, 1, -1 } };

        const std::vector<Normal3f> normals = { { 0, -1, 0 }, { 0, -1, 0 }, { 0, -1, 0 }, { 0, -1, 0 }, { 0, 1, 0 },
                                                { 0, 1, 0 },  { 0, 1, 0 },  { 0, 1, 0 },  { 1, 0, 0 },  { 1, 0, 0 },
                                                { 1, 0, 0 },  { 1, 0, 0 },  { 0, 0, 1 },  { 0, 0, 1 },  { 0, 0, 1 },
                                                { 0, 0, 1 },  { -1, 0, 0 }, { -1, 0, 0 }, { -1, 0, 0 }, { -1, 0, 0 },
                                                { 0, 0, -1 }, { 0, 0, -1 }, { 0, 0, -1 }, { 0, 0, -1 } };

        const std::vector<Point2f> tex_coords = { { 0, 1 }, { 1, 1 }, { 1, 0 }, { 0, 0 }, { 0, 1 }, { 1, 1 },
                                                  { 1, 0 }, { 0, 0 }, { 0, 1 }, { 1, 1 }, { 1, 0 }, { 0, 0 },
                                                  { 0, 1 }, { 1, 1 }, { 1, 0 }, { 0, 0 }, { 0, 1 }, { 1, 1 },
                                                  { 1, 0 }, { 0, 0 }, { 0, 1 }, { 1, 1 }, { 1, 0 }, { 0, 0 } };

        const std::vector<Point3i> faces = { { 0, 1, 2 },    { 3, 0, 2 },    { 4, 5, 6 },    { 7, 4, 6 },
                                             { 8, 9, 10 },   { 11, 8, 10 },  { 12, 13, 14 }, { 15, 12, 14 },
                                             { 16, 17, 18 }, { 19, 16, 18 }, { 20, 21, 22 }, { 23, 20, 22 } };

        this->vertices.resize(vertices.size());
        // this->normals.resize(normals.size());
        this->uvs.resize(tex_coords.size());
        this->faces.resize(faces.size());

        bool valid = true;
        for (size_t i = 0; i < vertices.size(); ++i) {
            this->vertices[i] = to_world * vertices[i];
            // this->normals[i]  = (to_world * normals[i]).norm(valid);
            this->uvs[i] = tex_coords[i];
            this->bounding_box.expand_by(this->vertices[i]);
        }

        for (size_t i = 0; i < faces.size(); ++i) {
            this->faces[i] = faces[i];
        }
    }

    [[nodiscard]] std::string to_string() const override {
        std::ostringstream oss;
        oss << "Cube[" << std::endl
            << "  vertices = " << vertices.size() << "," << std::endl
            << "  faces = " << faces.size() << "," << std::endl
            << "]";
        return oss.str();
    }
};

REGISTER_CLASS(Cube, "cube");

M_NAMESPACE_END

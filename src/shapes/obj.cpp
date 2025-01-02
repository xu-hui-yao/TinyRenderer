#include <components/mesh.h>
#include <core/array.h>
#include <filesystem/resolver.h>
#include <fstream>
#include <iostream>
#include <unordered_map>

M_NAMESPACE_BEGIN
class WavefrontOBJ : public Mesh {
public:
    explicit WavefrontOBJ(const PropertyList &properties) {
        typedef std::unordered_map<OBJVertex, uint32_t, OBJVertexHash> VertexMap;
        filesystem::path filename = get_file_resolver()->resolve(filesystem::path(properties.get_string("filename")));
        std::ifstream is(filename.str());
        if (!is)
            throw std::runtime_error("Unable to open OBJ file " + filename.str());

        Transform4f transform = properties.get_transform("to_world", Transform4f());
        std::vector<Point3f> positions;
        std::vector<Point2f> tex_coords;
        std::vector<Normal3f> normals;
        std::vector<uint32_t> indices;
        std::vector<OBJVertex> vertices;
        VertexMap vertex_map;

        std::string line_str;
        while (std::getline(is, line_str)) {
            std::istringstream line(line_str);

            std::string prefix;
            line >> prefix;

            if (prefix == "v") {
                Point3f p;
                line >> p.x() >> p.y() >> p.z();
                p = transform * p;
                bounding_box.expand_by(p);
                positions.emplace_back(p);
            } else if (prefix == "vt") {
                Point2f tc;
                line >> tc.x() >> tc.y();
                tex_coords.emplace_back(tc);
            } else if (prefix == "vn") {
                bool valid = true;
                Normal3f n;
                line >> n.x() >> n.y() >> n.z();
                normals.emplace_back((transform * n.norm(valid)).norm(valid));
                if (!valid) {
                    std::cerr << "Invalid normal vector provided!" << std::endl;
                }
            } else if (prefix == "f") {
                std::string v1, v2, v3, v4;
                line >> v1 >> v2 >> v3 >> v4;
                OBJVertex vertex[6];
                int n_vertices = 3;

                vertex[0] = OBJVertex(v1);
                vertex[1] = OBJVertex(v2);
                vertex[2] = OBJVertex(v3);

                if (!v4.empty()) {
                    /* This is a quad, split into two triangles */
                    vertex[3]  = OBJVertex(v4);
                    vertex[4]  = vertex[0];
                    vertex[5]  = vertex[2];
                    n_vertices = 6;
                }
                /* Convert to an indexed vertex list */
                for (int i = 0; i < n_vertices; ++i) {
                    const OBJVertex &v           = vertex[i];
                    VertexMap::const_iterator it = vertex_map.find(v);
                    if (it == vertex_map.end()) {
                        vertex_map[v] = static_cast<uint32_t>(vertices.size());
                        indices.push_back(static_cast<uint32_t>(vertices.size()));
                        vertices.push_back(v);
                    } else {
                        indices.push_back(it->second); // Record vertex index of the triangle
                    }
                }
            }
        }

        this->faces.resize(static_cast<long long>(indices.size() / 3));
        for (uint32_t i = 0; i < faces.size(); ++i) {
            this->faces[i] = Point3i(static_cast<int>(indices[i * 3]), static_cast<int>(indices[i * 3 + 1]),
                                     static_cast<int>(indices[i * 3 + 2]));
        }

        this->vertices.resize(static_cast<long long>(vertices.size()));
        for (uint32_t i = 0; i < vertices.size(); ++i)
            this->vertices[i] = Point3f(positions.at(vertices[i].p - 1));

        if (!normals.empty()) {
            this->normals.resize(static_cast<long long>(vertices.size()));
            for (uint32_t i = 0; i < vertices.size(); ++i)
                this->normals[i] = normals.at(vertices[i].n - 1);
        }

        if (!tex_coords.empty()) {
            this->uvs.resize(static_cast<long long>(vertices.size()));
            for (uint32_t i = 0; i < vertices.size(); ++i)
                this->uvs[i] = tex_coords.at(vertices[i].uv - 1);
        }

        this->m_name = properties.get_string("name", filename.str());
    }

private:
    // Vertex indices used by the OBJ format
    struct OBJVertex {
        uint32_t p  = static_cast<uint32_t>(-1);
        uint32_t n  = static_cast<uint32_t>(-1);
        uint32_t uv = static_cast<uint32_t>(-1);

        explicit OBJVertex() = default;

        explicit OBJVertex(const std::string &string) {
            std::vector<std::string> tokens = tokenize(string, "/", true);

            if (tokens.empty() || tokens.size() > 3)
                throw std::runtime_error("Invalid vertex data: " + string);

            p = to_uint(tokens[0]);

            if (tokens.size() >= 2 && !tokens[1].empty())
                uv = to_uint(tokens[1]);

            if (tokens.size() >= 3 && !tokens[2].empty())
                n = to_uint(tokens[2]);
        }

        bool operator==(const OBJVertex &v) const { return v.p == p && v.n == n && v.uv == uv; }
    };

    // Hash function for OBJVertex
    struct OBJVertexHash {
        std::size_t operator()(const OBJVertex &v) const {
            size_t hash = std::hash<uint32_t>()(v.p);
            hash        = hash * 37 + std::hash<uint32_t>()(v.uv);
            hash        = hash * 37 + std::hash<uint32_t>()(v.n);
            return hash;
        }
    };
};

REGISTER_CLASS(WavefrontOBJ, "obj");

M_NAMESPACE_END

#include <components/accelerate.h>
#include <components/mesh.h>
#include <core/frame.h>
#include <core/intersection.h>

M_NAMESPACE_BEGIN
class OctreeNode {
public:
    BoundingBox3f bbox;                                                 // The bounding box for the current node
    std::vector<std::pair<std::shared_ptr<Mesh>, uint32_t>> primitives; // Meshes in the leaf node
    std::shared_ptr<OctreeNode> children[8];                            // Children nodes (if any)
    bool is_leaf;

    explicit OctreeNode(const BoundingBox3f &bbox) : bbox(bbox), is_leaf(false) {
        std::fill(std::begin(children), std::end(children), nullptr);
    }
};

class OctreeAccel : public Accel {
public:
    explicit OctreeAccel(const PropertyList &properties) {
        max_depth              = properties.get_integer("max_depth", 8);
        max_triangles_per_leaf = properties.get_integer("leaf_max", 10);
        m_name                 = properties.get_string("name", "octree");
    }

    void add_mesh(const std::shared_ptr<Mesh> &mesh) override {
        for (uint32_t i = 0; i < mesh->get_triangle_count(); ++i) {
            primitives.emplace_back(mesh, i);
            bounding_box.expand_by(mesh->get_bounding_box(i));
        }
    }

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
        root = build_tree(bounding_box, primitives, 0);
    }

    [[nodiscard]] const BoundingBox3f &get_bounding_box() const override { return bounding_box; }

    bool ray_intersect(const Ray3f &ray, SurfaceIntersection3f &its, bool shadow_ray) const override {
        Ray3f ray_(ray);
        return ray_intersect_node(ray_, its, shadow_ray, root);
    }

    [[nodiscard]] bool ray_test(const Ray3f &ray) const override { return ray_test_node(ray, root); }

    [[nodiscard]] std::string to_string() const override { return "Octree Accelerate\n" + node_to_string(root, 0); }

private:
    int max_depth;
    int max_triangles_per_leaf;
    std::shared_ptr<OctreeNode> root;
    std::vector<std::pair<std::shared_ptr<Mesh>, uint32_t>> primitives;

    // Helper function to recursively traverse the octree and build the string representation
    [[nodiscard]] static std::string node_to_string(const std::shared_ptr<OctreeNode> &node, int depth) {
        if (!node)
            return "";

        std::string indent_str = std::string(depth * 4, ' '); // Indentation based on depth

        std::ostringstream oss;
        oss << indent_str << "Node: {\n";

        oss << indent_str << "  bounding box: " << indent(node->bbox.to_string(), 2 + depth * 4) << "\n";
        if (node->is_leaf) {
            oss << indent_str << "  Meshes: " << node->primitives.size() << " triangles\n";
        } else {
            oss << indent_str << "  Children: {\n";
            for (int i = 0; i < 8; ++i) {
                if (node->children[i]) {
                    oss << indent_str << "    Child " << depth << ", " << i << ": \n"
                        << node_to_string(node->children[i], depth + 1);
                }
            }
            oss << indent_str << "  }\n";
        }

        oss << indent_str << "}\n";

        return oss.str();
    }

    // Recursively builds the octree
    std::shared_ptr<OctreeNode>
    build_tree(const BoundingBox3f &bbox, const std::vector<std::pair<std::shared_ptr<Mesh>, uint32_t>> &primitive_list,
               int depth) {
        if (primitive_list.size() <= max_triangles_per_leaf || depth > max_depth) {
            auto node        = std::make_shared<OctreeNode>(bbox);
            node->is_leaf    = true;
            node->primitives = primitive_list; // Store meshes directly in the leaf node

            for (const auto &prim : primitive_list) {
                node->bbox.expand_by(prim.first->get_bounding_box(prim.second));
            }
            return node;
        }

        auto node = std::make_shared<OctreeNode>(bbox);

        // Recursively build child nodes
        for (int i = 0; i < 8; ++i) {
            std::vector<std::pair<std::shared_ptr<Mesh>, uint32_t>> child_primitive_list;
            auto child_bbox = get_child_bbox(bbox, i);

            // Use triangle's bounding box, not just a vertex check
            for (auto &primitive : primitive_list) {
                std::shared_ptr<Mesh> mesh = primitive.first;
                uint32_t triangle_index    = primitive.second;

                // Get the bounding box of the triangle
                BoundingBox3f triangle_bbox = mesh->get_bounding_box(triangle_index);
                // Check if the triangle's bounding box intersects the child node's bounding box
                if (child_bbox.overlaps(triangle_bbox)) {
                    child_primitive_list.emplace_back(mesh, triangle_index);
                }
            }

            // Only build child node if it contains primitives
            if (!child_primitive_list.empty()) {
                node->children[i] = build_tree(child_bbox, child_primitive_list, depth + 1);
            }
        }

        return node;
    }

    // Get the bounding box for the i-th child of the current node
    [[nodiscard]] static BoundingBox3f get_child_bbox(const BoundingBox3f &bbox, int child_idx) {
        Point3f min_p  = bbox.get_min();
        Point3f max_p  = bbox.get_max();
        Point3f center = bbox.get_center();

        // The child bounding box is defined by its center and its corner
        for (int i = 0; i < 3; ++i) {
            if (child_idx & 1 << i) {
                min_p(i) = center(i); // Set the max bound for this dimension
            } else {
                max_p(i) = center(i); // Set the min bound for this dimension
            }
        }

        return { min_p, max_p };
    }

    // Traverses the octree nodes to find the closest intersection with the ray
    bool static ray_intersect_node(Ray3f &ray, SurfaceIntersection3f &its, bool shadow_ray,
                                   const std::shared_ptr<OctreeNode> &node) {
        if (!node->bbox.ray_intersect(ray)) {
            return false; // No intersection with this node's bounding box
        }

        if (node->is_leaf) {
            bool hit = false;

            for (const auto &prim : node->primitives) {
                float u, v, t;
                if (prim.first->ray_intersect(prim.second, ray, u, v, t)) {
                    if (shadow_ray) {
                        return true; // Early exit for shadow ray
                    }
                    ray.max_t() = its.t = t; // Ensure it is the closet
                    its.uv              = Point2f(u, v);
                    its.mesh            = prim.first;
                    its.primitive_index = prim.second;
                    its.wi              = -ray.d();
                    hit                 = true;
                }
            }

            if (hit) {
                /* At this point, we now know that there is an intersection,
                   and we know the triangle index of the closest such intersection.

                   The following computes a number of additional properties which
                   characterize the intersection (normals, texture coordinates, etc..)
                */

                /* Find the barycentric coordinates */
                Vector3f bary(1 - its.uv.x() - its.uv.y(), its.uv.x(), its.uv.y());

                /* References to all relevant mesh buffers */
                std::shared_ptr mesh(its.mesh);
                const std::vector<Point3f> &v    = mesh->get_vertex_positions();
                const std::vector<Normal3f> &n   = mesh->get_vertex_normals();
                const std::vector<Point2f> &uv   = mesh->get_vertex_tex_coords();
                const std::vector<Point3i> &face = mesh->get_indices();

                /* Vertex indices of the triangle */
                uint32_t idx0 = face[its.primitive_index](0);
                uint32_t idx1 = face[its.primitive_index](1);
                uint32_t idx2 = face[its.primitive_index](2);

                Point3f p0 = v[idx0], p1 = v[idx1], p2 = v[idx2];

                /* Compute the intersection position accurately
                   using barycentric coordinates */
                its.p = p0 * bary.x() + p1 * bary.y() + p2 * bary.z();

                /* Compute proper texture coordinates if provided by the mesh */
                if (!uv.empty()) {
                    Point2f uv0 = uv[idx0], uv1 = uv[idx1], uv2 = uv[idx2];
                    its.uv = uv0 * bary.x() + uv1 * bary.y() + uv2 * bary.z();
                    // Compute UV deltas
                    Vector2f duv1 = uv1 - uv0;
                    Vector2f duv2 = uv2 - uv0;
                    // Compute edge vectors
                    Vector3f dp1 = p1 - p0;
                    Vector3f dp2 = p2 - p0;
                    float det    = duv1.x() * duv2.y() - duv1.y() * duv2.x();
                    if (abs(det) < M_EPSILON) {
                        its.dp_du = dp1;
                        its.dp_dv = dp2;
                    } else {
                        float inv_det = 1.0f / det;
                        its.dp_du     = (dp1 * duv2.y() - dp2 * duv1.y()) * inv_det;
                        its.dp_dv     = (-dp1 * duv2.x() + dp2 * duv1.x()) * inv_det;
                    }
                }

                /* Compute the geometry frame */
                bool valid          = true;
                its.n               = Normal3f((p1 - p0).cross(p2 - p0).norm(valid));
                its.geometric_frame = Frame3f(its.n);
#ifdef M_DEBUG
                if (!valid) {
                    std::cout << "Warning: geometric_frame normal is invalid" << std::endl;
                }
#endif

                if (!n.empty()) {
                    valid = true;
                    its.n = Normal3f((n[idx0] * bary.x() + n[idx1] * bary.y() + n[idx2] * bary.z()).norm(valid));
                    its.shading_frame = Frame3f(its.n);
#ifdef M_DEBUG
                    if (!valid) {
                        std::cout << "Warning: shading_frame normal is invalid" << std::endl;
                    }
#endif
                } else {
                    its.shading_frame = its.geometric_frame;
                }
                its.wi = its.to_local(its.wi);
            }

            return hit;
        } else {
            // Check all children recursively
            bool hit = false;
            for (const auto &i : node->children) {
                if (i && ray_intersect_node(ray, its, shadow_ray, i)) {
                    hit = true;
                }
            }
            return hit;
        }
    }

    bool static ray_test_node(const Ray3f &ray, const std::shared_ptr<OctreeNode> &node) {
        if (!node->bbox.ray_intersect(ray)) {
            return false;
        }

        if (node->is_leaf) {
            for (const auto &prim : node->primitives) {
                float u, v, t;
                if (prim.first->ray_intersect(prim.second, ray, u, v, t)) {
                    return true;
                }
            }
            return false;
        } else {
            bool hit = false;
            for (const auto &i : node->children) {
                if (i && ray_test_node(ray, i)) {
                    hit = true;
                }
            }
            return hit;
        }
    }
};

REGISTER_CLASS(OctreeAccel, "octree");

M_NAMESPACE_END

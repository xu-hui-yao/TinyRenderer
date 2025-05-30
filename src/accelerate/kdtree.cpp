#include <components/accelerate.h>
#include <components/mesh.h>
#include <core/frame.h>
#include <core/intersection.h>

M_NAMESPACE_BEGIN
class KDTreeNode {
public:
    BoundingBox3f bbox;                                                 // Bounding box for the node
    std::vector<std::pair<std::shared_ptr<Mesh>, uint32_t>> primitives; // Primitives in leaf nodes
    std::shared_ptr<KDTreeNode> left, right;                            // Children nodes
    bool is_leaf;                                                       // Indicates if the node is a leaf
    int split_axis;                                                     // Splitting axis (0=x, 1=y, 2=z)
    float split_position;                                               // Splitting position

    explicit KDTreeNode(const BoundingBox3f &bbox)
        : bbox(bbox), left(nullptr), right(nullptr), is_leaf(false), split_axis(-1), split_position(0.0f) {}
};

class KDTreeAccel : public Accel {
public:
    explicit KDTreeAccel(const PropertyList &properties) {
        max_depth               = properties.get_integer("max_depth", 16);
        max_primitives_per_leaf = properties.get_integer("leaf_max", 10);
        m_name                  = properties.get_string("name", "kdtree");
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
        root = build_tree(primitives, bounding_box, 0);
    }

    [[nodiscard]] const BoundingBox3f &get_bounding_box() const override { return bounding_box; }

    bool ray_intersect(const Ray3f &ray, SurfaceIntersection3f &its, bool shadow_ray) const override {
        Ray3f ray_(ray);
        return ray_intersect_node(ray_, its, shadow_ray, root);
    }

    [[nodiscard]] bool ray_test(const Ray3f &ray) const override { return ray_test_node(ray, root); }

    [[nodiscard]] std::string to_string() const override { return "KDTree Accelerate\n"; }

private:
    int max_depth;
    int max_primitives_per_leaf;
    std::shared_ptr<KDTreeNode> root;
    BoundingBox3f bounding_box;
    std::vector<std::pair<std::shared_ptr<Mesh>, uint32_t>> primitives;

    // Recursive KD-Tree construction
    std::shared_ptr<KDTreeNode> build_tree(const std::vector<std::pair<std::shared_ptr<Mesh>, uint32_t>> &primitives,
                                           const BoundingBox3f &bbox, int depth) {
        if (primitives.size() <= max_primitives_per_leaf || depth >= max_depth) {
            auto node        = std::make_shared<KDTreeNode>(bbox);
            node->is_leaf    = true;
            node->primitives = primitives;
            return node;
        }

        int axis             = bbox.get_major_axis();
        float split_position = bbox.get_center()(axis);

        auto left_primitives  = std::vector<std::pair<std::shared_ptr<Mesh>, uint32_t>>();
        auto right_primitives = std::vector<std::pair<std::shared_ptr<Mesh>, uint32_t>>();

        for (const auto &prim : primitives) {
            const auto &tri_bbox = prim.first->get_bounding_box(prim.second);
            if (tri_bbox.get_max()(axis) <= split_position) {
                left_primitives.push_back(prim);
            } else if (tri_bbox.get_min()(axis) >= split_position) {
                right_primitives.push_back(prim);
            } else {
                left_primitives.push_back(prim);
                right_primitives.push_back(prim);
            }
        }

        auto node            = std::make_shared<KDTreeNode>(bbox);
        node->split_axis     = axis;
        node->split_position = split_position;

        if (!left_primitives.empty()) {
            BoundingBox3f left_bbox   = bbox;
            left_bbox.get_max()(axis) = split_position;
            node->left                = build_tree(left_primitives, left_bbox, depth + 1);
        }

        if (!right_primitives.empty()) {
            BoundingBox3f right_bbox   = bbox;
            right_bbox.get_min()(axis) = split_position;
            node->right                = build_tree(right_primitives, right_bbox, depth + 1);
        }

        return node;
    }

    // Ray intersection with KD-Tree nodes
    bool static ray_intersect_node(Ray3f &ray, SurfaceIntersection3f &its, bool shadow_ray,
                                   const std::shared_ptr<KDTreeNode> &node) {
        if (!node->bbox.ray_intersect(ray)) {
            return false;
        }

        if (node->is_leaf) {
            bool hit = false;

            for (const auto &prim : node->primitives) {
                float u, v, t;
                if (prim.first->ray_intersect(prim.second, ray, u, v, t)) {
                    if (shadow_ray) {
                        return true; // Early exit for shadow ray
                    }
                    ray.max_t() = its.t = t;
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
            bool hit_left = false, hit_right = false;
            if (node->left) {
                hit_left = ray_intersect_node(ray, its, shadow_ray, node->left);
            }
            if (node->right) {
                hit_right = ray_intersect_node(ray, its, shadow_ray, node->right);
            }
            return hit_left || hit_right;
        }
    }

    bool static ray_test_node(const Ray3f &ray, const std::shared_ptr<KDTreeNode> &node) {
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
            bool hit_left = false, hit_right = false;
            if (node->left) {
                hit_left = ray_test_node(ray, node->left);
            }
            if (node->right) {
                hit_right = ray_test_node(ray, node->right);
            }
            return hit_left || hit_right;
        }
    }
};

REGISTER_CLASS(KDTreeAccel, "kdtree");

M_NAMESPACE_END

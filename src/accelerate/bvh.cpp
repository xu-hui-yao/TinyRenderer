#include <algorithm>
#include <components/accelerate.h>
#include <components/mesh.h>
#include <core/frame.h>
#include <core/intersection.h>
#include <memory>

M_NAMESPACE_BEGIN

namespace {
// Number of buckets used for binned SAH cost evaluation
constexpr int SAH_BUCKET_COUNT = 12;

// Per-triangle info gathered once before building; only alive during construct()
struct PrimitiveInfo {
    uint32_t mesh_id;
    uint32_t tri_idx;
    BoundingBox3f bbox;
    Point3f centroid;
};

// Reference stored in the reordered primitive array after build (8 bytes, POD)
struct PrimRef {
    uint32_t mesh_id;
    uint32_t tri_idx;
};

// Flat (array-based) BVH node - no pointers, no shared_ptr, directly relocatable
struct LinearBVHNode {
    BoundingBox3f bbox;
    uint32_t offset       = 0; // leaf: index into ordered_prims / interior: index of the right child
    uint16_t n_primitives = 0; // 0 => interior node
    uint8_t axis          = 0; // interior node: split axis used to order traversal
};

// Temporary pointer-based tree used only while running the SAH builder.
// This never exists at ray-tracing time; it is flattened away right after construct().
struct BuildNode {
    BoundingBox3f bbox;
    std::unique_ptr<BuildNode> left, right;
    uint32_t start = 0, count = 0; // valid for leaves: range in ordered_prims
    int axis       = 0;
};
} // namespace

class BVHAccel : public Accel {
public:
    explicit BVHAccel(const PropertyList &properties) {
        max_depth               = properties.get_integer("max_depth", 64);
        max_primitives_per_leaf = properties.get_integer("leaf_max", 4);
        m_name                  = properties.get_string("name", "bvh");
    }

    void add_mesh(const std::shared_ptr<Mesh> &mesh) override {
        m_meshes.emplace_back(mesh);
        bounding_box.expand_by(mesh->get_bounding_box());
    }

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
        // Gather per-triangle info once; this working array is only used during the build
        std::vector<PrimitiveInfo> prim_info;
        uint32_t total_triangles = 0;
        for (const auto &mesh : m_meshes) {
            total_triangles += mesh->get_triangle_count();
        }
        prim_info.reserve(total_triangles);

        for (uint32_t mesh_id = 0; mesh_id < m_meshes.size(); ++mesh_id) {
            const auto &mesh      = m_meshes[mesh_id];
            uint32_t triangle_count = mesh->get_triangle_count();
            for (uint32_t tri_idx = 0; tri_idx < triangle_count; ++tri_idx) {
                PrimitiveInfo info;
                info.mesh_id  = mesh_id;
                info.tri_idx  = tri_idx;
                info.bbox     = mesh->get_bounding_box(tri_idx);
                info.centroid = mesh->get_centroid(tri_idx);
                prim_info.push_back(info);
            }
        }

        nodes.clear();
        ordered_prims.clear();

        if (prim_info.empty()) {
            return;
        }

        ordered_prims.reserve(prim_info.size());
        std::unique_ptr<BuildNode> root =
            build_tree(prim_info, 0, static_cast<uint32_t>(prim_info.size()), ordered_prims, 0);

        uint32_t node_count = count_nodes(root.get());
        nodes.assign(node_count, LinearBVHNode{});
        uint32_t next_free = 0;
        flatten_tree(root.get(), next_free);
    }

    [[nodiscard]] const BoundingBox3f &get_bounding_box() const override { return bounding_box; }

    bool ray_intersect(const Ray3f &ray, SurfaceIntersection3f &its, bool shadow_ray) const override {
        if (nodes.empty()) {
            return false;
        }

        Ray3f ray_(ray);
        uint32_t to_visit[64];
        int stack_ptr    = 0;
        uint32_t current = 0;

        bool hit             = false;
        uint32_t best_mesh   = 0;
        uint32_t best_tri    = 0;
        float best_u = 0, best_v = 0;

        while (true) {
            const LinearBVHNode &node = nodes[current];
            float near_t, far_t;

            if (node.bbox.ray_intersect(ray_, near_t, far_t) && far_t >= ray_.min_t() && near_t <= ray_.max_t()) {
                if (node.n_primitives > 0) {
                    // Leaf node: test contained triangles
                    for (uint16_t i = 0; i < node.n_primitives; ++i) {
                        const PrimRef &pr = ordered_prims[node.offset + i];
                        const Mesh *mesh  = m_meshes[pr.mesh_id].get();

                        float u, v, t;
                        if (mesh->ray_intersect(pr.tri_idx, ray_, u, v, t)) {
                            if (shadow_ray) {
                                return true; // Early exit for shadow ray
                            }
                            ray_.max_t() = t;
                            best_mesh    = pr.mesh_id;
                            best_tri     = pr.tri_idx;
                            best_u       = u;
                            best_v       = v;
                            hit          = true;
                        }
                    }

                    if (stack_ptr == 0) {
                        break;
                    }
                    current = to_visit[--stack_ptr];
                } else {
                    // Interior node: visit the child nearer to the ray origin first
                    if (ray_.d()(node.axis) < 0.0f) {
                        to_visit[stack_ptr++] = current + 1;
                        current               = node.offset;
                    } else {
                        to_visit[stack_ptr++] = node.offset;
                        current               = current + 1;
                    }
                }
            } else {
                if (stack_ptr == 0) {
                    break;
                }
                current = to_visit[--stack_ptr];
            }
        }

        if (hit) {
            /* At this point, we now know that there is an intersection,
               and we know the triangle index of the closest such intersection.

               The following computes a number of additional properties which
               characterize the intersection (normals, texture coordinates, etc..)
            */
            its.t   = ray_.max_t();
            its.uv  = Point2f(best_u, best_v);

            /* Store the mesh as a plain index - zero atomic refcount operations,
               fully GPU-uploadable (matches TSurfaceIntersection::mesh_id) */
            const Mesh *mesh    = m_meshes[best_mesh].get();
            its.mesh_id         = best_mesh;
            its.primitive_index = best_tri;
            its.wi              = -ray_.d();

            /* Find the barycentric coordinates */
            Vector3f bary(1 - best_u - best_v, best_u, best_v);

            /* References to all relevant mesh buffers */
            const std::vector<Point3f> &v    = mesh->get_vertex_positions();
            const std::vector<Normal3f> &n   = mesh->get_vertex_normals();
            const std::vector<Point2f> &uv   = mesh->get_vertex_tex_coords();
            const std::vector<Point3i> &face = mesh->get_indices();

            /* Vertex indices of the triangle */
            uint32_t idx0 = face[best_tri](0);
            uint32_t idx1 = face[best_tri](1);
            uint32_t idx2 = face[best_tri](2);

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
    }

    [[nodiscard]] bool ray_test(const Ray3f &ray) const override {
        if (nodes.empty()) {
            return false;
        }

        uint32_t to_visit[64];
        int stack_ptr    = 0;
        uint32_t current = 0;

        while (true) {
            const LinearBVHNode &node = nodes[current];
            float near_t, far_t;

            if (node.bbox.ray_intersect(ray, near_t, far_t) && far_t >= ray.min_t() && near_t <= ray.max_t()) {
                if (node.n_primitives > 0) {
                    for (uint16_t i = 0; i < node.n_primitives; ++i) {
                        const PrimRef &pr = ordered_prims[node.offset + i];
                        const Mesh *mesh  = m_meshes[pr.mesh_id].get();

                        float u, v, t;
                        if (mesh->ray_intersect(pr.tri_idx, ray, u, v, t)) {
                            return true; // Any-hit: no need to keep searching
                        }
                    }

                    if (stack_ptr == 0) {
                        break;
                    }
                    current = to_visit[--stack_ptr];
                } else {
                    if (ray.d()(node.axis) < 0.0f) {
                        to_visit[stack_ptr++] = current + 1;
                        current               = node.offset;
                    } else {
                        to_visit[stack_ptr++] = node.offset;
                        current               = current + 1;
                    }
                }
            } else {
                if (stack_ptr == 0) {
                    break;
                }
                current = to_visit[--stack_ptr];
            }
        }

        return false;
    }

    [[nodiscard]] std::string to_string() const override {
        return std::string("BVH Accelerate[\n") + "  node count = " + std::to_string(nodes.size()) +
               ",\n  primitive count = " + std::to_string(ordered_prims.size()) +
               ",\n  leaf max = " + std::to_string(max_primitives_per_leaf) + "\n]";
    }

    [[nodiscard]] bool export_gpu_bvh(GPUScene &scene) const override {
        if (nodes.empty()) {
            return false;
        }

        scene.bvh_nodes.resize(nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i) {
            const LinearBVHNode &n = nodes[i];
            GPUBVHNode g;
            g.bmin   = { n.bbox.get_min().x(), n.bbox.get_min().y(), n.bbox.get_min().z() };
            g.bmax   = { n.bbox.get_max().x(), n.bbox.get_max().y(), n.bbox.get_max().z() };
            g.offset = n.offset;
            g.meta   = gpu_bvh_pack_meta(n.n_primitives, n.axis);
            scene.bvh_nodes[i] = g;
        }

        // Convert each (mesh_id, tri_idx) primitive reference into a global
        // triangle id (an index into scene.indices/3), matching the numbering
        // Scene::build_gpu_scene() already assigned to scene.meshes. Leaf
        // node offsets stay valid unchanged: this is a 1:1, order-preserving
        // re-typing of the same array, not a reordering.
        scene.bvh_primitives.resize(ordered_prims.size());
        for (size_t i = 0; i < ordered_prims.size(); ++i) {
            const PrimRef &pr             = ordered_prims[i];
            const GPUMeshInfo &mesh_info  = scene.meshes[pr.mesh_id];
            scene.bvh_primitives[i]       = mesh_info.base_index / 3 + pr.tri_idx;
        }

        return true;
    }

private:
    int max_depth;
    int max_primitives_per_leaf;
    BoundingBox3f bounding_box;

    std::vector<std::shared_ptr<Mesh>> m_meshes; // one entry per mesh (small), not per triangle
    std::vector<LinearBVHNode> nodes;            // flat node array, uploadable as-is
    std::vector<PrimRef> ordered_prims;          // triangle refs reordered to match leaves

    // ------------------------------------------------------------------
    // SAH builder: operates on a temporary pointer tree (build time only)
    // ------------------------------------------------------------------
    std::unique_ptr<BuildNode> build_tree(std::vector<PrimitiveInfo> &prims, uint32_t start, uint32_t end,
                                          std::vector<PrimRef> &out_prims, int depth) const {
        auto node = std::make_unique<BuildNode>();

        BoundingBox3f node_bbox;
        for (uint32_t i = start; i < end; ++i) {
            node_bbox.expand_by(prims[i].bbox);
        }
        node->bbox         = node_bbox;
        uint32_t prim_count = end - start;

        auto make_leaf = [&]() {
            node->start = static_cast<uint32_t>(out_prims.size());
            node->count = prim_count;
            for (uint32_t i = start; i < end; ++i) {
                out_prims.push_back({ prims[i].mesh_id, prims[i].tri_idx });
            }
        };

        if (prim_count <= static_cast<uint32_t>(max_primitives_per_leaf) || depth >= max_depth) {
            make_leaf();
            return node;
        }

        // Choose the split axis as the major axis of the centroid bounds
        BoundingBox3f centroid_bounds;
        for (uint32_t i = start; i < end; ++i) {
            centroid_bounds.expand_by(prims[i].centroid);
        }
        int axis     = centroid_bounds.get_major_axis();
        float c_min  = centroid_bounds.get_min()(axis);
        float c_max  = centroid_bounds.get_max()(axis);

        if (c_max - c_min < M_EPSILON) {
            // All centroids coincide on every axis - nothing meaningful to split
            make_leaf();
            return node;
        }

        // ---------------- Binned SAH split search ----------------
        struct Bucket {
            BoundingBox3f bbox;
            uint32_t count = 0;
        };
        Bucket buckets[SAH_BUCKET_COUNT];

        float inv_extent = static_cast<float>(SAH_BUCKET_COUNT) / (c_max - c_min);
        auto bucket_of    = [&](float c) {
            int b = static_cast<int>((c - c_min) * inv_extent);
            return M_MIN(M_MAX(b, 0), SAH_BUCKET_COUNT - 1);
        };

        for (uint32_t i = start; i < end; ++i) {
            int b = bucket_of(prims[i].centroid(axis));
            buckets[b].count += 1;
            buckets[b].bbox.expand_by(prims[i].bbox);
        }

        BoundingBox3f prefix_bbox[SAH_BUCKET_COUNT], suffix_bbox[SAH_BUCKET_COUNT];
        uint32_t prefix_count[SAH_BUCKET_COUNT], suffix_count[SAH_BUCKET_COUNT];

        BoundingBox3f running_bbox;
        uint32_t running_count = 0;
        for (int i = 0; i < SAH_BUCKET_COUNT; ++i) {
            running_bbox.expand_by(buckets[i].bbox);
            running_count += buckets[i].count;
            prefix_bbox[i]  = running_bbox;
            prefix_count[i] = running_count;
        }

        running_bbox  = BoundingBox3f();
        running_count = 0;
        for (int i = SAH_BUCKET_COUNT - 1; i >= 0; --i) {
            running_bbox.expand_by(buckets[i].bbox);
            running_count += buckets[i].count;
            suffix_bbox[i]  = running_bbox;
            suffix_count[i] = running_count;
        }

        float node_area = node_bbox.get_surface_area();
        float best_cost  = M_MAX_FLOAT;
        int best_split   = -1;
        for (int i = 0; i < SAH_BUCKET_COUNT - 1; ++i) {
            uint32_t nl = prefix_count[i], nr = suffix_count[i + 1];
            if (nl == 0 || nr == 0) {
                continue;
            }
            float area_l = prefix_bbox[i].get_surface_area();
            float area_r = suffix_bbox[i + 1].get_surface_area();
            float cost   = static_cast<float>(nl) * area_l + static_cast<float>(nr) * area_r;
            if (cost < best_cost) {
                best_cost  = cost;
                best_split = i;
            }
        }

        bool force_median_split = false;
        if (best_split < 0) {
            // Degenerate bucket distribution: fall back to an equal-count split if the
            // leaf would otherwise be pathologically large; a small leaf is fine as-is.
            if (prim_count > static_cast<uint32_t>(max_primitives_per_leaf) * 8) {
                force_median_split = true;
            } else {
                make_leaf();
                return node;
            }
        } else {
            float safe_area   = node_area > 0.0f ? node_area : M_EPSILON;
            float split_cost  = 0.5f + best_cost / safe_area; // relative to unit traverse+intersect cost
            bool should_split = split_cost < static_cast<float>(prim_count) ||
                                prim_count > static_cast<uint32_t>(max_primitives_per_leaf) * 8;
            if (!should_split) {
                make_leaf();
                return node;
            }
        }

        uint32_t mid;
        if (force_median_split) {
            mid = start + prim_count / 2;
            std::nth_element(prims.begin() + start, prims.begin() + mid, prims.begin() + end,
                             [axis](const PrimitiveInfo &a, const PrimitiveInfo &b) {
                                 return a.centroid(axis) < b.centroid(axis);
                             });
        } else {
            auto mid_it = std::partition(prims.begin() + start, prims.begin() + end,
                                         [&](const PrimitiveInfo &p) { return bucket_of(p.centroid(axis)) <= best_split; });
            mid = static_cast<uint32_t>(mid_it - prims.begin());

            if (mid == start || mid == end) {
                mid = start + prim_count / 2;
                std::nth_element(prims.begin() + start, prims.begin() + mid, prims.begin() + end,
                                 [axis](const PrimitiveInfo &a, const PrimitiveInfo &b) {
                                     return a.centroid(axis) < b.centroid(axis);
                                 });
            }
        }

        node->axis  = axis;
        node->left  = build_tree(prims, start, mid, out_prims, depth + 1);
        node->right = build_tree(prims, mid, end, out_prims, depth + 1);

        return node;
    }

    static uint32_t count_nodes(const BuildNode *node) {
        if (!node) {
            return 0;
        }
        return 1 + count_nodes(node->left.get()) + count_nodes(node->right.get());
    }

    // Depth-first flatten: writes into `nodes` by index only, so no reference to any
    // element is ever kept alive across a recursive call (the vector never reallocates
    // since it was already sized to its final length before this runs).
    uint32_t flatten_tree(const BuildNode *node, uint32_t &next_free) {
        uint32_t my_index = next_free++;
        LinearBVHNode linear;
        linear.bbox = node->bbox;

        if (!node->left && !node->right) {
            linear.offset       = node->start;
            linear.n_primitives = static_cast<uint16_t>(node->count);
            linear.axis          = 0;
            nodes[my_index]      = linear;
        } else {
            linear.n_primitives  = 0;
            linear.axis          = static_cast<uint8_t>(node->axis);
            nodes[my_index]      = linear;
            flatten_tree(node->left.get(), next_free);
            uint32_t right_index      = flatten_tree(node->right.get(), next_free);
            nodes[my_index].offset    = right_index;
        }

        return my_index;
    }
};

REGISTER_CLASS(BVHAccel, "bvh");

M_NAMESPACE_END

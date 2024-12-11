#include <iostream>
#include <components/accelerate.h>
#include <core/intersection.h>
#include <core/frame.h>

M_NAMESPACE_BEGIN
	class NaiveAccel : public Accel {
	public:
		explicit NaiveAccel(const PropertyList& properties) {
			m_name = properties.get_string("name", "naive");
		}

		void add_mesh(const std::shared_ptr<Mesh>& mesh) override {
			this->meshes.emplace_back(mesh);
			bounding_box.expand_by(mesh->get_bounding_box());
		}

		void construct() override {
#ifdef M_DEBUG
			std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
		}

		[[nodiscard]] const BoundingBox3f& get_bounding_box() const override {
			return bounding_box;
		}

		bool ray_intersect(const Ray3f& ray_, SurfaceIntersection3f& its, bool shadow_ray) const override {
			bool found_intersection = false; // Was an intersection found so far
			auto f = static_cast<uint32_t>(-1); // Triangle index of the closest intersection

			Ray3f ray(ray_); // Make a copy of the ray (we will need to update its 'max_t' value)

			/* Brute force search through all triangles */
			for (auto& mesh : this->meshes) {
				for (uint32_t idx = 0; idx < mesh->get_triangle_count(); ++idx) {
					float u, v, t;
					if (mesh->ray_intersect(idx, ray, u, v, t)) {
						/* An intersection was found! Can terminate
						   immediately if this is a shadow ray query */
						if (shadow_ray)
							return true;
						ray.max_t() = its.t = t;
						its.uv = Point2f(u, v);
						its.mesh = mesh;
						f = idx;
						found_intersection = true;
					}
				}
			}

			if (found_intersection) {
				/* At this point, we now know that there is an intersection,
				   and we know the triangle index of the closest such intersection.

				   The following computes a number of additional properties which
				   characterize the intersection (normals, texture coordinates, etc..)
				*/

				/* Find the barycentric coordinates */
				Vector3f bary(1 - its.uv.x() - its.uv.y(), its.uv.x(), its.uv.y());

				/* References to all relevant mesh buffers */
				std::shared_ptr mesh(its.mesh);
				const std::vector<Point3f>& v = mesh->get_vertex_positions();
				const std::vector<Normal3f>& n = mesh->get_vertex_normals();
				const std::vector<Point2f>& uv = mesh->get_vertex_tex_coords();
				const std::vector<Point3i>& face = mesh->get_indices();

				/* Vertex indices of the triangle */
				uint32_t idx0 = face[f](0), idx1 = face[f](1), idx2 = face[f](2);

				Point3f p0 = v[idx0], p1 = v[idx1], p2 = v[idx2];

				/* Compute the intersection position accurately
				   using barycentric coordinates */
				its.p = p0 * bary.x() + p1 * bary.y() + p2 * bary.z();

				/* Compute proper texture coordinates if provided by the mesh */
				if (!uv.empty())
					its.uv = uv[idx0] * bary.x() + uv[idx1] * bary.y() + uv[idx2] * bary.z();

				/* Compute the geometry frame */
				bool valid = true;
				its.n = Normal3f((p1 - p0).cross(p2 - p0).norm(valid));
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

			return found_intersection;
		}

		[[nodiscard]] bool ray_test(const Ray3f& ray) const override {
			for (auto& mesh : this->meshes) {
				for (uint32_t idx = 0; idx < mesh->get_triangle_count(); ++idx) {
					float u, v, t;
					if (mesh->ray_intersect(idx, ray, u, v, t)) {
						return true;
					}
				}
			}

			return false;
		}

		[[nodiscard]] std::string to_string() const override {
			return "Naive Accelerate\n";
		}

	private:
		std::vector<std::shared_ptr<Mesh>> meshes;
	};

	REGISTER_CLASS(NaiveAccel, "naive");

M_NAMESPACE_END

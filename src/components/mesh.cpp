#include <components/mesh.h>
#include <core/bounding_box.h>
#include <components/bsdf.h>
#include <components/emitter.h>
#include <core/warp.h>
#include <core/intersection.h>

M_NAMESPACE_BEGIN
	Mesh::Mesh() = default;

	void Mesh::construct() {
#ifdef M_DEBUG
		std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif

		if (emitter) {
			emitter->construct();
		}

		if (!bsdf) {
			throw std::runtime_error("Mesh::construct: bsdf is null");
		}
		bsdf->construct();

		// Construct pmf of triangle
		if (get_triangle_count() <= 0) {
			throw std::runtime_error("Mesh::construct: triangle count <= 0");
		}

		// Compute the PMF for the mesh triangles based on their surface area
		std::vector<float> face_areas(get_triangle_count());

		for (uint32_t i = 0; i < get_triangle_count(); ++i) {
			// Get vertex indices for the triangle
			const auto& face = faces[i];
			const Point3f& p0 = vertices[face.x()];
			const Point3f& p1 = vertices[face.y()];
			const Point3f& p2 = vertices[face.z()];

			// Compute the area of the triangle using cross product
			Vector3f edge1 = p1 - p0;
			Vector3f edge2 = p2 - p0;
			face_areas[i] = 0.5f * edge1.cross(edge2).magnitude();

			// emitter weight add
			if (this->emitter && this->emitter->is_spatial_varying() && !uvs.empty()) {
				const auto& uv0 = uvs[face.x()];
				const auto& uv1 = uvs[face.y()];
				const auto& uv2 = uvs[face.z()];
				SurfaceIntersection3f si0, si1, si2;
				si0.uv = uv0;
				si1.uv = uv1;
				si2.uv = uv2;
				bool active = true;
				const auto& color0 = this->emitter->eval(si0, active);
				const auto& color1 = this->emitter->eval(si1, active);
				const auto& color2 = this->emitter->eval(si2, active);
				float color = (color0.luminance() + color1.luminance() + color2.luminance()) / 3.0f;
				face_areas[i] *= color;
			}
		}

		// Initialize the discrete distribution with the computed areas
		m_area_pmf = std::make_shared<DiscreteDistribution1f>(DiscreteDistribution1f(face_areas));
	}

	uint32_t Mesh::get_triangle_count() const { return static_cast<uint32_t>(faces.size()); }

	uint32_t Mesh::get_vertex_count() const { return static_cast<uint32_t>(vertices.size()); }

	float Mesh::surface_area(uint32_t index) const {
		uint32_t i0 = faces[index].x(), i1 = faces[index].y(), i2 = faces[index].z();
		const Point3f p0 = Point3f(vertices[i0]), p1 = Point3f(vertices[i1]), p2 = Point3f(vertices[i2]);
		return 0.5f * (p1 - p0).cross(p2 - p0).magnitude();
	}

	const BoundingBox3f& Mesh::get_bounding_box() const { return bounding_box; }

	bool Mesh::ray_intersect(uint32_t index, const Ray3f& ray, float& u, float& v, float& t) const {
		uint32_t i0 = faces[index].x(), i1 = faces[index].y(), i2 = faces[index].z();
		const Point3f p0 = Point3f(vertices[i0]), p1 = Point3f(vertices[i1]), p2 = Point3f(vertices[i2]);
		auto edge1 = p1 - p0, edge2 = p2 - p0;

		/* Begin calculating determinant - also used to calculate U parameter */
		auto p_vec = ray.d().cross(edge2);

		/* If determinant is near zero, ray lies in plane of triangle */
		float det = edge1.dot(p_vec);

		if (det > -1e-6 && det < 1e-6)
			return false;
		float inv_det = 1.0f / det;

		/* Calculate distance from v[0] to ray origin */
		Vector3f t_vec = Vector3f(ray.o() - p0);

		/* Calculate U parameter and test bounds */
		u = t_vec.dot(p_vec) * inv_det;
		if (u < 0.0f || u > 1.0f)
			return false;

		/* Prepare to test V parameter */
		auto q_vec = t_vec.cross(edge1);

		/* Calculate V parameter and test bounds */
		v = ray.d().dot(q_vec) * inv_det;
		if (v < 0.0f || u + v > 1.0f)
			return false;

		/* Ray intersects triangle -> compute t */
		t = edge2.dot(q_vec) * inv_det;

		return t >= ray.min_t() && t <= ray.max_t();
	}

	BoundingBox3f Mesh::get_bounding_box(uint32_t index) const {
		auto temp = Point3f(vertices[faces[index].x()]);
		BoundingBox3f result(temp);
		result.expand_by(Point3f(vertices[faces[index].y()]));
		result.expand_by(Point3f(vertices[faces[index].z()]));
		return result;
	}

	Point3f Mesh::get_centroid(uint32_t index) const {
		return (vertices[faces[index].x()] + vertices[faces[index].y()] + vertices[faces[index].z()]) * (1.0f / 3.0f);
	}

	const std::vector<Point3f>& Mesh::get_vertex_positions() const { return vertices; }

	const std::vector<Normal3f>& Mesh::get_vertex_normals() const { return normals; }

	const std::vector<Point2f>& Mesh::get_vertex_tex_coords() const { return uvs; }

	const std::vector<Point3i>& Mesh::get_indices() const { return faces; }

	bool Mesh::is_emitter() const { return emitter != nullptr; }

	std::shared_ptr<Emitter>& Mesh::get_emitter() { return emitter; }

	const std::shared_ptr<Emitter>& Mesh::get_emitter() const { return emitter; }

	const std::shared_ptr<BSDF>& Mesh::get_bsdf() const { return bsdf; }

	void Mesh::add_child(const std::shared_ptr<Object>& obj) {
		switch (obj->get_class_type()) {
		case EBSDF:
			if (this->bsdf)
				throw std::runtime_error("Mesh: tried to register multiple BSDF instances!");
			bsdf = std::dynamic_pointer_cast<BSDF>(obj);
			break;

		case EEmitter: {
			if (this->emitter)
				throw std::runtime_error("Mesh: tried to register multiple Emitter instances!");
			this->emitter = std::dynamic_pointer_cast<Emitter>(obj);
		}
		break;

		default:
			throw std::runtime_error(
				"Mesh::add_child(<" + class_type_name(obj->get_class_type()) + ">) is not supported!");
		}
	}

	std::string Mesh::to_string() const {
		return std::string("Mesh[\n") +
			std::string("  name = ") + m_name + std::string(",\n") +
			std::string("  vertex count = ") + std::to_string(vertices.size()) + std::string(",\n") +
			std::string("  triangle count = ") + std::to_string(faces.size()) + std::string(",\n") +
			std::string("  BSDF = ") + (bsdf ? indent(bsdf->to_string()) : std::string("null")) + std::string(",\n") +
			std::string("  Emitter = ") + (emitter ? indent(emitter->to_string()) : std::string("null")) +
			std::string("\n") +
			std::string("]");
	}

	PositionSample3f Mesh::sample_position(const Point2f& sample_, bool& active) const {
		int face_idx;
		Point2f sample = sample_;

		std::tie(face_idx, sample.y()) = m_area_pmf->sample_reuse(sample.y());
		Point3i fi = this->faces[face_idx];

		Point3f p0 = vertices[fi.x()];
		Point3f p1 = vertices[fi.y()];
		Point3f p2 = vertices[fi.z()];

		Vector3f e0 = p1 - p0;
		Vector3f e1 = p2 - p0;
		Point2f b = square_to_uniform_triangle(sample);

		PositionSample3f ps;
		ps.p = e0 * b.x() + e1 * b.y() + p0;
		ps.pdf = m_area_pmf->normalization();
		ps.delta = false;
		ps.t = 0;

		if (!uvs.empty()) {
			Point2f uv0 = uvs[fi.x()];
			Point2f uv1 = uvs[fi.y()];
			Point2f uv2 = uvs[fi.z()];

			ps.uv = uv0 * (1.0f - b.x() - b.y()) + uv1 * b.x() + uv2 * b.y();
		} else {
			ps.uv = b;
		}

		if (!normals.empty()) {
			Normal3f n0 = normals[fi.x()];
			Normal3f n1 = normals[fi.y()];
			Normal3f n2 = normals[fi.z()];

			ps.n = n0 * (1.0f - b.x() - b.y()) + n1 * b.x() + n2 * b.y();
		} else {
			ps.n = e0.cross(e1);
		}
		ps.n.self_norm(active);

		return ps;
	}

	float Mesh::pdf_position(const PositionSample3f& ps, bool& active) const {
		return m_area_pmf->normalization();
	}

	DirectionSample3f Mesh::sample_direction(const Intersection3f& it, const Point2f& sample, bool& active) const {
		DirectionSample3f ds(sample_position(sample, active));
		ds.d = ds.p - it.p;

		float dist_squared = ds.d.square_magnitude();
		ds.dist = sqrt(dist_squared);
		ds.d /= ds.dist;

		float dp = abs(ds.d.dot(ds.n));
		ds.pdf *= dp != 0 ? dist_squared / dp : 0.0f;

		return ds;
	}

	float Mesh::pdf_direction(const Intersection3f& it, const DirectionSample3f& ds, bool& active) const {
		float pdf = pdf_position(ds, active);
		float dp = abs(ds.d.dot(ds.n));

		pdf *= dp != 0.0f ? ds.dist * ds.dist / dp : 0.0f;
		return pdf;
	}

M_NAMESPACE_END

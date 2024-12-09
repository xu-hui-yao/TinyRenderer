#pragma once

#include <core/array.h>
#include <core/frame.h>
#include <core/ray.h>

M_NAMESPACE_BEGIN
	/**
      * \brief Intersection data structure
      *
      * This data structure records local information about a ray-triangle intersection.
      * This includes the position, traveled ray distance, uv coordinates, as well
      * as two local coordinate frames (one that corresponds to the true
      * geometry, and one that is used for shading computations).
      */
	template <typename Scalar_, DeviceType Device_>
	class TIntersection {
	public:
		static constexpr DeviceType Device = Device_;

		typedef Scalar_ Scalar;
		typedef TArray<Scalar, 3, ArrayType::Vector, Device> VectorType;
		typedef TArray<Scalar, 3, ArrayType::Point, Device> PointType;
		typedef TArray<Scalar, 2, ArrayType::Point, Device> PointType2;
		typedef TArray<Scalar, 3, ArrayType::Normal, Device> NormalType;
		typedef TFrame<Scalar, Device> FrameType;
		typedef TRay<Scalar, 3, Device> RayType;

		PointType p; // Position of the surface intersection
		Scalar t; // Un-occluded distance along the ray
		NormalType n;

		M_HOST_DEVICE
		TIntersection(const PointType& _p, const Scalar& _t, const NormalType& _n) : p(_p), t(_t), n(_n) {}

		M_HOST_DEVICE TIntersection() : p(PointType(0)), t(M_MAX_FLOAT), n(1) {}

		M_HOST_DEVICE TIntersection& operator=(const TSurfaceIntersection<Scalar, Device>& other) {
			p = other.p;
			t = other.t;
			n = other.n;
			return *this;
		}

		M_HOST_DEVICE[[nodiscard]] bool is_valid() const {
			return t < M_MAX_FLOAT;
		}

		// Spawn a semi-infinite ray towards the given direction
		M_HOST_DEVICE RayType spawn_ray(const Vector3f& d) const {
			bool valid = true;
			return RayType(offset_p(d), d, 0, M_MAX_FLOAT, valid);
		}

		// Spawn a finite ray towards the given position
		M_HOST_DEVICE RayType spawn_ray_to(const Point3f& t) const {
			bool valid = true;
			PointType o = offset_p(t - p);
			VectorType d = t - o;
			Scalar dist = d.magnitude();
			d /= dist;
			return RayType(o, d, 0, dist * (1.f - M_EPSILON), valid);
		}

		// Return a human-readable string summary of this frame
		[[nodiscard]] std::string to_string() const {
			std::ostringstream oss;
			oss << "Intersect[\n"
				"  p = " << p.to_string() << ",\n"
				"  t = " << t << ",\n"
				"  n = " << n.to_string() << "\n"
				"]\n";
			return oss.str();
		}

	private:
		/**
		 * Compute an offset position, used when spawning a ray from this
		 * interaction. When the interaction is on the surface of a shape, the
		 * position is offset along the surface normal to prevent self intersection.
		 */
		M_HOST_DEVICE PointType offset_p(const VectorType& d) const {
			Scalar mag = (1.0f + p.abs().max_value()) * M_EPSILON;
			mag = n.dot(d) >= 0 ? mag : -mag;
			return n * mag + p;
		}
	};

	template <typename Scalar_, DeviceType Device_>
	class TSurfaceIntersection : public TIntersection<Scalar_, Device_> {
	public:
		static constexpr DeviceType Device = Device_;

		typedef Scalar_ Scalar;
		typedef TArray<Scalar, 3, ArrayType::Vector, Device> VectorType;
		typedef TArray<Scalar, 3, ArrayType::Point, Device> PointType;
		typedef TArray<Scalar, 2, ArrayType::Point, Device> PointType2;
		typedef TArray<Scalar, 3, ArrayType::Normal, Device> NormalType;
		typedef TFrame<Scalar, Device> FrameType;
		typedef TRay<Scalar, 3, Device> RayType;
		typedef TIntersection<Scalar, Device> Base;

		PointType2 uv;
		FrameType shading_frame;
		FrameType geometric_frame;
		VectorType wi;
		uint32_t primitive_index = 0;
		std::shared_ptr<Mesh> mesh = nullptr;

		M_HOST_DEVICE TSurfaceIntersection() : TIntersection<Scalar, Device>() {}

		M_HOST_DEVICE explicit TSurfaceIntersection(const TPositionSample<Scalar, Device>& ps)
			: Base(ps.p, 0, ps.n), uv(ps.uv), shading_frame(FrameType(ps.n)), wi(0) {}

		M_HOST_DEVICE VectorType to_world(const VectorType& v) const {
			return shading_frame.to_world(v);
		}

		M_HOST_DEVICE VectorType to_local(const VectorType& v) const {
			return shading_frame.to_local(v);
		}

		// Return a human-readable string summary of this frame
		[[nodiscard]] std::string to_string() const {
			std::ostringstream oss;
			oss << "Surface Intersect[\n"
				"  p = " << this->p.to_string() << ",\n"
				"  t = " << this->t << ",\n"
				"  n = " << this->n.to_string() << "\n"
				"  uv = " << uv.to_string() << "\n"
				"  shading frame = " << indent(shading_frame.to_string()) << "\n"
				"  wi = " << wi.to_string() << "\n"
				"  primitive index = " << primitive_index << "\n"
				"  mesh = " << (mesh ? "exists" : "null") << "\n"
				"]\n";
			return oss.str();
		}
	};

M_NAMESPACE_END

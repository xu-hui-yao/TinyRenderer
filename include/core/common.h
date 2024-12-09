#pragma once

#include <vector>
#include <filesystem/resolver.h>

#define M_NAMESPACE_BEGIN namespace tiny_renderer {
#define M_NAMESPACE_END }

M_NAMESPACE_BEGIN
#define M_EPSILON	   1e-4
#define M_PI           3.14159265358979323846
#define M_INV_PI       0.31830988618379067154
#define M_INV_TWOPI    0.15915494309189533577
#define M_INV_FOUR_PI  0.07957747154594766788
#define M_SQRT_TWO     1.41421356237309504880
#define M_INV_SQRT_TWO 0.70710678118654752440
#define M_MAX_FLOAT    3.402823466e+38

#define M_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define M_MAX(a,b) (((a) > (b)) ? (a) : (b))

#ifdef __CUDACC__
#define M_HOST_DEVICE __host__ __device__
#else
#define M_HOST_DEVICE
#endif

	enum class DeviceType { CPU, GPU };

	enum class ArrayType { Vector, Point, Normal };

	template <typename Scalar, int Rows, int Cols, DeviceType Device>
	class TMatrix;
	template <typename Scalar, DeviceType Device>
	class TTensor;
	template <typename Scalar, int Dimension, ArrayType ArrayType, DeviceType Device>
	class TArray;
	template <typename Scalar, int Dimension, DeviceType Device>
	class TRay;
	template <typename Scalar, int Dimension, DeviceType Device>
	class TBoundingBox;
	template <typename Scalar, DeviceType Device>
	class TTransform;
	template <typename Scalar, DeviceType Device>
	class TFrame;
	template <typename Scalar, int Dimension, DeviceType Device>
	class TRGBSpectrum;
	template <typename Scalar, DeviceType Device>
	class TIntersection;
	template <typename Scalar, DeviceType Device>
	class TSurfaceIntersection;
	template <typename Scalar, DeviceType Device>
	class TPositionSample;
	template <typename Scalar, DeviceType Device>
	class TDirectionSample;
	template <typename Scalar, DeviceType Device>
	class TBSDFSample;
	template <typename Scalar, typename Index>
	class TDiscreteDistribution;

	typedef TArray<float, 1, ArrayType::Vector, DeviceType::CPU> Vector1f;
	typedef TArray<float, 2, ArrayType::Vector, DeviceType::CPU> Vector2f;
	typedef TArray<float, 3, ArrayType::Vector, DeviceType::CPU> Vector3f;
	typedef TArray<float, 4, ArrayType::Vector, DeviceType::CPU> Vector4f;
	typedef TArray<int, 1, ArrayType::Vector, DeviceType::CPU> Vector1i;
	typedef TArray<int, 2, ArrayType::Vector, DeviceType::CPU> Vector2i;
	typedef TArray<int, 3, ArrayType::Vector, DeviceType::CPU> Vector3i;
	typedef TArray<int, 4, ArrayType::Vector, DeviceType::CPU> Vector4i;
	typedef TArray<float, 1, ArrayType::Point, DeviceType::CPU> Point1f;
	typedef TArray<float, 2, ArrayType::Point, DeviceType::CPU> Point2f;
	typedef TArray<float, 3, ArrayType::Point, DeviceType::CPU> Point3f;
	typedef TArray<float, 4, ArrayType::Point, DeviceType::CPU> Point4f;
	typedef TArray<int, 1, ArrayType::Point, DeviceType::CPU> Point1i;
	typedef TArray<int, 2, ArrayType::Point, DeviceType::CPU> Point2i;
	typedef TArray<int, 3, ArrayType::Point, DeviceType::CPU> Point3i;
	typedef TArray<int, 4, ArrayType::Point, DeviceType::CPU> Point4i;
	typedef TArray<float, 2, ArrayType::Normal, DeviceType::CPU> Normal2f;
	typedef TArray<float, 3, ArrayType::Normal, DeviceType::CPU> Normal3f;
	typedef TTransform<float, DeviceType::CPU> Transform4f;
	typedef TBoundingBox<float, 1, DeviceType::CPU> BoundingBox1f;
	typedef TBoundingBox<float, 2, DeviceType::CPU> BoundingBox2f;
	typedef TBoundingBox<float, 3, DeviceType::CPU> BoundingBox3f;
	typedef TBoundingBox<float, 4, DeviceType::CPU> BoundingBox4f;
	typedef TBoundingBox<int, 1, DeviceType::CPU> BoundingBox1i;
	typedef TBoundingBox<int, 2, DeviceType::CPU> BoundingBox2i;
	typedef TBoundingBox<int, 3, DeviceType::CPU> BoundingBox3i;
	typedef TBoundingBox<int, 4, DeviceType::CPU> BoundingBox4i;
	typedef TRay<float, 2, DeviceType::CPU> Ray2f;
	typedef TRay<float, 3, DeviceType::CPU> Ray3f;
	typedef TFrame<float, DeviceType::CPU> Frame3f;
	typedef TMatrix<float, 4, 4, DeviceType::CPU> Matrix4f;
	typedef TTensor<float, DeviceType::CPU> TensorXf;
	typedef TRGBSpectrum<float, 3, DeviceType::CPU> Color3f;
	typedef TRGBSpectrum<float, 4, DeviceType::CPU> Color4f;
	typedef TIntersection<float, DeviceType::CPU> Intersection3f;
	typedef TSurfaceIntersection<float, DeviceType::CPU> SurfaceIntersection3f;
	typedef TPositionSample<float, DeviceType::CPU> PositionSample3f;
	typedef TDirectionSample<float, DeviceType::CPU> DirectionSample3f;
	typedef TBSDFSample<float, DeviceType::CPU> BSDFSample3f;
	typedef TDiscreteDistribution<float, int> DiscreteDistribution1f;

	class BSDF;
	class BlockGenerator;
	class Camera;
	class Film;
	class Integrator;
	class Emitter;
	class Mesh;
	class Object;
	class ObjectFactory;
	class ReconstructionFilter;
	class Sampler;
	class Scene;
	class Accel;

#define M_FLOAT2_OP(op) float2 operator op (const float2& other) const {\
	float2 result{0, 0};\
	result.x = this->x op other.x;\
	result.y = this->y op other.y;\
	return result;\
}
#define M_FLOAT2_SELF_OP(op) float2& operator op (const float2& other) {\
	this->x op other.x;\
	this->y op other.y;\
	return *this;\
}

	struct float2 {
		float x, y;
		M_FLOAT2_OP(+)
		M_FLOAT2_OP(-)
		M_FLOAT2_OP(*)
		M_FLOAT2_OP(/)
		M_FLOAT2_SELF_OP(+=)
		M_FLOAT2_SELF_OP(-=)
		M_FLOAT2_SELF_OP(*=)
		M_FLOAT2_SELF_OP(/=)

		bool operator ==(const float2& other) const {
			return x == other.x && y == other.y;
		}

		bool operator !=(const float2& other) const {
			return x != other.x || y != other.y;
		}
	};

#define M_FLOAT3_OP(op) float3 operator op (const float3& other) const {\
	float3 result{0, 0, 0};\
	result.x = this->x op other.x;\
	result.y = this->y op other.y;\
	result.z = this->z op other.z;\
	return result;\
}
#define M_FLOAT3_SELF_OP(op) float3& operator op (const float3& other) {\
	this->x op other.x;\
	this->y op other.y;\
	this->z op other.z;\
	return *this;\
}

	struct float3 {
		float x, y, z;
		M_FLOAT3_OP(+)
		M_FLOAT3_OP(-)
		M_FLOAT3_OP(*)
		M_FLOAT3_OP(/)
		M_FLOAT3_SELF_OP(+=)
		M_FLOAT3_SELF_OP(-=)
		M_FLOAT3_SELF_OP(*=)
		M_FLOAT3_SELF_OP(/=)

		bool operator ==(const float3& other) const {
			return x == other.x && y == other.y && z == other.z;
		}

		bool operator !=(const float3& other) const {
			return x != other.x || y != other.y || z != other.z;
		}
	};

#define M_INT3_OP(op) int3 operator op (const int3& other) const {\
int3 result{0, 0, 0};\
result.x = this->x op other.x;\
result.y = this->y op other.y;\
result.z = this->z op other.z;\
return result;\
}
#define M_INT3_SELF_OP(op) int3& operator op (const int3& other) {\
this->x op other.x;\
this->y op other.y;\
this->z op other.z;\
return *this;\
}

	struct int3 {
		int x, y, z;
		M_INT3_OP(+)
		M_INT3_OP(-)
		M_INT3_OP(*)
		M_INT3_OP(/)
		M_INT3_SELF_OP(+=)
		M_INT3_SELF_OP(-=)
		M_INT3_SELF_OP(*=)
		M_INT3_SELF_OP(/=)

		bool operator ==(const int3& other) const {
			return x == other.x && y == other.y && z == other.z;
		}

		bool operator !=(const int3& other) const {
			return x != other.x || y != other.y || z != other.z;
		}
	};

	// Convert radians to degrees
	template <typename Scalar>
	M_HOST_DEVICE Scalar rad_to_deg(Scalar value) { return value * (180.0 / M_PI); }

	// Convert degrees to radians
	template <typename Scalar>
	M_HOST_DEVICE Scalar deg_to_rad(Scalar value) { return value * (M_PI / 180.0); }

	// Linearly interpolate between two values
	template <typename Scalar>
	M_HOST_DEVICE Scalar interpolate(Scalar t, Scalar v1, Scalar v2) {
		return (static_cast<Scalar>(1) - t) * v1 + t * v2;
	}

	template <typename Scalar>
	M_HOST_DEVICE Scalar sqrt(Scalar value) {
#ifdef __CUDA_ARCH__
	return sqrtf(value);
#else
		return std::sqrt(value);
#endif
	}

	template <typename Scalar>
	M_HOST_DEVICE Scalar safe_sqrt(Scalar value) {
#ifdef __CUDA_ARCH__
		value = M_MAX(0, value);
		return sqrtf(value);
#else
		value = M_MAX(0, value);
		return std::sqrt(value);
#endif
	}

	template <typename Scalar>
	M_HOST_DEVICE static Scalar clamp(Scalar value, Scalar min, Scalar max) {
		if (value < min)
			return min;
		else if (value > max)
			return max;
		else return value;
	}

	template <typename Scalar>
	M_HOST_DEVICE static Scalar pow(Scalar base, Scalar exponent) {
#ifdef __CUDA_ARCH__
		return powf(base, exponent);
#else
		return std::pow(base, exponent);
#endif
	}

	// Always-positive modulo operation
	template <typename Scalar>
	M_HOST_DEVICE int mod(Scalar a, Scalar b) {
		Scalar r = a % b;
		return r < 0 ? r + b : r;
	}

	template <typename Scalar>
	M_HOST_DEVICE Scalar cos(Scalar x) {
#ifdef __CUDA_ARCH__
		return cosf(x);
#else
		return std::cos(x);
#endif
	}

	template <typename Scalar>
	M_HOST_DEVICE Scalar sin(Scalar x) {
#ifdef __CUDA_ARCH__
		return sinf(x);
#else
		return std::sin(x);
#endif
	}

	template <typename Scalar>
	M_HOST_DEVICE Scalar tan(Scalar x) {
#ifdef __CUDA_ARCH__
		return tanf(x);
#else
		return std::tan(x);
#endif
	}

	template <typename Scalar>
	M_HOST_DEVICE Scalar atan2(Scalar y, Scalar x) {
#ifdef __CUDA_ARCH__
		return atan2f(y, x);  // CUDA device implementation
#else
		return std::atan2(y, x); // Standard C++ implementation
#endif
	}

	template <typename Scalar>
	M_HOST_DEVICE Scalar abs(Scalar x) {
		return x < 0 ? -x : x;
	}

	template <typename T>
	M_HOST_DEVICE T copysign(T value, T sign) {
		return value >= T(0) == sign >= T(0) ? value : -value;
	}

	template <typename T>
	M_HOST_DEVICE T mulsign(T value, T sign) {
		return sign >= 0 ? value : -value;
	}

	extern std::string indent(const std::string& string, int amount = 2);

	// Convert a string to lower case
	extern std::string to_lower(const std::string& value);

	// Convert a string into a boolean value
	extern bool to_bool(const std::string& str);

	// Convert a string into a signed integer value
	extern int to_int(const std::string& str);

	// Convert a string into an unsigned integer value
	extern unsigned int to_uint(const std::string& str);

	// Convert a string into a floating point value
	extern float to_float(const std::string& str);

	// Convert a string into a 3D vector
	extern Vector3f to_vector3f(const std::string& str);

	// Convert a string into a 3D point
	extern Point3f to_point3f(const std::string& str);

	// Convert a string into RGB color
	extern Color3f to_color3f(const std::string& str);

	// Tokenize a string into a list by splitting at 'delim'
	extern std::vector<std::string> tokenize(const std::string& s, const std::string& delim = ", ",
	                                         bool include_empty = false);

	// Check if a string ends with another string
	extern bool end_with(const std::string& value, const std::string& ending);

	// Convert a time value in milliseconds into a human-readable string
	extern std::string time_string(double time, bool precise = false);

	// Convert a memory amount in bytes into a human-readable string
	extern std::string mem_string(size_t size, bool precise = false);

	extern std::shared_ptr<filesystem::resolver> get_file_resolver();

M_NAMESPACE_END

#pragma once

#include <filesystem/resolver.h>
#include <vector>

#define M_NAMESPACE_BEGIN namespace tiny_renderer {
#define M_NAMESPACE_END }

M_NAMESPACE_BEGIN
#define M_EPSILON 1e-4
#define M_PI 3.14159265358979323846
#define M_INV_PI 0.31830988618379067154
#define M_INV_TWOPI 0.15915494309189533577
#define M_INV_FOUR_PI 0.07957747154594766788
#define M_SQRT_TWO 1.41421356237309504880
#define M_INV_SQRT_TWO 0.70710678118654752440
#define M_MAX_FLOAT 3.402823466e+38

#define M_MIN(a, b) (((a) < (b)) ? (a) : (b))
#define M_MAX(a, b) (((a) > (b)) ? (a) : (b))

enum class ArrayType { Vector, Point, Normal };

template <typename Scalar, int Rows, int Cols> class TMatrix;
template <typename Scalar> class TTensor;
template <typename Scalar, int Dimension, ArrayType ArrayType> class TArray;
template <typename Scalar, int Dimension> class TRay;
template <typename Scalar, int Dimension> class TBoundingBox;
template <typename Scalar> class TTransform;
template <typename Scalar> class TFrame;
template <typename Scalar, int Dimension> class TRGBSpectrum;
template <typename Scalar> class TIntersection;
template <typename Scalar> class TSurfaceIntersection;
template <typename Scalar> class TPositionSample;
template <typename Scalar> class TDirectionSample;
template <typename Scalar> class TBSDFSample;
template <typename Scalar, typename Index> class TDiscreteDistribution;
template <typename Scalar, typename Index> class TDiscreteDistribution2D;
template <typename Scalar> class TMicrofacetDistribution;

typedef TArray<float, 1, ArrayType::Vector> Vector1f;
typedef TArray<float, 2, ArrayType::Vector> Vector2f;
typedef TArray<float, 3, ArrayType::Vector> Vector3f;
typedef TArray<float, 4, ArrayType::Vector> Vector4f;
typedef TArray<int, 1, ArrayType::Vector> Vector1i;
typedef TArray<int, 2, ArrayType::Vector> Vector2i;
typedef TArray<int, 3, ArrayType::Vector> Vector3i;
typedef TArray<int, 4, ArrayType::Vector> Vector4i;
typedef TArray<float, 1, ArrayType::Point> Point1f;
typedef TArray<float, 2, ArrayType::Point> Point2f;
typedef TArray<float, 3, ArrayType::Point> Point3f;
typedef TArray<float, 4, ArrayType::Point> Point4f;
typedef TArray<int, 1, ArrayType::Point> Point1i;
typedef TArray<int, 2, ArrayType::Point> Point2i;
typedef TArray<int, 3, ArrayType::Point> Point3i;
typedef TArray<int, 4, ArrayType::Point> Point4i;
typedef TArray<float, 2, ArrayType::Normal> Normal2f;
typedef TArray<float, 3, ArrayType::Normal> Normal3f;
typedef TTransform<float> Transform4f;
typedef TBoundingBox<float, 1> BoundingBox1f;
typedef TBoundingBox<float, 2> BoundingBox2f;
typedef TBoundingBox<float, 3> BoundingBox3f;
typedef TBoundingBox<float, 4> BoundingBox4f;
typedef TBoundingBox<int, 1> BoundingBox1i;
typedef TBoundingBox<int, 2> BoundingBox2i;
typedef TBoundingBox<int, 3> BoundingBox3i;
typedef TBoundingBox<int, 4> BoundingBox4i;
typedef TRay<float, 2> Ray2f;
typedef TRay<float, 3> Ray3f;
typedef TFrame<float> Frame3f;
typedef TMatrix<float, 4, 4> Matrix4f;
typedef TTensor<float> TensorXf;
typedef TRGBSpectrum<float, 3> Color3f;
typedef TRGBSpectrum<float, 4> Color4f;
typedef TIntersection<float> Intersection3f;
typedef TSurfaceIntersection<float> SurfaceIntersection3f;
typedef TPositionSample<float> PositionSample3f;
typedef TDirectionSample<float> DirectionSample3f;
typedef TBSDFSample<float> BSDFSample3f;
typedef TDiscreteDistribution<float, int> DiscreteDistribution1f;
typedef TDiscreteDistribution2D<float, int> DiscreteDistribution2f;
typedef TMicrofacetDistribution<float> MicrofacetDistribution1f;

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

#define M_FLOAT2_OP(op)                                                                                                \
    float2 operator op(const float2 &other) const {                                                                    \
        float2 result{ 0, 0 };                                                                                         \
        result.x = this->x op other.x;                                                                                 \
        result.y = this->y op other.y;                                                                                 \
        return result;                                                                                                 \
    }
#define M_FLOAT2_SELF_OP(op)                                                                                           \
    float2 &operator op(const float2 & other) {                                                                        \
        this->x op other.x;                                                                                            \
        this->y op other.y;                                                                                            \
        return *this;                                                                                                  \
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

    bool operator==(const float2 &other) const { return x == other.x && y == other.y; }

    bool operator!=(const float2 &other) const { return x != other.x || y != other.y; }
};

#define M_FLOAT3_OP(op)                                                                                                \
    float3 operator op(const float3 &other) const {                                                                    \
        float3 result{ 0, 0, 0 };                                                                                      \
        result.x = this->x op other.x;                                                                                 \
        result.y = this->y op other.y;                                                                                 \
        result.z = this->z op other.z;                                                                                 \
        return result;                                                                                                 \
    }
#define M_FLOAT3_SELF_OP(op)                                                                                           \
    float3 &operator op(const float3 & other) {                                                                        \
        this->x op other.x;                                                                                            \
        this->y op other.y;                                                                                            \
        this->z op other.z;                                                                                            \
        return *this;                                                                                                  \
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

    bool operator==(const float3 &other) const { return x == other.x && y == other.y && z == other.z; }

    bool operator!=(const float3 &other) const { return x != other.x || y != other.y || z != other.z; }
};

#define M_INT3_OP(op)                                                                                                  \
    int3 operator op(const int3 &other) const {                                                                        \
        int3 result{ 0, 0, 0 };                                                                                        \
        result.x = this->x op other.x;                                                                                 \
        result.y = this->y op other.y;                                                                                 \
        result.z = this->z op other.z;                                                                                 \
        return result;                                                                                                 \
    }
#define M_INT3_SELF_OP(op)                                                                                             \
    int3 &operator op(const int3 & other) {                                                                            \
        this->x op other.x;                                                                                            \
        this->y op other.y;                                                                                            \
        this->z op other.z;                                                                                            \
        return *this;                                                                                                  \
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

    bool operator==(const int3 &other) const { return x == other.x && y == other.y && z == other.z; }

    bool operator!=(const int3 &other) const { return x != other.x || y != other.y || z != other.z; }
};

// Convert radians to degrees
template <typename Scalar> Scalar rad_to_deg(Scalar value) { return value * (180.0 / M_PI); }

// Convert degrees to radians
template <typename Scalar> Scalar deg_to_rad(Scalar value) { return value * (M_PI / 180.0); }

// Linearly interpolate between two values
template <typename Scalar> Scalar interpolate(Scalar t, Scalar v1, Scalar v2) {
    return (static_cast<Scalar>(1) - t) * v1 + t * v2;
}

template <typename Scalar> Scalar sqrt(Scalar value) {
    return std::sqrt(value);
}

template <typename Scalar> Scalar safe_sqrt(Scalar value) {
    value = M_MAX(0, value);
    return std::sqrt(value);
}

template <typename Scalar> static Scalar clamp(Scalar value, Scalar min, Scalar max) {
    if (value < min)
        return min;
    else if (value > max)
        return max;
    else
        return value;
}

template <typename Scalar> static Scalar pow(Scalar base, Scalar exponent) {
    return std::pow(base, exponent);
}

// Always-positive modulo operation
template <typename Scalar> int mod(Scalar a, Scalar b) {
    Scalar r = a % b;
    return r < 0 ? r + b : r;
}

template <typename Scalar> Scalar cos(Scalar x) {
    return std::cos(x);
}

template <typename Scalar> Scalar safe_acos(Scalar x) {
    x = clamp(x, Scalar(-1), Scalar(1));
    return std::acos(x);
}

template <typename Scalar> Scalar sin(Scalar x) {
    return std::sin(x);
}

template <typename Scalar> Scalar tan(Scalar x) {
    return std::tan(x);
}

template <typename Scalar> Scalar atan2(Scalar y, Scalar x) {
    return std::atan2(y, x);
}

template <typename Scalar> Scalar asin(Scalar x) {
    return std::asin(x);
}

template <typename Scalar> Scalar acos(Scalar x) {
    return std::acos(x);
}

template <typename Scalar> Scalar abs(Scalar x) { return x < 0 ? -x : x; }

template <typename Scalar> Scalar floor(Scalar x) {
    return std::floor(x);
}

template <typename Scalar> Scalar ceil(Scalar x) {
    return std::ceil(x);
}

template <typename Scalar> Scalar lerp(Scalar a, Scalar b, Scalar t) { return b * t + a * (Scalar(1) - t); }

template <typename T> T copysign(T value, T sign) { return value >= T(0) == sign >= T(0) ? value : -value; }

template <typename T, typename S> T mulsign(T value, S sign) { return sign >= 0 ? value : -value; }

template <typename Value, size_t n> Value estrin_impl(const Value &x, const Value (&coeff)[n]) {
    constexpr size_t n_rec = (n - 1) / 2, n_fma = n / 2;

    Value coeff_rec[n_rec + 1];

    for (size_t i = 0; i < n_fma; ++i)
        coeff_rec[i] = x * coeff[2 * i + 1] + coeff[2 * i];

    if constexpr (n_rec == n_fma)
        coeff_rec[n_rec] = coeff[n - 1];

    if constexpr (n_rec == 0)
        return coeff_rec[0];
    else
        return estrin_impl(sqr(x), coeff_rec);
}

template <typename Value, size_t n> Value horner_impl(const Value &x, const Value (&coeff)[n]) {
    Value accum = coeff[n - 1];

    for (size_t i = 1; i < n; ++i)
        accum = x * accum + coeff[n - 1 - i];

    return accum;
}

template <typename Value, typename... Ts> Value estrin(const Value &x, Ts... ts) {
    Value coeffs[]{ Value(ts)... };
    return estrin_impl(x, coeffs);
}

template <typename Value, typename... Ts> Value horner(const Value &x, Ts... ts) {
    Value coeffs[]{ Value(ts)... };
    return horner_impl(x, coeffs);
}

extern std::string indent(const std::string &string, int amount = 2);

// Convert a string to lower case
extern std::string to_lower(const std::string &value);

// Convert a string into a boolean value
extern bool to_bool(const std::string &str);

// Convert a string into a signed integer value
extern int to_int(const std::string &str);

// Convert a string into an unsigned integer value
extern unsigned int to_uint(const std::string &str);

// Convert a string into a floating point value
extern float to_float(const std::string &str);

// Convert a string into a 3D vector
extern Vector3f to_vector3f(const std::string &str);

// Convert a string into a 3D point
extern Point3f to_point3f(const std::string &str);

// Convert a string into RGB color
extern Color3f to_color3f(const std::string &str);

// Tokenize a string into a list by splitting at 'delim'
extern std::vector<std::string> tokenize(const std::string &s, const std::string &delim = ", ",
                                         bool include_empty = false);

// Check if a string ends with another string
extern bool end_with(const std::string &value, const std::string &ending);

// Convert a time value in milliseconds into a human-readable string
extern std::string time_string(double time, bool precise = false);

// Convert a memory amount in bytes into a human-readable string
extern std::string mem_string(size_t size, bool precise = false);

extern std::shared_ptr<filesystem::resolver> get_file_resolver();

M_NAMESPACE_END

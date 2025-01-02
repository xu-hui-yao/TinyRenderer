#include <core/warp.h>

M_NAMESPACE_BEGIN
// Low-distortion concentric square to disk mapping by Peter Shirley
Point2f square_to_uniform_disk_concentric(const Point2f &sample) {
    float x = 2.f * sample.x() - 1.f;
    float y = 2.f * sample.y() - 1.f;

    bool is_zero         = x == 0.f && y == 0.f;
    bool quadrant_1_or_3 = abs(x) < abs(y);

    float r  = quadrant_1_or_3 ? y : x;
    float rp = quadrant_1_or_3 ? x : y;

    float phi = 0.25f * static_cast<float>(M_PI) * rp / r; // [0, \pi/4]
    if (quadrant_1_or_3)
        phi = 0.5f * static_cast<float>(M_PI) - phi; // [\pi/4, \pi/2]
    if (is_zero)
        phi = 0.f;

    float s = sin(phi);
    float c = cos(phi);

    return { r * c, r * s };
}

// Inverse of the mapping square_to_uniform_disk_concentric
Point2f uniform_disk_to_square_concentric(const Point2f &p) {
    bool quadrant_0_or_2 = abs(p.x()) > abs(p.y());
    float r_sign         = quadrant_0_or_2 ? p.x() : p.y();
    float r              = copysign(sqrt(p.x() * p.x() + p.y() * p.y()), r_sign);

    float phi = atan2(copysign(p.y(), r_sign), copysign(p.x(), r_sign));

    float t = 4.0f / static_cast<float>(M_PI) * phi;
    t       = quadrant_0_or_2 ? t : 2.f - t;
    t *= r;

    float a = quadrant_0_or_2 ? r : t;
    float b = quadrant_0_or_2 ? t : r;

    return { (a + 1.0f) * 0.5f, (b + 1.0f) * 0.5f };
}

// Low-distortion warping technique based on concentric disk mapping for cosine-weighted hemisphere
Vector3f square_to_cosine_hemisphere(const Point2f &sample) {
    // Low-distortion warping technique based on concentric disk mapping
    Point2f p = square_to_uniform_disk_concentric(sample);

    // Guard against numerical precisions
    float z = sqrt(1.0f - p.square_magnitude());

    return { p.x(), p.y(), z };
}

// Inverse of the mapping square_to_cosine_hemisphere
Point2f cosine_hemisphere_to_square(const Vector3f &v) {
    return uniform_disk_to_square_concentric(Point2f(v.x(), v.y()));
}

// Density of square_to_cosine_hemisphere() with respect to solid angles
float square_to_cosine_hemisphere_pdf(const Vector3f &v) { return static_cast<float>(M_PI) * v.z(); }

// Uniformly sample a vector on the unit hemisphere with respect to solid angles
Vector3f square_to_uniform_hemisphere(const Point2f &sample) {
    // Low-distortion warping technique based on concentric disk mapping
    Point2f p = square_to_uniform_disk_concentric(sample);
    float z   = 1.0f - p.square_magnitude();
    p *= sqrt(z + 1.0f);
    return { p.x(), p.y(), z };
}

// Inverse of the mapping square_to_uniform_hemisphere
Point2f uniform_hemisphere_to_square(const Vector3f &v) {
    Point2f p(v.x(), v.y());
    return uniform_disk_to_square_concentric(p * 1.0f / sqrt(v.z() + 1.0f));
}

// Density of square_to_uniform_hemisphere() with respect to solid angles
float square_to_uniform_hemisphere_pdf(const Vector3f &v) { return M_PI * 0.5f; }

Point2f square_to_uniform_triangle(const Point2f &sample) {
    float t = safe_sqrt(1.0f - sample.x());
    return { 1.0f - t, t * sample.y() };
}

Point2f uniform_triangle_to_square(const Point2f &p) {
    float t = 1 - p.x();
    return { 1.0f - t * t, p.y() / t };
}

float square_to_uniform_triangle_pdf(const Point2f &p) { return 2.0f; }

Vector3f square_to_uniform_sphere(const Point2f &sample) {
    float z = 1.0f - 2.0f * sample.y();
    float r = safe_sqrt(1.0f - z * z);
    float s = sin(2.0f * static_cast<float>(M_PI) * sample.x());
    float c = cos(2.0f * static_cast<float>(M_PI) * sample.y());
    return Vector3f({ r * c, r * s, z });
}

Point2f uniform_sphere_to_square(const Vector3f &p) {
    float phi = atan2(p.y(), p.x()) * static_cast<float>(M_INV_TWOPI);
    return Point2f({ phi < 0.0f ? phi + 1.0f : phi, (1.0f - p.z()) * 0.5f });
}

float square_to_uniform_sphere_pdf(const Vector3f &v) { return M_INV_FOUR_PI; }

M_NAMESPACE_END

#pragma once

#include <core/array.h>
#include <core/common.h>

M_NAMESPACE_BEGIN
/**
 * Low-distortion concentric square to disk mapping
 * This function maps a point from a square (in the range [0, 1] x [0, 1]) to a
 * uniform disk. The mapping minimizes distortion compared to a simple polar
 * transformation. This is used in various graphics techniques, particularly in
 * sampling.
 *
 * Arguments:
 * \param sample: A point (x, y) in the square, where each coordinate is in the
 * range [0, 1].
 *
 * \return
 * - A point in 2D space (x, y) on a unit disk. The coordinates are also in the
 * range [-1, 1].
 */
Point2f square_to_uniform_disk_concentric(const Point2f &sample);

/**
 * Inverse of the mapping square_to_uniform_disk_concentric
 * This function takes a point in the disk and maps it back to a square,
 * reversing the transformation applied by square_to_uniform_disk_concentric.
 *
 * Arguments:
 * \param p: A point (x, y) on the disk, where both coordinates are in the range
 * [-1, 1].
 *
 * \return
 * - A point (x, y) in the square, where both coordinates are in the range [0,
 * 1].
 */
Point2f uniform_disk_to_square_concentric(const Point2f &p);

/**
 * Low-distortion warping technique based on concentric disk mapping for
 * cosine-weighted hemisphere This function maps a square point to a vector on
 * the hemisphere, following a cosine-weighted distribution. The method uses the
 * concentric disk mapping to avoid sampling bias and is typically used for
 * importance sampling.
 *
 * Arguments:
 * \param sample: A point (x, y) in the square, where each coordinate is in the
 * range [0, 1].
 *
 * \return
 * - A 3D vector (x, y, z) on the hemisphere with cosine-weighted distribution.
 */
Vector3f square_to_cosine_hemisphere(const Point2f &sample);

/**
 * Inverse of the mapping square_to_cosine_hemisphere
 * This function takes a 3D vector on the cosine-weighted hemisphere and maps it
 * back to a square. The mapping is the inverse of the
 * square_to_cosine_hemisphere mapping.
 *
 * Arguments:
 * \param v: A 3D vector (x, y, z) on the cosine-weighted hemisphere.
 *
 * \return
 * A point (x, y) in the square, where each coordinate is in the range [0, 1].
 */
Point2f cosine_hemisphere_to_square(const Vector3f &v);

/**
 * Density of square_to_cosine_hemisphere() with respect to solid angles
 * This function returns the probability density of sampling a point on the
 * cosine-weighted hemisphere. It is used to adjust for the non-uniform
 * distribution when sampling directions.
 *
 * Arguments:
 * \param v: A 3D vector (x, y, z) on the cosine-weighted hemisphere.
 *
 * \return
 * - The probability density of the sampled vector with respect to solid angles.
 */
float square_to_cosine_hemisphere_pdf(const Vector3f &v);

/**
 * Uniformly sample a vector on the unit hemisphere with respect to solid angles
 * This function maps a square point to a uniformly distributed vector on the
 * hemisphere. The mapping uses the concentric disk method and is ideal for
 * Monte Carlo sampling.
 *
 * Arguments:
 * \param sample: A point (x, y) in the square, where each coordinate is in the
 * range [0, 1].
 *
 * \return
 * - A 3D vector (x, y, z) on the unit hemisphere, with a uniform distribution.
 */
Vector3f square_to_uniform_hemisphere(const Point2f &sample);

/**
 * Inverse of the mapping square_to_uniform_hemisphere
 * This function takes a 3D vector on the unit hemisphere and maps it back to a
 * square. The mapping is the inverse of the square_to_uniform_hemisphere
 * mapping.
 *
 * Arguments:
 * \param v: A 3D vector (x, y, z) on the unit hemisphere.
 *
 * \return
 * - A point (x, y) in the square, where each coordinate is in the range [0, 1].
 */
Point2f uniform_hemisphere_to_square(const Vector3f &v);

/**
 * Density of square_to_uniform_hemisphere() with respect to solid angles
 * This function returns the probability density of sampling a point on the
 * uniform hemisphere. It adjusts for the non-uniform distribution of sampled
 * directions.
 *
 * Arguments:
 * \param v: A 3D vector (x, y, z) on the unit hemisphere.
 *
 * \return
 * - The probability density of the sampled vector with respect to solid angles.
 */
float square_to_uniform_hemisphere_pdf(const Vector3f &v);

Point2f square_to_uniform_triangle(const Point2f &sample);

Point2f uniform_triangle_to_square(const Point2f &p);

float square_to_uniform_triangle_pdf(const Point2f &p);

Vector3f square_to_uniform_sphere(const Point2f &sample);

Point2f uniform_sphere_to_square(const Vector3f &p);

float square_to_uniform_sphere_pdf(const Vector3f &v);

M_NAMESPACE_END

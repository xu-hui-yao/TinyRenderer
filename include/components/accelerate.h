#pragma once

#include <components/mesh.h>

M_NAMESPACE_BEGIN
/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel : public Object {
public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    virtual void add_mesh(const std::shared_ptr<Mesh> &mesh) = 0;

    // Build the acceleration data structure (currently a no-op)
    void construct() override = 0;

    // Return an axis-aligned box that bounds the scene
    [[nodiscard]] virtual const BoundingBox3f &get_bounding_box() const = 0;

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadow_ray
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    virtual bool ray_intersect(const Ray3f &ray, SurfaceIntersection3f &its, bool shadow_ray) const = 0;

    /**
     * \brief Intersect a ray with the shapes comprising the scene and return a
     * boolean specifying whether an intersection was found.
     *
     * Testing for the mere presence of intersections is considerably faster
     * than finding an actual intersection, hence this function should be
     * preferred over \ref ray_intersect() when geometric information about the
     * first visible intersection is not needed.
     *
     * \return \c true if an intersection was found
     */
    [[nodiscard]] virtual bool ray_test(const Ray3f &ray) const = 0;

    [[nodiscard]] EClassType get_class_type() const override { return EAccelerate; }

    [[nodiscard]] std::string to_string() const override = 0;

protected:
    BoundingBox3f bounding_box;
};

M_NAMESPACE_END

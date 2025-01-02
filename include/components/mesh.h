#pragma once

#include <components/object.h>
#include <core/array.h>
#include <core/bounding_box.h>
#include <core/distribution.h>
#include <core/record.h>

M_NAMESPACE_BEGIN
/**
 * \brief Triangle mesh
 *
 * This class stores a triangle mesh object and provides numerous functions
 * for querying the individual triangles. Subclasses of \c Mesh implement
 * the specifics of how to create its contents (e.g. by loading from an
 * external file)
 */
class Mesh : public Object {
public:
    // Create an empty mesh
    Mesh();

    // Release all memory
    ~Mesh() override = default;

    // Initialize internal data structures (called once by the XML parser)
    void construct() override;

    // Return the total number of triangles in this shape
    [[nodiscard]] uint32_t get_triangle_count() const;

    // Return the total number of vertices in this shape
    [[nodiscard]] uint32_t get_vertex_count() const;

    // Return the surface area of the given triangle
    [[nodiscard]] float surface_area(uint32_t index) const;

    // Return an axis-aligned bounding box of the entire mesh
    [[nodiscard]] const BoundingBox3f &get_bounding_box() const;

    // Return an axis-aligned bounding box containing the given triangle
    [[nodiscard]] BoundingBox3f get_bounding_box(uint32_t index) const;

    // Return the centroid of the given triangle
    [[nodiscard]] Point3f get_centroid(uint32_t index) const;

    /** \brief Ray-triangle intersection test
     *
     * Uses the algorithm by Moeller and Trumbore discussed at
     * <tt>http://www.acm.org/jgt/papers/MollerTrumbore97/code.html</tt>.
     *
     * Note that the test only applies to a single triangle in the mesh.
     * An acceleration data structure like \ref BVH is needed to search
     * for intersections against many triangles.
     *
     * \param index
     *    Index of the triangle that should be intersected
     * \param ray
     *    The ray segment to be used for the intersection query
     * \param t
     *    Upon success, \a t contains the distance from the ray origin to the
     *    intersection point,
     * \param u
     *   Upon success, \c u will contain the 'U' component of the intersection
     *   in barycentric coordinates
     * \param v
     *   Upon success, \c v will contain the 'V' component of the intersection
     *   in barycentric coordinates
     * \return
     *   \c true if an intersection has been detected
     */
    bool ray_intersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const;

    [[nodiscard]] const std::vector<Point3f> &get_vertex_positions() const;

    [[nodiscard]] const std::vector<Normal3f> &get_vertex_normals() const;

    [[nodiscard]] const std::vector<Point2f> &get_vertex_tex_coords() const;

    [[nodiscard]] const std::vector<Point3i> &get_indices() const;

    [[nodiscard]] bool is_emitter() const;

    std::shared_ptr<Emitter> &get_emitter();

    [[nodiscard]] const std::shared_ptr<Emitter> &get_emitter() const;

    [[nodiscard]] const std::shared_ptr<BSDF> &get_bsdf() const;

    void add_child(const std::shared_ptr<Object> &child) override;

    [[nodiscard]] std::string to_string() const override;

    [[nodiscard]] EClassType get_class_type() const override { return EMesh; }

    // =============================================================
    //! @{ \name Sampling routines
    // =============================================================

    /**
     * \brief Sample a point on the surface of this shape
     *
     * The sampling strategy is ideally uniform over the surface, though
     * implementations are allowed to deviate from a perfectly uniform
     * distribution as long as this is reflected in the returned probability
     * density.
     *
     * \param sample
     *     A uniformly distributed 2D point on the domain <tt>[0,1]^2</tt>
     *
     * \param active
     *
     * \return
     *     A \ref PositionSample instance describing the generated sample
     */
    [[nodiscard]] virtual PositionSample3f sample_position(const Point2f &sample, bool &active) const;

    /**
     * \brief Query the probability density of \ref sample_position() for
     * a particular point on the surface.
     *
     * \param ps
     *     A position record describing the sample in question
     *
     * \param active
     *
     * \return
     *     The probability density per unit area
     */
    [[nodiscard]] virtual float pdf_position(const PositionSample3f &ps, bool &active) const;

    /**
     * \brief Sample a direction towards this shape with respect to solid
     * angles measured at a reference position within the scene
     *
     * An ideal implementation of this interface would achieve a uniform solid
     * angle density within the surface region that is visible from the
     * reference position <tt>it.p</tt> (though such an ideal implementation
     * is usually neither feasible nor advisable due to poor efficiency).
     *
     * The function returns the sampled position and the inverse probability
     * per unit solid angle associated with the sample.
     *
     * When the Shape subclass does not supply a custom implementation of this
     * function, the \ref Shape class reverts to a fallback approach that
     * piggybacks on \ref sample_position(). This will generally lead to a
     * suboptimal sample placement and higher variance in Monte Carlo
     * estimators using the samples.
     *
     * \param it
     *    A reference position somewhere within the scene.
     *
     * \param sample
     *     A uniformly distributed 2D point on the domain <tt>[0,1]^2</tt>
     *
     * \param active
     *
     * \return
     *     A \ref DirectionSample instance describing the generated sample
     */
    [[nodiscard]] virtual DirectionSample3f sample_direction(const Intersection3f &it, const Point2f &sample,
                                                             bool &active) const;

    /**
     * \brief Query the probability density of \ref sample_direction()
     *
     * \param it
     *    A reference position somewhere within the scene.
     *
     * \param ds
     *     A position record describing the sample in question
     *
     * \param active
     *
     * \return
     *     The probability density per unit solid angle
     */
    [[nodiscard]] virtual float pdf_direction(const Intersection3f &it, const DirectionSample3f &ds,
                                              bool &active) const;

protected:
    std::vector<Point3f> vertices;              // Vertex positions
    std::vector<Normal3f> normals;              // Vertex normals
    std::vector<Point2f> uvs;                   // Vertex texture coordinates
    std::vector<Point3i> faces;                 // Faces
    std::shared_ptr<BSDF> bsdf       = nullptr; // BSDF of the surface
    std::shared_ptr<Emitter> emitter = nullptr; // Associated emitter, if any
    BoundingBox3f bounding_box;                 // Bounding box of the mesh
    std::shared_ptr<DiscreteDistribution1f> m_area_pmf = nullptr;
};

M_NAMESPACE_END

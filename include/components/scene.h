#pragma once

#include <components/object.h>
#include <components/accelerate.h>
#include <components/integrator.h>
#include <components/sampler.h>

M_NAMESPACE_BEGIN
	/**
	 * \brief Main scene data structure
	 *
	 * This class holds information on scene objects and is responsible for
	 * coordinating rendering jobs. It also provides useful query routines that
	 * are mostly used by the \ref Integrator implementations.
	 */
	class Scene : public Object {
	public:
		explicit Scene(const PropertyList& properties);

		~Scene() override = default;

		[[nodiscard]] const std::shared_ptr<Camera>& get_camera() const { return camera; }

		[[nodiscard]] const std::vector<std::shared_ptr<Mesh>>& get_meshes() const { return meshes; }

		[[nodiscard]] const std::shared_ptr<Accel>& get_accel() const { return accel; }

		[[nodiscard]] const std::shared_ptr<Integrator>& get_integrator() const { return integrator; }

		[[nodiscard]] const std::shared_ptr<Sampler>& get_sampler() const { return sampler; }

		void construct() override;

		void add_child(const std::shared_ptr<Object>& obj) override;

		[[nodiscard]] std::string to_string() const override;

		[[nodiscard]] EClassType get_class_type() const override { return EScene; }

		/**
		 * \brief Direct illumination sampling routine
		 *
		 * This method implements stochastic connections to emitters, which is
		 * variously known as <em>emitter sampling</em>, <em>direct illumination
		 * sampling</em>, or <em>next event estimation</em>.
		 *
		 * The function expects a 3D reference location \c ref as input, which may
		 * influence the sampling process. Normally, this would be the location of
		 * a surface position being shaded. Ideally, the implementation of this
		 * function should then draw samples proportional to the scene's emission
		 * profile and the inverse square distance between the reference point and
		 * the sampled emitter position. However, approximations are acceptable as
		 * long as these are reflected in the returned Monte Carlo sampling weight.
		 *
		 * \param its
		 *    A 3D reference location within the scene, which may influence the
		 *    sampling process.
		 *
		 * \param sample
		 *    A uniformly distributed 2D random variate
		 *
		 * \param test_visibility
		 *    When set to \c true, a shadow ray will be cast to ensure that the
		 *    sampled emitter position and the reference point are mutually visible.
		 *
		 * \param active
		 *
		 * \return
		 *    A tuple <tt>(ds, spec)</tt> where
		 *    <ul>
		 *      <li>\c ds is a fully populated \ref DirectionSample3f data
		 *          structure, which provides further detail about the sampled
		 *          emitter position (e.g. its surface normal, solid angle density,
		 *          whether Dirac delta distributions were involved, etc.)</li>
		 *      <li>
		 *      <li>\c spec is a Monte Carlo sampling weight specifying the ratio
		 *          of the radiance incident from the emitter and the sample
		 *          probability per unit solid angle.</li>
		 *    </ul>
		 */
		[[nodiscard]] std::pair<DirectionSample3f, Color3f> sample_emitter_direction(
			const SurfaceIntersection3f& its, const Point2f& sample, bool test_visibility, bool& active) const;

		/**
		 * \brief Evaluate the PDF of direct illumination sampling
		 *
		 * This function evaluates the probability density (per unit solid angle)
		 * of the sampling technique implemented by the \ref
		 * sample_emitter_direct() function. The returned probability will always
		 * be zero when the emission profile contains a Dirac delta term (e.g.
		 * point or directional emitters/sensors).
		 *
		 * \param it
		 *    A 3D reference location within the scene, which may influence the
		 *    sampling process.
		 *
		 * \param ds
		 *    A direction sampling record, which specifies the query location.
		 *
		 * \param active
		 *
		 * \return
		 *    The solid angle density of the sample
		 */
		[[nodiscard]] float pdf_emitter_direction(const Intersection3f& it, const DirectionSample3f& ds,
		                                          bool& active) const;

		/**
	     * \brief Sample one emitter in the scene and rescale the input sample
	     * for reuse.
	     *
	     * Currently, the sampling scheme implemented by the \ref Scene class is
	     * very simplistic (uniform).
	     *
	     * \param sample
	     *    A uniformly distributed number in [0, 1).
	     *
	     * \param active
	     *
	     * \return
	     *    The index of the chosen emitter along with the sampling weight (equal
	     *    to the inverse PDF), and the transformed random sample for reuse.
	     */
		std::tuple<uint32_t, float, float> sample_emitter(float sample, bool& active) const;

		/**
	     * \brief Evaluate the discrete probability of the \ref sample_emitter() technique for the given emitter index.
	     */
		float pdf_emitter(uint32_t index, bool& active) const;

	private:
		std::vector<std::shared_ptr<Mesh>> meshes;
		std::vector<std::shared_ptr<Emitter>> emitters;
		std::shared_ptr<Camera> camera;
		std::shared_ptr<Accel> accel;
		std::shared_ptr<Integrator> integrator;
		std::shared_ptr<Sampler> sampler;

		uint32_t m_num_emitters;
		float m_emitter_pmf;
	};

M_NAMESPACE_END

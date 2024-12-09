#pragma once

#include <components/object.h>

M_NAMESPACE_BEGIN
	/**
	 * \brief Abstract integrator (i.e. a rendering technique)
	 *
	 * In Nori, the different rendering techniques are collectively referred to as
	 * integrators, since they perform integration over a high-dimensional
	 * space. Each integrator represents a specific approach for solving
	 * the light transport equation---usually favored in certain scenarios, but
	 * at the same time affected by its own set of intrinsic limitations.
	 */
	class Integrator : public Object {
	public:
		/// Release all memory
		~Integrator() override = default;

		void construct() override = 0;

		void add_child(const std::shared_ptr<Object>& child) override = 0;

		/// Perform an (optional) preprocess step
		virtual void preprocess(const std::shared_ptr<Scene>& scene) {}

		/**
		 * \brief Sample the incident radiance along a ray
		 *
		 * \param scene
		 *    A pointer to the underlying scene
		 * \param sampler
		 *    A pointer to a sample generator
		 * \param ray
		 *    The ray in question
		 * \param valid
		 *	  The valid
		 * \return
		 *    A (usually) unbiased estimate of the radiance in this direction
		 */
		[[nodiscard]] virtual Color3f li(const std::shared_ptr<Scene>& scene, std::shared_ptr<Sampler> sampler,
		                                 const Ray3f& ray, bool &valid) const = 0;

		/**
		 * \brief Return the type of object (i.e. Mesh/BSDF/etc.)
		 * provided by this instance
		 * */
		[[nodiscard]] EClassType get_class_type() const override { return EIntegrator; }
	};

M_NAMESPACE_END

#pragma once

#include <components/object.h>
#include <components/rfilter.h>

M_NAMESPACE_BEGIN
	/**
	 * \brief Generic camera interface
	 *
	 * This class provides an abstract interface to cameras in Nori and
	 * exposes the ability to sample their response function. By default, only
	 * a perspective camera implementation exists, but you may choose to
	 * implement other types (e.g. an environment camera, or a physically-based
	 * camera model that simulates the behavior actual lenses)
	 */
	class Camera : public Object {
	public:
		/**
		 * \brief Importance sample a ray according to the camera's response function
		 *
		 * \param ray
		 *    A ray data structure to be filled with a position
		 *    and direction value
		 *
		 * \param sample_position
		 *    Denotes the desired sample position on the film
		 *    expressed in fractional pixel coordinates
		 *
		 * \param aperture_sample
		 *    A uniformly distributed 2D vector that is used to sample
		 *    a position on the aperture of the sensor if necessary.
		 *
		 * \param valid
		 *
		 * \return
		 *    An importance weight associated with the sampled ray.
		 *    This accounts for the difference in the camera response
		 *    function and the sampling density.
		 */
		virtual Color3f sample_ray(Ray3f& ray, const Point2f& sample_position,
		                           const Point2f& aperture_sample, bool& valid) const = 0;

		void construct() override = 0;

		/// Return the size of the output image in pixels
		[[nodiscard]] const Vector2i& get_output_size() const { return output_size; }

		void add_child(const std::shared_ptr<Object>& child) override = 0;

		[[nodiscard]] std::shared_ptr<ReconstructionFilter> get_reconstruction_filter() const {
			return m_reconstruction_filter;
		}

		/**
		 * \brief Return the type of object (i.e. Mesh/Camera/etc.)
		 * provided by this instance
		 * */
		[[nodiscard]] EClassType get_class_type() const override { return ECamera; }

		/// Return a human-readable summary of this instance
		[[nodiscard]] std::string to_string() const override;

	protected:
		Vector2i output_size;
		Transform4f camera_to_world;
		std::shared_ptr<ReconstructionFilter> m_reconstruction_filter;
	};

M_NAMESPACE_END

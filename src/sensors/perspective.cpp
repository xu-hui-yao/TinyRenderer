#include <components/camera.h>

M_NAMESPACE_BEGIN
	/**
	 * \brief Perspective camera with depth of field
	 *
	 * This class implements a simple perspective Camera model. It uses an
	 * infinitesimally small aperture, creating an infinite depth of field.
	 */
	class PerspectiveCamera : public Camera {
	public:
		explicit PerspectiveCamera(const PropertyList& property_list) {
			/* Width and height in pixels. Default: 720p */
			output_size.x() = property_list.get_integer("width", 1280);
			output_size.y() = property_list.get_integer("height", 720);
			inverse_output_size = Vector2f(1.0f / static_cast<float>(output_size.x()),
			                               1.0f / static_cast<float>(output_size.y()));

			/* Specifies an optional Camera-to-world transformation. Default: none */
			camera_to_world = property_list.get_transform("to_world", Transform4f());

			/* Horizontal field of view in degrees */
			fov = property_list.get_float("fov", 30.0f);

			/* Near and far clipping planes in world-space units */
			near_clip = property_list.get_float("near_clip", 1e-4f);
			far_clip = property_list.get_float("far_clip", 1e4f);

			this->m_name = property_list.get_string("name", "camera");
		}

		~PerspectiveCamera() override = default;

		void construct() override {
#ifdef M_DEBUG
			std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif

			auto aspect = static_cast<float>(output_size.x()) / static_cast<float>(output_size.y());

			/* Project vectors in Camera space onto a plane at z=1:
			 *
			 *  xProj = cot * x / z
			 *  yProj = cot * y / z
			 *  zProj = (far * (z - near)) / (z * (far-near))
			 *  The cotangent factor ensures that the field of view is
			 *  mapped to the interval [-1, 1].
			 */
			float recip = 1.0f / (far_clip - near_clip),
			      cot = 1.0f / tan(deg_to_rad(fov / 2.0f));

			Matrix4f perspective = Matrix4f::zero();
			perspective(0, 0) = cot;
			perspective(1, 1) = cot;
			perspective(2, 2) = far_clip * recip;
			perspective(2, 3) = -near_clip * far_clip * recip;
			perspective(3, 2) = 1;

			/**
			 * Translation and scaling to shift the clip coordinates into the
			 * range from zero to one. Also takes the aspect ratio into account.
			 * inverse() represents this is a camera to world matrix.
			 */
			sample_to_camera = Transform4f(Matrix4f::scale(Vector4f(-0.5f, -0.5f * aspect, 1.0f, 1.0f)) *
				Matrix4f::translate(Vector4f(-1.0f, -1.0f / aspect, 0.0f, 0.0f)) * perspective).inverse();

			m_reconstruction_filter->construct();
		}

		Color3f sample_ray(Ray3f& ray, const Point2f& sample_position, const Point2f& aperture_sample,
		                   bool& valid) const override {
			/* Compute the corresponding position on the
			   near plane (in local Camera space) */
			Point3f near_p = sample_to_camera * Point3f(sample_position.x() * inverse_output_size.x(),
			                                            sample_position.y() * inverse_output_size.y(), 0.0f);

			/* Turn into a normalized ray direction, and
			   adjust the ray interval accordingly */
			Vector3f d = Vector3f(near_p.x(), near_p.y(), near_p.z()).norm(valid);
			float inv_z = 1.0f / d.z();

			ray.o() = camera_to_world * Point3f(0, 0, 0);
			ray.d() = camera_to_world * d;
			ray.min_t() = near_clip * inv_z;
			ray.max_t() = far_clip * inv_z;
			ray.update(valid);

			return Color3f(1.0f);
		}

		void add_child(const std::shared_ptr<Object>& child) override {
			switch (child->get_class_type()) {
			case EReconstructionFilter:
				if (m_reconstruction_filter) {
					throw std::runtime_error("There can only be one Reconstruction Filter per Camera!");
				}
				m_reconstruction_filter = std::dynamic_pointer_cast<ReconstructionFilter>(child);
				break;

			default:
				throw std::runtime_error(
					"Camera::add_child(<" + class_type_name(child->get_class_type()) + ">) is not supported!");
			}
		}

		/// Return a human-readable summary
		[[nodiscard]] std::string to_string() const override {
			return std::string("PerspectiveCamera[\n") +
				std::string("  c2w = ") + indent(camera_to_world.to_string()) + std::string(",\n") +
				std::string("  output size = ") + indent(output_size.to_string()) + std::string(",\n") +
				std::string("  fov = ") + std::to_string(fov) + std::string(",\n") +
				std::string("  clip = [near = ") + std::to_string(near_clip) + std::string(", far = ") +
				std::to_string(far_clip) + std::string("],\n") +
				std::string("  reconstruction filter = ") + (m_reconstruction_filter
					                                             ? indent(m_reconstruction_filter->to_string())
					                                             : std::string("None")) +
				std::string("\n]");
		}

	private:
		Vector2f inverse_output_size;
		Transform4f sample_to_camera;
		float fov;
		float near_clip;
		float far_clip;
	};

	REGISTER_CLASS(PerspectiveCamera, "perspective");

M_NAMESPACE_END

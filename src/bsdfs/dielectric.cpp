#include <components/bsdf.h>
#include <components/texture.h>
#include <core/frame.h>
#include <core/record.h>
#include <core/intersection.h>
#include <core/microfacet.h>
#include <core/fresnel.h>

M_NAMESPACE_BEGIN
	class Dielectric : public BSDF {
	public:
		explicit Dielectric(const PropertyList& properties) {
			m_name = properties.get_string("name", "rough dielectric");
			float ext_ior = properties.get_float("ext_ior", 1.0);
			float int_ior = properties.get_float("int_ior", 1.33);
			has_reflection = properties.get_boolean("has_reflection", true);
			has_transmission = properties.get_boolean("has_transmission", true);
			m_eta = int_ior / ext_ior;
		}

		void construct() override {
#ifdef M_DEBUG
			std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
		}

		void add_child(const std::shared_ptr<Object>& child) override {}

		[[nodiscard]] std::pair<BSDFSample3f, Color3f> sample(const SurfaceIntersection3f& si,
		                                                      float sample1,
		                                                      const Point2f& sample2,
		                                                      bool active) const override {
			float cos_theta_i = Frame3f::cos_theta(si.wi, active);

			auto [r_i, cos_theta_t, eta_it, eta_ti] = fresnel(cos_theta_i, m_eta);
			float t_i = 1.0f - r_i;

			BSDFSample3f bs(Vector3f(0));
			bool selected_r;
			float weight = 0.0f;
			if (has_reflection && has_transmission) {
				selected_r = sample1 <= r_i && active;
				bs.pdf = selected_r ? r_i : t_i;
				weight = 1.0f;
			} else if (has_reflection || has_transmission) {
				selected_r = has_reflection && active;
				bs.pdf = 1.0f;
				weight = has_reflection ? r_i : t_i;
			} else {
				return {bs, Color3f(0.0f)};
			}

			bool selected_t = !selected_r && active;

			bs.wo = selected_r ? reflect(si.wi) : refract(si.wi, cos_theta_t, eta_ti);
			bs.eta = selected_r ? 1.0f : eta_it;

			if (selected_t) {
				float factor = eta_ti;
				weight *= factor * factor;
			}

			return {bs, active ? Color3f(weight) : Color3f(0.0f)};
		}

		[[nodiscard]] Color3f eval(const SurfaceIntersection3f& si,
		                           const Vector3f& wo, bool active) const override {
			return Color3f(0.0f);
		}

		[[nodiscard]] float pdf(const SurfaceIntersection3f& si, const Vector3f& wo, bool active) const override {
			return 0.0f;
		}

		// Return a human-readable summary
		[[nodiscard]] std::string to_string() const override {
			return "Dielectric[\n"
				"  eta = " + std::to_string(m_eta) + "\n"
				"]";
		}

		[[nodiscard]] EClassType get_class_type() const override { return EBSDF; }

	private:
		float m_eta;
		bool has_reflection, has_transmission;
	};

	REGISTER_CLASS(Dielectric, "dielectric");

M_NAMESPACE_END

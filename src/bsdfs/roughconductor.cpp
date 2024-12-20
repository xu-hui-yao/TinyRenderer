#include <components/bsdf.h>
#include <components/texture.h>
#include <core/frame.h>
#include <core/record.h>
#include <core/intersection.h>
#include <core/microfacet.h>
#include <core/fresnel.h>

M_NAMESPACE_BEGIN
	class RoughConductor : public BSDF {
	public:
		explicit RoughConductor(const PropertyList& properties) : m_eta(nullptr), m_k(nullptr), m_alpha(nullptr),
		                                                          m_specular_reflectance(nullptr) {
			m_name = properties.get_string("name", "rough conductor");
		}

		void construct() override {
#ifdef M_DEBUG
			std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif

			if (!m_eta) {
				throw std::runtime_error("Rough conductor: eta does not exist");
			}
			if (!m_k) {
				throw std::runtime_error("Rough conductor: k does not exist");
			}
			if (!m_alpha) {
				throw std::runtime_error("Rough conductor: alpha does not exist");
			}
		}

		void add_child(const std::shared_ptr<Object>& child) override {
			switch (child->get_class_type()) {
			case ETexture:
				if (!this->m_eta) {
					m_eta = std::dynamic_pointer_cast<Texture>(child);
				} else if (!this->m_k) {
					m_k = std::dynamic_pointer_cast<Texture>(child);
				} else if (!this->m_alpha) {
					m_alpha = std::dynamic_pointer_cast<Texture>(child);
				} else if (!this->m_specular_reflectance) {
					m_specular_reflectance = std::dynamic_pointer_cast<Texture>(child);
				} else {
					throw std::runtime_error(
						"Rough conductor only supports eta and k and alpha and specular_reflectance");
				}
				break;

			default:
				throw std::runtime_error(
					"BSDF::add_child(<" + class_type_name(child->get_class_type()) + ">) is not supported!");
			}
		}

		[[nodiscard]] std::pair<BSDFSample3f, Color3f> sample(const SurfaceIntersection3f& si,
		                                                      float /* sample1 */,
		                                                      const Point2f& sample2,
		                                                      bool active) const override {
			float cos_theta_i = Frame3f::cos_theta(si.wi, active);
			BSDFSample3f bs(Vector3f(0));
			active &= cos_theta_i > 0.f;

			if (!active) {
				return {bs, Color3f(0)};
			}

			MicrofacetDistribution1f distribution(m_alpha->eval(si, active)(0));

			// Sample the microfacet normal
			Normal3f m;
			std::tie(m, bs.pdf) = distribution.sample(si.wi, sample2, active);

			// Perfect specular reflection based on the microfacet normal
			bs.wo = reflect(si.wi, m);
			bs.eta = 1.f;
			bs.delta = false;

			// Ensure that this is a valid sample
			active &= bs.pdf != 0.0f && Frame3f::cos_theta(bs.wo, active) > 0.f;

			float weight = distribution.smith_g1(bs.wo, m);

			// Jacobian of the half-direction mapping
			bs.pdf /= 4.0f * bs.wo.dot(m);

			// Evaluate the Fresnel factor
			Color3f f({
				fresnel_conductor(si.wi.dot(m), m_eta->eval(si, active)(0), m_k->eval(si, active)(0)),
				fresnel_conductor(si.wi.dot(m), m_eta->eval(si, active)(1), m_k->eval(si, active)(1)),
				fresnel_conductor(si.wi.dot(m), m_eta->eval(si, active)(2), m_k->eval(si, active)(2))
			});

			auto result = f * weight;

			if (m_specular_reflectance) {
				result *= m_specular_reflectance->eval(si, active);
			}

			return {bs, active ? result : Color3f(0)};
		}

		[[nodiscard]] Color3f eval(const SurfaceIntersection3f& si,
		                           const Vector3f& wo, bool active) const override {
			float cos_theta_i = Frame3f::cos_theta(si.wi, active),
			      cos_theta_o = Frame3f::cos_theta(wo, active);

			active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

			if (!active) {
				return Color3f(0);
			}

			// Calculate the half-direction vector
			Normal3f h((wo + si.wi).norm(active));

			// Construct a microfacet distribution matching the roughness values at the current surface position
			MicrofacetDistribution1f distribution(m_alpha->eval(si, active)(0));

			// Evaluate the microfacet normal distribution
			float d = distribution.eval(h);

			active &= d != 0.0f;

			// Evaluate Smith's shadow-masking function
			float g = distribution.g(si.wi, wo, h);

			// Evaluate the full microfacet model (except Fresnel)
			float temp = d * g / (4.0f * Frame3f::cos_theta(si.wi, active));

			// Evaluate the Fresnel factor
			Color3f f({
				fresnel_conductor(si.wi.dot(h), m_eta->eval(si, active)(0), m_k->eval(si, active)(0)),
				fresnel_conductor(si.wi.dot(h), m_eta->eval(si, active)(1), m_k->eval(si, active)(1)),
				fresnel_conductor(si.wi.dot(h), m_eta->eval(si, active)(2), m_k->eval(si, active)(2))
			});

			auto result = f * temp;

			if (m_specular_reflectance) {
				result *= m_specular_reflectance->eval(si, active);
			}

			return active ? result : Color3f(0);
		}

		[[nodiscard]] float pdf(const SurfaceIntersection3f& si, const Vector3f& wo, bool active) const override {
			float cos_theta_i = Frame3f::cos_theta(si.wi, active),
			      cos_theta_o = Frame3f::cos_theta(wo, active);

			Normal3f m((wo + si.wi).norm(active));

			/**
			 * Filter cases where the micro/macro-surface don't agree on the side.
			 * This logic is evaluated in smith_g1() called as part of the eval()
			 * and sample() methods and needs to be replicated in the probability
			 * density computation as well.
			 */
			active &= cos_theta_i > 0.f && cos_theta_o > 0.f && si.wi.dot(m) > 0.f && wo.dot(m) > 0.f;

			if (!active) {
				return 0.0f;
			}

			MicrofacetDistribution1f distribution(m_alpha->eval(si, active)(0));

			float result = distribution.pdf(si.wi, m) / (4.0f * wo.dot(m));

			return active ? result : 0.0f;
		}

		// Return a human-readable summary
		[[nodiscard]] std::string to_string() const override {
			return "RoughConductor[\n"
				"  alpha = " + indent(m_alpha->to_string(), 2) + "\n"
				"  eta = " + indent(m_eta->to_string(), 2) + "\n"
				"  k = " + indent(m_k->to_string(), 2) + "\n"
				"  specular_reflectance = " + indent(m_specular_reflectance->to_string(), 2) + "\n"
				"]";
		}

		[[nodiscard]] EClassType get_class_type() const override { return EBSDF; }

	private:
		std::shared_ptr<Texture> m_eta;
		std::shared_ptr<Texture> m_k;
		std::shared_ptr<Texture> m_alpha;
		std::shared_ptr<Texture> m_specular_reflectance;
	};

	REGISTER_CLASS(RoughConductor, "roughconductor");

M_NAMESPACE_END

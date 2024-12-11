#include <components/bsdf.h>
#include <components/texture.h>
#include <core/frame.h>
#include <core/record.h>
#include <core/intersection.h>
#include <core/microfacet.h>
#include <core/fresnel.h>

M_NAMESPACE_BEGIN
	class Conductor : public BSDF {
	public:
		explicit Conductor(const PropertyList& properties) : m_eta(nullptr), m_k(nullptr) {
			m_name = properties.get_string("name", "conductor");
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
		}

		void add_child(const std::shared_ptr<Object>& child) override {
			switch (child->get_class_type()) {
			case ETexture:
				if (!this->m_eta) {
					m_eta = std::dynamic_pointer_cast<Texture>(child);
				} else if (!this->m_k) {
					m_k = std::dynamic_pointer_cast<Texture>(child);
				} else {
					throw std::runtime_error("Rough conductor only supports eta and k and alpha");
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
			active &= cos_theta_i > 0.f;

			BSDFSample3f bs(Vector3f(0));

			if (!active) {
				return {bs, Color3f(0.0f)};
			}

			// Perfect specular reflection based on the microfacet normal
			bs.wo = reflect(si.wi);
			bs.eta = 1.f;
			bs.pdf = 1.f;

			// Evaluate the Fresnel factor
			Color3f f({
				fresnel_conductor(cos_theta_i, m_eta->eval(si, active)(0), m_k->eval(si, active)(0)),
				fresnel_conductor(cos_theta_i, m_eta->eval(si, active)(1), m_k->eval(si, active)(1)),
				fresnel_conductor(cos_theta_i, m_eta->eval(si, active)(2), m_k->eval(si, active)(2))
			});

			return {bs, active ? f : Color3f(0)};
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
			return "Conductor[\n"
				"  eta = " + indent(m_eta->to_string(), 2) + "\n"
				"  k = " + indent(m_k->to_string(), 2) + "\n"
				"]";
		}

		[[nodiscard]] EClassType get_class_type() const override { return EBSDF; }

	private:
		std::shared_ptr<Texture> m_eta;
		std::shared_ptr<Texture> m_k;
	};

	REGISTER_CLASS(Conductor, "conductor");

M_NAMESPACE_END

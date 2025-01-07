#include <components/bsdf.h>
#include <core/frame.h>
#include <core/intersection.h>
#include <core/record.h>

M_NAMESPACE_BEGIN
/**
 * \brief Two-sided BSDF adapter
 *
 * Wraps a BSDF to make it behave as a two-sided material.
 */
class TwoSided : public BSDF {
public:
    explicit TwoSided(const PropertyList &properties) : m_front_bsdf(nullptr), m_back_bsdf(nullptr) {
        m_name = properties.get_string("name", "two_sided");
    }

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
        if (!m_front_bsdf) {
            throw std::runtime_error("TwoSided: Front BSDF is not provided.");
        }

        m_front_bsdf->construct();

        if (!m_back_bsdf) {
            m_back_bsdf = m_front_bsdf; // Use the same BSDF for both sides if not explicitly provided
        } else {
            m_back_bsdf->construct();
        }

        m_flags = static_cast<BSDFFlags>(static_cast<uint32_t>(m_front_bsdf->get_flag()) |
                                         static_cast<uint32_t>(m_back_bsdf->get_flag())); // TODO: front and back
    }

    void add_child(const std::shared_ptr<Object> &child) override {
        switch (child->get_class_type()) {
            case EBSDF:
                if (!m_front_bsdf) {
                    m_front_bsdf = std::dynamic_pointer_cast<BSDF>(child);
                } else if (!m_back_bsdf) {
                    m_back_bsdf = std::dynamic_pointer_cast<BSDF>(child);
                } else {
                    throw std::runtime_error("TwoSided: tried to register more than two BSDF instances!");
                }
                break;

            default:
                throw std::runtime_error("TwoSided::add_child(<" + class_type_name(child->get_class_type()) +
                                         ">) is not supported!");
        }
    }

    [[nodiscard]] std::pair<BSDFSample3f, Color3f> sample(const SurfaceIntersection3f &si, float sample1,
                                                          const Point2f &sample2, bool active) const override {
        BSDFSample3f bs;
        Color3f value(0);

        if (Frame3f::cos_theta(si.wi, active) > 0.f) {
            // Front side
            if (m_front_bsdf) {
                std::tie(bs, value) = m_front_bsdf->sample(si, sample1, sample2, active);
            }
        } else {
            // Back side
            SurfaceIntersection3f flipped_si = si;
            flipped_si.wi.z() *= -1.f;

            if (m_back_bsdf) {
                std::tie(bs, value) = m_back_bsdf->sample(flipped_si, sample1, sample2, active);
                bs.wo.z() *= -1.f; // Flip outgoing direction back to original orientation
            }
        }

        return { bs, value };
    }

    [[nodiscard]] Color3f eval(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        if (Frame3f::cos_theta(si.wi, active) > 0.f) {
            // Front side
            return m_front_bsdf ? m_front_bsdf->eval(si, wo, active) : Color3f(0);
        } else {
            // Back side
            SurfaceIntersection3f flipped_si = si;
            flipped_si.wi.z() *= -1.f;

            Vector3f flipped_wo = wo;
            flipped_wo.z() *= -1.f;

            return m_back_bsdf ? m_back_bsdf->eval(flipped_si, flipped_wo, active) : Color3f(0);
        }
    }

    [[nodiscard]] float pdf(const SurfaceIntersection3f &si, const Vector3f &wo, bool active) const override {
        if (Frame3f::cos_theta(si.wi, active) > 0.f) {
            // Front side
            return m_front_bsdf ? m_front_bsdf->pdf(si, wo, active) : 0.f;
        } else {
            // Back side
            SurfaceIntersection3f flipped_si = si;
            flipped_si.wi.z() *= -1.f;

            Vector3f flipped_wo = wo;
            flipped_wo.z() *= -1.f;

            return m_back_bsdf ? m_back_bsdf->pdf(flipped_si, flipped_wo, active) : 0.f;
        }
    }

    [[nodiscard]] std::string to_string() const override {
        std::ostringstream oss;
        oss << "TwoSided[\n"
            << "  front_bsdf = " << (m_front_bsdf ? m_front_bsdf->to_string() : "null") << ",\n"
            << "  back_bsdf = " << (m_back_bsdf ? m_back_bsdf->to_string() : "null") << "\n"
            << "]";
        return oss.str();
    }

    [[nodiscard]] EClassType get_class_type() const override { return EBSDF; }

private:
    std::shared_ptr<BSDF> m_front_bsdf; //< BSDF for the front-facing side
    std::shared_ptr<BSDF> m_back_bsdf;  //< BSDF for the back-facing side (defaults to front BSDF if not provided)
};

REGISTER_CLASS(TwoSided, "two_sided");

M_NAMESPACE_END

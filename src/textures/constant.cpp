#include <components/texture.h>

M_NAMESPACE_BEGIN
	class Constant : public Texture {
	public:
		explicit Constant(const PropertyList& properties) : Texture(false) {
			m_name = properties.get_string("name", "constant");
			m_color = properties.get_color("color", Color3f({1.0f, 1.0f, 1.0f}));
		}

		void add_child(const std::shared_ptr<Object>& child) override {}

		void construct() override {
#ifdef M_DEBUG
			std::cout << "Construct" << class_type_name(get_class_type()) << std::endl;
#endif
		}

		Color3f eval(const SurfaceIntersection3f& si, bool active) override {
			return m_color;
		}

		[[nodiscard]] std::string to_string() const override {
			return "Constant {color = " + m_color.to_string() + "}";
		}

	private:
		Color3f m_color;
	};

	REGISTER_CLASS(Constant, "constant")

M_NAMESPACE_END

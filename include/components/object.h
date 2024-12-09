#pragma once

#include <parse/property_list.h>
#include <memory>
#include <functional>

M_NAMESPACE_BEGIN
	class Object {
	public:
		enum EClassType {
			EScene = 0,
			EMesh,
			EBSDF,
			EEmitter,
			ECamera,
			ESampler,
			EIntegrator,
			EReconstructionFilter,
			ETexture,
			EAccelerate,
			EClassTypeCount
		};

		virtual ~Object() = 0;

		[[nodiscard]] virtual EClassType get_class_type() const = 0;

		virtual void add_child(const std::shared_ptr<Object>& child);

		void set_parent(const std::shared_ptr<Object>& parent);

		// Construct from top to bottom
		virtual void construct();

		[[nodiscard]] virtual std::string to_string() const = 0;

		[[nodiscard]] const std::string& get_name() const;

		static std::string class_type_name(EClassType type) {
			switch (type) {
			case EScene: return "scene";
			case EMesh: return "mesh";
			case EBSDF: return "bsdf";
			case EEmitter: return "emitter";
			case ECamera: return "camera";
			case ESampler: return "sampler";
			case EIntegrator: return "integrator";
			case EReconstructionFilter: return "reconstruction_filter";
			case ETexture: return "texture";
			case EAccelerate: return "accelerate";
			default: return "<unknown>";
			}
		}

	protected:
		std::shared_ptr<Object> parent = nullptr;
		std::string m_name = "Object";
	};

	class ObjectFactory {
	public:
		typedef std::function<std::shared_ptr<Object>(const PropertyList&)> Constructor;

		/**
		 * \brief Register an object constructor with the object factory
		 *
		 * This function is called by the macro \ref NORI_REGISTER_CLASS
		 *
		 * \param name
		 *     An internal name that is associated with this class. This is the
		 *     'type' field found in the scene description XML files
		 *
		 * \param constr
		 *     A function pointer to an anonymous function that is
		 *     able to call the constructor of the class.
		 */
		static void register_class(const std::string& name, const Constructor& constr) {
			if (!constructors)
				constructors = std::make_shared<std::map<std::string, Constructor>>();
			(*constructors)[name] = constr;
		}

		/**
		 * \brief Construct an instance from the class of the given name
		 *
		 * \param name
		 *     An internal name that is associated with this class. This is the
		 *     'type' field found in the scene description XML files
		 *
		 * \param property_list
		 *     A list of properties that will be passed to the constructor
		 *     of the class.
		 */
		static std::shared_ptr<Object> create_instance(const std::string& name, const PropertyList& property_list) {
			if (!constructors || constructors->find(name) == constructors->end())
				throw std::runtime_error("A constructor for class " + name + " could not be found!");
			return (*constructors)[name](property_list);
		}

	private:
		static std::shared_ptr<std::map<std::string, Constructor>> constructors;
	};

	// Macro for registering an object constructor with the ObjectFactory
#define REGISTER_CLASS(cls, name) \
std::shared_ptr<Object> cls ##_create(const PropertyList &list) { \
	return std::make_shared<cls>(list); \
} \
static struct cls ##_ { \
	cls ##_() { \
		ObjectFactory::register_class(name, cls ##_create); \
	} \
} cls ##__;

M_NAMESPACE_END

#include <components/object.h>

M_NAMESPACE_BEGIN
	Object::~Object() = default;

	void Object::add_child(const std::shared_ptr<Object>& child) {
		throw std::runtime_error(
			"Object::add_child() is not implemented for objects of type " + class_type_name(get_class_type()));
	}

	void Object::construct() {}

	void Object::set_parent(const std::shared_ptr<Object>& parent) {
		this->parent = parent;
	}

	const std::string& Object::get_name() const {
		return m_name;
	}


	std::shared_ptr<std::map<std::string, ObjectFactory::Constructor>> ObjectFactory::constructors = nullptr;

M_NAMESPACE_END

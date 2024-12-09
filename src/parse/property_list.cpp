#include <parse/property_list.h>
#include <iostream>

M_NAMESPACE_BEGIN

#define DEFINE_PROPERTY_ACCESSOR(Type, TypeName, XmlName) \
void PropertyList::set_##TypeName(const std::string &name, const Type &value) { \
	if (properties.find(name) != properties.end()) \
		std::cerr << "Property \"" << name <<  "\" was specified multiple times!" << std::endl; \
	Property prop; \
	prop.value = value; \
	prop.type = Property::E##XmlName; \
	properties[name] = prop; \
} \
\
Type PropertyList::get_##TypeName(const std::string &name) const { \
	auto it = properties.find(name); \
	if (it == properties.end()) \
		throw std::runtime_error("Property '" + name + "' is missing!"); \
	if (it->second.type != Property::E##XmlName) \
		throw std::runtime_error("Property '" + name + "' has the wrong type! " \
		"(expected <" #XmlName ">)!"); \
	return std::get<Type>(it->second.value); \
} \
\
Type PropertyList::get_##TypeName(const std::string &name, const Type &default_value) const { \
	auto it = properties.find(name); \
	if (it == properties.end()) \
		return default_value; \
	if (it->second.type != Property::E##XmlName) \
		throw std::runtime_error("Property '" + name + "' has the wrong type! " \
		"(expected <" #XmlName ">)!"); \
	return std::get<Type>(it->second.value); \
}

DEFINE_PROPERTY_ACCESSOR(bool, boolean, Boolean)
DEFINE_PROPERTY_ACCESSOR(int, integer, Integer)
DEFINE_PROPERTY_ACCESSOR(float, float, Float)
DEFINE_PROPERTY_ACCESSOR(Color3f, color, Color)
DEFINE_PROPERTY_ACCESSOR(Point3f, point, Point)
DEFINE_PROPERTY_ACCESSOR(Vector3f, vector, Vector)
DEFINE_PROPERTY_ACCESSOR(std::string, string, String)
DEFINE_PROPERTY_ACCESSOR(Transform4f, transform, Transform)

M_NAMESPACE_END

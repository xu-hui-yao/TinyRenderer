#pragma once

#include <core/spectrum.h>
#include <core/transform.h>
#include <core/array.h>
#include <map>
#include <variant>

M_NAMESPACE_BEGIN

class PropertyList {
public:
	PropertyList() = default;

	// Set a boolean property
	void set_boolean(const std::string& name, const bool& value);

	// Get a boolean property, and throw an exception if it does not exist
	[[nodiscard]] bool get_boolean(const std::string& name) const;

	// Get a boolean property, and use a default value if it does not exist
	[[nodiscard]] bool get_boolean(const std::string& name, const bool& default_value) const;

	// Set an integer property
	void set_integer(const std::string& name, const int& value);

	// Get an integer property, and throw an exception if it does not exist
	[[nodiscard]] int get_integer(const std::string& name) const;

	// Get am integer property, and use a default value if it does not exist
	[[nodiscard]] int get_integer(const std::string& name, const int& default_value) const;

	// Set a float property
	void set_float(const std::string& name, const float& value);

	// Get a float property, and throw an exception if it does not exist
	[[nodiscard]] float get_float(const std::string& name) const;

	// Get a float property, and use a default value if it does not exist
	[[nodiscard]] float get_float(const std::string& name, const float& default_value) const;

	// Set a string property
	void set_string(const std::string& name, const std::string& value);

	// Get a string property, and throw an exception if it does not exist
	[[nodiscard]] std::string get_string(const std::string& name) const;

	// Get a string property, and use a default value if it does not exist
	[[nodiscard]] std::string get_string(const std::string& name, const std::string& default_value) const;

	// Set a color property
	void set_color(const std::string& name, const Color3f& value);

	// Get a color property, and throw an exception if it does not exist
	[[nodiscard]] Color3f get_color(const std::string& name) const;

	// Get a color property, and use a default value if it does not exist
	[[nodiscard]] Color3f get_color(const std::string& name, const Color3f& default_value) const;

	// Set a point property
	void set_point(const std::string& name, const Point3f& value);

	// Get a point property, and throw an exception if it does not exist
	[[nodiscard]] Point3f get_point(const std::string& name) const;

	// Get a point property, and use a default value if it does not exist
	[[nodiscard]] Point3f get_point(const std::string& name, const Point3f& default_value) const;

	// Set a vector property
	void set_vector(const std::string& name, const Vector3f& value);

	// Get a vector property, and throw an exception if it does not exist
	[[nodiscard]] Vector3f get_vector(const std::string& name) const;

	// Get a vector property, and use a default value if it does not exist
	[[nodiscard]] Vector3f get_vector(const std::string& name, const Vector3f& default_value) const;

	// Set a transform property
	void set_transform(const std::string& name, const Transform4f& value);

	// Get a transform property, and throw an exception if it does not exist
	[[nodiscard]] Transform4f get_transform(const std::string& name) const;

	// Get a transform property, and use a default value if it does not exist
	[[nodiscard]] Transform4f get_transform(const std::string& name, const Transform4f& default_value) const;

private:
	/* Custom variant data type */
	struct Property {
		enum {
			EBoolean,
			EInteger,
			EFloat,
			EString,
			EColor,
			EPoint,
			EVector,
			ETransform
		} type {};
		std::variant<bool, int, float, std::string, Color3f, Point3f, Vector3f, Transform4f> value;
	};

	std::map<std::string, Property> properties;
};

M_NAMESPACE_END

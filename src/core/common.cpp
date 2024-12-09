#include <core/array.h>
#include <core/common.h>
#include <core/spectrum.h>
#include <iomanip>

M_NAMESPACE_BEGIN

std::string indent(const std::string &string, int amount) {
	std::istringstream iss(string);
	std::ostringstream oss;
	std::string spacer(amount, ' ');
	bool first_line = true;
	for (std::string line; std::getline(iss, line); ) {
		if (!first_line)
			oss << spacer;
		oss << line;
		if (!iss.eof())
			oss << std::endl;
		first_line = false;
	}
	return oss.str();
}

std::string to_lower(const std::string& value) {
	std::string result;
	result.resize(value.size());
	std::transform(value.begin(), value.end(), result.begin(), ::tolower);
	return result;
}

bool to_bool(const std::string& str) {
	std::string value = to_lower(str);
	if (value == "false")
		return false;
	else if (value == "true")
		return true;
	else
		throw std::runtime_error("Could not parse boolean value " + str);
}

int to_int(const std::string& str) {
	char* end_ptr = nullptr;
	auto result = static_cast<int>(strtol(str.c_str(), &end_ptr, 10));
	if (*end_ptr != '\0')
		throw std::runtime_error("Could not parse integer value " + str);
	return result;
}

unsigned int to_uint(const std::string& str) {
	char* end_ptr = nullptr;
	auto result = static_cast<unsigned int>(strtoul(str.c_str(), &end_ptr, 10));
	if (*end_ptr != '\0')
		throw std::runtime_error("Could not parse unsigned int value " + str);
	return result;
}

float to_float(const std::string& str) {
	char* end_ptr = nullptr;
	auto result = strtof(str.c_str(), &end_ptr);
	if (*end_ptr != '\0')
		throw std::runtime_error("Could not parse float value " + str);
	return result;
}

Vector3f to_vector3f(const std::string& str) {
	std::vector<std::string> tokens = tokenize(str);
	if (tokens.size() != 3)
		throw std::runtime_error("Expected 3 values");
	Vector3f result;
	for (int i = 0; i < 3; ++i)
		result(i) = to_float(tokens[i]);
	return result;
}

Point3f to_point3f(const std::string& str) {
	std::vector<std::string> tokens = tokenize(str);
	if (tokens.size() != 3)
		throw std::runtime_error("Expected 3 values");
	Point3f result;
	for (int i = 0; i < 3; ++i)
		result(i) = to_float(tokens[i]);
	return result;
}

Color3f to_color3f(const std::string& str) {
	std::vector<std::string> tokens = tokenize(str);
	if (tokens.size() != 3)
		throw std::runtime_error("Expected 3 values");
	Color3f result;
	for (int i = 0; i < 3; ++i)
		result(i) = to_float(tokens[i]);
	return result;
}

std::vector<std::string> tokenize(const std::string& string, const std::string& delim, bool include_empty) {
	std::string::size_type last_pos = 0, pos = string.find_first_of(delim, last_pos);
	std::vector<std::string> tokens;

	while (last_pos != std::string::npos) {
		if (pos != last_pos || include_empty)
			tokens.push_back(string.substr(last_pos, pos - last_pos));
		last_pos = pos;
		if (last_pos != std::string::npos) {
			last_pos += 1;
			pos = string.find_first_of(delim, last_pos);
		}
	}

	return tokens;
}

bool end_with(const std::string& value, const std::string& ending) {
	if (ending.size() > value.size())
		return false;
	return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

std::string time_string(double time, bool precise) {
	if (std::isnan(time) || std::isinf(time))
		return "inf";

	std::string suffix = "ms";
	if (time > 1000) {
		time /= 1000;
		suffix = "s";
		if (time > 60) {
			time /= 60;
			suffix = "m";
			if (time > 60) {
				time /= 60;
				suffix = "h";
				if (time > 12) {
					time /= 12;
					suffix = "d";
				}
			}
		}
	}

	std::ostringstream os;
	os << std::setprecision(precise ? 4 : 1)
		<< std::fixed << time << suffix;

	return os.str();
}

std::string mem_string(size_t size, bool precise) {
	auto value = static_cast<double>(size);
	const char* suffixes[] = {
		"B", "KiB", "MiB", "GiB", "TiB", "PiB"
	};
	int suffix = 0;
	while (suffix < 5 && value > 1024.0f) {
		value /= 1024.0f;
		++suffix;
	}

	std::ostringstream os;
	os << std::setprecision(suffix == 0 ? 0 : (precise ? 4 : 1))
		<< std::fixed << value << " " << suffixes[suffix];

	return os.str();
}

std::shared_ptr<filesystem::resolver> get_file_resolver() {
	static std::shared_ptr<filesystem::resolver> resolver = std::make_shared<filesystem::resolver>();
	return resolver;
}

M_NAMESPACE_END

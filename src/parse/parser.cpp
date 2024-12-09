#include <parse/parser.h>
#include <parse/property_list.h>
#include <src/pugixml.hpp>
#include <fstream>
#include <set>

M_NAMESPACE_BEGIN
	std::shared_ptr<Object> load_from_xml(const std::string& filename) {
		pugi::xml_document doc;
		filesystem::path file = get_file_resolver()->resolve(filesystem::path(filename));
		get_file_resolver()->append(file.parent_path());
		pugi::xml_parse_result result = doc.load_file(file.str().c_str());

		auto offset = [&](ptrdiff_t pos) -> std::string {
			std::fstream is(filename);
			int line = 0, line_start = 0, offset_ = 0;
			while (is.good()) {
				char buffer[1024];
				is.read(buffer, sizeof(buffer));
				for (int i = 0; i < is.gcount(); ++i) {
					if (buffer[i] == '\n') {
						if (offset_ + i >= pos)
							return "line " + std::to_string(line + 1) + ", col " + std::to_string(pos - line_start);
						++line;
						line_start = offset_ + i;
					}
				}
				offset_ += static_cast<int>(is.gcount());
			}
			return "byte offset " + std::to_string(pos);
		};

		if (!result) {
			throw std::runtime_error("Error while parsing " + filename + ": " + result.description() +
				" (at " + offset(result.offset) + ")");
		}

		/* Set of supported XML tags */
		enum TagType {
			/* Object classes */
			EScene = Object::EScene,
			EMesh = Object::EMesh,
			EBSDF = Object::EBSDF,
			EEmitter = Object::EEmitter,
			ECamera = Object::ECamera,
			ESampler = Object::ESampler,
			EIntegrator = Object::EIntegrator,
			EReconstructionFilter = Object::EReconstructionFilter,
			ETexture = Object::ETexture,
			EAccelerate = Object::EAccelerate,

			/* Properties */
			EBoolean = Object::EClassTypeCount,
			EInteger,
			EFloat,
			EString,
			EPoint,
			EVector,
			EColor,
			ETransform,
			ETranslate,
			EMatrix,
			ERotate,
			EScale,
			ELookAt,

			EInvalid
		};

		/* Create a mapping from tag names to tag IDs */
		std::map<std::string, TagType> tags;
		tags["scene"] = EScene;
		tags["mesh"] = EMesh;
		tags["bsdf"] = EBSDF;
		tags["emitter"] = EEmitter;
		tags["camera"] = ECamera;
		tags["sampler"] = ESampler;
		tags["integrator"] = EIntegrator;
		tags["rfilter"] = EReconstructionFilter;
		tags["texture"] = ETexture;
		tags["accelerate"] = EAccelerate;

		tags["boolean"] = EBoolean;
		tags["integer"] = EInteger;
		tags["float"] = EFloat;
		tags["string"] = EString;
		tags["point"] = EPoint;
		tags["vector"] = EVector;
		tags["color"] = EColor;
		tags["transform"] = ETransform;
		tags["translate"] = ETranslate;
		tags["matrix"] = EMatrix;
		tags["rotate"] = ERotate;
		tags["scale"] = EScale;
		tags["look_at"] = ELookAt;

		/* Helper function to check if attributes are fully specified */
		auto check_attributes = [&](const pugi::xml_node& node, std::set<std::string> attrs) {
			for (auto attr : node.attributes()) {
				auto it = attrs.find(attr.name());
				if (it == attrs.end()) {
					throw std::runtime_error("Error while parsing " + filename + ": unexpected attribute " +
						attr.name() + " in " + node.name() + " at " + offset(node.offset_debug()));
				}
				attrs.erase(it);
			}
			if (!attrs.empty()) {
				throw std::runtime_error("Error while parsing " + filename + ": missing attribute " +
					*attrs.begin() + " in " + node.name() + " at " + offset(node.offset_debug()));
			}
		};

		Matrix4f transform;

		/* Helper function to parse XML node (recursive) */
		std::function<std::shared_ptr<Object>(pugi::xml_node&, PropertyList&, int)> parse_tag = [&](
			pugi::xml_node& node, PropertyList& list, int parent_tag) -> std::shared_ptr<Object> {
			/* Skip over comments */
			if (node.type() == pugi::node_comment || node.type() == pugi::node_declaration) {
				return nullptr;
			}

			if (node.type() != pugi::node_element) {
				throw std::runtime_error(
					"Error while parsing " + filename + ": unexpected content at " + offset(node.offset_debug()));
			}

			/* Look up the name of the current element */
			auto it = tags.find(node.name());
			if (it == tags.end()) {
				throw std::runtime_error("Error while parsing " + filename + ": unexpected tag " +
					node.name() + " at " + offset(node.offset_debug()));
			}
			int tag = it->second;

			/* Perform some safety checks to make sure that the XML tree really makes sense */
			bool has_parent = parent_tag != EInvalid;
			bool parent_is_object = has_parent && parent_tag < Object::EClassTypeCount;
			bool current_is_object = tag < Object::EClassTypeCount;
			bool parent_is_transform = parent_tag == ETransform;
			bool current_is_transform_op = tag == ETranslate || tag == ERotate || tag == EScale ||
				tag == ELookAt || tag == EMatrix;

			if (!has_parent && !current_is_object) {
				throw std::runtime_error("Error while parsing " + filename + ": root element " + node.name() +
					" must be a object (at " + offset(node.offset_debug()) + ")");
			}

			if (parent_is_transform != current_is_transform_op) {
				throw std::runtime_error("Error while parsing " + filename + ": transform nodes "
					"can only contain transform operations (at " + offset(node.offset_debug()) + ")");
			}

			if (has_parent && !parent_is_object && !(parent_is_transform && current_is_transform_op)) {
				throw std::runtime_error(
					"Error while parsing " + filename + ": node " + node.name() +
					" requires a object as parent (at " + offset(node.offset_debug()) + ")");
			}

			if (tag == EScene) {
				node.append_attribute("type") = "scene";
			} else if (tag == ETransform) {
				transform.set_identity();
			}

			PropertyList property_list;
			std::vector<std::shared_ptr<Object>> children;
			for (pugi::xml_node& ch : node.children()) {
				if (std::shared_ptr<Object> child = parse_tag(ch, property_list, tag))
					children.push_back(child);
			}

			std::shared_ptr<Object> result_ = nullptr;
			try {
				if (current_is_object) {
					check_attributes(node, {"type"});

					/* This is an object, first instantiate it */
					result_ = ObjectFactory::create_instance(node.attribute("type").value(), property_list);

#ifdef M_DEBUG
					std::cout << "Create instance " << Object::class_type_name(result_->get_class_type()) << std::endl;
#endif

					if (result_->get_class_type() != tag) {
						throw std::runtime_error(
							"Unexpectedly constructed an object "
							"of type <" + Object::class_type_name(result_->get_class_type()) + "> (expected type <" +
							Object::class_type_name(static_cast<Object::EClassType>(tag)) + ">): " + result_->
							to_string());
					}

					/* Add all children */
					for (const auto& ch : children) {
						result_->add_child(ch);
						ch->set_parent(result_);
#ifdef M_DEBUG
						std::cout << "Instance " << Object::class_type_name(result_->get_class_type()) << " add child "
							<< Object::class_type_name(ch->get_class_type()) << std::endl;
#endif
					}
				} else {
					/* This is a property */
					switch (tag) {
					case EString: {
						check_attributes(node, {"name", "value"});
						list.set_string(node.attribute("name").value(), node.attribute("value").value());
					}
					break;
					case EFloat: {
						check_attributes(node, {"name", "value"});
						list.set_float(node.attribute("name").value(), to_float(node.attribute("value").value()));
					}
					break;
					case EInteger: {
						check_attributes(node, {"name", "value"});
						list.set_integer(node.attribute("name").value(), to_int(node.attribute("value").value()));
					}
					break;
					case EBoolean: {
						check_attributes(node, {"name", "value"});
						list.set_boolean(node.attribute("name").value(), to_bool(node.attribute("value").value()));
					}
					break;
					case EPoint: {
						check_attributes(node, {"name", "value"});
						list.set_point(node.attribute("name").value(),
						               Point3f(to_point3f(node.attribute("value").value())));
					}
					break;
					case EVector: {
						check_attributes(node, {"name", "value"});
						list.set_vector(node.attribute("name").value(),
						                Vector3f(to_vector3f(node.attribute("value").value())));
					}
					break;
					case EColor: {
						check_attributes(node, {"name", "value"});
						list.set_color(node.attribute("name").value(),
						               Color3f(to_color3f(node.attribute("value").value())));
					}
					break;
					case ETransform: {
						check_attributes(node, {"name"});
						list.set_transform(node.attribute("name").value(), Transform4f(transform));
					}
					break;
					case ETranslate: {
						check_attributes(node, {"value"});
						Vector3f v = to_vector3f(node.attribute("value").value());
						transform = Transform4f::translate(v.x(), v.y(), v.z()).get_transform() * transform;
					}
					break;
					case EMatrix: {
						check_attributes(node, {"value"});
						std::vector<std::string> tokens = tokenize(node.attribute("value").value());
						if (tokens.size() != 16)
							throw std::runtime_error("Expected 16 values");
						Matrix4f matrix;
						for (int i = 0; i < 4; ++i)
							for (int j = 0; j < 4; ++j)
								matrix(i, j) = to_float(tokens[i * 4 + j]);
						transform = matrix * transform;
					}
					break;
					case EScale: {
						bool valid = true;
						check_attributes(node, {"value"});
						Vector3f v = to_vector3f(node.attribute("value").value());
						transform = Transform4f::scale(v, valid).get_transform() * transform;
						if (!valid) {
							throw std::runtime_error("Invalid scale");
						}
					}
					break;
					case ERotate: {
						bool valid = true;
						check_attributes(node, {"angle", "axis"});
						float angle = deg_to_rad(to_float(node.attribute("angle").value()));
						Vector3f axis = to_vector3f(node.attribute("axis").value());
						transform = Transform4f::rotate(angle, axis, valid).get_transform() * transform;
						if (!valid) {
							throw std::runtime_error("Rotate axis can not be normalized");
						}
					}
					break;
					case ELookAt: {
						check_attributes(node, {"origin", "target", "up"});
						Vector3f origin = to_vector3f(node.attribute("origin").value());
						Vector3f target = to_vector3f(node.attribute("target").value());
						Vector3f up = to_vector3f(node.attribute("up").value());

						bool valid = true;
						transform = Transform4f::look_at(origin, target, up, valid).get_transform() * transform;
						if (!valid) {
							throw std::runtime_error("Look at matrix initialized failed");
						}
					}
					break;

					default: throw std::runtime_error("Unhandled element " + std::string(node.name()));
					};
				}
			} catch (const std::exception& e) {
				throw std::runtime_error(
					"Error while parsing " + filename + ": " + e.what() + " (at " + offset(node.offset_debug()) + ")");
			}

			return result_;
		};

		PropertyList list;
		return parse_tag(*doc.begin(), list, EInvalid);
	}

M_NAMESPACE_END

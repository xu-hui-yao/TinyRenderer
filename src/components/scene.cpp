#include <components/scene.h>
#include <components/camera.h>
#include <iostream>
#include <core/intersection.h>

M_NAMESPACE_BEGIN
	Scene::Scene(const PropertyList&) {
		m_num_emitters = 0;
		m_emitter_pmf = 0;
	}

	void Scene::construct() {
		if (!camera) {
			throw std::runtime_error("No Camera was specified!");
		}

		if (meshes.empty()) {
			throw std::runtime_error("No mesh was specified!");
		}

		if (!integrator) {
			throw std::runtime_error("No Integrator was specified!");
		}

		if (!sampler) {
			throw std::runtime_error("No Sampler was specified!");
		}

		for (const auto& mesh : meshes) {
			mesh->construct();
			if (mesh->is_emitter()) {
				m_num_emitters += 1;
				emitters.emplace_back(mesh->get_emitter());
			}
			accel->add_mesh(mesh);
		}

		accel->construct();
		camera->construct();
		integrator->construct();
		sampler->construct();

		if (m_num_emitters == 0) {
			throw std::runtime_error("No emitter was specified!");
		}

		m_emitter_pmf = 1.0f / static_cast<float>(m_num_emitters);

#ifdef M_DEBUG
		std::cout << std::endl;
		std::cout << "Configuration: " << to_string() << std::endl;
		std::cout << std::endl;
#endif
	}

	void Scene::add_child(const std::shared_ptr<Object>& obj) {
		switch (obj->get_class_type()) {
		case EMesh: {
			auto mesh = std::dynamic_pointer_cast<Mesh>(obj);
			meshes.push_back(mesh);
		}
		break;

		case ECamera:
			if (camera) {
				throw std::runtime_error("There can only be one Camera per scene!");
			}
			camera = std::dynamic_pointer_cast<Camera>(obj);
			break;

		case EIntegrator:
			if (integrator) {
				throw std::runtime_error("There can only be one Integrator per scene!");
			}
			integrator = std::dynamic_pointer_cast<Integrator>(obj);
			break;

		case ESampler:
			if (sampler) {
				throw std::runtime_error("There can only be one Sampler per scene!");
			}
			sampler = std::dynamic_pointer_cast<Sampler>(obj);
			break;

		case EAccelerate:
			if (accel) {
				throw std::runtime_error("There can only be one Accel per scene!");
			}
			accel = std::dynamic_pointer_cast<Accel>(obj);
			break;

		default:
			throw std::runtime_error(
				"Scene::add_child(<" + class_type_name(obj->get_class_type()) + ">) is not supported!");
		}
	}

	std::string Scene::to_string() const {
		std::string mesh_string;
		for (size_t i = 0; i < meshes.size(); ++i) {
			mesh_string += std::string("  ") + indent(meshes[i]->to_string(), 2);
			if (i + 1 < meshes.size())
				mesh_string += ",";
			mesh_string += "\n";
		}

		return std::string("Scene[\n") +
			std::string("  Camera = ") + indent(camera->to_string()) + std::string(",\n") +
			std::string("  Integrator = ") + indent(integrator->to_string()) + std::string(",\n") +
			std::string("  Sampler = ") + indent(sampler->to_string()) + std::string("\n") +
			std::string("  Meshes = {\n") + std::string("  ") + indent(mesh_string) + std::string(" }\n") +
			std::string("]");
	}

	std::pair<DirectionSample3f, Color3f> Scene::sample_emitter_direction(
		const SurfaceIntersection3f& its, const Point2f& sample_, bool test_visibility, bool& active) const {
		Point2f sample(sample_);
		DirectionSample3f ds;
		Color3f spec;
		if (m_num_emitters > 1) {
			// Randomly pick an emitter
			auto [index, emitter_weight, sample_x_re] = sample_emitter(sample.x(), active);
			sample.x() = sample_x_re;

			// Sample a direction towards the emitter
			std::shared_ptr<Emitter> emitter = emitters[index];
			std::tie(ds, spec) = emitter->sample_direction(its, sample, active);

			// Account for the discrete probability of sampling this emitter
			ds.pdf *= pdf_emitter(index, active);
			spec *= emitter_weight;

			active &= ds.pdf != 0.0f;

			// Mark occluded samples as invalid if requested by the user
			if (test_visibility && active) {
				if (accel->ray_test(its.spawn_ray_to(ds.p))) {
					spec = Color3f(0.0f);
					ds.pdf = 0.0f;
				}
			}
		} else if (m_num_emitters == 1) {
			std::tie(ds, spec) = emitters[0]->sample_direction(its, sample, active);
			active &= ds.pdf != 0.0f;

			// Mark occluded samples as invalid if requested by the user
			if (test_visibility && active) {
				if (accel->ray_test(its.spawn_ray_to(ds.p))) {
					spec = Color3f(0.0f);
					ds.pdf = 0.0f;
				}
			}
		} else {
			spec = Color3f(0.0f);
		}

		return {ds, spec};
	}

	float Scene::pdf_emitter_direction(const Intersection3f& it, const DirectionSample3f& ds,
	                                   bool& active) const {
		return ds.emitter->pdf_direction(it, ds, active) * m_emitter_pmf;
	}

	std::tuple<uint32_t, float, float> Scene::sample_emitter(float sample, bool& active) const {
		if (m_num_emitters == 0) {
			return {static_cast<uint32_t>(-1), 0.0f, sample};
		} else if (m_num_emitters == 1) {
			return {static_cast<uint32_t>(0), 1.0f, sample};
		} else {
			auto num_emitters_f = static_cast<float>(m_num_emitters);
			float index_sample_scaled = sample * num_emitters_f;
			uint32_t index = M_MIN(static_cast<uint32_t>(index_sample_scaled), m_num_emitters - 1);
			return {index, num_emitters_f, index_sample_scaled - static_cast<float>(index)};
		}
	}

	float Scene::pdf_emitter(uint32_t index, bool& active) const {
		return m_emitter_pmf;
	}

	REGISTER_CLASS(Scene, "scene");

M_NAMESPACE_END

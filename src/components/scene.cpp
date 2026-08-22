#include <components/camera.h>
#include <components/scene.h>
#include <core/intersection.h>

M_NAMESPACE_BEGIN
Scene::Scene(const PropertyList &) {
    m_num_emitters = 0;
    m_emitter_pmf  = 0;
}

void Scene::construct() {
    if (!m_camera) {
        throw std::runtime_error("No Camera was specified!");
    }

    if (m_meshes.empty()) {
        throw std::runtime_error("No mesh was specified!");
    }

    if (!m_integrator) {
        throw std::runtime_error("No Integrator was specified!");
    }

    if (!m_sampler) {
        throw std::runtime_error("No Sampler was specified!");
    }

    for (const auto &mesh : m_meshes) {
        mesh->construct();
        if (mesh->is_emitter()) {
            m_num_emitters += 1;
            mesh->get_emitter()->construct();
            m_emitters.emplace_back(mesh->get_emitter());
        }
        m_accel->add_mesh(mesh);
    }

    if (m_environment) {
        m_emitters.emplace_back(m_environment);
        m_environment->construct();
        m_num_emitters += 1;
    }

    if (m_num_emitters == 0) {
        throw std::runtime_error("No emitter was specified!");
    }

    m_emitter_pmf = 1.0f / static_cast<float>(m_num_emitters);

    m_accel->construct();
    m_camera->construct();
    m_integrator->construct();
    m_sampler->construct();

    // Optional - a scene without a <denoiser> simply gets its raw render
    // written out, exactly as before this feature existed.
    if (m_denoiser) {
        m_denoiser->construct();
    }

    for (const auto &emitter : m_emitters) {
        emitter->set_scene(std::dynamic_pointer_cast<Scene>(shared_from_this()));
    }

#ifdef M_DEBUG
    std::cout << std::endl;
    std::cout << "Configuration: " << to_string() << std::endl;
    std::cout << std::endl;
#endif
}

void Scene::add_child(const std::shared_ptr<Object> &obj) {
    switch (obj->get_class_type()) {
        case EMesh: {
            auto mesh = std::dynamic_pointer_cast<Mesh>(obj);
            m_meshes.push_back(mesh);
        } break;

        case ECamera:
            if (m_camera) {
                throw std::runtime_error("There can only be one Camera per scene!");
            }
            m_camera = std::dynamic_pointer_cast<Camera>(obj);
            break;

        case EIntegrator:
            if (m_integrator) {
                throw std::runtime_error("There can only be one Integrator per scene!");
            }
            m_integrator = std::dynamic_pointer_cast<Integrator>(obj);
            break;

        case ESampler:
            if (m_sampler) {
                throw std::runtime_error("There can only be one Sampler per scene!");
            }
            m_sampler = std::dynamic_pointer_cast<Sampler>(obj);
            break;

        case EAccelerate:
            if (m_accel) {
                throw std::runtime_error("There can only be one Accel per scene!");
            }
            m_accel = std::dynamic_pointer_cast<Accel>(obj);
            break;

        case EEmitter:
            if (m_environment) {
                throw std::runtime_error("There can only be one Environment map per scene!");
            }
            m_environment = std::dynamic_pointer_cast<Emitter>(obj);
            break;

        case EDenoiser:
            if (m_denoiser) {
                throw std::runtime_error("There can only be one Denoiser per scene!");
            }
            m_denoiser = std::dynamic_pointer_cast<Denoiser>(obj);
            break;

        default:
            throw std::runtime_error("Scene::add_child(<" + class_type_name(obj->get_class_type()) +
                                     ">) is not supported!");
    }
}

bool Scene::ray_intersect(const Ray3f &ray, SurfaceIntersection3f &its, bool shadow_ray) const {
    if (m_accel->ray_intersect(ray, its, shadow_ray)) {
        return true;
    }

    if (m_environment) {
        its.t             = M_MAX_FLOAT;
        its.p             = ray.o() + ray.d() * its.t;
        its.n             = -ray.d();
        its.shading_frame = its.geometric_frame = Frame3f(its.n);
        its.dp_du = its.dp_dv = Vector3f({ 0, 0, 0 });
        its.mesh_id           = M_INVALID_INDEX;
        its.primitive_index   = -1;
        its.wi                = -ray.d();
        its.uv                = Point2f({ 0, 0 });
        return true;
    }

    return false;
}

std::string Scene::to_string() const {
    std::string mesh_string;
    for (size_t i = 0; i < m_meshes.size(); ++i) {
        mesh_string += std::string("  ") + indent(m_meshes[i]->to_string(), 2);
        if (i + 1 < m_meshes.size())
            mesh_string += ",";
        mesh_string += "\n";
    }

    return std::string("Scene[\n") + std::string("  Camera = ") + indent(m_camera->to_string()) + std::string(",\n") +
           std::string("  Integrator = ") + indent(m_integrator->to_string()) + std::string(",\n") +
           std::string("  Sampler = ") + indent(m_sampler->to_string()) + std::string("\n") +
           std::string("  Meshes = {\n") + std::string("  ") + indent(mesh_string) + std::string(" }\n") +
           std::string("]");
}

std::pair<DirectionSample3f, Color3f> Scene::sample_emitter_direction(const SurfaceIntersection3f &its,
                                                                      const Point2f &sample_, bool test_visibility,
                                                                      bool active) const {
    Point2f sample(sample_);
    DirectionSample3f ds;
    Color3f spec;
    if (m_num_emitters > 1) {
        // Randomly pick an emitter
        auto [index, emitter_weight, sample_x_re] = sample_emitter(sample.x(), active);
        sample.x()                                = sample_x_re;

        // Sample a direction towards the emitter
        std::shared_ptr<Emitter> emitter = m_emitters[index];
        std::tie(ds, spec)               = emitter->sample_direction(its, sample, active);

        // Account for the discrete probability of sampling this emitter
        ds.pdf *= pdf_emitter(index, active);
        spec *= emitter_weight;

        active &= ds.pdf != 0.0f;

        // Mark occluded samples as invalid if requested by the user
        if (test_visibility && active) {
            if (m_accel->ray_test(its.spawn_ray_to(ds.p))) {
                spec   = Color3f(0.0f);
                ds.pdf = 0.0f;
            }
        }
    } else if (m_num_emitters == 1) {
        std::tie(ds, spec) = m_emitters[0]->sample_direction(its, sample, active);
        active &= ds.pdf != 0.0f;

        // Mark occluded samples as invalid if requested by the user
        if (test_visibility && active) {
            if (m_accel->ray_test(its.spawn_ray_to(ds.p))) {
                spec   = Color3f(0.0f);
                ds.pdf = 0.0f;
            }
        }
    } else {
        spec = Color3f(0.0f);
    }

    return { ds, spec };
}

float Scene::pdf_emitter_direction(const Intersection3f &it, const DirectionSample3f &ds, bool active) const {
    return ds.emitter->pdf_direction(it, ds, active) * m_emitter_pmf;
}

std::tuple<uint32_t, float, float> Scene::sample_emitter(float sample, bool active) const {
    if (m_num_emitters == 0) {
        return { static_cast<uint32_t>(-1), 0.0f, sample };
    } else if (m_num_emitters == 1) {
        return { static_cast<uint32_t>(0), 1.0f, sample };
    } else {
        auto num_emitters_f       = static_cast<float>(m_num_emitters);
        float index_sample_scaled = sample * num_emitters_f;
        uint32_t index            = M_MIN(static_cast<uint32_t>(index_sample_scaled), m_num_emitters - 1);
        return { index, num_emitters_f, index_sample_scaled - static_cast<float>(index) };
    }
}

float Scene::pdf_emitter(uint32_t index, bool active) const { return m_emitter_pmf; }

GPUScene Scene::build_gpu_scene() const {
    GPUSceneBuilder builder;
    GPUScene &result = builder.result;

    result.meshes.reserve(m_meshes.size());

    for (uint32_t mesh_id = 0; mesh_id < m_meshes.size(); ++mesh_id) {
        const auto &mesh = m_meshes[mesh_id];

        const auto &positions = mesh->get_vertex_positions();
        const auto &normals   = mesh->get_vertex_normals();
        const auto &uvs       = mesh->get_vertex_tex_coords();
        const auto &faces     = mesh->get_indices();

        GPUMeshInfo info;
        info.base_vertex    = static_cast<uint32_t>(result.vertex_positions.size());
        info.vertex_count   = static_cast<uint32_t>(positions.size());
        info.base_index     = static_cast<uint32_t>(result.indices.size());
        info.triangle_count = static_cast<uint32_t>(faces.size());
        info.has_normals    = !normals.empty();
        info.has_uvs        = !uvs.empty();

        // Concatenate this mesh's vertex attributes into the global buffers.
        // Meshes without normals/UVs get filled with placeholders so that
        // vertex_positions/vertex_normals/vertex_uvs stay a fixed-format,
        // uniformly-indexable vertex buffer (GPUMeshInfo::has_normals/
        // has_uvs record whether a given mesh's range is "real" data).
        for (uint32_t v = 0; v < info.vertex_count; ++v) {
            result.vertex_positions.push_back(positions[v]);
            result.vertex_normals.push_back(info.has_normals ? normals[v] : Normal3f(0.0f));
            result.vertex_uvs.push_back(info.has_uvs ? uvs[v] : Point2f(0.0f));
        }

        // Convert mesh-local face indices to GLOBAL vertex indices by adding
        // this mesh's base_vertex, exactly like a typical GPU index buffer.
        for (uint32_t f = 0; f < info.triangle_count; ++f) {
            result.indices.push_back(info.base_vertex + static_cast<uint32_t>(faces[f](0)));
            result.indices.push_back(info.base_vertex + static_cast<uint32_t>(faces[f](1)));
            result.indices.push_back(info.base_vertex + static_cast<uint32_t>(faces[f](2)));
        }

        info.material_id = builder.add_material(mesh->get_bsdf());
        if (mesh->is_emitter()) {
            info.light_id = builder.add_light(mesh->get_emitter(), mesh_id);

            // Populate the area-weighted triangle-sampling CDF for this
            // light, exactly mirroring Mesh::m_area_pmf (DiscreteDistribution1f
            // built from per-triangle surface_area()) so GPU-side area
            // sampling matches the CPU renderer's distribution exactly (see
            // GPULight::cdf_offset/cdf_count/inv_total_area in gpu_scene.h).
            GPULight &light   = result.lights[info.light_id];
            light.cdf_offset  = static_cast<uint32_t>(result.light_triangle_cdf.size());
            light.cdf_count   = info.triangle_count;
            light.base_triangle = info.base_index / 3;

            float running_sum = 0.0f;
            for (uint32_t t = 0; t < info.triangle_count; ++t) {
                running_sum += mesh->surface_area(t);
                result.light_triangle_cdf.push_back(running_sum);
            }
            light.inv_total_area = running_sum > 0.0f ? 1.0f / running_sum : 0.0f;
            // Normalize the CDF to [0, 1] so the GPU-side binary search can
            // work directly against a canonical [0,1) sample u (matching
            // DiscreteDistribution1f::eval_cdf_normalized()'s convention).
            if (running_sum > 0.0f) {
                for (uint32_t t = 0; t < info.triangle_count; ++t) {
                    result.light_triangle_cdf[light.cdf_offset + t] *= light.inv_total_area;
                }
            }
        }

        // Denormalize this mesh's material/light onto each of its triangles,
        // so a GPU shading kernel can resolve a hit triangle in O(1) without
        // searching which mesh's [base_index/3, base_index/3+triangle_count)
        // range it falls into.
        result.triangle_material_id.insert(result.triangle_material_id.end(), info.triangle_count, info.material_id);
        result.triangle_light_id.insert(result.triangle_light_id.end(), info.triangle_count, info.light_id);

        result.meshes.push_back(info);
    }

    if (m_environment) {
        result.environment_light_id = builder.add_light(m_environment, M_INVALID_INDEX);
    }

    if (m_camera) {
        result.camera = m_camera->to_gpu_camera();
    }

    // Export the flat BVH for GPU-side traversal, if the configured Accel
    // implementation has one (currently only BVHAccel does). This must run
    // AFTER the mesh loop above, since it converts (mesh_id, local triangle
    // index) primitive references into global triangle ids using
    // result.meshes[mesh_id].base_index, which the loop just populated.
    if (m_accel) {
        (void)m_accel->export_gpu_bvh(result);
    }

    return result;
}

REGISTER_CLASS(Scene, "scene");

M_NAMESPACE_END

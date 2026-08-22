#include <components/bsdf.h>
#include <components/camera.h>
#include <components/integrator.h>
#include <components/scene.h>
#include <core/intersection.h>

M_NAMESPACE_BEGIN
class Path : public Integrator {
public:
    explicit Path(const PropertyList &properties) {
        m_max_depth = properties.get_integer("max_depth", 10);
        m_rr_depth  = properties.get_integer("rr_depth", 5);
        m_name      = properties.get_string("name", "path");
    }

    ~Path() override = default;

    void preprocess(const std::shared_ptr<Scene> &scene) override {}

    void construct() override {
#ifdef M_DEBUG
        std::cout << "Construct " << class_type_name(get_class_type()) << std::endl;
#endif
    }

    void add_child(const std::shared_ptr<Object> &child) override {}

    [[nodiscard]] Color3f li(const std::shared_ptr<Scene> &scene, std::shared_ptr<Sampler> sampler, const Ray3f &ray_,
                             bool &valid) const override {
        return li_impl(scene, sampler, ray_, valid, nullptr);
    }

    [[nodiscard]] Color3f li_aov(const std::shared_ptr<Scene> &scene, const std::shared_ptr<Sampler> &sampler,
                                 const Ray3f &ray_, bool &valid, AOVSample &aov) const override {
        return li_impl(scene, sampler, ray_, valid, &aov);
    }

private:
    /**
     * The path tracer proper. `aov` is null for the plain \ref li() entry
     * point, in which case not a single extra operation is performed and the
     * radiance estimate is bit-for-bit what it was before feature output
     * existed. When non-null, the primary hit's shading normal / depth and the
     * first usable surface albedo are recorded along the way for the denoisers
     * (see include/render/aov.h).
     */
    [[nodiscard]] Color3f li_impl(const std::shared_ptr<Scene> &scene, const std::shared_ptr<Sampler> &sampler,
                                  const Ray3f &ray_, bool &valid, AOVSample *aov) const {
        // Configure loop state
        Ray3f ray(ray_);
        Color3f throughput(1.0f);
        Color3f result(0.0f);
        float eta      = 1.0f;
        uint32_t depth = 0;
        bool valid_ray = false;

        // Variables caching information from the previous bounce
        Intersection3f prev_si;
        float prev_bsdf_pdf  = 1.0f;
        bool prev_bsdf_delta = true;

        // Path tracing loop
        for (uint32_t i = 0; i < m_max_depth && valid; i++) {
            // Ray intersect
            SurfaceIntersection3f its;
            bool is_intersect = scene->ray_intersect(ray, its, false);

            // ---------------------- Direct emission ----------------------

            bool hit_mesh_valid = its.mesh_id != M_INVALID_INDEX;
            // Resolves to a reference into Scene::m_meshes - no shared_ptr copy /
            // atomic refcount here, unlike the previous its.mesh-based version.
            const std::shared_ptr<Mesh> *hit_mesh = hit_mesh_valid ? &scene->get_mesh(its.mesh_id) : nullptr;

            // ------------------- Geometric AOV features -------------------
            // Recorded at the PRIMARY hit only: these describe the surface
            // this pixel actually shows, which is what the denoisers need in
            // order to know where the image's real edges are. A ray that
            // escapes into the environment leaves has_feature false, and the
            // denoisers then fall back to color-only weights for that pixel.
            if (aov && i == 0 && is_intersect && hit_mesh_valid) {
                aov->normal      = its.shading_frame.n;
                aov->depth       = its.t;
                aov->has_feature = true;
            }

            // If intersect an emitter
            if (is_intersect && (!hit_mesh_valid || (*hit_mesh)->is_emitter())) {
                std::shared_ptr<Emitter> hit_emitter = hit_mesh_valid ? (*hit_mesh)->get_emitter() : nullptr;
                if (!hit_mesh_valid) {
                    hit_emitter = scene->get_environment();
                }
                DirectionSample3f ds(its, prev_si, hit_emitter);
                float em_pdf = 0.0f;

                if (!prev_bsdf_delta) {
                    em_pdf = scene->pdf_emitter_direction(prev_si, ds, valid);
                }

                float mis_bsdf = mis_weight(prev_bsdf_pdf, em_pdf);

                result += throughput * ds.emitter->eval(its, valid) * mis_bsdf;
            }

            // Continue tracing the path at this point?
            bool active_next = depth + 1 < m_max_depth && is_intersect && hit_mesh_valid;

            if (!active_next) {
                break;
            }

            const std::shared_ptr<BSDF> &bsdf = (*hit_mesh)->get_bsdf();

            // ---------------------- Albedo AOV feature ----------------------
            // Take the albedo off the first surface that actually HAS a
            // meaningful one, i.e. the first non-delta (ESmooth) BSDF. Delta
            // surfaces (mirrors, glass) are transparent to this search for up
            // to M_AOV_MAX_DELTA_BOUNCES bounces, so that a texture seen in a
            // mirror still gets demodulated against its own albedo rather than
            // against a meaningless 1.
            if (aov && !aov->has_albedo && (bsdf->has_flag(ESmooth) || depth >= M_AOV_MAX_DELTA_BOUNCES)) {
                aov->albedo     = bsdf->albedo(its, active_next);
                aov->has_albedo = true;
            }

            // ---------------------- Emitter sampling ----------------------
            bool active_em = bsdf->has_flag(ESmooth);

            DirectionSample3f ds;
            Color3f em_weight;
            Vector3f wo;

            if (active_em) {
                std::tie(ds, em_weight) = scene->sample_emitter_direction(its, sampler->next2d(), true, active_em);
                active_em &= ds.pdf != 0.0f;
                wo = its.to_local(ds.d);
            }

            // ------ Evaluate BSDF * cos(theta) and sample direction -------
            float sample1   = sampler->next1d();
            Point2f sample2 = sampler->next2d();

            auto bsdf_val                   = bsdf->eval(its, wo, active_next);
            auto bsdf_pdf                   = bsdf->pdf(its, wo, active_next);
            auto [bsdf_sample, bsdf_weight] = bsdf->sample(its, sample1, sample2, active_next);

            // --------------- Emitter sampling contribution ----------------
            if (active_em) {
                float mis_em = ds.delta ? 1.0f : mis_weight(ds.pdf, bsdf_pdf);
                result += throughput * bsdf_val * em_weight * mis_em;
            }

            // ---------------------- BSDF sampling ----------------------
            ray = its.spawn_ray(its.to_world(bsdf_sample.wo));

            // ------ Update loop variables based on current interaction ------
            throughput *= bsdf_weight;
            eta *= bsdf_sample.eta;
            valid_ray |= valid && its.is_valid();

            // Information about the current vertex needed by the next iteration
            prev_si         = its;
            prev_bsdf_pdf   = bsdf_sample.pdf;
            prev_bsdf_delta = bsdf_sample.delta;

            // -------------------- Stopping criterion ---------------------
            depth += 1;
            float throughput_max = throughput.max_value();
            float rr_prob        = M_MIN(throughput_max * eta * eta, 0.95f);
            bool rr_active       = depth >= m_rr_depth;
            bool rr_continue     = sampler->next1d() < rr_prob;

            if (rr_active) {
                throughput *= 1.0f / rr_prob;
            }

            valid = (!rr_active || rr_continue) && throughput_max != 0.0f;
        }

        return result;
    }

public:
    [[nodiscard]] std::string to_string() const override {
        return std::string("Path[\n  max_depth=") + std::to_string(m_max_depth) + std::string("\n  rr_depth") +
               std::to_string(m_rr_depth) + std::string("\n]");
    }

    [[nodiscard]] int get_max_depth() const override { return m_max_depth; }
    [[nodiscard]] int get_rr_depth() const override { return m_rr_depth; }

private:
    int m_max_depth;
    int m_rr_depth;

    [[nodiscard]] static float mis_weight(float pdf_a, float pdf_b) {
        pdf_a *= pdf_a;
        pdf_b *= pdf_b;
        float w = pdf_a / (pdf_b + pdf_a);
        return pdf_a == 0 && pdf_b == 0 ? 0.0f : w;
    }
};

REGISTER_CLASS(Path, "path")

M_NAMESPACE_END

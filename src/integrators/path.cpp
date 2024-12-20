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
        std::cout << "Construct " << class_type_name(get_class_type())
                  << std::endl;
#endif
    }

    void add_child(const std::shared_ptr<Object> &child) override {}

    [[nodiscard]] Color3f li(const std::shared_ptr<Scene> &scene,
                             std::shared_ptr<Sampler> sampler,
                             const Ray3f &ray_, bool &valid) const override {
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
            bool is_intersect =
                scene->get_accel()->ray_intersect(ray, its, false);

            // ---------------------- Direct emission ----------------------

            // If intersect an emitter
            if (is_intersect && its.mesh->is_emitter()) {
                DirectionSample3f ds(its, prev_si);
                float em_pdf = 0.0f;

                if (!prev_bsdf_delta) {
                    em_pdf = scene->pdf_emitter_direction(prev_si, ds, valid);
                }

                float mis_bsdf = mis_weight(prev_bsdf_pdf, em_pdf);

                result += throughput * ds.emitter->eval(its, valid) * mis_bsdf;
            }

            // Continue tracing the path at this point?
            bool active_next = depth + 1 < m_max_depth && is_intersect;

            if (!active_next) {
                break;
            }

            std::shared_ptr<BSDF> bsdf = its.mesh->get_bsdf();

            // ---------------------- Emitter sampling ----------------------
            bool active_em = active_next;

            DirectionSample3f ds;
            Color3f em_weight;
            Vector3f wo;

            if (active_em) {
                std::tie(ds, em_weight) = scene->sample_emitter_direction(
                    its, sampler->next2d(), true, active_em);
                active_em &= ds.pdf != 0.0f;
                wo = its.to_local(ds.d);
            }

            // ------ Evaluate BSDF * cos(theta) and sample direction -------
            float sample1   = sampler->next1d();
            Point2f sample2 = sampler->next2d();

            auto bsdf_val = bsdf->eval(its, wo, active_next);
            auto bsdf_pdf = bsdf->pdf(its, wo, active_next);
            auto [bsdf_sample, bsdf_weight] =
                bsdf->sample(its, sample1, sample2, active_next);

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
            prev_bsdf_pdf   = bsdf_pdf;
            prev_bsdf_delta = bsdf_sample.delta;

            // -------------------- Stopping criterion ---------------------
            depth += 1;
            float throughput_max = throughput.max_value();
            float rr_prob        = M_MIN(throughput_max * eta * eta, 0.95f);
            bool rr_active       = depth >= m_rr_depth;
            bool rr_continue     = sampler->next1d() < rr_prob;

            valid = (!rr_active || rr_continue) && throughput_max != 0.0f;
        }

        return result;
    }

    [[nodiscard]] std::string to_string() const override {
        return std::string("Path[\n  max_depth=") +
               std::to_string(m_max_depth) + std::string("\n  rr_depth") +
               std::to_string(m_rr_depth) + std::string("\n]");
    }

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

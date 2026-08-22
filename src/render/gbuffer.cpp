#include <algorithm>
#include <components/bsdf.h>
#include <components/camera.h>
#include <components/mesh.h>
#include <core/intersection.h>
#include <render/aov.h>
#include <render/gbuffer.h>
#include <thread>
#include <vector>

M_NAMESPACE_BEGIN

void compute_gbuffer(const std::shared_ptr<Scene> &scene, FrameBufferSet &fbs, int thread_count) {
    const auto &camera   = scene->get_camera();
    const Vector2i size  = camera->get_output_size();
    const int rows = size.y(), cols = size.x();

    fbs.albedo = std::make_shared<Bitmap>(rows, cols, 3);
    fbs.normal = std::make_shared<Bitmap>(rows, cols, 3);
    fbs.depth  = std::make_shared<Bitmap>(rows, cols, 1);

    if (thread_count < 1)
        thread_count = static_cast<int>(std::thread::hardware_concurrency());
    thread_count = std::max(1, std::min(thread_count, rows));

    auto worker = [&](int y_begin, int y_end) {
        for (int y = y_begin; y < y_end; ++y) {
            for (int x = 0; x < cols; ++x) {
                // Pixel center, no jitter: the features are deterministic, so
                // there is nothing to average and one ray per pixel suffices.
                Point2f pixel_sample(static_cast<float>(x) + 0.5f, static_cast<float>(y) + 0.5f);

                Ray3f ray;
                bool valid = true;
                camera->sample_ray(ray, pixel_sample, Point2f(0.5f, 0.5f), valid);

                Color3f albedo(1.f);
                Normal3f normal(0.f);
                float depth = 0.f;

                if (valid) {
                    // Walk through delta (specular / refractive) surfaces just
                    // like Path::li_aov() does, so that a texture seen through
                    // glass or in a mirror is demodulated against its own
                    // albedo rather than against a meaningless 1.
                    Ray3f current(ray);
                    for (uint32_t bounce = 0; bounce <= M_AOV_MAX_DELTA_BOUNCES; ++bounce) {
                        SurfaceIntersection3f its;
                        if (!scene->ray_intersect(current, its, false) || its.mesh_id == M_INVALID_INDEX)
                            break;

                        if (bounce == 0) {
                            // Geometry always comes from the PRIMARY hit - it
                            // has to describe what this pixel actually shows,
                            // otherwise the edge-stopping weights would compare
                            // surfaces the viewer never sees.
                            normal = its.shading_frame.n;
                            depth  = its.t;
                        }

                        const std::shared_ptr<BSDF> &bsdf = scene->get_mesh(its.mesh_id)->get_bsdf();
                        if (!bsdf)
                            break;

                        if (bsdf->has_flag(ESmooth) || bounce == M_AOV_MAX_DELTA_BOUNCES) {
                            albedo = bsdf->albedo(its, true);
                            break;
                        }

                        // Delta surface: follow the (single, deterministic)
                        // specular direction. Sampling with fixed 0.5 variates
                        // is exact for a delta lobe, since it has only one
                        // outgoing direction to begin with.
                        auto [bs, weight] = bsdf->sample(its, 0.5f, Point2f(0.5f, 0.5f), true);
                        if (bs.pdf <= 0.f)
                            break;
                        current = its.spawn_ray(its.to_world(bs.wo));
                    }
                }

                (*fbs.albedo)(y, x, 0) = albedo(0);
                (*fbs.albedo)(y, x, 1) = albedo(1);
                (*fbs.albedo)(y, x, 2) = albedo(2);
                (*fbs.normal)(y, x, 0) = normal.x();
                (*fbs.normal)(y, x, 1) = normal.y();
                (*fbs.normal)(y, x, 2) = normal.z();
                (*fbs.depth)(y, x, 0)  = depth;
            }
        }
    };

    if (thread_count == 1) {
        worker(0, rows);
        return;
    }

    std::vector<std::thread> threads;
    threads.reserve(thread_count);
    int base = rows / thread_count, extra = rows % thread_count, y = 0;
    for (int t = 0; t < thread_count; ++t) {
        int count = base + (t < extra ? 1 : 0);
        threads.emplace_back([&worker, y, count] { worker(y, y + count); });
        y += count;
    }
    for (auto &th : threads)
        th.join();
}

M_NAMESPACE_END

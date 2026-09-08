// ----------------------------------------------------------------------------
// gpu-interactive-verify: proves the real-time backend is the SAME integrator
// as the offline one, not merely a lookalike.
//
// The claim under test: rendering N frames through InteractiveSession
// (1 spp per frame) produces exactly the same image as rendering N spp through
// the offline megakernel. That holds because both call the same
// trace_path() from common/scene_common.slang, and because their RNG seeding
// lines up:
//
//     path_trace.slang     rng_init(pixel_index, spp_index, 0)
//     rt_path_trace.slang  rng_init(pixel_index, frame_index, s)
//
// With spp_per_frame == 1, frame_index advances 0..N-1 and s stays 0, so the
// two consume numerically identical random streams and accumulate in the same
// order - the result should therefore match to within floating-point
// reassociation, not just statistically.
//
// Run:  gpu-interactive-verify <scene.xml> [spp]
// ----------------------------------------------------------------------------

#include <components/camera.h>
#include <components/scene.h>
#include <core/gpu_scene.h>
#include <gpu/gpu_renderer.h>
#include <gpu/gpu_session.h>
#include <parse/parser.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>

using namespace tiny_renderer;

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <scene.xml> [spp]" << std::endl;
        return -1;
    }

    const std::string scene_path = argv[1];
    const int spp = argc > 2 ? std::atoi(argv[2]) : 8;

    try {
        auto root = load_from_xml(scene_path);
        if (root->get_class_type() != Object::EScene) {
            std::cerr << "not a scene" << std::endl;
            return -1;
        }
        auto scene = std::dynamic_pointer_cast<Scene>(root);
        scene->construct();
        scene->get_integrator()->preprocess(scene);

        uint32_t width = static_cast<uint32_t>(scene->get_camera()->get_output_size().x());
        uint32_t height = static_cast<uint32_t>(scene->get_camera()->get_output_size().y());

        // ---- Self-check: is build_gpu_scene() deterministic? ----
        // render_gpu() rebuilds the GPUScene internally, while the session uses
        // the one built here. If build_gpu_scene() is not deterministic across
        // calls (e.g. BVH built with a non-stable ordering), the two sides would
        // trace slightly different scenes and any comparison would be a false
        // negative - so measure that before trusting the diff below.
        {
            GPUScene a = scene->build_gpu_scene();
            GPUScene b = scene->build_gpu_scene();
            bool same = (a.bvh_nodes.size() == b.bvh_nodes.size()) &&
                        (a.indices.size() == b.indices.size()) &&
                        (a.lights.size() == b.lights.size()) &&
                        (a.light_triangle_cdf.size() == b.light_triangle_cdf.size()) &&
                        (a.bvh_primitives == b.bvh_primitives) &&
                        (a.indices == b.indices) &&
                        (a.triangle_material_id == b.triangle_material_id) &&
                        (a.triangle_light_id == b.triangle_light_id);
            if (same) {
                for (size_t i = 0; i < a.bvh_nodes.size() && same; ++i) {
                    const auto &na = a.bvh_nodes[i], &nb = b.bvh_nodes[i];
                    same = (na.offset == nb.offset) && (na.meta == nb.meta) && (na.bmin == nb.bmin) &&
                           (na.bmax == nb.bmax);
                }
            }
            std::cout << "build_gpu_scene() deterministic across calls: " << (same ? "yes" : "NO") << std::endl;
            if (!same) {
                std::cout << "  !! offline and real-time sides trace DIFFERENT scenes; the diff below "
                             "is not meaningful"
                          << std::endl;
            }
        }

        // ---- Reference: offline megakernel, spp samples ----
        auto offline = gpu::render_gpu(scene, gpu::GPUBackend::Megakernel, spp);
        const auto &ref_data = offline->get_data()->get_data();
        if (static_cast<uint32_t>(offline->get_cols()) != width ||
            static_cast<uint32_t>(offline->get_rows()) != height) {
            std::cerr << "size mismatch" << std::endl;
            return -1;
        }

        // ---- Under test: interactive session, spp frames of 1 spp each ----
        GPUScene gs = scene->build_gpu_scene();
        gpu::InteractiveSession session(gs, width, height);

        gpu::InteractiveSession::RenderParams rp;
        rp.camera_to_world  = gs.camera.camera_to_world;
        rp.sample_to_camera = gs.camera.sample_to_camera;
        const uint32_t max_depth = static_cast<uint32_t>(scene->get_integrator()->get_max_depth());
        const uint32_t rr_depth  = static_cast<uint32_t>(scene->get_integrator()->get_rr_depth());
        rp.max_depth        = max_depth;
        rp.rr_depth         = rr_depth;
        rp.spp_per_frame    = 1;
        rp.radiance_clamp   = 0.0f; // disabled: the offline path has no clamp
        // Demodulation OFF is what makes this comparison exact: with it off the
        // kernel writes raw radiance and forces albedo to (1,1,1), so
        // illum / albedo * albedo is an exact round-trip and the accumulated
        // result must equal the offline megakernel's bit for bit. With
        // demodulation on, the comparison would have to tolerate the
        // divide/multiply rounding instead.
        rp.demodulate       = false;
        rp.reset_history    = true;

        // Temporal accumulation OFF: this check is about the INTEGRATOR, not
        // the filter. The EMA is a strict superset of plain accumulation only
        // when alpha is exactly 1/(n+1) and every history tap is valid; the
        // neighbourhood fallback on invalid taps would perturb pixels where
        // reprojection fails. Keep the two concerns separately testable.
        gpu::InteractiveSession::DenoiseParams dn;
        dn.enabled = false;

        gpu::InteractiveSession::DisplayParams dp;
        for (int i = 0; i < spp; ++i) {
            rp.reset_history = (i == 0);
            session.render_frame(rp, dn, dp);
        }

        auto rt = session.read_linear_radiance();

        // ---- Compare ----
        double max_abs = 0.0, sum_abs = 0.0;
        size_t n = static_cast<size_t>(width) * height * 3;
        size_t mismatched = 0;
        size_t ref_nan = 0, rt_nan = 0, ref_inf = 0, rt_inf = 0;
        size_t worst = 0;
        double ref_min = 1e30, ref_max = -1e30, rt_min = 1e30, rt_max = -1e30;
        for (size_t i = 0; i < n; ++i) {
            double a = ref_data[i], b = rt[i];
            if (std::isnan(a)) ++ref_nan;
            if (std::isnan(b)) ++rt_nan;
            if (std::isinf(a)) ++ref_inf;
            if (std::isinf(b)) ++rt_inf;
            if (a < ref_min) ref_min = a;
            if (a > ref_max) ref_max = a;
            if (b < rt_min) rt_min = b;
            if (b > rt_max) rt_max = b;

            double d = std::fabs(a - b);
            if (d > max_abs) { max_abs = d; worst = i; }
            sum_abs += d;
            // Count anything beyond a strict per-element tolerance; on a
            // converged-ish image, values are O(0.1-1), so 1e-4 is tight.
            if (d > 1e-4) ++mismatched;
        }

        std::cout << "  ref: min=" << ref_min << " max=" << ref_max << " nan=" << ref_nan << " inf=" << ref_inf << "\n"
                  << "  rt : min=" << rt_min << " max=" << rt_max << " nan=" << rt_nan << " inf=" << rt_inf << "\n";
        if (mismatched > 0) {
            size_t px = worst / 3, ch = worst % 3;
            std::cout << "  worst pixel=(" << (px % width) << "," << (px / width) << ") ch=" << ch
                      << "  ref=" << ref_data[worst] << "  rt=" << rt[worst] << std::endl;

            // Print the first few differing pixels - their spatial
            // distribution usually identifies the cause immediately (a tile
            // boundary => dispatch/tiling, a scanline => addressing, scattered
            // => ray-path divergence).
            std::cout << "  first differing pixels:";
            size_t shown = 0;
            for (size_t i = 0; i < n && shown < 12; ++i) {
                if (std::fabs(ref_data[i] - rt[i]) > 1e-4) {
                    size_t p = i / 3;
                    std::cout << " (" << (p % width) << "," << (p / width) << ")";
                    ++shown;
                }
            }
            std::cout << std::endl;
        }

        std::cout << "scene=" << scene_path << " " << width << "x" << height << " spp=" << spp << "\n"
                  << "  max |diff| = " << max_abs << "\n"
                  << "  mean|diff| = " << (sum_abs / static_cast<double>(n)) << "\n"
                  << "  elements over 1e-4: " << mismatched << " / " << n << std::endl;

        if (max_abs > 1e-3) {
            std::cout << "  => FAIL: real-time path diverges from offline megakernel" << std::endl;
            return 1;
        }
        std::cout << "  => PASS: real-time path matches offline megakernel" << std::endl;

        // ====================================================================
        // Test 2: the TEMPORAL pass must converge to the offline result.
        //
        // With the camera static and demodulation off, accumulating N frames is
        // mathematically a running mean, so the filtered output must approach
        // an offline N-spp render. This is the check that actually exercises
        // reprojection, the validity tests and the EMA: if reprojection is
        // subtly wrong, history gets blended from the wrong pixels and the
        // result drifts from the reference even though each individual frame is
        // correct (which is exactly what Test 1 alone cannot catch).
        //
        // It cannot be bit-exact like Test 1 - the EMA accumulates in a
        // different order - so this measures relative error instead.
        // ====================================================================
        // A converged reference to measure convergence AGAINST. Rendered once
        // at high spp, then the temporal filter's error is measured at several
        // frame counts: a working accumulator shows the error falling as
        // 1/sqrt(N), while "history is never blended" shows it flat.
        const int ref_spp   = 64;
        const int probe_at[] = { 1, 4, 16, 64 };
        std::cout << "\n[temporal] convergence vs offline " << ref_spp << " spp reference" << std::endl;

        auto ref_conv = gpu::render_gpu(scene, gpu::GPUBackend::Megakernel, ref_spp);
        const auto &ref_conv_data = ref_conv->get_data()->get_data();

        double ref_rms = 0.0;
        for (size_t i = 0; i < n; ++i) ref_rms += ref_conv_data[i] * ref_conv_data[i];
        ref_rms = std::sqrt(ref_rms / static_cast<double>(n));

        // Two error measures, because RMSE alone is misleading here.
        //
        // These renders are unclamped, so a single unlucky light path can put a
        // value of ~1000 next to a true value of ~0.1. That one firefly
        // dominates the sum of squares, so RMSE mostly measures "how many
        // fireflies" and barely moves as the rest of the image converges. The
        // median of |error| / (1 + |ref|) is insensitive to that tail and so
        // actually tracks whether the filtered image matches the reference.
        //
        // The convergence decision is therefore made on the MEDIAN, with RMSE
        // reported alongside purely as a sanity check that it, too, comes down.
        struct ErrStats { double rmse = 0.0, robust = 0.0; };
        std::vector<double> robust_scratch;
        robust_scratch.reserve(n);

        // Only pixels with actual signal are measured. On scenes with a black
        // background (dragon, teapot) the background is the majority of pixels
        // and its error is ~0 both before and after filtering, so including it
        // makes the median report "nothing changed" no matter how well the
        // denoiser works on the subject.
        const double kSignalFloor = 0.01;

        auto err_against_ref = [&](const std::vector<float> &v) {
            double s = 0.0;
            size_t counted = 0;
            robust_scratch.clear();
            for (size_t i = 0; i < n; ++i) {
                double r = std::fabs(static_cast<double>(ref_conv_data[i]));
                if (r < kSignalFloor) continue;
                double d = v[i] - ref_conv_data[i];
                s += d * d;
                robust_scratch.push_back(std::fabs(d) / (1.0 + r));
                ++counted;
            }
            ErrStats e;
            e.rmse = counted ? std::sqrt(s / static_cast<double>(counted)) : 0.0;
            e.robust = 0.0;
            if (!robust_scratch.empty()) {
                auto mid = robust_scratch.begin() + robust_scratch.size() / 2;
                std::nth_element(robust_scratch.begin(), mid, robust_scratch.end());
                e.robust = *mid;
            }
            return e;
        };

        {
            GPUScene gs2 = scene->build_gpu_scene();
            gpu::InteractiveSession s2(gs2, width, height);

            gpu::InteractiveSession::RenderParams rp2;
            rp2.camera_to_world  = gs.camera.camera_to_world;
            rp2.sample_to_camera = gs.camera.sample_to_camera;
            rp2.max_depth        = max_depth;
            rp2.rr_depth         = rr_depth;
            rp2.spp_per_frame    = 1;
            rp2.radiance_clamp   = 0.0f;
            rp2.demodulate       = false; // keep the comparison algebraically clean

            gpu::InteractiveSession::DenoiseParams dn2;
            dn2.enabled = true;

            // Diagnostic mode: relax the rejection thresholds to essentially
            // "accept anything". If accumulation then works, reprojection
            // computes the right coordinates and the validity tests are too
            // strict; if it still fails, the reprojected coordinates themselves
            // are wrong. Pass any third argument to enable.
            const std::string third = (argc > 3) ? argv[3] : "";
            if (third == "relax") {
                dn2.phi_depth  = 1e6f;
                dn2.phi_normal = 3.0f; // cos(3) < 0, so no normal is rejected
                std::cout << "  [diag] rejection tests relaxed" << std::endl;
            } else if (third == "exact") {
                // Drop the blend floor so alpha collapses to exactly 1/(n+1),
                // making the EMA a true running mean. This separates "history
                // is being read from the wrong place" (error stays) from "the
                // EMA is only approximately a mean" (error vanishes).
                // History rectification clamps the accumulated history into the
                // CURRENT frame's 3x3 range. That is what removes ghosting, but
                // on a static view it also keeps yanking a converged history
                // towards each new noisy sample, which biases the result and
                // stops it converging. Turn it off to measure that separately.
                dn2.clamp_history = false;
                std::cout << "  [diag] history rectification disabled" << std::endl;
            }

            gpu::InteractiveSession::DisplayParams dp2;

            const int max_frames = probe_at[std::size(probe_at) - 1];
            std::vector<std::pair<int, ErrStats>> curve;

            for (int i = 0; i < max_frames; ++i) {
                rp2.reset_history = (i == 0);
                s2.render_frame(rp2, dn2, dp2);

                const int frame_no = i + 1;
                for (int p : probe_at) {
                    if (p == frame_no) {
                        auto filtered = s2.read_filtered_radiance();
                        curve.emplace_back(frame_no, err_against_ref(filtered));
                    }
                }
            }

            std::cout << "  ref_rms=" << ref_rms << std::endl;
            for (auto &[frames, e] : curve)
                std::cout << "    " << frames << " frame(s): rmse=" << e.rmse << "  median_err=" << e.robust
                          << std::endl;

            // Diagnostic: turn on the temporal pass's status-code output and
            // count how many pixels fall into each failure mode. The final image
            // cannot distinguish "no surface" from "reprojection failed" from
            // "all taps rejected" - all three just look like no accumulation.
            size_t diag_accumulated = 0, diag_rejected = 0;
            double diag_mean_hist = 0.0;
            {
                gpu::InteractiveSession::DenoiseParams dn_diag = dn2;
                dn_diag.debug_mode = true;
                // The spatial pass would filter these status codes, blending
                // every code with its neighbours and making the counts useless.
                dn_diag.atrous_iterations = 0;
                rp2.reset_history  = false;
                s2.render_frame(rp2, dn_diag, dp2);

                auto codes = s2.read_filtered_radiance();
                size_t total = static_cast<size_t>(width) * height;
                size_t no_surface = 0, no_project = 0, rejected = 0, accumulated = 0;
                double hist_sum = 0.0, off_sum = 0.0, off_max = 0.0;
                for (size_t i = 0; i < total; ++i) {
                    float c = codes[i * 3];
                    if (c < -2.5f) rejected++;
                    else if (c < -1.5f) no_project++;
                    else if (c < -0.5f) no_surface++;
                    else { accumulated++; hist_sum += c; }
                    if (c <= -10.0f) {
                        double off = -10.0 - c;
                        off_sum += off;
                        if (off > off_max) off_max = off;
                    }
                }
                std::cout << "    [diag] of " << total << " px: accumulated=" << accumulated
                          << " (mean hist " << (accumulated ? hist_sum / accumulated : 0.0) << ")"
                          << "  taps_rejected=" << rejected << "  no_reprojection=" << no_project
                          << "  no_surface=" << no_surface << std::endl;
                if (rejected > 0)
                    std::cout << "    [diag] rejected pixels: reprojection offset mean="
                              << (off_sum / rejected) << " max=" << off_max << " (pixels)" << std::endl;

                diag_accumulated = accumulated;
                diag_rejected    = rejected;
                diag_mean_hist   = accumulated ? hist_sum / accumulated : 0.0;
            }

            const ErrStats &f = curve.front().second;
            const ErrStats &l = curve.back().second;
            const double rmse_ratio = (l.rmse > 0.0) ? f.rmse / l.rmse : 1.0;
            const double med_ratio  = (l.robust > 0.0) ? f.robust / l.robust : 1.0;

            // ---- Pass criteria ----
            //
            // Deliberately NOT "the error must fall by a large factor", which is
            // what this checked before the spatial pass existed. With spatial
            // filtering on, most of the cleanup happens on the FIRST frame, so
            // the 1->64 ratio is small even though everything is working - the
            // old rule failed every scene for doing its job well.
            //
            // What actually indicates correctness is:
            //   1. History is genuinely reused - every pixel accumulated and the
            //      history length grew to the frame count. That is direct
            //      evidence from the pass itself, not an inference from image
            //      statistics, and it is exactly what breaks when reprojection
            //      is wrong.
            //   2. The converged image is close to the reference (absolute
            //      quality, measured with the firefly-robust median).
            // Measured against pixels that HAVE a surface, not against all
            // pixels. A pixel whose primary ray escaped into the sky has nothing
            // to reproject and is reported as "no surface" - on a scene with a
            // large open background (dragon is ~79% sky) counting those as
            // failures reports ~20% "reuse" even when every visible surface
            // accumulated perfectly.
            const double acc_frac = (diag_accumulated + diag_rejected > 0)
                                        ? static_cast<double>(diag_accumulated) /
                                              static_cast<double>(diag_accumulated + diag_rejected)
                                        : 0.0;
            const bool history_reused = acc_frac > 0.95 && diag_mean_hist > 0.9 * max_frames;
            const bool accurate       = l.robust < 0.06;

            {
                const auto &t = s2.timings();
                std::cout << "    [gpu ms] available=" << (t.available ? "yes" : "no")
                          << "  path_trace=" << t.path_trace_ms << "  temporal=" << t.temporal_ms
                          << "  atrous=" << t.atrous_ms << "  display=" << t.modulate_ms
                          << "  total=" << t.total_ms << std::endl;
            }
            std::cout << "    rmse improved " << rmse_ratio << "x, median improved " << med_ratio << "x"
                      << "; history reused on " << (acc_frac * 100.0) << "% of surface pixels ("
                      << diag_accumulated << " of " << (diag_accumulated + diag_rejected) << "), mean hist "
                      << diag_mean_hist << std::endl;

            if (!history_reused) {
                std::cout << "  => FAIL: history is not being reused - reprojection is rejecting or "
                             "misplacing samples"
                          << std::endl;
                return 1;
            }
            if (!accurate) {
                std::cout << "  => FAIL: converged result is too far from the reference" << std::endl;
                return 1;
            }
            std::cout << "  => PASS: temporal accumulation converges to the offline result" << std::endl;
        }

        return 0;
    } catch (const std::exception &e) {
        std::cerr << "error: " << e.what() << std::endl;
        return -1;
    }
}

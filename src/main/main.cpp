#include <atomic>
#include <chrono>
#include <components/bitmap.h>
#include <components/camera.h>
#include <components/denoiser.h>
#include <components/scene.h>
#include <core/timer.h>
#include <filesystem/resolver.h>
#include <iostream>
#include <parse/parser.h>
#include <render/block.h>
#include <render/gbuffer.h>
#include <render/progress.h>
#include <thread>
#include <vector>

#ifdef M_ENABLE_GPU_BACKEND
#include <gpu/gpu_renderer.h>
#endif

using namespace tiny_renderer;

static void render_block(const std::shared_ptr<Scene> &scene, const std::shared_ptr<Sampler> &sampler,
                         ImageBlock &block) {
    const auto &camera     = scene->get_camera();
    const auto &integrator = scene->get_integrator();

    Point2i offset = block.get_offset();
    Vector2i size  = block.get_size();

    // Only pay for the AOV code path when a denoiser actually asked for it -
    // see ImageBlock's `with_aov` constructor flag.
    const bool want_aov = block.has_aov();

    /* Clear the block contents */
    block.clear();

    /* For each pixel and pixel sample sample */
    for (int y = 0; y < size.y(); ++y) {
        for (int x = 0; x < size.x(); ++x) {
            for (uint32_t i = 0; i < sampler->get_sample_count(); ++i) {
                Point2f pixel_sample =
                    Point2f(static_cast<float>(x + offset.x()), static_cast<float>(y + offset.y())) + sampler->next2d();
                Point2f aperture_sample = sampler->next2d();

                /* Sample a ray from the camera */
                Ray3f ray;
                bool valid    = true;
                Color3f value = camera->sample_ray(ray, pixel_sample, aperture_sample, valid);

                if (want_aov) {
                    /* Compute the incident radiance along with the primary
                       hit's denoising features (albedo / normal / depth) */
                    AOVSample aov;
                    value *= integrator->li_aov(scene, sampler, ray, valid, aov);

                    /* Store in the image block, routing this sample into
                       half-buffer A or B by its parity so that the two halves
                       remain statistically independent */
                    block.put(pixel_sample, value, aov, i);
                } else {
                    /* Compute the incident radiance */
                    value *= integrator->li(scene, sampler, ray, valid);

                    /* Store in the image block */
                    block.put(pixel_sample, value);
                }
            }
        }
    }
}

// Resolves `filename`'s output path (same convention used by both the CPU
// and GPU render paths: strip the scene file's extension, save as PNG next
// to it) and saves `bitmap` there. `tonemap` is forwarded to
// Bitmap::save_png() unchanged - see include/components/bitmap.h's
// ToneMapMode doc comment for what it does and why it defaults to None.
//
// `suffix` is appended to the stem, which is how the denoising path writes the
// unfiltered render alongside the filtered one (`scene_noisy.png` next to
// `scene.png`) for side-by-side comparison.
static void save_output(const std::shared_ptr<Bitmap> &bitmap, const std::string &filename,
                        ToneMapMode tonemap = ToneMapMode::None, const std::string &suffix = "") {
    filesystem::path file   = get_file_resolver()->resolve(filesystem::path(filename));
    std::string output_name = file.make_absolute().str();
    size_t lastdot          = output_name.find_last_of('.');

    if (lastdot != std::string::npos)
        output_name.erase(lastdot, std::string::npos);

    /* Save using the OpenEXR format */
    // bitmap->save_exr(output_name + ".exr");

    /* Save tone mapped (sRGB) output using the PNG format */
    bitmap->save_png(output_name + suffix + ".png", tonemap);
}

/**
 * \brief Run the scene's denoiser (if any) on a finished render and save both
 * the raw and the filtered result.
 *
 * Denoising deliberately happens here rather than inside Bitmap::save_png():
 * it must operate on LINEAR radiance, before any tonemapping or sRGB encoding,
 * which also means this single function serves both the CPU and the GPU render
 * path.
 *
 * Without a denoiser this is exactly equivalent to the previous
 * `save_output(result.to_bitmap(), ...)` call - same file name, same pixels.
 */
static void save_render(const FrameBufferSet &fbs, const std::shared_ptr<Scene> &scene, const std::string &filename,
                        ToneMapMode tonemap, bool dump_aov) {
    const auto &denoiser = scene->get_denoiser();

    if (dump_aov) {
        // Diagnostic dump of the denoiser's inputs. Useful when tuning the
        // sigmas: a blurry normal buffer or an all-white albedo buffer
        // immediately explains a disappointing result.
        if (fbs.albedo)
            save_output(fbs.albedo, filename, ToneMapMode::None, "_albedo");
        if (fbs.normal) {
            // Normals live in [-1, 1]; remap to [0, 1] so the PNG is viewable.
            const int h = fbs.normal->get_rows(), w = fbs.normal->get_cols();
            auto vis    = std::make_shared<Bitmap>(h, w, 3);
            for (int y = 0; y < h; ++y)
                for (int x = 0; x < w; ++x)
                    for (int c = 0; c < 3; ++c)
                        (*vis)(y, x, c) = 0.5f * ((*fbs.normal)(y, x, c) + 1.0f);
            save_output(vis, filename, ToneMapMode::None, "_normal");
        }
        if (fbs.variance)
            save_output(fbs.variance, filename, ToneMapMode::None, "_variance");
    }

    if (!denoiser) {
        save_output(fbs.color, filename, tonemap);
        return;
    }

    // Keep the unfiltered image around under a `_noisy` suffix: denoising is
    // lossy by nature and being able to A/B it against the raw render is
    // essential when judging whether the parameters are right.
    save_output(fbs.color, filename, tonemap, "_noisy");

    std::cout << "Denoising (" << denoiser->get_name() << ") .. " << std::flush;
    Timer timer;
    auto denoised = denoiser->denoise(DenoiseInput::from(fbs));
    std::cout << "done. (took " << timer.elapsed_string() << ")" << std::endl;

    save_output(denoised, filename, tonemap);
}

#ifdef M_ENABLE_GPU_BACKEND
// Renders `scene` on the GPU (see include/gpu/gpu_renderer.h) instead of
// spinning up the CPU's threaded ImageBlock loop below. Reuses the exact
// same Scene::construct()/preprocess()/save_png() calls as the CPU path, so
// the only thing that differs between `--gpu` and the default is which
// function actually produces the linear-radiance Bitmap.
static void render_on_gpu(const std::shared_ptr<Scene> &scene, const std::string &filename,
                          gpu::GPUBackend backend, ToneMapMode tonemap, bool show_progress, bool dump_aov) {
    scene->construct();

    scene->get_integrator()->preprocess(scene);

    const char *backend_name = backend == gpu::GPUBackend::Megakernel ? "megakernel" : "wavefront";
    std::cout << "Rendering on GPU (" << backend_name << ") .. " << std::endl;
    Timer timer;

    Vector2i output_size = scene->get_camera()->get_output_size();
    ProgressReporter progress(show_progress, output_size.x(), output_size.y(),
                              std::string("GPU ") + backend_name);

    // Called once per completed spp dispatch (see render_megakernel()'s/
    // render_wavefront()'s spp loops in src/gpu/gpu_renderer.cpp): this is
    // the natural "noise decreasing over time" progress unit for the GPU
    // path, since every pass accumulates over the SAME full image.
    auto on_progress = [&](int done, int total, const float *rgb, uint32_t w, uint32_t h) {
        progress.update(done, total, rgb, static_cast<int>(w), static_cast<int>(h));
    };

    // Ask for the half-buffers only when something will consume them, so a
    // plain --gpu render keeps its previous behavior and cost exactly.
    const bool want_aov = scene->get_denoiser() != nullptr || dump_aov;

    auto callback = show_progress ? gpu::ProgressCallback(on_progress) : nullptr;

    // The GPU backends do not export the geometric feature buffers, but those
    // are deterministic G-buffer quantities that do not need Monte Carlo
    // sampling at all - so rather than plumbing extra outputs through every
    // shader, they are recomputed on the CPU with one primary ray per pixel
    // (see compute_gbuffer). Combined with the real half-buffer variance above,
    // the denoisers then run with the SAME information they get on the CPU
    // path, instead of having to guess noise levels from local contrast and
    // blurring texture detail away in the process.
    FrameBufferSet fbs = want_aov ? gpu::render_gpu_with_aov(scene, backend, 0, callback)
                                  : FrameBufferSet{ gpu::render_gpu(scene, backend, 0, callback) };
    progress.finish();

    std::cout << "done. (took " << timer.elapsed_string() << ")" << std::endl;

    if (want_aov) {
        std::cout << "Computing feature buffers .. " << std::flush;
        Timer gbuffer_timer;
        compute_gbuffer(scene, fbs);
        std::cout << "done. (took " << gbuffer_timer.elapsed_string() << ")" << std::endl;
    }

    save_render(fbs, scene, filename, tonemap, dump_aov);
}
#endif

static void render(const std::shared_ptr<Scene> &scene, const std::string &filename, int thread_count,
                   ToneMapMode tonemap, bool show_progress, bool dump_aov) {
    scene->construct();

#ifdef M_DEBUG
    // Smoke-test the Stage 0 GPU flattening export layer: build it and print
    // a summary. This is read-only and does not affect the CPU render below;
    // it exists purely to catch regressions in the export code at run time.
    {
        GPUScene gpu_scene = scene->build_gpu_scene();
        std::cout << "GPUScene summary:" << std::endl;
        std::cout << "  vertices  = " << gpu_scene.vertex_positions.size() << std::endl;
        std::cout << "  triangles = " << gpu_scene.indices.size() / 3 << std::endl;
        std::cout << "  meshes    = " << gpu_scene.meshes.size() << std::endl;
        std::cout << "  materials = " << gpu_scene.materials.size() << std::endl;
        std::cout << "  textures  = " << gpu_scene.textures.size() << std::endl;
        std::cout << "  lights    = " << gpu_scene.lights.size() << std::endl;
        std::cout << "  env_light = " << (gpu_scene.environment_light_id != M_INVALID_INDEX ? "yes" : "no")
                  << std::endl;
        std::cout << "  bvh_nodes = " << gpu_scene.bvh_nodes.size() << std::endl;
        std::cout << "  bvh_prims = " << gpu_scene.bvh_primitives.size() << std::endl;
        std::cout << "  camera    = " << gpu_scene.camera.width << "x" << gpu_scene.camera.height << std::endl;
    }
#endif

    auto camera          = scene->get_camera();
    Vector2i output_size = camera->get_output_size();
    scene->get_integrator()->preprocess(scene);

    /* Create a block generator (i.e. a work scheduler) */
    BlockGenerator block_generator(output_size, M_BLOCK_SIZE);

    // Only accumulate the auxiliary denoising buffers when something is
    // actually going to consume them: a scene with no <denoiser> (and no
    // --denoise / --dump-aov flag) renders through the exact same code path,
    // with the same memory footprint and the same output, as before denoising
    // existed.
    const bool want_aov = scene->get_denoiser() != nullptr || dump_aov;

    /* Allocate memory for the entire output image and clear it */
    ImageBlock result(output_size, camera->get_reconstruction_filter(), want_aov);
    result.clear();

    /* Create threads */
    std::vector<std::thread> threads;
    int total_blocks      = block_generator.get_block_count(); // must be read BEFORE next() below decrements it
    int blocks_per_thread = total_blocks / thread_count;
    int remaining_blocks  = total_blocks % thread_count;

    std::cout << "Rendering .. " << std::endl;
    Timer timer;

    // Block-completion progress: each finished block is one natural unit of
    // "the image gradually filling in" (see include/render/progress.h) -
    // unlike the GPU path's per-spp noise-reduction progress, the CPU path
    // renders every block to its full sample count in one shot, so blocks
    // completed/total is the meaningful unit here instead.
    std::atomic<int> blocks_done{ 0 };
    ProgressReporter progress(show_progress, output_size.x(), output_size.y(), "CPU");

    threads.reserve(thread_count);
    for (int i = 0; i < thread_count; ++i) {
        threads.emplace_back([&, i] {
            // Each thread works on its own block range
            int start_block = i * blocks_per_thread;
            int end_block   = i == thread_count - 1 ? start_block + blocks_per_thread + remaining_blocks
                                                    : start_block + blocks_per_thread;

            ImageBlock block(Vector2i(M_BLOCK_SIZE), camera->get_reconstruction_filter(), want_aov);
            std::shared_ptr sampler(scene->get_sampler()->clone());

            for (int j = start_block; j < end_block; ++j) {
                block_generator.next(block);
                sampler->prepare(block); // The sampler of different block has different seeds.

                /* Render all contained pixels */
                render_block(scene, sampler, block);

                /* The image block has been processed. Now add it to
                   the "big" block that represents the entire image */
                result.put(block);
                blocks_done.fetch_add(1, std::memory_order_relaxed);
            }
        });
    }

    if (show_progress) {
        // Poll from the main thread while the workers above run, taking a
        // (mutex-protected, see ImageBlock::to_bitmap()) snapshot of the
        // full image so far for the console bar / live preview window.
        int done;
        while ((done = blocks_done.load(std::memory_order_relaxed)) < total_blocks) {
            auto snapshot = result.to_bitmap();
            progress.update(done, total_blocks, snapshot->get_data()->get_data(), output_size.x(), output_size.y());
            std::this_thread::sleep_for(std::chrono::milliseconds(150));
        }
    }

    // Join all threads
    for (auto &t : threads) {
        t.join();
    }

    if (show_progress) {
        auto snapshot = result.to_bitmap();
        progress.update(total_blocks, total_blocks, snapshot->get_data()->get_data(), output_size.x(),
                        output_size.y());
        progress.finish();
    }

    std::cout << "done. (took " << timer.elapsed_string() << ")" << std::endl;

    /* Now turn the rendered image block into
       a properly normalized bitmap */
    save_render(result.to_framebuffers(), scene, filename, tonemap, dump_aov);
}

int main(int argc, char **argv) {
    int thread_count = 1;
    bool use_gpu     = false;
    bool show_progress = false;
    std::string gpu_backend_name = "megakernel";
    ToneMapMode tonemap = ToneMapMode::None;
    // Empty = leave whatever the scene XML configured (possibly nothing).
    std::string denoiser_name;
    bool dump_aov = false;

    if (argc < 2) {
        std::cerr << "Syntax: " << argv[0]
                  << " <scene.xml> [--no-gui] [--threads N] [--gpu[=megakernel|wavefront]] [--tonemap=none|aces] "
                     "[--progress] [--denoise[=outlier|atrous|nlm]] [--dump-aov]"
                  << std::endl;
        return -1;
    }

    std::string scene_name;

    for (int i = 1; i < argc; ++i) {
        std::string token(argv[i]);
        if (token == "-t" || token == "--threads") {
            if (i + 1 >= argc) {
                std::cerr << "\"--threads\" argument expects a positive integer following it." << std::endl;
                return -1;
            }
            char *end_ptr;
            thread_count = strtol(argv[i + 1], &end_ptr, 10);
            i++;
            if (thread_count <= 0) {
                std::cerr << "\"--threads\" argument expects a positive integer following it." << std::endl;
                return -1;
            }

            continue;
        }
        if (token == "--gpu") {
            use_gpu = true;
            continue;
        }
        if (token.rfind("--gpu=", 0) == 0) {
            use_gpu          = true;
            gpu_backend_name = token.substr(std::string("--gpu=").size());
            continue;
        }
        if (token == "--progress") {
            show_progress = true;
            continue;
        }
        // Install a denoiser with its default parameters, overriding the
        // scene's <denoiser> block if it has one. Lets the same scene file be
        // rendered with and without denoising (and with different denoisers)
        // for comparison, without editing the XML.
        if (token == "--denoise") {
            denoiser_name = "atrous";
            continue;
        }
        if (token.rfind("--denoise=", 0) == 0) {
            denoiser_name = token.substr(std::string("--denoise=").size());
            continue;
        }
        // Write the denoiser input buffers out as PNGs next to the render.
        // Primarily a parameter-tuning aid - see save_render().
        if (token == "--dump-aov") {
            dump_aov = true;
            continue;
        }
        if (token.rfind("--tonemap=", 0) == 0) {
            std::string mode = token.substr(std::string("--tonemap=").size());
            if (mode == "none") {
                tonemap = ToneMapMode::None;
            } else if (mode == "aces") {
                tonemap = ToneMapMode::ACES;
            } else {
                std::cerr << "Unknown --tonemap mode \"" << mode << "\" (expected \"none\" or \"aces\")" << std::endl;
                return -1;
            }
            continue;
        }
#ifdef M_DEBUG
        thread_count = 1;
#endif

        filesystem::path path(argv[i]);

        if (path.extension() == "xml") {
            scene_name = argv[i];

            /* Add the parent directory of the scene file to the
               file resolver. That way, the XML file can reference
               resources (OBJ files, textures) using relative paths */
            get_file_resolver()->prepend(path.parent_path());
        }
    }

    if (scene_name.empty()) {
        std::cerr << "Please provide xml" << std::endl;
        return -1;
    }

    try {
        auto root(load_from_xml(scene_name));
        if (root->get_class_type() == Object::EScene) {
            auto scene = std::dynamic_pointer_cast<Scene>(root);

            if (!denoiser_name.empty()) {
                // Built through the same ObjectFactory the XML parser uses, so
                // --denoise=<type> accepts exactly the set of types registered
                // by src/denoisers/*.cpp and nothing else.
                auto obj = ObjectFactory::create_instance(denoiser_name, PropertyList());
                if (obj->get_class_type() != Object::EDenoiser) {
                    std::cerr << "\"" << denoiser_name << "\" is not a denoiser." << std::endl;
                    return -1;
                }
                scene->set_denoiser(std::dynamic_pointer_cast<Denoiser>(obj));
            }

            if (use_gpu) {
#ifdef M_ENABLE_GPU_BACKEND
                gpu::GPUBackend backend;
                if (gpu_backend_name == "megakernel") {
                    backend = gpu::GPUBackend::Megakernel;
                } else if (gpu_backend_name == "wavefront") {
                    backend = gpu::GPUBackend::Wavefront;
                } else {
                    std::cerr << "Unknown --gpu backend \"" << gpu_backend_name
                              << "\" (expected \"megakernel\" or \"wavefront\")" << std::endl;
                    return -1;
                }
                render_on_gpu(scene, scene_name, backend, tonemap, show_progress, dump_aov);
#else
                std::cerr << "This build was compiled with M_ENABLE_GPU=OFF; --gpu is unavailable." << std::endl;
                return -1;
#endif
            } else {
                render(scene, scene_name, thread_count, tonemap, show_progress, dump_aov);
            }
        }
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    return 0;
}
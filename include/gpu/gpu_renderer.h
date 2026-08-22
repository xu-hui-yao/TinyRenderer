#pragma once

// ============================================================================
// GPU backend entry point for the MAIN tiny-renderer executable.
//
// Unlike the standalone gpu-path-trace/gpu-wavefront debug tools (which each
// parse their own scene file and write their own tonemapped PNG for
// side-by-side comparison against the CPU renderer), this function is meant
// to be called from src/main/main.cpp on an already-constructed Scene, and
// it returns a Bitmap holding the averaged LINEAR radiance per pixel - not
// an already-tonemapped image. The caller then uses the EXACT SAME
// Bitmap::save_png()/save_exr() code path as the CPU renderer, so there is
// only ONE tonemap implementation in the whole codebase (see the bathroom2
// wavefront investigation: a second, divergent tonemap implementation is
// exactly how that bug was introduced).
// ============================================================================

#include <components/bitmap.h>
#include <components/scene.h>
#include <cstdint>
#include <functional>
#include <memory>
#include <render/framebuffer.h>

namespace tiny_renderer::gpu {

enum class GPUBackend {
    Megakernel, // src/gpu/shaders/path_trace.slang: one dispatch per spp, full bounce loop in-shader.
    Wavefront,  // src/gpu/shaders/wavefront_*.slang: 6 kernels/bounce, compacted work queues.
};

// Invoked after every completed spp dispatch (see src/gpu/gpu_renderer.cpp's
// render_megakernel()/render_wavefront() spp loops) with `done`/`total` spp
// passes and a row-major height*width*3 LINEAR-radiance snapshot of the
// accumulation buffer so far (`rgb`, sized `width`*`height`*3 floats - valid
// only for the duration of the call). Used by include/render/progress.h to
// drive --progress's console bar and optional live preview window; pass
// `nullptr` (the default) to skip progress reporting entirely, which also
// skips the (otherwise per-spp) accumulation-buffer readback and Bitmap
// allocation this requires.
using ProgressCallback = std::function<void(int done, int total, const float *rgb, uint32_t width, uint32_t height)>;

// Renders `scene` on the GPU using `backend`. `scene->construct()` must
// already have been called. `spp` samples per pixel are rendered; if
// `spp <= 0`, `scene->get_sampler()->get_sample_count()` is used instead,
// matching the CPU renderer's convention of reading sample count from the
// scene's <sampler> element. Throws std::runtime_error on any Vulkan/setup
// failure (e.g. no GPU available, or the scene's Accel implementation not
// exporting a flat BVH - only BVHAccel currently does).
[[nodiscard]] std::shared_ptr<Bitmap> render_gpu(const std::shared_ptr<Scene> &scene, GPUBackend backend,
                                                  int spp = 0, const ProgressCallback &progress = nullptr);

// As render_gpu(), but additionally fills in the denoiser input buffers that
// can be obtained without any shader modification (see FrameBufferSet):
//
//   - `color`            the same averaged linear radiance render_gpu() returns
//   - `color_a`/`color_b` the two independent half-sample estimates
//   - `variance`          per-pixel variance derived from those halves
//
// The half-buffers come for free from the way the GPU backends accumulate: the
// device-side accumulation buffer already holds a running (sum, count) pair per
// pixel and is host-visible, so snapshotting it once at the halfway point of
// the spp loop yields the first half's mean directly, and the second half's
// mean follows by subtracting that snapshot's sums from the final ones. No
// extra rendering work, no extra device memory, and - crucially - no
// divergence between the CPU and GPU sampling code.
//
// The geometric feature buffers (albedo / normal / depth) are NOT produced
// here, as those genuinely require the shaders to write additional per-pixel
// outputs. Denoisers therefore still run without feature guidance on GPU
// renders, but they do get a real, filter-weight-independent variance estimate
// instead of the crude spatial fallback - which is what made GPU denoising
// blur out texture detail.
[[nodiscard]] FrameBufferSet render_gpu_with_aov(const std::shared_ptr<Scene> &scene, GPUBackend backend, int spp = 0,
                                                 const ProgressCallback &progress = nullptr);

} // namespace tiny_renderer::gpu


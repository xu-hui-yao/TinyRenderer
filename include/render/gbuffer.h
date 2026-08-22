#pragma once

#include <components/scene.h>
#include <render/framebuffer.h>

M_NAMESPACE_BEGIN

/**
 * \brief Compute the denoiser feature buffers (albedo / normal / depth) by
 * tracing one camera ray per pixel on the CPU.
 *
 * These buffers describe the FIRST surface each pixel sees. That is a purely
 * deterministic, geometric quantity - it does not depend on light transport and
 * therefore needs no Monte Carlo sampling at all. One primary ray per pixel is
 * enough to compute it exactly.
 *
 * That property is what makes this function useful: the GPU backends accumulate
 * radiance entirely on the device and do not export any feature buffers, so a
 * GPU render would otherwise have to be denoised without feature guidance -
 * which means the filter cannot distinguish a texture edge or a silhouette from
 * noise, and visibly blurs both away. Rather than plumbing extra outputs
 * through every shader, the features are simply recomputed here against the
 * same Scene the GPU render was built from, at the cost of one shadow-ray-free
 * BVH traversal per pixel (a few percent of a typical render).
 *
 * The values produced are identical in meaning (and very nearly in value - the
 * only difference being that this uses one pixel-centered ray where the CPU
 * integrator averages jittered samples) to what ImageBlock accumulates during a
 * CPU render, so both paths hand the denoisers equivalent inputs.
 *
 * \param scene
 *    Must already be constructed (`scene->construct()` called).
 * \param fbs
 *    Filled in: `albedo`, `normal` and `depth` are allocated and written.
 *    Existing members (`color`, half-buffers, `variance`) are left untouched.
 * \param thread_count
 *    Number of worker threads; values < 1 mean hardware concurrency.
 */
void compute_gbuffer(const std::shared_ptr<Scene> &scene, FrameBufferSet &fbs, int thread_count = 0);

M_NAMESPACE_END

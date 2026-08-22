#pragma once

#include <core/common.h>
#include <core/spectrum.h>

// Number of floats per pixel in ImageBlock's auxiliary feature buffer:
//   [0..2] albedo RGB, [3..5] shading normal XYZ, [6] depth, [7] weight
#define M_AOV_FEATURE_CHANNELS 8

// How many consecutive delta (perfectly specular / refractive) bounces
// Path::li_aov() is willing to follow while looking for a surface to read the
// demodulation albedo off of. Without this, a mirror or a glass pane would
// yield albedo == 1 and the textured surface seen THROUGH it would lose the
// benefit of demodulation; with it, the reflected/refracted surface's own
// albedo is used instead. Kept small so that a long chain of specular bounces
// does not end up attributing a far-away surface's texture to this pixel.
#define M_AOV_MAX_DELTA_BOUNCES 3

M_NAMESPACE_BEGIN

/**
 * \brief Per-sample auxiliary (G-buffer-like) features produced alongside the
 * radiance estimate by \ref Integrator::li_aov().
 *
 * These are the "feature buffers" that every modern Monte Carlo denoiser is
 * built on: they are essentially noise-free (they come from the FIRST,
 * deterministic-per-pixel-position camera hit rather than from the stochastic
 * light transport that follows), so a denoiser can use them to decide which
 * neighboring pixels are actually part of the same surface and may therefore
 * be averaged together, and which must be kept separate to preserve an edge.
 *
 * `albedo` additionally enables *demodulation*: filtering
 * `radiance / albedo` (i.e. the smooth, low-frequency irradiance) and
 * re-multiplying afterwards keeps high-frequency texture detail perfectly
 * intact instead of blurring it away. This is standard practice in SVGF /
 * Intel OIDN and is usually worth more than switching to a fancier filter
 * kernel.
 */
struct AOVSample {
    // Reflectance of the first non-delta (ESmooth) BSDF along the path - see
    // BSDF::albedo(). Stays at 1 for paths that never reach a smooth surface
    // (e.g. a camera ray escaping straight into the environment map), which
    // conveniently makes demodulation a no-op for those pixels.
    Color3f albedo{ 1.f };

    // World-space shading normal at the first camera hit. Zero when the ray
    // escaped the scene.
    Normal3f normal{ 0.f };

    // Distance along the camera ray to the first hit
    // (`SurfaceIntersection::t`). Zero when the ray escaped the scene.
    float depth = 0.f;

    // False if the camera ray hit no geometry at all (environment / void), in
    // which case `normal` and `depth` carry no information.
    bool has_feature = false;

    // Internal bookkeeping used by Path::li_aov() while walking the path:
    // whether `albedo` has already been committed.
    bool has_albedo = false;
};

M_NAMESPACE_END

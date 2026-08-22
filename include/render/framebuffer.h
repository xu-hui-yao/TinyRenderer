#pragma once

#include <components/bitmap.h>
#include <memory>
#include <render/aov.h>

M_NAMESPACE_BEGIN

/**
 * \brief A complete set of denoiser inputs extracted from a finished render.
 *
 * Every buffer stores LINEAR radiance / feature data (no tonemapping, no sRGB
 * encoding), already normalized by the reconstruction filter weight - i.e.
 * exactly what \ref Bitmap::save_png() would be handed. Any of the optional
 * members may be null, and every Denoiser implementation must degrade
 * gracefully when they are (see \ref DenoiseInput).
 */
struct FrameBufferSet {
    // Full-sample-count radiance estimate, 3 channels. Never null.
    std::shared_ptr<Bitmap> color;

    // The two half-buffers: samples with an even index went into `color_a`,
    // odd ones into `color_b`, each normalized independently. Two
    // statistically independent estimates of the same image, which is what
    // makes robust variance estimation and cross-validation possible.
    std::shared_ptr<Bitmap> color_a;
    std::shared_ptr<Bitmap> color_b;

    // Per-pixel, per-channel variance OF THE MEAN (i.e. of `color`),
    // estimated from the half-buffers as 1/4 (A - B)^2 and then spatially
    // prefiltered. 3 channels.
    //
    // Why not the textbook sum-of-squares estimator? Because samples are
    // splatted through a reconstruction filter with non-uniform weights (see
    // ImageBlock::put()), so the effective sample count per pixel is not an
    // integer and the weighted-variance formula gets awkward. The half-buffer
    // difference is immune to that: A and B are filtered identically, so the
    // filter weights cancel out.
    std::shared_ptr<Bitmap> variance;

    // Feature buffers, 3 / 3 / 1 channels. Accumulated with a BOX filter (one
    // pixel per sample) rather than the reconstruction filter, so that
    // geometric edges in them stay perfectly sharp - a blurred normal buffer
    // would defeat the whole point of edge-stopping.
    std::shared_ptr<Bitmap> albedo;
    std::shared_ptr<Bitmap> normal;
    std::shared_ptr<Bitmap> depth;
};

/**
 * \brief Non-owning view of a \ref FrameBufferSet handed to
 * \ref Denoiser::denoise().
 *
 * `color` is mandatory; everything else is optional and may be null. This
 * keeps simple denoisers (e.g. the pure image-space outlier rejection filter)
 * usable on renders that carry no auxiliary buffers at all - including the
 * GPU backend's output.
 */
struct DenoiseInput {
    const Bitmap *color    = nullptr;
    const Bitmap *color_a  = nullptr;
    const Bitmap *color_b  = nullptr;
    const Bitmap *variance = nullptr;
    const Bitmap *albedo   = nullptr;
    const Bitmap *normal   = nullptr;
    const Bitmap *depth    = nullptr;

    [[nodiscard]] int rows() const { return color->get_rows(); }
    [[nodiscard]] int cols() const { return color->get_cols(); }

    [[nodiscard]] bool has_features() const { return albedo && normal && depth; }
    [[nodiscard]] bool has_half_buffers() const { return color_a && color_b; }

    static DenoiseInput from(const FrameBufferSet &fbs) {
        DenoiseInput in;
        in.color    = fbs.color.get();
        in.color_a  = fbs.color_a.get();
        in.color_b  = fbs.color_b.get();
        in.variance = fbs.variance.get();
        in.albedo   = fbs.albedo.get();
        in.normal   = fbs.normal.get();
        in.depth    = fbs.depth.get();
        return in;
    }
};

M_NAMESPACE_END

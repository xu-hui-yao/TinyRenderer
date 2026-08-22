#pragma once

#include <components/bitmap.h>
#include <core/array.h>
#include <core/common.h>
#include <core/spectrum.h>
#include <memory>
#include <mutex>
#include <render/framebuffer.h>
#include <vector>

#define M_BLOCK_SIZE 32 /* Block size used for parallelization */

M_NAMESPACE_BEGIN
class ImageBlock {
public:
    /**
     * \param with_aov
     *    When true, additionally allocate and accumulate the auxiliary buffers
     *    needed by the denoisers (see include/render/framebuffer.h): the two
     *    half-buffers used for variance estimation / cross-validation, plus
     *    the albedo / normal / depth feature buffers.
     *
     *    Defaults to false, in which case NOTHING about this class's behavior
     *    changes: no extra memory is allocated, put() executes the exact same
     *    sequence of floating point operations as before, and to_bitmap()
     *    produces bit-for-bit identical output. Denoising is strictly opt-in.
     */
    ImageBlock(const Vector2i &size, const std::shared_ptr<ReconstructionFilter> &filter, bool with_aov = false);
    ~ImageBlock();

    void set_offset(const Point2i &offset);
    [[nodiscard]] const Point2i &get_offset() const;

    void set_size(const Point2i &size);
    [[nodiscard]] const Vector2i &get_size() const;

    [[nodiscard]] int get_border_size() const;

    [[nodiscard]] bool has_aov() const { return m_with_aov; }

    void clear(const Color3f &color = Color3f(0)) const;
    void put(const Point2f &pos, const Color3f &value);

    /**
     * \brief AOV-aware variant of put().
     *
     * Splats `value` exactly like put() above (same filter, same arithmetic),
     * and additionally routes it into half-buffer A or B depending on the
     * parity of `sample_index`, and accumulates `aov`'s features with a box
     * filter. A no-op w.r.t. the auxiliary buffers if this block was
     * constructed with `with_aov == false`.
     */
    void put(const Point2f &pos, const Color3f &value, const AOVSample &aov, uint32_t sample_index);

    void put(const ImageBlock &block);

    std::shared_ptr<Bitmap> to_bitmap() const;
    void from_bitmap(const Bitmap &bitmap);

    /**
     * \brief Extract every buffer this block holds as a denoiser-ready set of
     * linear-radiance bitmaps.
     *
     * `color` is always produced (identical to to_bitmap()). The remaining
     * members are only populated when this block was constructed with
     * `with_aov == true`; `variance` is derived from the half-buffers as
     * 1/4 (A - B)^2 and spatially prefiltered with a 3x3 box (a raw
     * half-buffer difference has only one degree of freedom and is far too
     * noisy to drive a filter directly).
     */
    [[nodiscard]] FrameBufferSet to_framebuffers() const;

    [[nodiscard]] std::string to_string() const;

private:
    // Shared implementation of the two put() overloads above. `aov` may be
    // null (plain radiance splat).
    void put_impl(const Point2f &pos, const Color3f &value, const AOVSample *aov, uint32_t sample_index);

    Point2i m_offset;
    Vector2i m_size;
    int m_border_size = 0;

    std::unique_ptr<float[]> m_filter;
    std::vector<float> m_weights_x, m_weights_y;
    float m_filter_radius = 0;
    float m_lookup_factor = 0;

    mutable std::mutex m_mutex;
    std::unique_ptr<Color4f[]> m_color;
    int m_rows, m_cols;

    // ---------------------------------------------------------------------
    // Optional denoiser inputs (only allocated when m_with_aov is set).
    // ---------------------------------------------------------------------
    bool m_with_aov = false;

    // Same layout/stride as m_color. Even-indexed samples go into A, odd ones
    // into B, using the SAME reconstruction filter weights as m_color, so
    // that (A + B) reproduces m_color exactly and the filter weights cancel
    // out of the (A - B) variance estimate.
    std::unique_ptr<Color4f[]> m_color_a;
    std::unique_ptr<Color4f[]> m_color_b;

    // m_rows * m_cols * M_AOV_FEATURE_CHANNELS floats, box-filtered.
    std::unique_ptr<float[]> m_feature;
};

class BlockGenerator {
public:
    BlockGenerator(const Vector2i &size, int block_size);

    bool next(ImageBlock &block);

    bool has_next() const;

    [[nodiscard]] int get_block_count() const;

private:
    enum EDirection { ERight = 0, EDown, ELeft, EUp };

    Point2i m_block;
    Vector2i m_num_blocks, m_size;
    int m_block_size, m_num_steps, m_blocks_left, m_steps_left, m_direction;
    mutable std::mutex m_mutex;
};

M_NAMESPACE_END

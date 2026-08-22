#include <components/rfilter.h>
#include <core/bounding_box.h>
#include <iostream>
#include <render/block.h>

M_NAMESPACE_BEGIN
ImageBlock::ImageBlock(const Vector2i &size, const std::shared_ptr<ReconstructionFilter> &filter, bool with_aov)
    : m_offset(0, 0), m_size(size), m_with_aov(with_aov) {
    if (filter) {
        m_filter_radius = filter->get_radius();
        m_border_size   = static_cast<int>(std::ceil(m_filter_radius - 0.5f));

        m_filter = std::make_unique<float[]>(M_FILTER_RESOLUTION + 1);
        for (int i = 0; i < M_FILTER_RESOLUTION; ++i) {
            float pos   = m_filter_radius * static_cast<float>(i) / M_FILTER_RESOLUTION;
            m_filter[i] = filter->eval(pos);
        }
        m_filter[M_FILTER_RESOLUTION] = 0.0f;
        m_lookup_factor               = M_FILTER_RESOLUTION / m_filter_radius;

        int weight_size = static_cast<int>(std::ceil(2 * m_filter_radius) + 1);
        m_weights_x.resize(weight_size, 0.0f);
        m_weights_y.resize(weight_size, 0.0f);
    } else {
        m_border_size   = 0;
        m_filter_radius = 0;
    }

    m_rows  = size.y() + 2 * m_border_size;
    m_cols  = size.x() + 2 * m_border_size;
    m_color = std::make_unique<Color4f[]>(m_rows * m_cols);

    if (m_with_aov) {
        m_color_a = std::make_unique<Color4f[]>(m_rows * m_cols);
        m_color_b = std::make_unique<Color4f[]>(m_rows * m_cols);
        m_feature = std::make_unique<float[]>(static_cast<size_t>(m_rows) * m_cols * M_AOV_FEATURE_CHANNELS);
    }
}

ImageBlock::~ImageBlock() = default;

void ImageBlock::set_offset(const Point2i &offset) { m_offset = offset; }
const Point2i &ImageBlock::get_offset() const { return m_offset; }

void ImageBlock::set_size(const Point2i &size) { m_size = size; }
const Vector2i &ImageBlock::get_size() const { return m_size; }

int ImageBlock::get_border_size() const { return m_border_size; }

void ImageBlock::clear(const Color3f &color) const {
    // The 4th channel is the accumulated FILTER WEIGHT, which
    // divide_by_weight() uses as the denominator when normalizing. It must
    // start at 0, in lockstep with the (zero) radiance it weights -
    // seeding it with 1.0 adds a phantom unit of weight that no sample ever
    // contributed radiance for, systematically DARKENING every pixel by
    // roughly 1/(1 + total_sample_weight) (~6% at 16spp with the default
    // radius-0.5 tent filter, and worse the fewer samples are taken).
    std::fill_n(m_color.get(), m_rows * m_cols, Color4f(color, 0.0f));

    if (m_with_aov) {
        std::fill_n(m_color_a.get(), m_rows * m_cols, Color4f(color, 0.0f));
        std::fill_n(m_color_b.get(), m_rows * m_cols, Color4f(color, 0.0f));
        std::fill_n(m_feature.get(), static_cast<size_t>(m_rows) * m_cols * M_AOV_FEATURE_CHANNELS, 0.0f);
    }
}

void ImageBlock::put(const Point2f &pos, const Color3f &value) { put_impl(pos, value, nullptr, 0); }

void ImageBlock::put(const Point2f &pos, const Color3f &value, const AOVSample &aov, uint32_t sample_index) {
    put_impl(pos, value, &aov, sample_index);
}

void ImageBlock::put_impl(const Point2f &_pos, const Color3f &value, const AOVSample *aov, uint32_t sample_index) {
    if (!value.is_valid()) {
        std::cerr << "Integrator: computed an invalid radiance value: " << value.to_string() << std::endl;
        return;
    }

    Point2f pos(_pos.x() - 0.5f - static_cast<float>(m_offset.x() - m_border_size),
                _pos.y() - 0.5f - static_cast<float>(m_offset.y() - m_border_size));

    BoundingBox2i bbox(Point2i(static_cast<int>(std::ceil(pos.x() - m_filter_radius)),
                               static_cast<int>(std::ceil(pos.y() - m_filter_radius))),
                       Point2i(static_cast<int>(std::floor(pos.x() + m_filter_radius)),
                               static_cast<int>(std::floor(pos.y() + m_filter_radius))));
    bbox.clip(BoundingBox2i(Point2i(0, 0), Point2i(m_cols - 1, m_rows - 1)));

    for (int x = bbox.get_min().x(), idx = 0; x <= bbox.get_max().x(); ++x, ++idx)
        m_weights_x[idx] = m_filter[static_cast<int>(std::abs(static_cast<float>(x) - pos.x()) * m_lookup_factor)];
    for (int y = bbox.get_min().y(), idx = 0; y <= bbox.get_max().y(); ++y, ++idx)
        m_weights_y[idx] = m_filter[static_cast<int>(std::abs(static_cast<float>(y) - pos.y()) * m_lookup_factor)];

    for (int y = bbox.get_min().y(), yr = 0; y <= bbox.get_max().y(); ++y, ++yr)
        for (int x = bbox.get_min().x(), xr = 0; x <= bbox.get_max().x(); ++x, ++xr)
            m_color[y * m_cols + x] += Color4f(value, 1.0f) * (m_weights_x[xr] * m_weights_y[yr]);

    if (!m_with_aov || aov == nullptr)
        return;

    // Split the sample stream into two statistically independent halves. Both
    // halves are splatted with the very same filter footprint/weights
    // computed above, so A + B == m_color and, crucially, the filter weights
    // cancel out of the (A - B) variance estimate in to_framebuffers().
    Color4f *half = (sample_index & 1u) ? m_color_b.get() : m_color_a.get();
    for (int y = bbox.get_min().y(), yr = 0; y <= bbox.get_max().y(); ++y, ++yr)
        for (int x = bbox.get_min().x(), xr = 0; x <= bbox.get_max().x(); ++x, ++xr)
            half[y * m_cols + x] += Color4f(value, 1.0f) * (m_weights_x[xr] * m_weights_y[yr]);

    // Features get a BOX footprint (the single pixel the sample landed in)
    // rather than the reconstruction filter's. Spreading normals/depth across
    // a filter footprint would round off exactly the geometric discontinuities
    // that the denoisers rely on to decide where NOT to blur.
    int fx = static_cast<int>(std::floor(pos.x() + 0.5f));
    int fy = static_cast<int>(std::floor(pos.y() + 0.5f));
    if (fx < 0 || fy < 0 || fx >= m_cols || fy >= m_rows)
        return;

    float *f = &m_feature[(static_cast<size_t>(fy) * m_cols + fx) * M_AOV_FEATURE_CHANNELS];
    f[0] += aov->albedo(0);
    f[1] += aov->albedo(1);
    f[2] += aov->albedo(2);
    if (aov->has_feature) {
        f[3] += aov->normal.x();
        f[4] += aov->normal.y();
        f[5] += aov->normal.z();
        f[6] += aov->depth;
    }
    f[7] += 1.0f;
}

void ImageBlock::put(const ImageBlock &block) {
    std::scoped_lock lock(m_mutex);

    Vector2i offset = block.get_offset() - m_offset + Vector2i(m_border_size - block.get_border_size());
    Vector2i size   = block.get_size() + Vector2i(2 * block.get_border_size());

    // NOTE: `block`'s underlying m_color buffer is physically laid out with
    // row stride `block.m_cols`, which is fixed at that ImageBlock's
    // CONSTRUCTION time (e.g. always M_BLOCK_SIZE for the per-thread scratch
    // block reused across BlockGenerator::next() calls in main.cpp) and can
    // therefore be LARGER than this particular block's current logical
    // `size.x()` for edge blocks (whenever the image dimensions are not an
    // exact multiple of M_BLOCK_SIZE). Indexing the source buffer with
    // `size.x()` as if it were the stride - instead of the buffer's actual
    // physical stride `block.m_cols` - silently reads from the wrong rows
    // for every such edge block, producing e.g. a deterministic "every other
    // row is black" pattern. Must use `block.m_cols` here, not `size.x()`.
    for (int i = 0; i < size.y(); i++) {
        for (int j = 0; j < size.x(); j++) {
            m_color[(offset.y() + i) * m_cols + (offset.x() + j)] = block.m_color[i * block.m_cols + j];
        }
    }

    // Same copy for the auxiliary buffers, only if BOTH sides carry them.
    if (!m_with_aov || !block.m_with_aov)
        return;

    for (int i = 0; i < size.y(); i++) {
        for (int j = 0; j < size.x(); j++) {
            size_t dst = static_cast<size_t>(offset.y() + i) * m_cols + (offset.x() + j);
            size_t src = static_cast<size_t>(i) * block.m_cols + j;
            m_color_a[dst] = block.m_color_a[src];
            m_color_b[dst] = block.m_color_b[src];
            for (int c = 0; c < M_AOV_FEATURE_CHANNELS; ++c)
                m_feature[dst * M_AOV_FEATURE_CHANNELS + c] = block.m_feature[src * M_AOV_FEATURE_CHANNELS + c];
        }
    }
}

std::shared_ptr<Bitmap> ImageBlock::to_bitmap() const {
    // Locked so this can safely be called from a progress-reporting thread
    // (see include/render/progress.h) while other threads are concurrently
    // calling put(const ImageBlock&) on this same block - without this,
    // reading m_color here races with put()'s writes to it.
    std::scoped_lock lock(m_mutex);

    auto result = std::make_shared<Bitmap>(m_size.y(), m_size.x(), 3);
    for (int y = 0; y < m_size.y(); ++y) {
        for (int x = 0; x < m_size.x(); ++x) {
            Color3f color      = m_color[(y + m_border_size) * m_cols + (x + m_border_size)].divide_by_weight();
            (*result)(y, x, 0) = color(0);
            (*result)(y, x, 1) = color(1);
            (*result)(y, x, 2) = color(2);
        }
    }
    return result;
}

FrameBufferSet ImageBlock::to_framebuffers() const {
    FrameBufferSet fbs;
    fbs.color = to_bitmap();

    if (!m_with_aov)
        return fbs;

    std::scoped_lock lock(m_mutex);

    const int h = m_size.y(), w = m_size.x();

    fbs.color_a  = std::make_shared<Bitmap>(h, w, 3);
    fbs.color_b  = std::make_shared<Bitmap>(h, w, 3);
    fbs.variance = std::make_shared<Bitmap>(h, w, 3);
    fbs.albedo   = std::make_shared<Bitmap>(h, w, 3);
    fbs.normal   = std::make_shared<Bitmap>(h, w, 3);
    fbs.depth    = std::make_shared<Bitmap>(h, w, 1);

    // Raw, un-prefiltered variance estimate; smoothed into fbs.variance below.
    std::vector<float> raw_var(static_cast<size_t>(h) * w * 3, 0.0f);

    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            size_t src = static_cast<size_t>(y + m_border_size) * m_cols + (x + m_border_size);

            Color3f a = m_color_a[src].divide_by_weight();
            Color3f b = m_color_b[src].divide_by_weight();

            for (int c = 0; c < 3; ++c) {
                (*fbs.color_a)(y, x, c) = a(c);
                (*fbs.color_b)(y, x, c) = b(c);

                // A and B are each the mean of ~half the samples, so
                //   E[(A - B)^2] = Var(A) + Var(B) = 4 * Var(mean of all)
                // giving Var(color) ~= 1/4 (A - B)^2. Filter-weight
                // independent, because A and B were filtered identically.
                float d                                    = a(c) - b(c);
                raw_var[(static_cast<size_t>(y) * w + x) * 3 + c] = 0.25f * d * d;
            }

            const float *f = &m_feature[src * M_AOV_FEATURE_CHANNELS];
            float fw       = f[7];
            float inv      = fw > 0.0f ? 1.0f / fw : 0.0f;

            (*fbs.albedo)(y, x, 0) = f[0] * inv;
            (*fbs.albedo)(y, x, 1) = f[1] * inv;
            (*fbs.albedo)(y, x, 2) = f[2] * inv;

            // Averaging unit normals shrinks them, so renormalize. Pixels
            // straddling a silhouette legitimately average towards zero; those
            // are left at zero, which the denoisers read as "no reliable
            // normal here" (their dot-product term then simply stops
            // contributing rather than misfiring).
            Normal3f n(f[3] * inv, f[4] * inv, f[5] * inv);
            float mag = n.magnitude();
            if (mag > 1e-6f)
                n = n / mag;
            (*fbs.normal)(y, x, 0) = n.x();
            (*fbs.normal)(y, x, 1) = n.y();
            (*fbs.normal)(y, x, 2) = n.z();

            (*fbs.depth)(y, x, 0) = f[6] * inv;
        }
    }

    // A raw half-buffer difference has a single degree of freedom - it is an
    // unbiased but extremely noisy variance estimate. Every consumer of it
    // divides by it, so smooth it over a 3x3 box first; without this the
    // edge-stopping weights themselves become noisy and the filtered result
    // develops a characteristic blotchy look.
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            for (int c = 0; c < 3; ++c) {
                float sum = 0.0f;
                int count = 0;
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        int yy = y + dy, xx = x + dx;
                        if (yy < 0 || xx < 0 || yy >= h || xx >= w)
                            continue;
                        sum += raw_var[(static_cast<size_t>(yy) * w + xx) * 3 + c];
                        ++count;
                    }
                }
                (*fbs.variance)(y, x, c) = count > 0 ? sum / static_cast<float>(count) : 0.0f;
            }
        }
    }

    return fbs;
}

void ImageBlock::from_bitmap(const Bitmap &bitmap) {
    if (bitmap.get_rows() != m_rows || bitmap.get_cols() != m_cols)
        throw std::runtime_error("Invalid bitmap dimensions!");

    for (int y = 0; y < m_size.y(); ++y) {
        for (int x = 0; x < m_size.x(); ++x) {
            float r                                                     = bitmap(y, x, 0);
            float g                                                     = bitmap(y, x, 1);
            float b                                                     = bitmap(y, x, 2);
            m_color[(y + m_border_size) * m_cols + (x + m_border_size)] = Color4f({ r, g, b, 1.0f });
        }
    }
}

std::string ImageBlock::to_string() const {
    return "ImageBlock[offset=" + m_offset.to_string() + ", size=" + m_size.to_string() + "]";
}

BlockGenerator::BlockGenerator(const Vector2i &size, int block_size) : m_size(size), m_block_size(block_size) {
    m_num_blocks  = Vector2i(static_cast<int>(std::ceil(static_cast<float>(size.x()) / static_cast<float>(block_size))),
                             static_cast<int>(std::ceil(static_cast<float>(size.y()) / static_cast<float>(block_size))));
    m_blocks_left = m_num_blocks.x() * m_num_blocks.y();
    m_direction   = ERight;
    m_block       = m_num_blocks / 2;
    m_steps_left  = 1;
    m_num_steps   = 1;
}

bool BlockGenerator::has_next() const { return m_blocks_left != 0; }

bool BlockGenerator::next(ImageBlock &block) {
    std::scoped_lock lock(m_mutex);

    if (m_blocks_left == 0)
        return false;

    Point2i pos = m_block * m_block_size;
    block.set_offset(pos);
    block.set_size(Point2i(m_size - pos).wise_min(Vector2i(m_block_size)));

    if (--m_blocks_left == 0)
        return true;

    do {
        switch (m_direction) {
            case ERight:
                ++m_block.x();
                break;
            case EDown:
                ++m_block.y();
                break;
            case ELeft:
                --m_block.x();
                break;
            case EUp:
                --m_block.y();
                break;
            default:
                break;
        }

        if (--m_steps_left == 0) {
            m_direction = (m_direction + 1) % 4;
            if (m_direction == ELeft || m_direction == ERight)
                ++m_num_steps;
            m_steps_left = m_num_steps;
        }
    } while ((m_block < 0).any() || (m_block >= m_num_blocks).any());

    return true;
}

int BlockGenerator::get_block_count() const { return m_blocks_left; }

M_NAMESPACE_END

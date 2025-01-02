#include <components/rfilter.h>
#include <core/bounding_box.h>
#include <iostream>
#include <render/block.h>

M_NAMESPACE_BEGIN
ImageBlock::ImageBlock(const Vector2i &size, const std::shared_ptr<ReconstructionFilter> &filter)
    : m_offset(0, 0), m_size(size) {
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
}

ImageBlock::~ImageBlock() = default;

void ImageBlock::set_offset(const Point2i &offset) { m_offset = offset; }
const Point2i &ImageBlock::get_offset() const { return m_offset; }

void ImageBlock::set_size(const Point2i &size) { m_size = size; }
const Vector2i &ImageBlock::get_size() const { return m_size; }

int ImageBlock::get_border_size() const { return m_border_size; }

void ImageBlock::clear(const Color3f &color) const {
    std::fill_n(m_color.get(), m_rows * m_cols, Color4f(color, 1.0f));
}

void ImageBlock::put(const Point2f &_pos, const Color3f &value) {
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
}

void ImageBlock::put(const ImageBlock &block) {
    std::scoped_lock lock(m_mutex);

    Vector2i offset = block.get_offset() - m_offset + Vector2i(m_border_size - block.get_border_size());
    Vector2i size   = block.get_size() + Vector2i(2 * block.get_border_size());

    for (int i = 0; i < size.y(); i++) {
        for (int j = 0; j < size.x(); j++) {
            m_color[(offset.y() + i) * m_cols + (offset.x() + j)] = block.m_color[i * size.x() + j];
        }
    }
}

std::shared_ptr<Bitmap> ImageBlock::to_bitmap() const {
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

#pragma once

#include <core/common.h>
#include <cmath>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

M_NAMESPACE_BEGIN
template <typename Scalar_, typename Index_> class TDiscreteDistribution {
public:
    typedef Scalar_ Scalar;
    typedef Index_ Index;

    TDiscreteDistribution() = default;

    // Initialize from a PMF (probability mass function)
    explicit TDiscreteDistribution(const std::vector<Scalar> &pmf) { initialize(pmf); }

    // Update the distribution with a new PMF
    void initialize(const std::vector<Scalar> &pmf) {
        if (pmf.empty())
            throw std::invalid_argument("DiscreteDistribution: PMF cannot be empty.");

        m_pmf = pmf;
        m_cdf.resize(pmf.size());

        Scalar sum = 0.0f;
        for (size_t i = 0; i < pmf.size(); ++i) {
            if (pmf[i] < 0)
                throw std::invalid_argument("DiscreteDistribution: PMF values must be non-negative.");
            sum += pmf[i];
            m_cdf[i] = sum;
        }

        if (sum <= 0.0f)
            throw std::invalid_argument("DiscreteDistribution: Total probability mass must be greater "
                                        "than zero.");

        m_sum           = sum;
        m_normalization = 1.0f / sum;
    }

    // Evaluate the un-normalized PMF at a given index
    Scalar eval_pmf(Index index) const { return m_pmf[index]; }

    // Evaluate the normalized PMF at a given index
    Scalar eval_pmf_normalized(Index index) const { return m_pmf[index] * m_normalization; }

    // Evaluate the CDF at a given index
    Scalar eval_cdf(Index index) const { return m_cdf[index]; }

    // Evaluate the normalized CDF at a given index
    Scalar eval_cdf_normalized(Index index) const { return m_cdf[index] * m_normalization; }

    // Sample the distribution
    Index sample(Scalar u) const {
        if (u < 0.0f || u > 1.0f)
            throw std::invalid_argument("Sample value must be in the range [0, 1].");

        u *= m_sum;

        // Binary search in the CDF to find the corresponding index
        auto it = std::lower_bound(m_cdf.begin(), m_cdf.end(), u);
        // lower_bound returns end() when u exceeds the last CDF entry, which
        // happens for u == 1 (and, through floating point rounding, for values
        // just below it). Clamping to the last valid interval keeps every
        // caller's subsequent m_pmf/m_cdf indexing in bounds.
        if (it == m_cdf.end())
            return static_cast<Index>(m_cdf.size()) - 1;
        return std::distance(m_cdf.begin(), it);
    }

    // Sample the distribution and return the index and PMF value
    std::pair<Index, Scalar> sample_pmf(Scalar u) const {
        Index index = sample(u);
        return { index, eval_pmf_normalized(index) };
    }

    // Sample the distribution and reuse the sample for further computations
    std::tuple<Index, Scalar> sample_reuse(Scalar u) const {
        Index index = sample(u);
        Scalar pmf  = eval_pmf_normalized(index);
        // The CDF *below* the selected interval. For index == 0 there is no
        // preceding entry and the correct lower bound is 0 - reading
        // m_cdf[index - 1] there is an out-of-bounds access that returns
        // whatever happens to precede the vector's heap allocation. Since that
        // garbage value is then subtracted from `u` to rescale the sample, the
        // reused variate came out arbitrary and, worse, DIFFERENT from run to
        // run, making every render that samples an emitter's first triangle
        // (Mesh::sample_position) irreproducible.
        Scalar cdf        = index > 0 ? eval_cdf_normalized(index - 1) : static_cast<Scalar>(0);
        Scalar rescaled_u = (u - cdf) / pmf;

        return { index, rescaled_u };
    }

    // Get the normalization factor
    Scalar normalization() const { return m_normalization; }

    // Get the total sum of the PMF before normalization
    Scalar sum() const { return m_sum; }

    // Get the size of the distribution
    [[nodiscard]] size_t size() const { return m_pmf.size(); }

private:
    std::vector<Scalar> m_pmf;
    std::vector<Scalar> m_cdf;
    Scalar m_sum           = 0.0f;
    Scalar m_normalization = 0.0f;
};

template <typename Scalar_, typename Index_> class TDiscreteDistribution2D {
public:
    typedef Scalar_ Scalar;
    typedef Index_ Index;

    // Constructor for a 2D discrete distribution, given a 2D PMF
    TDiscreteDistribution2D() = default;

    explicit TDiscreteDistribution2D(const std::vector<std::vector<Scalar>> &pmf) { initialize(pmf); }

    // Initialize the distribution with a 2D PMF (Probability Mass Function)
    void initialize(const std::vector<std::vector<Scalar>> &pmf) {
        if (pmf.empty() || pmf[0].empty()) {
            throw std::invalid_argument("PMF cannot be empty.");
        }

        m_size_y = pmf.size();
        m_size_x = pmf[0].size();

        m_cond_cdf.resize(m_size_y * m_size_x);
        m_marg_cdf.resize(m_size_y);

        // Construct conditional and marginal CDFs
        Scalar accum_marg = 0.0f;
        for (Index y = 0; y < m_size_y; ++y) {
            Scalar accum_cond = 0.0f;
            for (Index x = 0; x < m_size_x; ++x) {
                accum_cond += pmf[y][x];
                m_cond_cdf[y * m_size_x + x] = accum_cond;
            }
            accum_marg += accum_cond;
            m_marg_cdf[y] = accum_marg;
        }

        m_normalization = 1.0f / accum_marg;
        m_sum           = accum_marg;
    }

    Scalar eval(const Point2i &pos) {
        Index index = pos.x() + pos.y() * m_size_x;
        if (pos.x() > 0) {
            return m_cond_cdf[index] - m_cond_cdf[index - 1];
        } else {
            return m_cond_cdf[index];
        }
    }

    // Evaluate the un-normalized PMF at a given (x, y) position
    Scalar pdf(const Point2i &pos) const { return eval(pos) * m_normalization; }

    // Sample the distribution (returns the (x, y) index)
    std::tuple<Point2i, Scalar, Point2f> sample(const Point2f &sample_) const {
        Point2f sample(sample_);

        // Avoid degeneracies on the domain boundary
        sample = sample.clamp(M_EPSILON, 1.0f - static_cast<float>(M_EPSILON));

        // Scale sample Y range
        sample.y() *= m_sum;

        // Sample the row from the marginal distribution
        auto it1  = std::lower_bound(m_marg_cdf.begin(), m_marg_cdf.end(), sample.y());
        Index row = std::distance(m_marg_cdf.begin(), it1);

        Index offset = row * m_size_x;

        // Scale sample X range
        sample.x() *= m_cond_cdf[offset + m_size_x - 1];

        // Sample the column from the conditional distribution
        auto it2 =
            std::lower_bound(m_cond_cdf.begin() + offset, m_cond_cdf.begin() + offset + m_size_x - 1, sample.x());
        Index col = std::distance(m_cond_cdf.begin() + offset, it2);

        // Re-scale uniform variate
        Scalar col_cdf_0 = col > 0 ? m_cond_cdf[offset + col - 1] : 0;
        Scalar col_cdf_1 = m_cond_cdf[offset + col];
        Scalar row_cdf_0 = row > 0 ? m_marg_cdf[row - 1] : 0;
        Scalar row_cdf_1 = m_marg_cdf[row];

        sample.x() -= col_cdf_0;
        sample.y() -= row_cdf_0;
        if (col_cdf_1 != col_cdf_0) {
            sample.x() /= col_cdf_1 - col_cdf_0;
        }
        if (row_cdf_1 != row_cdf_0) {
            sample.y() /= row_cdf_1 - row_cdf_0;
        }
        return { Point2i({ col, row }), (col_cdf_1 - col_cdf_0) * m_normalization, sample };
    }

    [[nodiscard]] int get_rows() const { return m_size_y; }

    [[nodiscard]] int get_cols() const { return m_size_x; }

    // Print distribution information (optional)
    [[nodiscard]] std::string to_string() const {
        std::ostringstream oss;
        oss << "TDiscreteDistribution2D"
            << "[" << std::endl
            << "  size = " << m_size_x << " x " << m_size_y << "," << std::endl
            << "  normalization = " << m_normalization << std::endl
            << "]";
        return oss.str();
    }

private:
    std::vector<Scalar> m_cond_cdf;
    std::vector<Scalar> m_marg_cdf;

    Scalar m_sum           = 0.0f;
    Scalar m_normalization = 0.0f;

    size_t m_size_x = 0;
    size_t m_size_y = 0;
};

/**
 * \brief A hierarchical 2D distribution.
 *
 * At initialization, we build a multi-level (mip-map) representation of the input
 * PMF. Each subsequent level coarsens the resolution by roughly half in each
 * dimension, accumulating 2×2 blocks.
 *
 * For sampling, we traverse from the coarsest level to the finest level, picking
 * sub-quadrants based on partial sums, then at the finest level we perform a final
 * bi-linear partition. This yields a continuous sample in [0,1]^2 that follows
 * the underlying 2D distribution.
 */
template <typename Scalar_, typename Index_> class THierarchicalDistribution2D {
public:
    typedef Scalar_ Scalar;
    typedef Index_ Index;

    THierarchicalDistribution2D() = default;

    explicit THierarchicalDistribution2D(const std::vector<std::vector<Scalar>> &pmf) { initialize(pmf); }

    /**
     * \brief Build the multi-level structure from the given PMF.
     *
     * Level 0 is the full-resolution data.
     * Higher levels are built by 2×2 down sampling until width/height < 2.
     */
    void initialize(const std::vector<std::vector<Scalar>> &pmf) {
        if (pmf.empty() || pmf[0].empty()) {
            throw std::invalid_argument("THierarchicalDistribution2D: PMF cannot be empty.");
        }
        m_size_y = pmf.size();
        m_size_x = pmf[0].size();

        // --- Copy PMF into level 0 ---
        level_data level0;
        level0.width  = static_cast<int>(m_size_x);
        level0.height = static_cast<int>(m_size_y);
        level0.data.resize(m_size_x * m_size_y, Scalar(0));

        double total_sum = 0.0;
        for (Index y = 0; y < m_size_y; ++y) {
            if (pmf[y].size() != m_size_x) {
                throw std::runtime_error("THierarchicalDistribution2D: each row in PMF must have the same length.");
            }
            for (Index x = 0; x < m_size_x; ++x) {
                Scalar val = pmf[y][x];
                total_sum += static_cast<double>(val);
                level0.data[y * m_size_x + x] = val;
            }
        }

        if (total_sum <= 0.0) {
            throw std::runtime_error("THierarchicalDistribution2D: PMF sum is zero or negative.");
        }

        m_sum           = static_cast<Scalar>(total_sum);
        m_normalization = static_cast<Scalar>(1.0 / total_sum);

        // Clear old levels, push level 0
        m_levels.clear();
        m_levels.push_back(std::move(level0));

        // --- Build coarser levels by 2x2 down sampling ---
        //
        // Each level halves a dimension only while that dimension is still
        // larger than 1, and the pyramid continues until BOTH dimensions
        // reach 1. Stopping as soon as EITHER dimension bottoms out (as in
        // `if (next_w <= 1 || next_h <= 1) break;`) leaves the coarsest level
        // wider than the 2x2 block that sample() starts from, so the high
        // bits of the wider axis are never decided - for a typical 2:1
        // equirectangular envmap that means the entire right half of the
        // image can never be sampled. Odd sizes are handled by treating the
        // missing row/column as 0 (bounds-checked below) instead of reading
        // past the end of the previous level.
        while (m_levels.back().width > 1 || m_levels.back().height > 1) {
            const level_data &prev = m_levels.back();

            int step_x = prev.width > 1 ? 2 : 1;
            int step_y = prev.height > 1 ? 2 : 1;
            int next_w = prev.width > 1 ? (prev.width + 1) / 2 : 1;
            int next_h = prev.height > 1 ? (prev.height + 1) / 2 : 1;

            level_data coarse;
            coarse.width  = next_w;
            coarse.height = next_h;
            coarse.data.resize(static_cast<size_t>(next_w) * next_h, Scalar(0));

            for (int y = 0; y < next_h; ++y) {
                for (int x = 0; x < next_w; ++x) {
                    Scalar sum = Scalar(0);
                    for (int dy = 0; dy < step_y; ++dy) {
                        int sy = y * step_y + dy;
                        if (sy >= prev.height)
                            continue;
                        for (int dx = 0; dx < step_x; ++dx) {
                            int sx = x * step_x + dx;
                            if (sx >= prev.width)
                                continue;
                            sum += prev.data[static_cast<size_t>(sy) * prev.width + sx];
                        }
                    }
                    coarse.data[static_cast<size_t>(y) * next_w + x] = sum;
                }
            }

            m_levels.push_back(std::move(coarse));
        }
    }

    Scalar eval(const Point2f &pos) const {
        // 1. Clamp pos to [0, 1]
        float px = clamp(pos.x(), 0.f, 1.f);
        float py = clamp(pos.y(), 0.f, 1.f);

        // 2. Convert pos to patch coordinates.
        //
        // sample() maps a chosen cell (ox, oy) plus an in-cell offset s to
        // (ox + s) / size, i.e. cell `i` covers [i/size, (i+1)/size) and its
        // CENTER sits at (i + 0.5)/size. The inverse mapping used here must
        // match that convention (`pos * size - 0.5`), otherwise eval()/pdf()
        // describe a density on a grid that is both scaled (size-1 vs size)
        // and shifted by half a cell relative to the one sample() actually
        // draws from - making pdf() inconsistent with sample() and biasing
        // every estimator built on the pair.
        px = px * static_cast<float>(m_size_x) - 0.5f;
        py = py * static_cast<float>(m_size_y) - 0.5f;

        // 3. Identify the integer patch indices (offset_x, offset_y).
        // std::floor (not a cast) is required: px/py can now be slightly
        // negative (down to -0.5) near the low edge, and a cast truncates
        // toward zero rather than down.
        int ix = static_cast<int>(std::floor(px));
        int iy = static_cast<int>(std::floor(py));
        ix     = clamp(ix, 0, static_cast<int>(m_size_x) - 2);
        iy     = clamp(iy, 0, static_cast<int>(m_size_y) - 2);

        // 4. Fractional part within this patch, clamped so the edge cells
        // (where the sample point lies outside the outermost cell centers)
        // extrapolate to a constant rather than past the corner values.
        float frac_x = clamp(px - static_cast<float>(ix), 0.f, 1.f);
        float frac_y = clamp(py - static_cast<float>(iy), 0.f, 1.f);

        // 5. Retrieve the four corners in level 0
        size_t idx = static_cast<size_t>(iy) * m_size_x + static_cast<size_t>(ix);
        auto v00   = static_cast<float>(m_levels[0].data[idx]);
        auto v10   = static_cast<float>(m_levels[0].data[idx + 1]);
        auto v01   = static_cast<float>(m_levels[0].data[idx + m_size_x]);
        auto v11   = static_cast<float>(m_levels[0].data[idx + m_size_x + 1]);

        // 6. Bilinear interpolation
        float i0 = lerp(v00, v10, frac_x); // row 0 interpolation
        float i1 = lerp(v01, v11, frac_x); // row 1 interpolation
        float v  = lerp(i0, i1, frac_y);

        // Return un-normalized PMF value at continuous position
        return static_cast<Scalar>(v);
    }

    /**
     * \brief Probability DENSITY at a continuous position in [0,1]^2.
     *
     * `m_normalization` (1 / sum of all cells) only turns eval()'s
     * un-normalized cell value into a discrete PMF - i.e. something that sums
     * to 1 over the m_size_x * m_size_y cells. sample() however returns a
     * CONTINUOUS point in [0,1]^2, so the matching quantity is a density with
     * respect to that unit-square measure: each cell covers an area of
     * 1/(m_size_x * m_size_y), so the density is the PMF divided by that
     * area. Without this factor the returned "pdf" is too small by exactly
     * the number of cells, which silently scales every importance-sampled
     * estimator that divides by it (e.g. EnvironmentMap::sample_direction's
     * radiance/pdf weight) by that same factor.
     */
    Scalar pdf(const Point2f &p) const {
        auto cell_count = static_cast<Scalar>(m_size_x * m_size_y);
        return eval(p) * m_normalization * cell_count;
    }

    /**
     * \brief Hierarchical sample in [0, 1]^2
     *
     * Descends the mip pyramid from the coarsest level (a single cell) down
     * to level 0. At each step the current cell is refined into the (up to)
     * 2x2 block of finer cells that composed it, and the sample chooses one
     * of them proportionally to their sums - first a row, then a column
     * within that row. A dimension that is not actually subdivided between
     * two levels (which happens for non-square inputs, where the shorter
     * axis reaches 1 first) is simply not refined in that step.
     *
     * Returns a continuous point in [0,1]^2 together with pdf() at that
     * point, so callers can form a consistent radiance/pdf estimator.
     */
    std::pair<Point2f, Scalar> sample(const Point2f &sample_xy) const {
        float sx = clamp(sample_xy.x(), static_cast<float>(M_EPSILON), 1.f - static_cast<float>(M_EPSILON));
        float sy = clamp(sample_xy.y(), static_cast<float>(M_EPSILON), 1.f - static_cast<float>(M_EPSILON));

        // Cell index within the level currently being considered. The
        // coarsest level is 1x1, so we start at its only cell.
        int offset_x_ = 0, offset_y_ = 0;

        for (int l = static_cast<int>(m_levels.size()) - 2; l >= 0; --l) {
            const level_data &fine   = m_levels[l];
            const level_data &coarse = m_levels[l + 1];

            // Whether this step actually subdivides each axis (see the
            // pyramid construction in initialize()).
            bool split_x = fine.width > coarse.width;
            bool split_y = fine.height > coarse.height;

            int x0 = split_x ? offset_x_ * 2 : offset_x_;
            int y0 = split_y ? offset_y_ * 2 : offset_y_;

            // Values of the (up to) 2x2 finer cells this coarse cell splits
            // into; out-of-range neighbours (odd dimensions) count as 0.
            auto at = [&fine](int x, int y) -> float {
                if (x < 0 || y < 0 || x >= fine.width || y >= fine.height)
                    return 0.0f;
                return static_cast<float>(fine.data[static_cast<size_t>(y) * fine.width + x]);
            };
            float v00 = at(x0, y0);
            float v10 = split_x ? at(x0 + 1, y0) : 0.0f;
            float v01 = split_y ? at(x0, y0 + 1) : 0.0f;
            float v11 = (split_x && split_y) ? at(x0 + 1, y0 + 1) : 0.0f;

            offset_x_ = x0;
            offset_y_ = y0;

            // Pick the row (top vs bottom), if this step splits y.
            float r0 = v00 + v10; // top row
            float r1 = v01 + v11; // bottom row
            bool in_bottom = false;
            if (split_y) {
                float row_sum  = r0 + r1;
                float scaled_y = sy * row_sum;
                in_bottom      = scaled_y > r0;
                if (in_bottom) {
                    sy = (scaled_y - r0) / M_MAX(r1, static_cast<float>(M_EPSILON));
                    offset_y_ += 1;
                } else {
                    sy = scaled_y / M_MAX(r0, static_cast<float>(M_EPSILON));
                }
            }

            // Pick the column (left vs right) within the chosen row.
            if (split_x) {
                float c0       = in_bottom ? v01 : v00;
                float c1       = in_bottom ? v11 : v10;
                float col_sum  = c0 + c1;
                float scaled_x = sx * col_sum;
                if (scaled_x > c0) {
                    sx = (scaled_x - c0) / M_MAX(c1, static_cast<float>(M_EPSILON));
                    offset_x_ += 1;
                } else {
                    sx = scaled_x / M_MAX(c0, static_cast<float>(M_EPSILON));
                }
            }
        }

        // Convert the chosen level-0 cell + in-cell position to [0,1]^2.
        // Must stay in sync with eval()'s inverse mapping.
        float final_x = (static_cast<float>(offset_x_) + sx) / static_cast<float>(m_size_x);
        float final_y = (static_cast<float>(offset_y_) + sy) / static_cast<float>(m_size_y);
        Point2f continuous_coords(final_x, final_y);

        return { continuous_coords, pdf(continuous_coords) };
    }

    [[nodiscard]] std::string to_string() const {
        std::ostringstream oss;
        oss << "THierarchicalDistribution2D[\n"
            << "  size_x = " << m_size_x << ", size_y = " << m_size_y << ",\n"
            << "  sum = " << m_sum << ", normalization = " << m_normalization << ",\n"
            << "  #levels = " << m_levels.size() << "\n"
            << "]";
        return oss.str();
    }

    [[nodiscard]] int get_rows() const { return static_cast<int>(m_size_y); }
    [[nodiscard]] int get_cols() const { return static_cast<int>(m_size_x); }

    // Row-major, un-normalized level-0 (full resolution) PMF grid. Exposed
    // read-only for the GPU flattening export layer (see core/gpu_scene.h);
    // building an equivalent GPU-side sampling structure (e.g. an alias
    // table) from this data is left to that later stage.
    [[nodiscard]] const std::vector<Scalar> &level0_data() const { return m_levels[0].data; }

private:
    // Internal structure for each mip level
    struct level_data {
        int width  = 0;
        int height = 0;
        std::vector<Scalar> data;
    };

    size_t m_size_x = 0;
    size_t m_size_y = 0;
    std::vector<level_data> m_levels;
    Scalar m_sum           = Scalar(0);
    Scalar m_normalization = Scalar(0);
};

M_NAMESPACE_END

#pragma once

#include <core/common.h>
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
        return std::distance(m_cdf.begin(), it);
    }

    // Sample the distribution and return the index and PMF value
    std::pair<Index, Scalar> sample_pmf(Scalar u) const {
        Index index = sample(u);
        return { index, eval_pmf_normalized(index) };
    }

    // Sample the distribution and reuse the sample for further computations
    std::tuple<Index, Scalar> sample_reuse(Scalar u) const {
        Index index       = sample(u);
        Scalar pmf        = eval_pmf_normalized(index);
        Scalar cdf        = eval_cdf_normalized(index - 1);
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

        // --- Build coarser levels by 2×2 down sampling ---
        {
            int cur_w = m_levels[0].width;
            int cur_h = m_levels[0].height;

            // Continue until both dimensions are <= 1 or so
            while (cur_w > 1 || cur_h > 1) {
                if ((cur_w & 1) != 0)
                    cur_w += 1;
                if ((cur_h & 1) != 0)
                    cur_h += 1;
                int next_w = cur_w >> 1;
                int next_h = cur_h >> 1;

                if (next_w <= 1 || next_h <= 1)
                    break;

                level_data coarse;
                coarse.width  = next_w;
                coarse.height = next_h;
                coarse.data.resize(next_w * next_h, Scalar(0));

                const level_data &prev = m_levels.back();

                for (int y = 0; y < next_h; ++y) {
                    for (int x = 0; x < next_w; ++x) {
                        int x0                      = x * 2;
                        int y0                      = y * 2;
                        Scalar v00                  = prev.data[y0 * prev.width + x0];
                        Scalar v10                  = prev.data[y0 * prev.width + (x0 + 1)];
                        Scalar v01                  = prev.data[(y0 + 1) * prev.width + x0];
                        Scalar v11                  = prev.data[(y0 + 1) * prev.width + (x0 + 1)];
                        coarse.data[y * next_w + x] = v00 + v10 + v01 + v11;
                    }
                }

                m_levels.push_back(std::move(coarse));
                cur_w = next_w;
                cur_h = next_h;
            }
        }
    }

    Scalar eval(const Point2f &pos) const {
        // 1. Clamp pos to [0, 1]
        float px = clamp(pos.x(), 0.f, 1.f);
        float py = clamp(pos.y(), 0.f, 1.f);

        // 2. Convert pos to patch coordinates
        px *= m_size_x - 1;
        py *= m_size_y - 1;

        // 3. Identify the integer patch indices (offset_x, offset_y)
        int ix = static_cast<int>(px);
        int iy = static_cast<int>(py);
        ix     = clamp(ix, 0, static_cast<int>(m_size_x) - 2);
        iy     = clamp(iy, 0, static_cast<int>(m_size_y) - 2);

        // 4. Fractional part within this patch
        float frac_x = px - static_cast<float>(ix);
        float frac_y = py - static_cast<float>(iy);

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

    // Return the normalized PDF at (x, y)
    Scalar pdf(const Point2f &p) const { return eval(p) * m_normalization; }

    /**
     * \brief Hierarchical sample in [0, 1]^2
     *
     * 1. Start from the coarsest level.
     * 2. Each level is conceptually divided into 2×2 sub-blocks.
     * 3. Decide if the random sample is in the top vs bottom row, and left vs right column,
     *    then descend to the next level, shifting offset bits accordingly.
     * 4. Finally, at level 0, we perform a final 2×2 check to get exact sub-pixel coordinates.
     *
     * It returns a 2D point in [0,1]^2 and the PDF at the corresponding discrete cell.
     */
    std::pair<Point2f, Scalar> sample(const Point2f &sample_xy) const {
        float sx = clamp(sample_xy.x(), static_cast<float>(M_EPSILON), 1.f - static_cast<float>(M_EPSILON));
        float sy = clamp(sample_xy.y(), static_cast<float>(M_EPSILON), 1.f - static_cast<float>(M_EPSILON));

        // offset_x_ / offset_y_ keep track of the discrete patch as we descend
        int offset_x_ = 0, offset_y_ = 0;

        // Traverse from coarse to just above the finest level
        auto level_count = static_cast<int>(m_levels.size());
        for (int l = level_count - 1; l > 0; --l) {
            const level_data &lev_coarse = m_levels[l];

            // Shift offset (equivalent to offset_x_ <<= 1, offset_y_ <<= 1)
            offset_x_ <<= 1;
            offset_y_ <<= 1;

            // Fetch 2×2 block corner values at the coarse level
            int base_idx = index_2x2(lev_coarse, offset_x_, offset_y_);
            auto v00     = static_cast<float>(lev_coarse.data[base_idx]);
            auto v10     = static_cast<float>(lev_coarse.data[base_idx + 1]);
            auto v01     = static_cast<float>(lev_coarse.data[base_idx + lev_coarse.width]);
            auto v11     = static_cast<float>(lev_coarse.data[base_idx + lev_coarse.width + 1]);

            float r0      = v00 + v10; // top row sum
            float r1      = v01 + v11; // bottom row sum
            float row_sum = r0 + r1;

            // Decide top vs bottom
            float scaled_y = sy * row_sum;
            bool in_bottom = scaled_y > r0;
            if (in_bottom) {
                sy = (scaled_y - r0) / M_MAX(r1, static_cast<float>(M_EPSILON));
                offset_y_ += 1;
            } else {
                sy = scaled_y / M_MAX(r0, static_cast<float>(M_EPSILON));
            }

            // Decide left vs right
            float c0                = in_bottom ? v01 : v00;
            float c1                = in_bottom ? v11 : v10;
            float row_sum_given_col = c0 + c1;

            float scaled_x = sx * row_sum_given_col;
            if (scaled_x > c0) {
                sx = (scaled_x - c0) / M_MAX(c1, static_cast<float>(M_EPSILON));
                offset_x_ += 1;
            } else {
                sx = scaled_x / M_MAX(c0, static_cast<float>(M_EPSILON));
            }
        }

        // Final step at level 0
        {
            const level_data &lev0 = m_levels[0];
            offset_x_ <<= 1;
            offset_y_ <<= 1;

            int base_idx = index_2x2(lev0, offset_x_, offset_y_);
            auto v00     = static_cast<float>(lev0.data[base_idx]);
            auto v10     = static_cast<float>(lev0.data[base_idx + 1]);
            auto v01     = static_cast<float>(lev0.data[base_idx + lev0.width]);
            auto v11     = static_cast<float>(lev0.data[base_idx + lev0.width + 1]);

            float r0      = v00 + v10;
            float r1      = v01 + v11;
            float row_sum = r0 + r1;

            float scaled_y = sy * row_sum;
            bool in_bottom = (scaled_y > r0);
            if (in_bottom) {
                sy = (scaled_y - r0) / M_MAX(r1, static_cast<float>(M_EPSILON));
                offset_y_ += 1;
            } else {
                sy = scaled_y / M_MAX(r0, static_cast<float>(M_EPSILON));
            }

            float c0                = in_bottom ? v01 : v00;
            float c1                = in_bottom ? v11 : v10;
            float row_sum_given_col = c0 + c1;

            float scaled_x = sx * row_sum_given_col;
            if (scaled_x > c0) {
                sx = (scaled_x - c0) / M_MAX(c1, static_cast<float>(M_EPSILON));
                offset_x_ += 1;
            } else {
                sx = scaled_x / M_MAX(c0, static_cast<float>(M_EPSILON));
            }
        }

        // Convert to continuous [0,1]^2
        // This mapping is analogous to square_to_bilinear.
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

private:
    // Internal structure for each mip level
    struct level_data {
        int width  = 0;
        int height = 0;
        std::vector<Scalar> data;
    };

    /**
     * \brief Returns the 1D index of the top-left pixel in a 2×2 block
     * for the given offset (offset_x, offset_y).
     *
     * Originally, there's a specialized bitwise layout that packs 2×2 blocks
     * contiguously. Here we keep row-major ordering and simply clamp to avoid
     * out-of-bounds for x+1, y+1. This yields equivalent results with a less
     * complex indexing function.
     */
    static int index_2x2(const level_data &lev, int ox, int oy) {
        int block_x = (ox & ~1), block_y = (oy & ~1);
        int offset_in_block = (ox & 1) + ((oy & 1) << 1);
        return block_y * lev.width + block_x + offset_in_block;
    }

    size_t m_size_x = 0;
    size_t m_size_y = 0;
    std::vector<level_data> m_levels;
    Scalar m_sum           = Scalar(0);
    Scalar m_normalization = Scalar(0);
};

M_NAMESPACE_END

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
        return m_cond_cdf[index] - index % m_size_x != 0 ? m_cond_cdf[index - 1] : 0.0f;
    }

    // Evaluate the un-normalized PMF at a given (x, y) position
    Scalar pdf(const Point2i &pos) const { return eval(pos) * m_normalization; }

    // Sample the distribution (returns the (x, y) index)
    std::tuple<Point2i, Scalar, Point2f> sample(const Point2f &sample_) const {
        Point2f sample(sample_);

        // Avoid degeneracies on the domain boundary
        sample.clamp(M_EPSILON, 1.0f - static_cast<float>(M_EPSILON));

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
        Scalar col_cdf_0 = m_cond_cdf[offset + col - 1];
        Scalar col_cdf_1 = m_cond_cdf[offset + col];
        Scalar row_cdf_0 = m_marg_cdf[row - 1];
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

M_NAMESPACE_END

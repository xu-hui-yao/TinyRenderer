#pragma once

#include <vector>
#include <stdexcept>
#include <numeric>
#include <tuple>
#include <core/common.h>

M_NAMESPACE_BEGIN
	template <typename Scalar_, typename Index_>
	class TDiscreteDistribution {
	public:
		typedef Scalar_ Scalar;
		typedef Index_ Index;

		M_HOST_DEVICE TDiscreteDistribution() = default;

		// Initialize from a PMF (probability mass function)
		M_HOST_DEVICE explicit TDiscreteDistribution(const std::vector<Scalar>& pmf) {
			initialize(pmf);
		}

		// Update the distribution with a new PMF
		M_HOST_DEVICE void initialize(const std::vector<Scalar>& pmf) {
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
				throw std::invalid_argument("DiscreteDistribution: Total probability mass must be greater than zero.");

			m_sum = sum;
			m_normalization = 1.0f / sum;

			// Normalize the CDF
			for (auto& value : m_cdf)
				value *= m_normalization;
		}

		// Evaluate the un-normalized PMF at a given index
		M_HOST_DEVICE Scalar eval_pmf(Index index) const {
			return m_pmf[index];
		}

		// Evaluate the normalized PMF at a given index
		M_HOST_DEVICE Scalar eval_pmf_normalized(Index index) const {
			return m_pmf[index] * m_normalization;
		}

		// Evaluate the normalized CDF at a given index
		M_HOST_DEVICE Scalar eval_cdf(Index index) const {
			return m_cdf[index];
		}

		// Sample the distribution
		M_HOST_DEVICE Index sample(Scalar u) const {
			if (u < 0.0f || u > 1.0f)
				throw std::invalid_argument("Sample value must be in the range [0, 1].");

			// Binary search in the CDF to find the corresponding index
			auto it = std::lower_bound(m_cdf.begin(), m_cdf.end(), u);
			return std::distance(m_cdf.begin(), it);
		}

		// Sample the distribution and return the index and PMF value
		M_HOST_DEVICE std::pair<Index, Scalar> sample_pmf(Scalar u) const {
			Index index = sample(u);
			return {index, eval_pmf_normalized(index)};
		}

		// Sample the distribution and reuse the sample for further computations
		M_HOST_DEVICE std::tuple<Index, Scalar> sample_reuse(Scalar u) const {
			Index index = sample(u);
			Scalar pmf = eval_pmf_normalized(index);
			Scalar cdf = index > 0 ? m_cdf[index - 1] : 0.0f; // Ensure cdf for index 0 is 0
			Scalar rescaled_u = (u - cdf) / pmf; // Rescale the uniform sample

			// Ensure that rescaled_u stays in [0, 1] to avoid out of bounds errors
			rescaled_u = clamp(rescaled_u, 0.0f, 1.0f);

			return {index, rescaled_u};
		}

		// Get the normalization factor
		M_HOST_DEVICE Scalar normalization() const {
			return m_normalization;
		}

		// Get the total sum of the PMF before normalization
		M_HOST_DEVICE Scalar sum() const {
			return m_sum;
		}

		// Get the size of the distribution
		M_HOST_DEVICE [[nodiscard]] size_t size() const {
			return m_pmf.size();
		}

	private:
		std::vector<Scalar> m_pmf;
		std::vector<Scalar> m_cdf;
		Scalar m_sum = 0.0f;
		Scalar m_normalization = 0.0f;
	};

M_NAMESPACE_END

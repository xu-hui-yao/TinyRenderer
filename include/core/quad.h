#pragma once

#include <cmath>
#include <core/common.h>
#include <map>
#include <mutex>
#include <stdexcept>
#include <utility>
#include <vector>

M_NAMESPACE_BEGIN

template <typename Scalar> std::pair<Scalar, Scalar> legendre_pd(int l, Scalar x) {
    Scalar l_cur = Scalar(0), d_cur = Scalar(0);

    if (l > 1) {
        Scalar l_p_pred = Scalar(1), l_pred = x;
        Scalar d_p_pred = Scalar(0), d_pred = Scalar(1);
        Scalar k0 = Scalar(3), k1 = Scalar(2), k2 = Scalar(1);

        for (int ki = 2; ki <= l; ++ki) {
            l_cur = (k0 * x * l_pred - k2 * l_p_pred) / k1;
            d_cur = d_p_pred + k0 * l_pred;

            l_p_pred = l_pred;
            l_pred   = l_cur;
            d_p_pred = d_pred;
            d_pred   = d_cur;

            k2 = k1;
            k0 += Scalar(2);
            k1 += Scalar(1);
        }
    } else {
        if (l == 0) {
            l_cur = Scalar(1);
            d_cur = Scalar(0);
        } else {
            l_cur = x;
            d_cur = Scalar(1);
        }
    }

    return { l_cur, d_cur };
}

/**
 * \brief Computes the nodes and weights of a Gauss-Legendre quadrature
 * (aka "Gaussian quadrature") rule with the given number of evaluations.
 *
 * Integration is over the interval \f$[-1, 1]\f$. Gauss-Legendre quadrature
 * maximizes the order of exactly integrable polynomials achieves this up to
 * degree \f$2n-1\f$ (where \f$n\f$ is the number of function evaluations).
 *
 * This method is numerically well-behaved until about \f$n=200\f$
 * and then becomes progressively less accurate. It is generally not a
 * good idea to go much higher---in any case, a composite or
 * adaptive integration scheme will be superior for large \f$n\f$.
 *
 * \param n
 *     Desired number of evaluation points
 *
 * \return
 *     A tuple (nodes, weights) storing the nodes and weights of the
 *     quadrature rule.
 *
 * \remark
 *     This is the raw, uncached implementation. Prefer gauss_legendre()
 *     below, which memoizes the result - the rule depends on nothing but
 *     `n`, yet computing it costs a Newton solve (up to 20 iterations of an
 *     O(n) Legendre recurrence) per node plus two heap allocations. Call
 *     sites like TMicrofacetDistribution::eval_reflectance() ask for the
 *     same handful of `n` values thousands of times per scene, so
 *     recomputing it was pure waste.
 */
template <typename Scalar> std::pair<std::vector<Scalar>, std::vector<Scalar>> gauss_legendre_compute(int n) {
    if (n < 1) {
        throw std::invalid_argument("gauss_legendre: n must be >= 1");
    }

    std::vector<Scalar> nodes(n), weights(n);

    n--;

    if (n == 0) {
        nodes[0]   = Scalar(0);
        weights[0] = Scalar(2);
    } else if (n == 1) {
        nodes[0]   = -std::sqrt(Scalar(1.0) / Scalar(3.0));
        nodes[1]   = -nodes[0];
        weights[0] = weights[1] = Scalar(1);
    } else {
        int m = (n + 1) / 2;
        for (int i = 0; i < m; ++i) {
            Scalar x = -std::cos(Scalar(2 * i + 1) / Scalar(2 * n + 2) * Scalar(M_PI));
            int it   = 0;

            while (true) {
                if (++it > 20) {
                    throw std::runtime_error("gauss_legendre: did not converge after 20 iterations");
                }

                auto [l_val, l_der] = legendre_pd(n + 1, x);
                Scalar step         = l_val / l_der;
                x -= step;

                if (std::abs(step) <= Scalar(4) * abs(x) * M_EPSILON) {
                    break;
                }
            }

            auto [l_val, l_der] = legendre_pd(n + 1, x);
            weights[i] = weights[n - i] = Scalar(2) / ((Scalar(1) - x * x) * (l_der * l_der));
            nodes[i]                    = x;
            nodes[n - i]                = -x;
        }

        if (n % 2 == 0) {
            auto [l_val, l_der] = legendre_pd(n + 1, Scalar(0));
            weights[n / 2]      = Scalar(2) / (l_der * l_der);
            nodes[n / 2]        = Scalar(0);
        }
    }

    return { nodes, weights };
}

/**
 * \brief Memoized front end for gauss_legendre_compute().
 *
 * Returns a reference into a process-wide, mutex-guarded cache keyed by `n`,
 * so a given rule is computed at most once per (Scalar, n) pair. A std::map
 * is used deliberately: it is node-based, so previously handed-out
 * references stay valid no matter how many further entries other threads
 * insert afterwards.
 *
 * Bind the result with `const auto &[nodes, weights]` to avoid copying the
 * two vectors out of the cache.
 */
template <typename Scalar> const std::pair<std::vector<Scalar>, std::vector<Scalar>> &gauss_legendre(int n) {
    static std::mutex mutex;
    static std::map<int, std::pair<std::vector<Scalar>, std::vector<Scalar>>> cache;

    std::lock_guard<std::mutex> lock(mutex);
    auto it = cache.find(n);
    if (it == cache.end()) {
        it = cache.emplace(n, gauss_legendre_compute<Scalar>(n)).first;
    }
    return it->second;
}

M_NAMESPACE_END

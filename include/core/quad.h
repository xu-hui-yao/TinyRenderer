#pragma once

#include <cmath>
#include <core/common.h>
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
 */
template <typename Scalar> std::pair<std::vector<Scalar>, std::vector<Scalar>> gauss_legendre(int n) {
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

M_NAMESPACE_END

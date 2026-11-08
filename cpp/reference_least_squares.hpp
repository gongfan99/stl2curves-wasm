// Unbounded linear-loss TRF path corresponding to scipy.optimize.least_squares.
// Derived from SciPy's optimize/_lsq/trf.py and common.py (BSD-3-Clause).
// Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
// License: ../docs/licenses/SCIPY_BSD.txt
#pragma once
#include "internal.hpp"
#include <Eigen/SVD>
#include <limits>
namespace stl2curves::detail::reference_features {
inline Fit trust_region_fit(const std::function<Eigen::VectorXd(const Eigen::VectorXd &)> &fun,
                            Eigen::VectorXd x) {
    auto f = fun(x);
    const unsigned n = unsigned(x.size()), m = unsigned(f.size()), maximum = 100 * n;
    auto jacobian = [&](const Eigen::VectorXd &at, const Eigen::VectorXd &value) {
        Eigen::MatrixXd J(m, n);
        for (unsigned j = 0; j < n; ++j) {
            double h = std::sqrt(std::numeric_limits<double>::epsilon()) * (at[j] >= 0 ? 1 : -1) *
                       std::max(1., std::abs(at[j]));
            auto shifted = at;
            shifted[j] += h;
            double dx = shifted[j] - at[j];
            J.col(j) = (fun(shifted) - value) / dx;
        }
        return J;
    };
    auto J = jacobian(x, f);
    double cost = .5 * f.squaredNorm(), delta = x.norm(), alpha = 0;
    if (delta == 0)
        delta = 1;
    unsigned evaluations = 1;
    bool terminated = false;
    while (evaluations < maximum) {
        Eigen::VectorXd gradient = J.transpose() * f;
        if (gradient.lpNorm<Eigen::Infinity>() < 1e-8) {
            terminated = true;
            break;
        }
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(J, Eigen::ComputeThinU | Eigen::ComputeThinV);
        const auto singular = svd.singularValues();
        Eigen::VectorXd uf = svd.matrixU().transpose() * f;
        Eigen::VectorXd suf = singular.array() * uf.array();
        bool full_rank = m >= n && singular[singular.size() - 1] >
                                       std::numeric_limits<double>::epsilon() * m * singular[0];
        auto phi_at = [&](double a) {
            Eigen::ArrayXd denom = singular.array().square() + a;
            double length = (suf.array() / denom).matrix().norm();
            double prime = -(suf.array().square() / denom.cube()).sum() / length;
            return std::pair<double, double>{length - delta, prime};
        };
        auto trust_step = [&]() -> Eigen::VectorXd {
            if (full_rank) {
                Eigen::VectorXd step = -svd.matrixV() * (uf.array() / singular.array()).matrix();
                if (step.norm() <= delta) {
                    alpha = 0;
                    return step;
                }
            }
            double upper = suf.norm() / delta, lower = 0;
            if (full_rank) {
                auto [phi, prime] = phi_at(0);
                lower = -phi / prime;
            }
            if (!full_rank && alpha == 0)
                alpha = std::max(.001 * upper, std::sqrt(lower * upper));
            for (unsigned iteration = 0; iteration < 10; ++iteration) {
                if (alpha < lower || alpha > upper)
                    alpha = std::max(.001 * upper, std::sqrt(lower * upper));
                auto [phi, prime] = phi_at(alpha);
                if (phi < 0)
                    upper = alpha;
                double ratio = phi / prime;
                lower = std::max(lower, alpha - ratio);
                alpha -= (phi + delta) * ratio / delta;
                if (std::abs(phi) < .01 * delta)
                    break;
            }
            Eigen::VectorXd step =
                -svd.matrixV() * (suf.array() / (singular.array().square() + alpha)).matrix();
            step *= delta / step.norm();
            return step;
        };
        double actual = -1, new_cost = cost;
        Eigen::VectorXd next = x, new_f = f;
        while (actual <= 0 && evaluations < maximum) {
            Eigen::VectorXd step = trust_step();
            double predicted = -(.5 * (J * step).squaredNorm() + gradient.dot(step));
            next = x + step;
            new_f = fun(next);
            ++evaluations;
            double length = step.norm();
            if (!new_f.allFinite()) {
                delta = .25 * length;
                continue;
            }
            new_cost = .5 * new_f.squaredNorm();
            actual = cost - new_cost;
            double ratio =
                predicted > 0 ? actual / predicted : (predicted == 0 && actual == 0 ? 1 : 0);
            double new_delta = delta;
            if (ratio < .25)
                new_delta = .25 * length;
            else if (ratio > .75 && length > .95 * delta)
                new_delta *= 2;
            if ((actual < 1e-8 * cost && ratio > .25) || length < 1e-8 * (1e-8 + x.norm())) {
                terminated = true;
                break;
            }
            alpha *= delta / new_delta;
            delta = new_delta;
        }
        if (actual > 0) {
            x = next;
            f = new_f;
            cost = new_cost;
            J = jacobian(x, f);
        }
        if (terminated)
            break;
    }
    return {x, cost, terminated};
}
} // namespace stl2curves::detail::reference_features

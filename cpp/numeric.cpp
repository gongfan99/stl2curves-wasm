#include "internal.hpp"
#include <Eigen/Cholesky>
#include <atomic>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <mutex>
#include <numeric>
#include <set>
#include <thread>

namespace stl2curves::detail {
void trace(const char *stage) {
#ifndef __EMSCRIPTEN__
    if (std::getenv("STL2CURVES_TRACE"))
        std::cerr << "[stl2curves] " << stage << std::endl;
#else
    (void)stage;
#endif
}
unsigned worker_count() {
#if STL2CURVES_THREADS
    const auto n = std::thread::hardware_concurrency();
    return std::min(8u, std::max(1u, n > 2 ? n - 2 : 1));
#else
    return 1;
#endif
}
unsigned parallel_for(std::size_t count, const std::function<void(std::size_t)> &task) {
    const unsigned n = std::min<std::size_t>(count, worker_count());
    if (n <= 1 || count < 4) {
        for (std::size_t i = 0; i < count; ++i)
            task(i);
        return 0;
    }
#if STL2CURVES_THREADS
    std::atomic<std::size_t> next{0};
    std::atomic<unsigned> background{0};
    std::atomic<bool> stop{false};
    std::exception_ptr error;
    std::mutex mutex;
    auto run = [&](bool worker) {
        try {
            while (!stop) {
                auto i = next.fetch_add(1);
                if (i >= count)
                    break;
                if (worker)
                    ++background;
                task(i);
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(mutex);
            if (!error)
                error = std::current_exception();
            stop = true;
        }
    };
    std::vector<std::thread> threads;
    // Joining on creation failure prevents std::terminate during stack unwinding.
    try {
        for (unsigned i = 1; i < n; ++i)
            threads.emplace_back([&] { run(true); });
    } catch (...) {
        stop = true;
        for (auto &t : threads)
            t.join();
        throw;
    }
    run(false);
    for (auto &t : threads)
        t.join();
    if (error)
        std::rethrow_exception(error);
    return background.load();
#else
    return 0;
#endif
}
std::vector<unsigned> vertices(const Mesh &mesh, const std::vector<unsigned> &triangles) {
    std::vector<unsigned> ids;
    ids.reserve(triangles.size() * 3);
    for (auto t : triangles)
        for (auto v : mesh.triangles[t])
            ids.push_back(v);
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return ids;
}
std::vector<std::vector<unsigned>> components(const std::vector<std::vector<unsigned>> &adjacency) {
    std::vector<std::vector<unsigned>> out;
    std::vector<bool> seen(adjacency.size(), false);
    for (unsigned i = 0; i < adjacency.size(); ++i) {
        if (seen[i])
            continue;
        out.push_back({i});
        seen[i] = true;
        for (std::size_t k = 0; k < out.back().size(); ++k)
            for (auto j : adjacency[out.back()[k]]) {
                if (!seen[j]) {
                    seen[j] = true;
                    out.back().push_back(j);
                }
            }
        std::sort(out.back().begin(), out.back().end());
    }
    return out;
}
Fit least_squares(const std::function<Eigen::VectorXd(const Eigen::VectorXd &)> &f,
                  Eigen::VectorXd x, unsigned iterations, double scale) {
    auto cost = [&](const Eigen::VectorXd &r) {
        return (2 * scale * scale * (((r.array() / scale).square() + 1).sqrt() - 1)).sum();
    };
    auto r = f(x);
    double c = cost(r), damping = 1e-3;
    if (!r.allFinite() || !std::isfinite(c))
        return {x, c, false};
    for (unsigned it = 0; it < iterations; ++it) {
        Eigen::MatrixXd J(r.size(), x.size());
        for (Eigen::Index j = 0; j < x.size(); ++j) {
            auto shifted = x;
            double h = 1.4901161193847656e-8 * std::max(1.0, std::abs(x[j]));
            shifted[j] += h;
            J.col(j) = (f(shifted) - r) / h;
        }
        Eigen::VectorXd w = (1 + (r.array() / scale).square()).sqrt().inverse().matrix();
        Eigen::MatrixXd A = J.transpose() * w.asDiagonal() * J;
        Eigen::VectorXd g = J.transpose() * w.asDiagonal() * r;
        if (g.lpNorm<Eigen::Infinity>() < 1e-10)
            return {x, c, true};
        A.diagonal().array() += damping * (A.diagonal().array().abs() + 1e-12);
        Eigen::VectorXd step = A.ldlt().solve(-g);
        if (!step.allFinite())
            break;
        Eigen::VectorXd candidate = x + step, rr = f(candidate);
        double cc = cost(rr);
        if (std::isfinite(cc) && cc < c) {
            x = candidate;
            r = rr;
            if (std::abs(c - cc) < 1e-13 * (1 + c) || step.norm() < 1e-10 * (1 + x.norm()))
                return {x, cc, true};
            c = cc;
            damping = std::max(1e-12, damping / 3);
        } else
            damping = std::min(1e12, damping * 10);
    }
    return {x, c, false};
}
std::string safe_name(std::string name) {
    for (auto &c : name)
        if (std::string("<>:\"/\\|?*").find(c) != std::string::npos || (unsigned char)c < 32)
            c = '_';
    while (!name.empty() && (name.back() == ' ' || name.back() == '.'))
        name.pop_back();
    return name.empty() ? "object" : name;
}
} // namespace stl2curves::detail

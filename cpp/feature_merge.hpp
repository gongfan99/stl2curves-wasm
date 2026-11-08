// Derived from stl2curves, Copyright (c) 2026 Duvdivan. MIT license.
#pragma once
#include "feature_geometry.hpp"
namespace stl2curves::detail::reference_features {
inline bool strip_of(const Model &a, const Model &b) {
    if (a.kind == "sphere" && !a.sphere_cap && b.kind == "sphere" && !b.sphere_cap)
        return std::abs(a.r - b.r) <= .15 * a.r && (a.a - b.a).norm() <= .15 * a.r;
    if ((a.kind == "sphere" && !a.sphere_cap) || (b.kind == "sphere" && !b.sphere_cap))
        return false;
    if (std::abs(a.d.dot(b.d)) < std::cos(10 * pi / 180))
        return false;
    if (a.kind == "cylinder" && b.kind == "cylinder")
        return std::abs(a.r - b.r) <= .15 * a.r;
    if ((a.kind == "torus" || a.sphere_cap) && (b.kind == "torus" || b.sphere_cap)) {
        double ar = a.sphere_cap ? a.r : a.tube, br = b.sphere_cap ? b.r : b.tube;
        double ac = a.sphere_cap ? 0 : a.r, bc = b.sphere_cap ? 0 : b.r;
        return std::abs(ar - br) <= .15 * ar && std::abs(ac - bc) <= .15 * std::max(ac, ar);
    }
    return false;
}
inline bool holds(const Pass &pass, const Model &model, const std::vector<unsigned> &facets,
                  double tol, bool convex) {
    if (!fits(model, points(pass.m, pass.fs, facets), tol))
        return false;
    for (auto f : facets)
        if ((convex ? 1 : -1) * pass.fs.normals[f].dot(model.normal(pass.fs.centers[f])) <
            std::cos(pi / 6))
            return false;
    return true;
}
inline bool simple_outline(const Pass &pass, const std::vector<unsigned> &facets) {
    std::map<Edge, unsigned> edges;
    for (auto f : facets)
        for (auto t : pass.fs.triangles[f]) {
            auto tri = pass.m.triangles[t];
            for (unsigned j = 0; j < 3; ++j)
                ++edges[Edge(tri[j], tri[(j + 1) % 3])];
        }
    std::map<unsigned, unsigned> visits;
    for (auto &e : edges)
        if (e.second == 1) {
            if (++visits[e.first.a] > 2 || ++visits[e.first.b] > 2)
                return false;
        }
    return true;
}
std::vector<unsigned> vertex_ids(const Pass &pass, const Found &feature) {
    std::set<unsigned> v;
    for (auto f : feature.facets)
        v.insert(pass.fs.verts[f].begin(), pass.fs.verts[f].end());
    return {v.begin(), v.end()};
}
inline double range(std::vector<double> v) {
    if (v.empty())
        return 0;
    auto a = std::minmax_element(v.begin(), v.end());
    return *a.second - *a.first;
}
inline double unwrapped_range(std::vector<double> angles) {
    std::sort(angles.begin(), angles.end());
    if (angles.empty())
        return 0;
    double previous = angles[0], at = previous, lo = at, hi = at;
    for (unsigned i = 1; i < angles.size(); ++i) {
        double step = std::remainder(angles[i] - previous, 2 * pi);
        at += step;
        lo = std::min(lo, at);
        hi = std::max(hi, at);
        previous = angles[i];
    }
    return hi - lo;
}
inline Found joint_surface(const Pass &pass, const Found &a, const Found &b) {
    if (a.convex != b.convex || (a.model.kind == "sphere" && !a.model.sphere_cap) ||
        (b.model.kind == "sphere" && !b.model.sphere_cap) || !strip_of(a.model, b.model))
        return {};
    const auto &big = a.facets.size() >= b.facets.size() ? a : b;
    std::vector<unsigned> fids = a.facets;
    fids.insert(fids.end(), b.facets.begin(), b.facets.end());
    unsigned stray = std::max(1u, unsigned(.25 * std::min(a.facets.size(), b.facets.size())));
    auto fit = [&](const std::vector<unsigned> &ids) {
        return tube_fit(points(pass.m, pass.fs, ids), big.model);
    };
    auto model = fit(fids);
    if (!model || !strip_of(big.model, model))
        return {};
    std::vector<unsigned> keep;
    for (auto f : fids)
        if (fits(model, points(pass.m, pass.fs, {f}), pass.fs.tolerance))
            keep.push_back(f);
    if (fids.size() - keep.size() > stray)
        return {};
    if (keep.size() < fids.size()) {
        model = fit(keep);
        if (!model || !strip_of(big.model, model))
            return {};
        for (auto f : keep)
            if (!fits(model, points(pass.m, pass.fs, {f}), pass.fs.tolerance))
                return {};
    }
    if (!holds(pass, model, keep, pass.fs.tolerance, a.convex))
        return {};
    auto P = points(pass.m, pass.fs, keep);
    std::vector<double> rho, z, angles, profile;
    double sum_rho = 0;
    for (auto &p : P) {
        Vec q = p - model.a;
        double zz = q.dot(model.d);
        Vec w = q - zz * model.d;
        double rr = w.norm();
        rho.push_back(rr);
        z.push_back(zz);
        sum_rho += rr;
        angles.push_back(std::atan2(w.dot(model.e2), w.dot(model.e1)));
        profile.push_back(std::atan2(zz - model.zc, rr - (model.sphere_cap ? 0 : model.r)));
    }
    double radius = model.sphere_cap ? model.r : model.tube;
    double reach;
    if (model.kind == "cylinder" || model.kind == "cone") {
        radius = model.r;
        reach = std::max(range(z), sum_rho / P.size() * std::min(unwrapped_range(angles), 2 * pi));
    } else
        reach =
            std::max((model.sphere_cap ? 0 : model.r) * std::min(unwrapped_range(angles), 2 * pi),
                     radius * unwrapped_range(profile));
    if (reach < radius)
        return {};
    if (!simple_outline(pass, keep))
        return {};
    auto result = pass.describe(model, keep, true);
    if (result)
        result.convex = a.convex;
    return result;
}
inline std::vector<Found> merge_same_surface(const Pass &pass, std::vector<Found> features) {
    std::vector<std::set<unsigned>> verts;
    for (auto &f : features) {
        auto vs = vertex_ids(pass, f);
        verts.emplace_back(vs.begin(), vs.end());
    }
    auto near = [&](unsigned i, bool before) {
        std::vector<unsigned> out;
        for (unsigned j = 0; j < features.size(); ++j) {
            if (!features[j] || ((j < i) != before) || j == i)
                continue;
            bool share = false;
            for (auto v : verts[j])
                if (verts[i].count(v)) {
                    share = true;
                    break;
                }
            if (share)
                out.push_back(j);
        }
        return out;
    };
    auto sure = [](const Found &f) { return f.facets.size() >= 30 || f.rolled; };
    auto merge = [&](unsigned i, unsigned j) {
        const auto a = features[i], b = features[j];
        std::vector<std::pair<const Found *, const Found *>> order{{&a, &b}, {&b, &a}};
        if (b.facets.size() > a.facets.size())
            std::swap(order[0], order[1]);
        const Found *holder = nullptr;
        for (auto pair : order) {
            const auto &x = *pair.first;
            const auto &y = *pair.second;
            bool ok =
                sure(x) && strip_of(x.model, y.model)
                    ? holds(pass, x.model, y.facets, std::max(.015, pass.fs.tolerance), x.convex)
                    : x.convex == y.convex &&
                          holds(pass, x.model, y.facets, pass.fs.tolerance, x.convex);
            if (ok) {
                holder = &x;
                break;
            }
        }
        Found both;
        if (holder) {
            auto ids = a.facets;
            ids.insert(ids.end(), b.facets.begin(), b.facets.end());
            if (simple_outline(pass, ids)) {
                const auto &other = holder == &a ? b : a;
                bool anchored = sure(*holder) || (other.convex == holder->convex &&
                                                  holds(pass, other.model, holder->facets,
                                                        pass.fs.tolerance, other.convex));
                both = pass.describe(holder->model, ids, anchored, !anchored);
                if (both)
                    both.convex = holder->convex;
            }
        }
        if (!both)
            both = joint_surface(pass, a, b);
        if (!both)
            return false;
        auto vs = vertex_ids(pass, both);
        verts[i] = {vs.begin(), vs.end()};
        verts[j].clear();
        features[i] = std::move(both);
        features[j] = {};
        return true;
    };
    unsigned i = 0;
    while (i < features.size()) {
        if (!features[i]) {
            ++i;
            continue;
        }
        bool joined = false;
        for (auto j : near(i, false))
            if (merge(i, j)) {
                bool again = true;
                while (again) {
                    again = false;
                    for (auto k : near(i, true))
                        if (merge(k, i)) {
                            i = k;
                            again = true;
                            break;
                        }
                }
                joined = true;
                break;
            }
        if (!joined)
            ++i;
    }
    features.erase(
        std::remove_if(features.begin(), features.end(), [](const Found &f) { return !f; }),
        features.end());
    return features;
}
} // namespace stl2curves::detail::reference_features

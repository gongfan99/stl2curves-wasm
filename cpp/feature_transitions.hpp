// Derived from stl2curves, Copyright (c) 2026 Duvdivan. MIT license.
#pragma once
#include "feature_cleanup.hpp"
namespace stl2curves::detail::reference_features {
struct TransitionArea {
    std::vector<unsigned> facets, parts;
};
inline double transition_radius(const Pass &pass, const Found &feature) {
    const auto &model = feature.model;
    if (model.kind == "cone") {
        double r = 0;
        for (auto f : feature.facets) {
            Vec q = pass.fs.centers[f] - model.a;
            q -= q.dot(model.d) * model.d;
            r += norm3(q);
        }
        return r / feature.facets.size();
    }
    return radius(model);
}
inline std::vector<TransitionArea>
transition_areas(const Pass &pass, const std::vector<Found> &features, double fit_dev = .01) {
    auto owner = owners(pass, features);
    const auto &mesh = pass.m;
    const auto &fs = pass.fs;
    const double small = .02 * (*std::max_element(fs.areas.begin(), fs.areas.end()));
    std::vector<bool> member(owner.size(), false), seen(owner.size(), false);
    auto unit = [&](unsigned f) {
        return owner[f] >= 0 ? features[owner[f]].facets : std::vector<unsigned>{f};
    };
    auto bounds = [&](const std::vector<unsigned> &facets) {
        auto P = points(mesh, fs, facets);
        Vec lo = P.front(), hi = lo;
        for (auto &p : P) {
            lo = lo.cwiseMin(p);
            hi = hi.cwiseMax(p);
        }
        return Vec(hi - lo);
    };
    auto minor = [&](unsigned k) {
        const auto &f = features[k];
        if ((f.model.kind == "sphere" && !f.model.sphere_cap && !f.trimmed) || f.vouched ||
            f.has_parts)
            return false;
        double r = transition_radius(pass, f);
        if (r <= 0 || r > 10)
            return false;
        if (f.facets.size() <= 8)
            return true;
        return norm3(bounds(f.facets)) <= 3 * r;
    };
    for (unsigned k = 0; k < features.size(); ++k)
        if (minor(k))
            for (auto f : features[k].facets)
                member[f] = true;
    for (unsigned f = 0; f < owner.size(); ++f)
        if (owner[f] < 0 && fs.areas[f] <= small && !fs.neighbors[f].empty())
            member[f] = true;
    auto parts_of = [&](const std::vector<unsigned> &facets) {
        std::set<unsigned> parts;
        for (auto f : facets)
            if (owner[f] >= 0)
                parts.insert(owner[f]);
        return std::vector<unsigned>(parts.begin(), parts.end());
    };
    std::vector<TransitionArea> out;
    for (unsigned seed = 0; seed < owner.size(); ++seed) {
        if (!member[seed] || owner[seed] < 0 || seen[seed])
            continue;
        auto area = unit(seed);
        for (auto f : area)
            seen[f] = true;
        auto stack = area;
        while (!stack.empty()) {
            auto a = stack.back();
            stack.pop_back();
            for (auto b : fs.neighbors[a]) {
                if (seen[b] || !member[b] ||
                    fs.normals[a].dot(fs.normals[b]) < std::cos(20 * pi / 180))
                    continue;
                auto u = unit(b);
                for (auto f : u)
                    seen[f] = true;
                area.insert(area.end(), u.begin(), u.end());
                stack.insert(stack.end(), u.begin(), u.end());
            }
        }
        std::sort(area.begin(), area.end());
        area.erase(std::unique(area.begin(), area.end()), area.end());
        auto parts = parts_of(area);
        if (parts.size() < 2)
            continue;
        std::set<unsigned> border;
        for (auto a : area)
            for (auto b : fs.neighbors[a])
                if (owner[b] >= 0 && !member[b])
                    border.insert(owner[b]);
        std::vector<double> radii;
        for (auto k : border) {
            const auto &m = features[k].model;
            double r = (m.kind == "sphere" && !m.sphere_cap) ? 0 : radius(m);
            if (r > 0 && r <= 10)
                radii.push_back(r);
        }
        if (radii.size() < 2)
            continue;
        auto rr = std::minmax_element(radii.begin(), radii.end());
        double lo = *rr.first, hi = *rr.second;
        std::set<unsigned> long_facets;
        for (auto a : area)
            if (owner[a] < 0 && bounds({a}).maxCoeff() > 2 * hi)
                long_facets.insert(a);
        if (!long_facets.empty()) {
            std::set<unsigned> left;
            for (auto a : area)
                if (!long_facets.count(a))
                    left.insert(a);
            std::vector<std::vector<unsigned>> groups;
            while (!left.empty()) {
                std::set<unsigned> group;
                for (auto f : unit(*left.begin()))
                    if (left.count(f))
                        group.insert(f);
                for (auto f : group)
                    left.erase(f);
                std::vector<unsigned> front(group.begin(), group.end());
                while (!front.empty()) {
                    auto a = front.back();
                    front.pop_back();
                    for (auto b : fs.neighbors[a])
                        if (left.count(b))
                            for (auto f : unit(b))
                                if (left.erase(f)) {
                                    group.insert(f);
                                    front.push_back(f);
                                }
                }
                groups.emplace_back(group.begin(), group.end());
            }
            if (groups.empty())
                continue;
            auto score = [&](const std::vector<unsigned> &g) {
                std::size_t owned = 0;
                for (auto f : g)
                    owned += owner[f] >= 0;
                return std::make_pair(owned, g.size());
            };
            area =
                *std::max_element(groups.begin(), groups.end(), [&](const auto &a, const auto &b) {
                    return score(a) < score(b);
                });
            parts = parts_of(area);
            if (parts.size() < 2)
                continue;
        }
        std::size_t weight = 0, odd = 0;
        for (auto k : parts) {
            weight += features[k].facets.size();
            double r = transition_radius(pass, features[k]);
            if (r < lo / 2 || r > 2 * hi)
                odd += features[k].facets.size();
        }
        if (odd > .2 * weight)
            continue;
        double size = norm3(bounds(area));
        if (size < lo || size > 6 * hi)
            continue;
        auto P = points(mesh, fs, area);
        bool exact = false;
        for (auto k : parts)
            if (fits(features[k].model, P, fit_dev))
                exact = true;
        if (!exact)
            out.push_back({std::move(area), std::move(parts)});
    }
    return out;
}
} // namespace stl2curves::detail::reference_features

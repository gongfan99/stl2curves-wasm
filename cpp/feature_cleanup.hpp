// Derived from stl2curves, Copyright (c) 2026 Duvdivan. MIT license.
#pragma once
#include "feature_merge.hpp"
namespace stl2curves::detail::reference_features {
inline std::vector<int> owners(const Pass &pass, const std::vector<Found> &features) {
    std::vector<int> owner(pass.fs.triangles.size(), -1);
    for (unsigned k = 0; k < features.size(); ++k)
        for (auto f : features[k].facets)
            owner[f] = int(k);
    return owner;
}
inline std::vector<Found> fill_pinholes(const Pass &pass, std::vector<Found> features) {
    auto owner = owners(pass, features);
    std::vector<bool> seen(owner.size(), false);
    std::map<unsigned, std::vector<unsigned>> grown;
    for (unsigned seed = 0; seed < owner.size(); ++seed) {
        if (owner[seed] >= 0 || seen[seed])
            continue;
        std::vector<unsigned> group, stack{seed};
        std::set<unsigned> around;
        seen[seed] = true;
        while (!stack.empty() && group.size() <= 4) {
            auto g = stack.back();
            stack.pop_back();
            group.push_back(g);
            for (auto h : pass.fs.neighbors[g])
                if (owner[h] >= 0)
                    around.insert(owner[h]);
                else if (!seen[h]) {
                    seen[h] = true;
                    stack.push_back(h);
                }
        }
        if (!stack.empty() || group.size() > 4 || around.size() != 1)
            continue;
        unsigned k = *around.begin();
        auto vs = vertex_ids(pass, features[k]);
        std::set<unsigned> corners(vs.begin(), vs.end());
        bool on = true;
        for (auto g : group) {
            for (auto v : pass.fs.verts[g])
                if (!corners.count(v))
                    on = false;
            if (std::abs(features[k].model.signed_at(pass.fs.centers[g])) >
                std::max(features[k].worst, pass.fs.tolerance))
                on = false;
        }
        if (on)
            grown[k].insert(grown[k].end(), group.begin(), group.end());
    }
    for (auto &entry : grown) {
        auto ids = features[entry.first].facets;
        ids.insert(ids.end(), entry.second.begin(), entry.second.end());
        auto f = pass.describe(features[entry.first].model, ids, true);
        if (f) {
            f.convex = features[entry.first].convex;
            features[entry.first] = std::move(f);
        }
    }
    return features;
}
inline std::vector<Found> absorb_straddlers(const Pass &pass, std::vector<Found> features) {
    auto owner = owners(pass, features);
    bool joined = true, again = false, first = true;
    while (joined) {
        joined = false;
        again = !first;
        first = false;
        for (unsigned g = 0; g < owner.size(); ++g) {
            if (owner[g] >= 0)
                continue;
            auto P = points(pass.m, pass.fs, {g});
            Vec lo = P.front(), hi = lo;
            for (auto &p : P) {
                lo = lo.cwiseMin(p);
                hi = hi.cwiseMax(p);
            }
            bool needle = pass.fs.areas[g] < .05 * std::pow((hi - lo).maxCoeff(), 2);
            if (pass.fs.triangles[g].size() < 2 && !needle)
                continue;
            std::set<unsigned> beside;
            for (auto h : pass.fs.neighbors[g])
                if (owner[h] >= 0)
                    beside.insert(owner[h]);
            for (auto k : beside) {
                const auto &f = features[k];
                if (!fits(f.model, P, 5 * pass.fs.tolerance))
                    continue;
                if (again) {
                    auto vs = vertex_ids(pass, f);
                    std::set<unsigned> corners(vs.begin(), vs.end());
                    bool inside = true;
                    for (auto v : pass.fs.verts[g])
                        if (!corners.count(v))
                            inside = false;
                    if (!inside)
                        continue;
                }
                if (!needle && std::abs(pass.fs.normals[g].dot(
                                   f.model.normal(pass.fs.centers[g]))) < std::cos(pi / 6))
                    continue;
                bool near = true;
                double limit = std::max(1.5 * f.worst, 5 * pass.fs.tolerance);
                for (auto t : pass.fs.triangles[g]) {
                    auto tri = pass.m.triangles[t];
                    if (std::abs(f.model.signed_at(pass.top.centers[t])) > limit)
                        near = false;
                    for (unsigned j = 0; j < 3; ++j)
                        if (std::abs(f.model.signed_at(
                                (pass.m.points[tri[j]] + pass.m.points[tri[(j + 1) % 3]]) * .5)) >
                            limit)
                            near = false;
                }
                if (!near)
                    continue;
                auto ids = f.facets;
                ids.push_back(g);
                auto wider = pass.describe(f.model, ids, false, true);
                if (!wider)
                    wider = pass.describe(f.model, ids, true);
                if (wider) {
                    wider.convex = f.convex;
                    features[k] = std::move(wider);
                    owner[g] = k;
                    joined = true;
                    break;
                }
            }
        }
    }
    return features;
}
inline Eigen::Vector2d circle_fixed_r(const std::vector<Eigen::Vector2d> &p, double radius,
                                      Eigen::Vector2d c) {
    for (unsigned it = 0; it < 30; ++it) {
        Eigen::MatrixXd J(p.size(), 2);
        Eigen::VectorXd b(p.size());
        for (unsigned i = 0; i < p.size(); ++i) {
            auto q = p[i] - c;
            double d = std::max(q.norm(), 1e-12);
            J.row(i) = -q / d;
            b[i] = -(d - radius);
        }
        Eigen::Vector2d step = lstsq(J, b);
        c += step;
        if (step.cwiseAbs().maxCoeff() < 1e-13)
            break;
    }
    return c;
}
inline Model with_radius(Model model, const std::vector<Vec> &p, double radius) {
    if (model.kind == "sphere" && !model.sphere_cap) {
        for (unsigned it = 0; it < 30; ++it) {
            Eigen::MatrixXd J(p.size(), 3);
            Eigen::VectorXd b(p.size());
            for (unsigned i = 0; i < p.size(); ++i) {
                Vec q = p[i] - model.a;
                double d = std::max(q.norm(), 1e-12);
                J.row(i) = -q / d;
                b[i] = -(d - radius);
            }
            Vec step = lstsq(J, b);
            model.a += step;
            if (step.cwiseAbs().maxCoeff() < 1e-13)
                break;
        }
        model.r = radius;
        return model;
    }
    if (model.kind == "cylinder") {
        std::vector<Eigen::Vector2d> xy;
        for (auto &q : p)
            xy.emplace_back(q.dot(model.e1), q.dot(model.e2));
        auto c = circle_fixed_r(xy, radius, {model.a.dot(model.e1), model.a.dot(model.e2)});
        model.at(c.x() * model.e1 + c.y() * model.e2);
        model.r = radius;
        return model;
    }
    if (model.kind == "torus" || model.sphere_cap) {
        Vec base = model.a;
        std::vector<Eigen::Vector2d> rz;
        for (auto &q : p) {
            Vec d = q - base;
            double z = d.dot(model.d);
            rz.emplace_back((d - z * model.d).norm(), z);
        }
        double zc = model.zc, rc = model.sphere_cap ? 0 : model.r;
        if (rc == 0) {
            for (unsigned it = 0; it < 30; ++it) {
                Eigen::MatrixXd J(p.size(), 1);
                Eigen::VectorXd b(p.size());
                for (unsigned i = 0; i < p.size(); ++i) {
                    double dz = rz[i].y() - zc, d = std::max(std::hypot(rz[i].x(), dz), 1e-12);
                    J(i, 0) = -dz / d;
                    b[i] = -(d - radius);
                }
                double step = lstsq(J, b)[0];
                zc += step;
                if (std::abs(step) < 1e-13)
                    break;
            }
            model.at(base);
            model.zc = zc;
            model.r = radius;
            model.tube = radius;
            return model;
        }
        auto c = circle_fixed_r(rz, radius, {rc, zc});
        model.at(base);
        model.zc = c.y();
        model.r = c.x();
        model.tube = radius;
        return model;
    }
    return {};
}
inline double radius(const Model &model) {
    if (model.kind == "sphere" || model.kind == "cylinder")
        return model.r;
    if (model.kind == "torus")
        return model.tube;
    return 0;
}
inline Found rebuilt(const Pass &pass, const Found &f, const Model &model) {
    auto result = pass.describe(model, f.facets, false, true);
    if (!result)
        result = pass.describe(model, f.facets);
    if (result)
        result.convex = f.convex;
    return result;
}
inline std::vector<Found> snap(const Pass &pass, std::vector<Found> features) {
    double tolerance = 2 * pass.fs.tolerance;
    for (auto &f : features) {
        const auto &model = f.model;
        if (model.kind == "sphere" && !model.sphere_cap)
            continue;
        Eigen::Index k;
        model.d.cwiseAbs().maxCoeff(&k);
        Vec axis = Vec::Zero();
        axis[k] = model.d[k] > 0 ? 1 : -1;
        double dot = std::abs(model.d.dot(axis));
        if (dot >= 1 - 1e-12 || dot < std::cos(.1 * pi / 180))
            continue;
        auto P = points(pass.m, pass.fs, f.facets);
        Model new_model;
        if (model.kind == "cylinder") {
            new_model = with_radius(cylinder(axis, P), P, model.r);
        } else {
            new_model = model;
            Vec base = model.a;
            new_model.frame(axis);
            new_model.at(base);
        }
        if (fits(new_model, P, tolerance)) {
            auto changed = rebuilt(pass, f, new_model);
            if (changed)
                f = std::move(changed);
        }
    }
    std::vector<std::pair<double, unsigned>> radii;
    for (unsigned i = 0; i < features.size(); ++i) {
        double r = radius(features[i].model);
        if (r > 0)
            radii.emplace_back(r, i);
    }
    std::stable_sort(radii.begin(), radii.end(),
                     [](const auto &a, const auto &b) { return a.first < b.first; });
    std::vector<std::vector<std::pair<double, unsigned>>> groups;
    for (auto &entry : radii) {
        if (groups.empty() || entry.first - groups.back().back().first > .002 * entry.first)
            groups.push_back({});
        groups.back().push_back(entry);
    }
    for (auto &group : groups) {
        double total = 0, weighted = 0;
        for (auto &entry : group) {
            double area = 0;
            for (auto f : features[entry.second].facets)
                area += pass.fs.areas[f];
            weighted += entry.first * area;
            total += area;
        }
        double target = weighted / total;
        for (auto &entry : group) {
            if (entry.first == target)
                continue;
            auto &f = features[entry.second];
            auto P = points(pass.m, pass.fs, f.facets);
            auto model = with_radius(f.model, P, target);
            if (fits(model, P, tolerance)) {
                auto changed = rebuilt(pass, f, model);
                if (changed)
                    f = std::move(changed);
            }
        }
    }
    return absorb_straddlers(pass, std::move(features));
}
} // namespace stl2curves::detail::reference_features

// Derived from stl2curves, Copyright (c) 2026 Duvdivan. MIT license.
#pragma once
#include "feature_transitions.hpp"
#include <functional>
namespace stl2curves::detail::reference_features {
using V2 = Eigen::Vector2d;
using Ids = std::vector<unsigned>;
inline double percentile(std::vector<double> v, double q) {
    std::sort(v.begin(), v.end());
    double at = (v.size() - 1) * q;
    unsigned lo = unsigned(at), hi = std::min<unsigned>(lo + 1, v.size() - 1);
    return v[lo] + (v[hi] - v[lo]) * (at - lo);
}
inline double median(std::vector<double> v) {
    return percentile(std::move(v), .5);
}
inline double arc_span(std::vector<double> v) {
    if (v.empty())
        return 0;
    for (auto &x : v) {
        x = std::fmod(x, 2 * pi);
        if (x < 0)
            x += 2 * pi;
    }
    std::sort(v.begin(), v.end());
    double gap = v[0] + 2 * pi - v.back();
    for (unsigned i = 1; i < v.size(); ++i)
        gap = std::max(gap, v[i] - v[i - 1]);
    return 2 * pi - gap;
}
struct RollConfig {
    V2 n1, n2, p, w;
    double sign, span;
    RollConfig(V2 a, V2 b, V2 origin, double s)
        : n1(a), n2(b), p(origin), w(-s * (a + b) / (1 + a.dot(b))), sign(s),
          span(std::acos(std::clamp(a.dot(b), -1., 1.))) {}
    std::pair<double, bool> radii(V2 x) const {
        V2 q = x - p;
        double qa = w.squaredNorm() - 1, qb = -2 * q.dot(w), qc = q.squaredNorm(),
               disc = qb * qb - 4 * qa * qc, twice = std::abs(qb / (2 * qa)),
               root = std::sqrt(std::max(disc, 0.));
        double r = std::max((-qb - root) / (2 * qa), (-qb + root) / (2 * qa));
        return {r, disc < -8 * twice * std::abs(qa) * .01 || !(r > 0) || !std::isfinite(r)};
    }
    std::vector<double> both(const std::vector<V2> &X) const {
        double qa = w.squaredNorm() - 1;
        if (std::abs(qa) < 1e-12)
            return {};
        std::vector<double> lo, hi;
        for (auto x : X) {
            V2 q = x - p;
            double qb = -2 * q.dot(w), qc = q.squaredNorm(), disc = qb * qb - 4 * qa * qc;
            if (disc >= 0) {
                double root = std::sqrt(disc);
                lo.push_back((-qb - root) / (2 * qa));
                hi.push_back((-qb + root) / (2 * qa));
            }
        }
        if (lo.size() < 3)
            return {};
        lo.insert(lo.end(), hi.begin(), hi.end());
        return lo;
    }
    V2 center(double r) const { return p + r * w; }
    std::pair<double, double> heights(V2 x) const { return {(x - p).dot(n1), (x - p).dot(n2)}; }
};
struct RollSection {
    const Pass &pass;
    unsigned A, B = 0;
    bool flat = true;
    Found wall;
    Ids group_A;
    std::vector<bool> a_mask, b_mask;
    Vec origin, d, e1, e2, normal_A, normal_B;
    V2 n1, n2, p;
    double t0, t1, width_A, width_B, area_A, area_B, rough = 0, rho = 0, mid = 0;
    RollSection(const Pass &P, unsigned a, unsigned b)
        : pass(P), A(a), B(b), group_A{a}, a_mask(P.fs.areas.size(), false), b_mask(a_mask) {
        const auto &fs = P.fs;
        normal_A = fs.normals[a];
        normal_B = fs.normals[b];
        d = normalized(normal_A.cross(normal_B));
        Eigen::Matrix3d M;
        M.row(0) = normal_A;
        M.row(1) = normal_B;
        M.row(2) = d;
        origin = M.fullPivLu().solve(Vec(normal_A.dot(P.m.points[fs.verts[a][0]]),
                                         normal_B.dot(P.m.points[fs.verts[b][0]]), 0));
        e1 = normal_A;
        e2 = d.cross(e1);
        n1 = {normal_A.dot(e1), normal_A.dot(e2)};
        n2 = {normal_B.dot(e1), normal_B.dot(e2)};
        p = V2::Zero();
        auto pa = points(P.m, fs, {a}), pb = points(P.m, fs, {b});
        std::vector<double> ta, tb, wa, wb;
        for (auto q : pa) {
            ta.push_back(q.dot(d));
            wa.push_back(q.dot(normal_A.cross(d)));
        }
        for (auto q : pb) {
            tb.push_back(q.dot(d));
            wb.push_back(q.dot(normal_B.cross(d)));
        }
        t0 = std::max(*std::min_element(ta.begin(), ta.end()),
                      *std::min_element(tb.begin(), tb.end()));
        t1 = std::min(*std::max_element(ta.begin(), ta.end()),
                      *std::max_element(tb.begin(), tb.end()));
        width_A = range(wa);
        width_B = range(wb);
        area_A = fs.areas[a];
        area_B = fs.areas[b];
        a_mask[a] = b_mask[b] = true;
    }
    RollSection(const Pass &P, unsigned a, const Found &X, Ids group)
        : pass(P), A(a), flat(false), wall(X), group_A(std::move(group)),
          a_mask(P.fs.areas.size(), false), b_mask(a_mask) {
        auto face = P.m.points[P.fs.verts[a][0]];
        d = P.fs.normals[a];
        normal_A = d;
        origin = X.model.a + (face - X.model.a).dot(X.model.d) * X.model.d;
        origin -= origin.dot(d) * d;
        Model frame;
        frame.frame(d);
        e1 = frame.e1;
        e2 = frame.e2;
        auto V = points(P.m, P.fs, X.facets);
        std::vector<double> z, rr, angles;
        double sum_r = 0;
        for (auto q : V) {
            Vec v = q - origin;
            double zz = v.dot(d);
            z.push_back(zz);
            Vec radial = v - zz * d;
            double r = radial.norm();
            rr.push_back(r);
            sum_r += r;
            angles.push_back(std::atan2(radial.dot(e2), radial.dot(e1)));
        }
        double k = 0, c0 = sum_r / V.size();
        if (range(z) > 1e-9) {
            Eigen::MatrixXd M(z.size(), 2);
            Eigen::VectorXd y(z.size());
            for (unsigned i = 0; i < z.size(); ++i) {
                M.row(i) << z[i], 1;
                y[i] = rr[i];
            }
            auto fit = polyfit(M, y);
            k = fit[0];
            c0 = fit[1];
        }
        n1 = {0, 1};
        double g = (X.convex ? 1 : -1) / std::hypot(1., k);
        n2 = {g, -g * k};
        double za = (face - origin).dot(d);
        p = {c0 + k * za, za};
        rho = std::max(p[0], 1e-6);
        for (auto &x : angles) {
            x = std::fmod(x, 2 * pi);
            if (x < 0)
                x += 2 * pi;
        }
        std::sort(angles.begin(), angles.end());
        double gap = -1;
        unsigned at = 0;
        for (unsigned i = 0; i < angles.size(); ++i) {
            double diff = (i + 1 < angles.size() ? angles[i + 1] : angles[0] + 2 * pi) - angles[i];
            if (diff > gap) {
                gap = diff;
                at = i;
            }
        }
        double span = gap > 5 * pi / 180 ? 2 * pi - gap : 2 * pi;
        mid = angles[(at + 1) % angles.size()] + span / 2;
        t0 = -span / 2 * rho;
        t1 = span / 2 * rho;
        width_A = 0;
        for (auto q : points(P.m, P.fs, group_A))
            width_A = std::max(width_A, std::abs(coords(q)[0] - p[0]));
        width_B = range(z) * std::hypot(1., k);
        area_A = area_B = 0;
        for (auto id : group_A) {
            a_mask[id] = true;
            area_A += P.fs.areas[id];
        }
        for (auto id : X.facets) {
            b_mask[id] = true;
            area_B += P.fs.areas[id];
        }
        for (auto q : V)
            rough = std::max(rough, std::abs((coords(q) - p).dot(n2)));
    }
    V2 coords(const Vec &q) const {
        Vec v = q - origin;
        if (flat)
            return {v.dot(e1), v.dot(e2)};
        double z = v.dot(d);
        return {(v - z * d).norm(), z};
    }
    double along(const Vec &q) const {
        if (flat)
            return q.dot(d);
        Vec v = q - origin;
        double angle = std::atan2(v.dot(e2), v.dot(e1));
        double u = std::fmod(angle - mid + pi, 2 * pi);
        if (u < 0)
            u += 2 * pi;
        return (u - pi) * rho;
    }
    Model model(const V2 &o, double r) const {
        Model m;
        m.frame(d);
        m.kind = flat ? "cylinder" : "torus";
        m.at(flat ? Vec(origin + o[0] * e1 + o[1] * e2) : origin);
        m.zc = flat ? 0 : o[1];
        m.r = flat ? r : o[0];
        m.tube = flat ? 0 : r;
        return m;
    }
    double profile_span(const Ids &facets) const {
        std::vector<double> angles;
        for (auto id : facets) {
            Vec v = pass.fs.centers[id] - origin;
            Vec radial = v - v.dot(d) * d;
            radial = divided(radial, std::max(radial.norm(), 1e-12));
            angles.push_back(
                std::atan2(pass.fs.normals[id].dot(d), pass.fs.normals[id].dot(radial)));
        }
        return arc_span(angles);
    }
};
struct FilletPort {
    const Pass &pass;
    const Mesh &m;
    const FeatureFacets &fs;
    std::vector<Ids> nbrs;
    std::vector<bool> used;
    std::vector<double> vouched;
    double fillet_tol = .1, grow_tol = .25;
    FilletPort(const Pass &P) : pass(P), m(P.m), fs(P.fs), used(fs.areas.size(), false) {
        for (auto &list : fs.beside) {
            FeatureNeighborOrder order;
            for (auto id : list)
                order.add(id);
            nbrs.push_back(order.values());
        }
    }
    Ids around(const Ids &group) const {
        std::set<unsigned> s(group.begin(), group.end()), out;
        for (auto f : group)
            for (auto h : nbrs[f])
                if (!s.count(h))
                    out.insert(h);
        return {out.begin(), out.end()};
    }
    Ids spread(const Ids &start, const std::function<bool(unsigned)> &accept) const {
        std::set<unsigned> group(start.begin(), start.end());
        Ids frontier = start;
        while (!frontier.empty()) {
            std::set<unsigned> cand;
            for (auto g : frontier)
                for (auto h : nbrs[g])
                    if (!group.count(h))
                        cand.insert(h);
            frontier.clear();
            for (auto h : cand)
                if (accept(h)) {
                    frontier.push_back(h);
                    group.insert(h);
                }
        }
        return {group.begin(), group.end()};
    }
    std::optional<std::array<double, 3>> radius(const std::vector<V2> &core,
                                                const RollConfig &cfg) const {
        auto radii = cfg.both(core);
        radii.erase(
            std::remove_if(radii.begin(), radii.end(), [](double r) { return !(r > 1e-3); }),
            radii.end());
        if (radii.empty())
            return {};
        std::sort(radii.begin(), radii.end());
        double step = std::max(1e-3, .01 * median(radii));
        unsigned first = 0, best_at = 0, best_n = 0;
        for (unsigned i = 1; i <= radii.size(); ++i)
            if (i == radii.size() || radii[i] - radii[i - 1] > step) {
                if (i - first > best_n) {
                    best_at = first;
                    best_n = i - first;
                }
                first = i;
            }
        double r =
            median(std::vector<double>(radii.begin() + best_at, radii.begin() + best_at + best_n));
        auto gap = [&](double radius) {
            auto o = cfg.center(radius);
            std::vector<double> ds;
            for (auto x : core)
                ds.push_back(std::abs((x - o).norm() - radius));
            return percentile(ds, .75);
        };
        double floor = 2 * fs.tolerance + 2 * fs.noise;
        auto cap = [&](double limit, double rr) {
            return std::min(limit, .04 * rr) + 2 * fs.noise;
        };
        double err = gap(r);
        auto priors = vouched;
        std::stable_sort(priors.begin(), priors.end(),
                         [&](double a, double b) { return std::abs(a - r) < std::abs(b - r); });
        for (auto rv : priors) {
            if (std::abs(rv - r) > .05 * rv)
                continue;
            double e = gap(rv);
            if (e <= cap(.05, rv) &&
                (err > cap(fillet_tol, r) || e <= std::max(2 * err, err + fillet_tol / 2)))
                return std::array<double, 3>{
                    rv, e, std::min(cap(std::max(grow_tol, .05), rv), std::max(3 * e, floor))};
        }
        if (err <= cap(fillet_tol, r))
            return std::array<double, 3>{r, err,
                                         std::min(cap(grow_tol, r), std::max(3 * err, floor))};
        return {};
    }
    Found fillet(const RollSection &sec) {
        double small = std::min(fs.areas[sec.A], sec.flat ? fs.areas[sec.B] : sec.area_B) / 3,
               exact = 2 * fs.tolerance + 2 * fs.noise, tol = exact + fillet_tol,
               flat = std::cos(pi / 180);
        if (sec.t1 <= sec.t0)
            return {};
        auto reaches_B = [&](const Ids &facets) {
            auto near = around(facets);
            for (auto id : near)
                if (sec.b_mask[id])
                    return true;
            if (sec.flat)
                return false;
            for (auto id : near)
                if (!used[id] && !sec.a_mask[id])
                    for (auto v : fs.verts[id])
                        if (std::abs((sec.coords(m.points[v]) - sec.p).dot(sec.n2)) <=
                            tol / 2 + exact)
                            return true;
            return false;
        };
        for (double sign : {1., -1.}) {
            RollConfig cfg(sec.n1, sec.n2, sec.p, sign);
            std::vector<double> guess(fs.areas.size(), -1);
            auto judge = [&](unsigned id) -> double {
                if (used[id] || sec.a_mask[id] || sec.b_mask[id] || fs.areas[id] > small)
                    return NAN;
                double absA = 0, absB = 0, maxA = -1e300, maxB = -1e300, lo = 1e300, hi = -1e300,
                       sum = 0;
                bool bad = false;
                for (auto v : fs.verts[id]) {
                    auto x = sec.coords(m.points[v]);
                    auto heights = cfg.heights(x);
                    absA = std::max(absA, std::abs(heights.first));
                    absB = std::max(absB, std::abs(heights.second));
                    maxA = std::max(maxA, sign * heights.first);
                    maxB = std::max(maxB, sign * heights.second);
                    auto rr = cfg.radii(x);
                    bad |= rr.second;
                    double r = rr.second ? 1 : rr.first;
                    lo = std::min(lo, r);
                    hi = std::max(hi, r);
                    sum += r;
                }
                bool alongA = fs.normals[id].dot(sec.normal_A) > flat,
                     alongB = sec.flat ? fs.normals[id].dot(sec.normal_B) > flat
                                       : std::abs(fs.normals[id].dot(
                                             sec.wall.model.normal(fs.centers[id]))) > flat;
                if ((alongA && absA <= tol) || (alongB && absB <= tol) || absA <= exact ||
                    absB <= exact || maxA > tol || maxB > tol || bad ||
                    hi > std::exp(.2) * std::max(lo, 1e-9))
                    return NAN;
                return sum / fs.verts[id].size();
            };
            auto guesses = [&](unsigned id) {
                if (guess[id] == -1)
                    guess[id] = judge(id);
                return guess[id];
            };
            Ids seeds;
            for (auto id : around(sec.group_A))
                if (std::isfinite(guesses(id)))
                    seeds.push_back(id);
            if (seeds.empty())
                continue;
            std::stable_sort(seeds.begin(), seeds.end(),
                             [&](unsigned a, unsigned b) { return guess[a] < guess[b]; });
            std::vector<Ids> groups;
            unsigned first = 0;
            for (unsigned i = 1; i <= seeds.size(); ++i)
                if (i == seeds.size() ||
                    std::log(guess[seeds[i]]) - std::log(guess[seeds[first]]) > .3) {
                    groups.emplace_back(seeds.begin() + first, seeds.begin() + i);
                    first = i;
                }
            std::stable_sort(groups.begin(), groups.end(),
                             [](const Ids &a, const Ids &b) { return a.size() > b.size(); });
            Ids band;
            double middle = 0;
            for (auto group : groups) {
                std::vector<double> values;
                for (auto id : group)
                    values.push_back(guess[id]);
                double mid = median(values), rough = std::exp(-.3) * mid;
                if (sec.t1 - sec.t0 < rough || std::min(sec.width_A, sec.width_B) < .3 * rough)
                    continue;
                auto grown = spread(group, [&](unsigned id) {
                    return std::abs(std::log(guesses(id) / mid)) <= .3;
                });
                if (reaches_B(grown)) {
                    band = grown;
                    seeds = group;
                    middle = mid;
                    break;
                }
            }
            if (band.empty())
                continue;
            auto V = points(m, fs, band);
            std::vector<V2> core;
            for (auto q : V) {
                double t = sec.along(q);
                if (t >= sec.t0 - exact && t <= sec.t1 + exact)
                    core.push_back(sec.coords(q));
            }
            if (core.size() < 4) {
                double reach = std::exp(.3) * middle + exact;
                Ids strips;
                for (auto id : band) {
                    double lo = 1e300, hi = -1e300;
                    for (auto v : fs.verts[id]) {
                        double t = sec.along(m.points[v]);
                        lo = std::min(lo, t);
                        hi = std::max(hi, t);
                    }
                    if (lo < sec.t1 - exact && hi > sec.t0 + exact && lo >= sec.t0 - reach &&
                        hi <= sec.t1 + reach)
                        strips.push_back(id);
                }
                if (!strips.empty()) {
                    core.clear();
                    for (auto q : points(m, fs, strips))
                        core.push_back(sec.coords(q));
                }
            }
            if (core.size() < 4)
                continue;
            auto fit = radius(core, cfg);
            if (!fit)
                continue;
            double r = (*fit)[0], err = (*fit)[1],
                   allowed = std::max((*fit)[2], std::min(sec.rough, grow_tol + 2 * fs.noise));
            if (sec.t1 - sec.t0 < r || std::min(sec.width_A, sec.width_B) < .3 * r)
                continue;
            auto o = cfg.center(r);
            if (!o.allFinite())
                continue;
            double run_on = r + exact;
            std::vector<bool> in_band(fs.areas.size(), false);
            for (auto id : band)
                in_band[id] = true;
            auto off = [&](unsigned id) {
                double result = 0;
                for (auto v : fs.verts[id])
                    result = std::max(result, std::abs((sec.coords(m.points[v]) - o).norm() - r));
                return result;
            };
            auto on = [&](unsigned id, bool band_only, double slack) {
                double lo = 1e300, hi = -1e300;
                for (auto v : fs.verts[id]) {
                    double t = sec.along(m.points[v]);
                    lo = std::min(lo, t);
                    hi = std::max(hi, t);
                }
                bool mine = band_only ? in_band[id]
                                      : (!used[id] && !sec.a_mask[id] && !sec.b_mask[id] &&
                                         fs.areas[id] <= small);
                return mine && lo < sec.t1 - exact && hi > sec.t0 + exact &&
                       lo >= sec.t0 - run_on && hi <= sec.t1 + run_on && off(id) <= allowed + slack;
            };
            Ids start;
            for (auto id : seeds)
                if (on(id, true, 0))
                    start.push_back(id);
            Ids facets;
            if (!start.empty())
                facets = spread(start, [&](unsigned id) { return on(id, true, 0); });
            if (!facets.empty())
                facets = spread(facets, [&](unsigned id) { return on(id, false, fillet_tol / 4); });
            if (facets.empty() || !reaches_B(facets))
                continue;
            std::set<std::array<long long, 3>> normals;
            for (auto id : facets) {
                auto n = fs.normals[id];
                normals.insert(
                    {std::llrint(n[0] * 100), std::llrint(n[1] * 100), std::llrint(n[2] * 100)});
            }
            if (normals.size() < 2 || (!sec.flat && sec.profile_span(facets) < 5 * pi / 180))
                continue;
            bool between = false;
            std::vector<double> along;
            for (auto q : points(m, fs, facets)) {
                auto h = cfg.heights(sec.coords(q));
                if (-sign * h.first > tol && -sign * h.second > tol)
                    between = true;
                along.push_back(sec.along(q));
            }
            if (!between || range(along) < r)
                continue;
            auto span = [&](const Ids &ids) {
                std::vector<double> angles;
                for (auto q : points(m, fs, ids)) {
                    V2 x = sec.coords(q) - o;
                    angles.push_back(std::atan2(x[1], x[0]));
                }
                return arc_span(angles);
            };
            if (span(facets) > cfg.span + 8 * pi / 180)
                continue;
            auto whole = spread(facets, [&](unsigned id) {
                return !used[id] && !sec.a_mask[id] && !sec.b_mask[id] && off(id) <= allowed;
            });
            if (whole.size() > facets.size() && span(whole) > cfg.span + 8 * pi / 180)
                continue;
            double largest = 0;
            for (auto id : facets)
                largest = std::max(largest, fs.areas[id]);
            if (std::min(sec.area_A, sec.area_B) < 3 * largest)
                continue;
            auto model = sec.model(o, r);
            auto out = pass.describe(model, facets, true);
            if (out) {
                out.convex = sign > 0;
                out.rolled = true;
                last_error = err;
                return out;
            }
        }
        return {};
    }
    double last_error = 0;
    bool could_join(unsigned a, unsigned b) const {
        RollSection sec(pass, a, b);
        double near[2], width[2], lo[2], hi[2];
        unsigned i = 0;
        for (auto id : {a, b}) {
            auto n = fs.normals[id];
            std::vector<double> across, along;
            for (auto q : points(m, fs, {id})) {
                q -= sec.origin;
                across.push_back(q.dot(n.cross(sec.d)));
                along.push_back(q.dot(sec.d));
            }
            auto range_a = std::minmax_element(across.begin(), across.end());
            near[i] = *range_a.first > 0    ? *range_a.first
                      : *range_a.second < 0 ? -*range_a.second
                                            : 0;
            width[i] = range(across);
            auto range_t = std::minmax_element(along.begin(), along.end());
            lo[i] = *range_t.first;
            hi[i] = *range_t.second;
            ++i;
        }
        if (std::min(near[0], near[1]) <= 0)
            return true;
        if (std::abs(std::log(near[0] / near[1])) > .7)
            return false;
        double theta = std::acos(std::clamp(fs.normals[a].dot(fs.normals[b]), -1., 1.));
        double r = std::sqrt(near[0] * near[1]) / std::tan(theta / 2), rough = std::exp(-.7) * r;
        return std::min(width[0], width[1]) >= .3 * rough &&
               std::min(hi[0], hi[1]) - std::max(lo[0], lo[1]) >= rough;
    }
    struct PairOrder {
        using Key = std::pair<int, int>;
        struct Slot {
            Key value;
            std::uint64_t hash = 0;
            bool used = false;
        };
        std::vector<Slot> slots = std::vector<Slot>(8);
        std::size_t used = 0;
        static std::uint64_t hash(Key key) {
            std::uint64_t a = 2870177450012600261ULL;
            for (auto k : {key.first, key.second}) {
                auto h = std::uint64_t(k == -1 ? -2 : k);
                a += h * 14029467366897019727ULL;
                a = (a << 31) | (a >> 33);
                a *= 11400714785074694791ULL;
            }
            a += 2 ^ (2870177450012600261ULL ^ 3527539ULL);
            return a == UINT64_MAX ? 1546275796 : a;
        }
        bool insert(Key key, std::uint64_t h) {
            const std::uint64_t mask = slots.size() - 1;
            std::uint64_t at = h & mask, perturb = h;
            for (;;) {
                unsigned count = at + 9 <= mask ? 10 : 1;
                for (unsigned j = 0; j < count; ++j) {
                    auto &slot = slots[at + j];
                    if (!slot.used) {
                        slot = {key, h, true};
                        return true;
                    }
                    if (slot.hash == h && slot.value == key)
                        return false;
                }
                perturb >>= 5;
                at = (5 * at + 1 + perturb) & mask;
            }
        }
        void add(Key key) {
            if (!insert(key, hash(key)))
                return;
            ++used;
            if (used * 5 < (slots.size() - 1) * 3)
                return;
            std::size_t n = 8, minimum = used * (used > 50000 ? 2 : 4);
            while (n <= minimum)
                n *= 2;
            auto old = std::move(slots);
            slots.assign(n, {});
            for (auto &slot : old)
                if (slot.used)
                    insert(slot.value, slot.hash);
        }
        std::vector<Key> values() const {
            std::vector<Key> out;
            for (auto &slot : slots)
                if (slot.used)
                    out.push_back(slot.value);
            return out;
        }
    };
    std::vector<Found> find() {
        Ids faces;
        for (unsigned id = 0; id < fs.areas.size(); ++id)
            if (fs.areas[id] >= 2 && !used[id])
                faces.push_back(id);
        std::stable_sort(faces.begin(), faces.end(),
                         [&](auto a, auto b) { return fs.areas[a] > fs.areas[b]; });
        std::vector<bool> is_face(fs.areas.size(), false);
        for (auto id : faces)
            is_face[id] = true;
        using Key = PairOrder::Key;
        auto linked = [&](const std::vector<Found> &walls) {
            PairOrder pairs;
            std::vector<int> wall_at(fs.areas.size(), -1);
            for (unsigned k = 0; k < walls.size(); ++k)
                for (auto id : walls[k].facets)
                    wall_at[id] = k;
            double lo = std::cos(170 * pi / 180), hi = std::cos(30 * pi / 180),
                   square = std::cos(2 * pi / 180), beside = std::sin(2 * pi / 180);
            for (auto A : faces) {
                if (used[A])
                    continue;
                double small = fs.areas[A] / 3;
                std::set<unsigned> seen{A};
                Ids frontier{A};
                for (unsigned ring = 0; ring < 48; ++ring) {
                    Ids step;
                    for (auto g : frontier)
                        for (auto h : nbrs[g]) {
                            if (!seen.insert(h).second)
                                continue;
                            if (walls.empty()) {
                                double dot = fs.normals[A].dot(fs.normals[h]);
                                if (is_face[h] && lo <= dot && dot <= hi)
                                    pairs.add({std::min(A, h), std::max(A, h)});
                            } else if (wall_at[h] >= 0) {
                                auto &X = walls[wall_at[h]];
                                double along = std::abs(fs.normals[A].dot(X.model.d));
                                if (along >= square || (along <= beside && X.model.k == 0))
                                    pairs.add({A, -1 - wall_at[h]});
                            }
                            if (!used[h] && fs.areas[h] <= small)
                                step.push_back(h);
                        }
                    frontier = std::move(step);
                    if (frontier.empty())
                        break;
                }
            }
            auto out = pairs.values();
            std::stable_sort(out.begin(), out.end(), [&](Key a, Key b) {
                auto area = [&](Key k) {
                    return walls.empty() ? std::min(fs.areas[k.first], fs.areas[k.second])
                                         : fs.areas[k.first];
                };
                return area(a) > area(b);
            });
            return out;
        };
        double exact = 2 * fs.tolerance + 2 * fs.noise;
        std::map<unsigned, Ids> groups;
        auto coplanar = [&](unsigned a) {
            if (!groups.count(a)) {
                Ids ids{a};
                auto point = m.points[fs.verts[a][0]];
                for (unsigned i = 0; i < fs.areas.size(); ++i)
                    if (i != a && !used[i] &&
                        fs.normals[i].dot(fs.normals[a]) > std::cos(.25 * pi / 180) &&
                        std::abs((fs.centers[i] - point).dot(fs.normals[a])) <= exact)
                        ids.push_back(i);
                groups[a] = std::move(ids);
            }
            Ids out;
            for (auto id : groups[a])
                if (!used[id])
                    out.push_back(id);
            return out;
        };
        std::vector<Found> found;
        auto take = [&](const std::vector<Key> &pairs, const std::vector<Found> &walls) {
            std::vector<Key> waiting;
            for (unsigned round = 0; round < 2; ++round) {
                const auto current = round == 0 ? pairs : waiting;
                for (auto item : current) {
                    unsigned A = item.first;
                    int B = item.second;
                    Found got;
                    if (B < 0) {
                        if (used[A])
                            continue;
                        auto &X = walls[-1 - B];
                        if (std::abs(fs.normals[A].dot(X.model.d)) >= std::cos(2 * pi / 180))
                            got = fillet(RollSection(pass, A, X, coplanar(A)));
                    } else {
                        if (used[A] || used[B])
                            continue;
                        if (round == 0 && !could_join(A, B))
                            continue;
                        got = fillet(RollSection(pass, A, B));
                        if (!got)
                            got = fillet(RollSection(pass, B, A));
                    }
                    if (!got || (round == 0 && last_error > .005 + fs.noise)) {
                        if (round == 0)
                            waiting.push_back(item);
                        continue;
                    }
                    for (auto id : got.facets)
                        used[id] = true;
                    found.push_back(got);
                    if (last_error <= .005 + fs.noise) {
                        double r = std::nearbyint(radius_of(got.model) * 1e6) / 1e6;
                        vouched.push_back(r);
                        std::sort(vouched.begin(), vouched.end());
                        vouched.erase(std::unique(vouched.begin(), vouched.end()), vouched.end());
                    }
                }
            }
        };
        take(linked({}), {});
        std::vector<Found> walls;
        for (auto &f : found)
            if (f.model.kind == "cylinder" || f.model.kind == "cone")
                walls.push_back(f);
        std::vector<unsigned> parent(walls.size());
        std::iota(parent.begin(), parent.end(), 0);
        auto root = [&](unsigned a) {
            while (parent[a] != a) {
                parent[a] = parent[parent[a]];
                a = parent[a];
            }
            return a;
        };
        double wall_tol = exact + fillet_tol / 2;
        for (unsigned a = 0; a < walls.size(); ++a)
            for (unsigned b = a + 1; b < walls.size(); ++b) {
                auto &X = walls[a];
                auto &Y = walls[b];
                if (root(a) == root(b) || X.convex != Y.convex ||
                    std::abs(X.model.d.dot(Y.model.d)) < std::cos(.5 * pi / 180) ||
                    (Y.model.a - X.model.a).cross(X.model.d).norm() > 10 * wall_tol)
                    continue;
                if (fits(X.model, points(m, fs, Y.facets), wall_tol) &&
                    fits(Y.model, points(m, fs, X.facets), wall_tol))
                    parent[root(b)] = root(a);
            }
        std::map<unsigned, std::vector<unsigned>> joined;
        std::vector<unsigned> order;
        for (unsigned k = 0; k < walls.size(); ++k) {
            auto id = root(k);
            if (!joined.count(id))
                order.push_back(id);
            joined[id].push_back(k);
        }
        std::vector<Found> whole;
        for (auto id : order) {
            auto &g = joined[id];
            auto best = g[0];
            for (auto k : g)
                if (walls[k].facets.size() > walls[best].facets.size())
                    best = k;
            auto f = walls[best];
            f.convex = walls[g[0]].convex;
            f.facets.clear();
            for (auto k : g)
                f.facets.insert(f.facets.end(), walls[k].facets.begin(), walls[k].facets.end());
            whole.push_back(f);
        }
        if (!whole.empty())
            take(linked(whole), whole);
        return found;
    }
    static double radius_of(const Model &m) { return m.kind == "torus" ? m.tube : m.r; }
};
} // namespace stl2curves::detail::reference_features

// Derived from stl2curves, Copyright (c) 2026 Duvdivan. MIT license.
#pragma once
#include "feature_facets.hpp"
#include <Eigen/QR>
#include <Eigen/SVD>
#include <optional>
namespace stl2curves::detail::reference_features {
using namespace stl2curves::detail;
constexpr double pi = 3.14159265358979323846;
inline double norm3(const Vec &v) {
    return std::sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
}
inline Vec divided(const Vec &v, double d) {
    return Vec(v.x() / d, v.y() / d, v.z() / d);
}
inline Vec normalized(const Vec &v) {
    return divided(v, norm3(v));
}
inline Eigen::VectorXd lstsq(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
    Eigen::JacobiSVD<Eigen::MatrixXd> solve(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    solve.setThreshold(std::numeric_limits<double>::epsilon() * std::max(A.rows(), A.cols()));
    return solve.solve(b);
}
inline Eigen::Vector2d polyfit(Eigen::MatrixXd A, const Eigen::VectorXd &b) {
    Eigen::Vector2d scale(A.col(0).norm(), A.col(1).norm());
    A.col(0) /= scale[0];
    A.col(1) /= scale[1];
    Eigen::Vector2d x = lstsq(A, b);
    return x.cwiseQuotient(scale);
}
struct Model {
    std::string kind;
    Vec a = Vec::Zero(), d = Vec::UnitZ(), e1 = Vec::UnitX(), e2 = Vec::UnitY();
    double r = 0, k = 0, tube = 0, zc = 0;
    Vec center() const { return (kind == "torus" || sphere_cap) ? Vec(a + zc * d) : a; }
    void at(const Vec &p) { a = p - p.dot(d) * d; }
    bool sphere_cap = false;
    explicit operator bool() const { return !kind.empty(); }
    void frame(Vec axis) {
        d = normalized(axis);
        e1 = std::abs(d.x()) < .9 ? Vec(0, d.z(), -d.y()) : Vec(-d.z(), 0, d.x());
        e1 = normalized(e1);
        e2 = d.cross(e1);
    }
    double signed_at(const Vec &p) const {
        Vec q = p - a;
        if (kind == "sphere" && !sphere_cap)
            return norm3(q) - r;
        double z = q.dot(d), rho = norm3(q - z * d);
        if (kind == "torus" || sphere_cap)
            return std::hypot(rho - (sphere_cap ? 0 : r), z - zc) - (sphere_cap ? r : tube);
        return (rho - r - k * z) / std::hypot(1., k);
    }
    Vec normal(const Vec &p) const {
        Vec q = p - a;
        if (kind == "sphere" && !sphere_cap)
            return divided(q, std::max(norm3(q), 1e-12));
        double z = q.dot(d);
        Vec w = q - z * d;
        double rho = norm3(w);
        Vec radial = rho < 1e-9 ? Vec::Zero() : divided(w, std::max(rho, 1e-12));
        double nr, nz;
        if (kind == "torus" || sphere_cap) {
            nr = rho - (sphere_cap ? 0 : r);
            double length = std::max(std::hypot(nr, z - zc), 1e-12);
            nr /= length;
            nz = (z - zc) / length;
        } else {
            double length = std::hypot(1., k);
            nr = 1 / length;
            nz = -k / length;
        }
        return radial * nr + nz * d;
    }
    double scale() const {
        return kind == "sphere" ? r : kind == "torus" ? tube : std::max(std::abs(r), 1.);
    }
};
inline std::vector<Vec> points(const Mesh &m, const FeatureFacets &fs,
                               const std::vector<unsigned> &fids) {
    std::set<unsigned> vs;
    for (auto f : fids)
        vs.insert(fs.verts[f].begin(), fs.verts[f].end());
    std::vector<Vec> p;
    for (auto v : vs)
        p.push_back(m.points[v]);
    return p;
}
inline std::optional<Vec> circle(const std::vector<Eigen::Vector2d> &p) {
    std::set<std::pair<long long, long long>> unique;
    for (auto &q : p)
        unique.emplace(std::nearbyint(q.x() * 1e5), std::nearbyint(q.y() * 1e5));
    if (unique.size() < 3)
        return {};
    Eigen::MatrixXd A(p.size(), 3);
    Eigen::VectorXd b(p.size());
    for (unsigned i = 0; i < p.size(); ++i) {
        A.row(i) << p[i].x(), p[i].y(), 1;
        b[i] = p[i].squaredNorm();
    }
    Vec sol = lstsq(A, b);
    double cx = sol.x() / 2, cy = sol.y() / 2, r2 = sol.z() + cx * cx + cy * cy;
    if (r2 <= 0)
        return {};
    return Vec(cx, cy, std::sqrt(r2));
}
inline Model cylinder(Vec axis, const std::vector<Vec> &p) {
    Model m;
    m.frame(axis);
    std::vector<Eigen::Vector2d> xy;
    for (auto &q : p)
        xy.emplace_back(q.dot(m.e1), q.dot(m.e2));
    auto c = circle(xy);
    if (!c || c->z() > 1e4)
        return {};
    m.kind = "cylinder";
    m.at(c->x() * m.e1 + c->y() * m.e2);
    m.r = c->z();
    return m;
}
inline Model sphere(const std::vector<Vec> &p) {
    if (p.size() < 4)
        return {};
    Vec mean = Vec::Zero();
    for (auto &q : p)
        mean += q;
    mean /= double(p.size());
    Eigen::MatrixXd A(p.size(), 4), X(p.size(), 3);
    Eigen::VectorXd b(p.size());
    for (unsigned i = 0; i < p.size(); ++i) {
        A.row(i) << 2 * p[i].x(), 2 * p[i].y(), 2 * p[i].z(), 1;
        b[i] = p[i].squaredNorm();
        X.row(i) = p[i] - mean;
    }
    if (X.jacobiSvd().singularValues()[2] < 1e-6)
        return {};
    Eigen::Vector4d sol = lstsq(A, b);
    double r2 = sol[3] + sol.head<3>().squaredNorm();
    if (r2 <= 0 || r2 >= 1e8)
        return {};
    Model m;
    m.kind = "sphere";
    m.a = sol.head<3>();
    m.r = std::sqrt(r2);
    return m;
}
inline Model line_on_direction(Vec d, const std::vector<Vec> &p, double tol) {
    Model m;
    m.frame(d);
    Eigen::MatrixXd A(p.size(), 5);
    Eigen::VectorXd b(p.size());
    std::vector<double> x, y, z;
    for (unsigned i = 0; i < p.size(); ++i) {
        double xx = p[i].dot(m.e1), yy = p[i].dot(m.e2), zz = p[i].dot(m.d);
        x.push_back(xx);
        y.push_back(yy);
        z.push_back(zz);
        A.row(i) << 2 * xx, 2 * yy, 1, zz, zz * zz;
        b[i] = xx * xx + yy * yy;
    }
    auto zr = std::minmax_element(z.begin(), z.end());
    if (*zr.second - *zr.first < 1e-6)
        return {};
    Eigen::VectorXd sol = lstsq(A, b);
    Eigen::MatrixXd Z(p.size(), 2);
    Eigen::VectorXd rho(p.size());
    for (unsigned i = 0; i < p.size(); ++i) {
        Z.row(i) << z[i], 1;
        rho[i] = std::hypot(x[i] - sol[0], y[i] - sol[1]);
    }
    Eigen::Vector2d profile = polyfit(Z, rho);
    m.k = profile[0];
    m.r = profile[1];
    if (std::abs(m.k) > std::tan(80 * pi / 180))
        return {};
    if (std::abs(m.k) < 1e-4 || std::abs(m.k) * (*zr.second - *zr.first) / 2 <= tol) {
        m.k = 0;
        m.r = rho.mean();
    }
    m.kind = m.k == 0 ? "cylinder" : "cone";
    m.at(sol[0] * m.e1 + sol[1] * m.e2);
    return m;
}
inline bool fits(const Model &m, const std::vector<Vec> &p, double tol) {
    if (!m)
        return false;
    for (auto &q : p)
        if (std::abs(m.signed_at(q)) > tol)
            return false;
    return true;
}
inline Model on_axis(Vec a, Vec d, const std::vector<Vec> &p, double tol) {
    Model m;
    m.frame(d);
    m.a = a - a.dot(m.d) * m.d;
    std::vector<Eigen::Vector2d> rz;
    Eigen::MatrixXd A(p.size(), 2);
    Eigen::VectorXd rho(p.size());
    double lo = 1e100, hi = -1e100;
    for (unsigned i = 0; i < p.size(); ++i) {
        Vec q = p[i] - m.a;
        double z = q.dot(m.d), r = (q - z * m.d).norm();
        rz.emplace_back(r, z);
        lo = std::min(lo, z);
        hi = std::max(hi, z);
        A.row(i) << z, 1;
        rho[i] = r;
    }
    if (hi - lo > 1e-6) {
        Eigen::Vector2d profile = polyfit(A, rho);
        m.k = profile[0];
        m.r = profile[1];
        if (std::abs(m.k) < 1e-4 || std::abs(m.k) * (hi - lo) / 2 <= tol) {
            m.k = 0;
            m.r = rho.mean();
        }
        m.kind = m.k == 0 ? "cylinder" : "cone";
        if (std::abs(m.k) < std::tan(80 * pi / 180) && fits(m, p, tol))
            return m;
    }
    auto c = circle(rz);
    if (!c || c->z() > 1e4)
        return {};
    m.r = std::abs(c->x()) < .001 ? 0 : c->x();
    m.zc = c->y();
    m.tube = c->z();
    m.k = 0;
    if (m.r == 0) {
        m.kind = "sphere";
        m.sphere_cap = true;
        m.r = m.tube;
    } else
        m.kind = "torus";
    return m;
}
struct Found {
    Model model;
    std::vector<unsigned> facets;
    bool convex = false;
    bool trimmed = false;
    bool rolled = false;
    bool vouched = false, has_parts = false;
    double worst = 0, area = 0;
    double u0 = 0, span = 2 * pi, lo = 0, hi = 0;
    explicit operator bool() const { return bool(model); }
};
inline Model tube_fit(const std::vector<Vec> &points, const Model &model) {
    if (!model)
        return {};
    bool torus = model.kind == "torus" || model.sphere_cap, cone = model.kind == "cone";
    const Vec d0 = model.d, e1 = model.e1, e2 = model.e2, c0 = model.center();
    unsigned n = torus ? 7 : cone ? 6 : 5;
    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    if (torus) {
        x[5] = model.sphere_cap ? 0 : model.r;
        x[6] = model.sphere_cap ? model.r : model.tube;
    } else {
        x[4] = model.r;
        if (cone)
            x[5] = model.k;
    }
    auto residual = [&](const Eigen::VectorXd &v) -> Eigen::VectorXd {
        Vec d = normalized(d0 + v[0] * e1 + v[1] * e2);
        Vec c = c0 + v[2] * e1 + v[3] * e2;
        if (torus)
            c += v[4] * d0;
        Eigen::VectorXd f(points.size());
        for (unsigned i = 0; i < points.size(); ++i) {
            Vec q = points[i] - c;
            double z = q.dot(d), rho = (q - z * d).norm();
            f[i] = torus  ? std::hypot(rho - v[5], z) - v[6]
                   : cone ? (rho - v[4] - v[5] * z) / std::hypot(1., v[5])
                          : rho - v[4];
        }
        return f;
    };
    auto f = residual(x);
    for (unsigned it = 0; it < 30; ++it) {
        Eigen::MatrixXd J(points.size(), n);
        for (unsigned k = 0; k < n; ++k) {
            double h = 1e-7 * std::max(1., std::abs(x[k]));
            auto nudged = x;
            nudged[k] += h;
            J.col(k) = (residual(nudged) - f) / h;
        }
        Eigen::JacobiSVD<Eigen::MatrixXd> solve(J, Eigen::ComputeThinU | Eigen::ComputeThinV);
        solve.setThreshold(std::numeric_limits<double>::epsilon() * std::max(J.rows(), J.cols()));
        Eigen::VectorXd step = solve.solve(-f), next = x + step;
        auto f2 = residual(next);
        while (f2.squaredNorm() > f.squaredNorm() && step.cwiseAbs().maxCoeff() > 1e-14) {
            step /= 2;
            next = x + step;
            f2 = residual(next);
        }
        if (f2.squaredNorm() > f.squaredNorm())
            break;
        x = next;
        f = f2;
        if (step.cwiseAbs().maxCoeff() < 1e-12)
            break;
    }
    if (!x.allFinite())
        return {};
    Model result;
    Vec final_axis = normalized(d0 + x[0] * e1 + x[1] * e2);
    result.frame(final_axis);
    Vec c = c0 + x[2] * e1 + x[3] * e2;
    if (torus) {
        c += x[4] * d0;
        result.at(c);
        result.zc = c.dot(final_axis);
        result.r = std::abs(x[5]);
        result.tube = std::abs(x[6]);
        result.kind = "torus";
        if (result.r == 0) {
            result.kind = "sphere";
            result.sphere_cap = true;
            result.r = result.tube;
        }
    } else {
        result.a = c - c.dot(result.d) * result.d;
        result.r = cone ? x[4] - x[5] * c.dot(result.d) : std::abs(x[4]);
        result.k = cone ? x[5] : 0;
        result.kind = cone ? "cone" : "cylinder";
    }
    return result;
}
struct Pass {
    const Mesh &m;
    const Topology &top;
    FeatureFacets fs;
    std::vector<bool> free;
    mutable std::vector<bool> explored;
    bool loose = true;
    std::vector<std::pair<Model, double>> axes;
    Pass(const Mesh &mesh, const Topology &t, const std::vector<Patch> &known)
        : m(mesh), top(t), fs(mesh, t), free(fs.triangles.size(), true),
          explored(free.size(), false) {
        std::vector<int> o(m.triangles.size(), -1);
        for (unsigned k = 0; k < known.size(); ++k) {
            const auto &patch = known[k];
            if (patch.kind == "plane")
                continue;
            Model model;
            model.kind = patch.kind;
            model.frame(patch.axis);
            model.a = patch.origin;
            model.r = patch.radius;
            model.k = patch.kind == "cone" ? patch.secondary : 0;
            model.tube = patch.kind == "torus" ? patch.secondary : 0;
            if (patch.kind == "torus") {
                model.zc = model.a.dot(model.d);
                model.at(model.a);
            }
            double size = patch.kind == "torus" ? model.r + model.tube : model.r;
            bool duplicate = false;
            if (patch.kind != "sphere")
                for (auto &entry : axes)
                    if (std::abs(model.d.dot(entry.first.d)) > std::cos(.5 * pi / 180) &&
                        (model.a - entry.first.a).cross(entry.first.d).norm() < .004) {
                        entry.second = std::max(entry.second, 2 * size + 3);
                        duplicate = true;
                        break;
                    }
            if (!duplicate && patch.kind != "sphere")
                axes.emplace_back(model, 2 * size + 3);
            for (auto t : known[k].triangles)
                o[t] = k;
        }
        for (unsigned f = 0; f < free.size(); ++f)
            for (auto t : fs.triangles[f])
                if (o[t] >= 0)
                    free[f] = false;
    }
    double size(unsigned f) const {
        Vec lo = m.points[fs.verts[f].front()], hi = lo;
        for (auto v : fs.verts[f]) {
            lo = lo.cwiseMin(m.points[v]);
            hi = hi.cwiseMax(m.points[v]);
        }
        return (hi - lo).maxCoeff();
    }
    bool sliver(unsigned f) const {
        double s = size(f);
        return fs.areas[f] < .05 * s * s;
    }
    bool test(const Model &model, unsigned f) const {
        if (!free[f])
            return false;
        double worst = 0;
        std::vector<Vec> ns;
        for (auto v : fs.verts[f]) {
            worst = std::max(worst, std::abs(model.signed_at(m.points[v])));
            ns.push_back(model.normal(m.points[v]));
        }
        double facing = std::abs(fs.normals[f].dot(model.normal(fs.centers[f])));
        bool aligned;
        if (loose) {
            double least = 1;
            for (auto &a : ns)
                for (auto &b : ns)
                    least = std::min(least, a.dot(b));
            double spread = std::acos(std::clamp(least, -1., 1.));
            aligned = facing >= std::cos(std::min(8 * pi / 180 + spread, pi / 3));
        } else
            aligned = facing >= std::cos(pi / 6) || sliver(f);
        return worst <= fs.tolerance &&
               std::abs(model.signed_at(fs.centers[f])) <= (loose ? .08 : .25) && aligned;
    }
    bool bend(const Model &model, unsigned a, unsigned b) const {
        if ((size(a) < .1 * model.scale() && size(b) < .1 * model.scale()) || sliver(a) ||
            sliver(b))
            return true;
        double actual = std::acos(std::clamp(fs.normals[a].dot(fs.normals[b]), -1., 1.)),
               predicted = std::acos(std::clamp(
                   model.normal(fs.centers[a]).dot(model.normal(fs.centers[b])), -1., 1.));
        return std::abs(actual - predicted) <= 8 * pi / 180;
    }
    std::vector<unsigned> grow(const Model &model, const std::vector<unsigned> &seed) const {
        if (!model)
            return {};
        bool convex = fs.normals[seed.front()].dot(model.normal(fs.centers[seed.front()])) > 0;
        std::set<unsigned> got;
        for (auto f : seed) {
            if (!test(model, f) || (fs.normals[f].dot(model.normal(fs.centers[f])) > 0) != convex ||
                !bend(model, seed.front(), f))
                return {};
            got.insert(f);
        }
        std::vector<unsigned> front = seed;
        while (!front.empty()) {
            std::set<unsigned> candidates;
            for (auto f : front)
                for (auto g : fs.neighbors[f])
                    if (!got.count(g))
                        candidates.insert(g);
            front.clear();
            for (auto f : candidates)
                if (test(model, f) &&
                    (fs.normals[f].dot(model.normal(fs.centers[f])) > 0) == convex) {
                    bool right = false;
                    for (auto g : fs.neighbors[f])
                        if (got.count(g) && bend(model, g, f))
                            right = true;
                    if (right) {
                        got.insert(f);
                        front.push_back(f);
                    }
                }
        }
        return {got.begin(), got.end()};
    }
    Model refit(const Model &model, const std::vector<unsigned> &facets, bool fixed) const {
        auto p = points(m, fs, facets);
        if (model.kind == "sphere" && !model.sphere_cap)
            return sphere(p);
        if (fixed)
            return on_axis(model.a, model.d, p, fs.tolerance);
        if (model.kind != "cylinder" && model.kind != "cone")
            return {};
        std::map<long long, std::vector<Vec>> rings;
        for (auto &q : p)
            rings[std::nearbyint((q - model.a).dot(model.d) / .01)].push_back(q);
        std::vector<Vec> centers;
        for (auto &ring : rings) {
            if (ring.second.size() < 3)
                continue;
            Vec center = Vec::Zero();
            for (auto &q : ring.second)
                center += q;
            center /= double(ring.second.size());
            Eigen::MatrixXd centered(ring.second.size(), 3);
            for (unsigned i = 0; i < ring.second.size(); ++i)
                centered.row(i) = ring.second[i] - center;
            Model frame;
            frame.frame(centered.jacobiSvd(Eigen::ComputeThinV).matrixV().col(2));
            std::vector<Eigen::Vector2d> xy;
            for (auto &q : ring.second)
                xy.emplace_back((q - center).dot(frame.e1), (q - center).dot(frame.e2));
            auto c = circle(xy);
            if (c)
                centers.push_back(center + c->x() * frame.e1 + c->y() * frame.e2);
        }
        if (centers.size() < 2 || (centers.back() - centers.front()).norm() < 1e-6)
            return {};
        Vec axis = normalized(centers.back() - centers.front());
        return model.kind == "cylinder" ? cylinder(axis, p)
                                        : line_on_direction(axis, p, fs.tolerance);
    }
    Found describe(const Model &model, const std::vector<unsigned> &ids, bool anchored = false,
                   bool strict = false) const {
        Found out;
        if (ids.size() < 3 || !model)
            return out;
        if (model.kind != "sphere" || model.sphere_cap) {
            auto arc = [](std::vector<double> values) -> std::optional<Eigen::Vector2d> {
                for (auto &v : values) {
                    v = std::fmod(v, 2 * pi);
                    if (v < 0)
                        v += 2 * pi;
                }
                std::sort(values.begin(), values.end());
                std::vector<double> unique;
                for (unsigned i = 0; i < values.size(); ++i)
                    if (i == 0 || values[i] - values[i - 1] > 1e-6)
                        unique.push_back(values[i]);
                if (unique.size() < 2)
                    return {};
                double gap = -1;
                unsigned where = 0;
                for (unsigned i = 0; i < unique.size(); ++i) {
                    double next = i + 1 < unique.size() ? unique[i + 1] : unique[0] + 2 * pi;
                    double g = next - unique[i];
                    if (g > gap) {
                        gap = g;
                        where = i;
                    }
                }
                if (unique.size() >= 6 && gap <= pi / 4)
                    return Eigen::Vector2d(0, 2 * pi);
                return Eigen::Vector2d(unique[(where + 1) % unique.size()], 2 * pi - gap);
            };
            std::vector<double> around, profile;
            for (auto &p : points(m, fs, ids)) {
                Vec q = p - model.a;
                double z = q.dot(model.d);
                Vec w = q - z * model.d;
                double rho = w.norm();
                if (rho > fs.tolerance)
                    around.push_back(std::atan2(w.dot(model.e2), w.dot(model.e1)));
                if (model.kind == "torus" || model.sphere_cap)
                    profile.push_back(
                        std::atan2(z - model.zc, rho - (model.sphere_cap ? 0 : model.r)));
            }
            auto angle = arc(around);
            if (!angle)
                return {};
            if (model.kind == "torus" || model.sphere_cap) {
                auto v = arc(profile);
                if (!v)
                    return {};
                double rc = model.sphere_cap ? 0 : model.r,
                       r = model.sphere_cap ? model.r : model.tube;
                if (v->y() >= 2 * pi - 1e-9 && (angle->y() >= 2 * pi || rc - r <= fs.tolerance))
                    return {};
                for (unsigned i = 0; i <= 8; ++i)
                    if (rc + r * std::cos(v->x() + v->y() * i / 8) < -fs.tolerance)
                        return {};
            }
        }
        std::map<Edge, unsigned> edge_counts;
        for (auto f : ids)
            for (auto t : fs.triangles[f]) {
                auto tri = m.triangles[t];
                for (unsigned i = 0; i < 3; ++i)
                    ++edge_counts[Edge(tri[i], tri[(i + 1) % 3])];
            }
        std::vector<Edge> rim;
        for (auto &e : edge_counts)
            if (e.second == 1)
                rim.push_back(e.first);
        auto P = points(m, fs, ids);
        auto free_edges = rim;
        double angular_span = 2 * pi, profile_span = 2 * pi, profile_lo = 0, profile_hi = 0;
        if (model.kind == "sphere" && !model.sphere_cap) {
            std::set<unsigned> rv;
            for (auto &e : rim) {
                rv.insert(e.a);
                rv.insert(e.b);
            }
            bool flat = false;
            Vec mid = Vec::Zero(), flat_normal = Vec::UnitZ();
            if (rv.size() >= 3) {
                for (auto v : rv)
                    mid += m.points[v];
                mid /= double(rv.size());
                Eigen::MatrixXd centered(rv.size(), 3);
                unsigned row = 0;
                for (auto v : rv)
                    centered.row(row++) = m.points[v] - mid;
                flat_normal = centered.jacobiSvd(Eigen::ComputeThinV).matrixV().col(2);
                flat = true;
                for (auto v : rv)
                    if (std::abs((m.points[v] - mid).dot(flat_normal)) > fs.tolerance)
                        flat = false;
            }
            if (flat) {
                Vec delta = Vec::Zero();
                for (auto &p : P)
                    delta += p - mid;
                if (delta.dot(flat_normal) < 0)
                    flat_normal = -flat_normal;
                Model cap = model;
                cap.frame(flat_normal);
                cap.zc = model.a.dot(flat_normal);
                cap.at(model.a);
                cap.sphere_cap = true;
                cap.tube = cap.r;
                return describe(cap, ids, anchored, strict);
            }
            if (!flat) {
                std::vector<Vec> ns;
                for (auto &e : rim) {
                    Vec n = (m.points[e.a] - model.a).cross(m.points[e.b] - model.a);
                    if (n.norm() < 1e-9)
                        continue;
                    n = normalized(n);
                    bool same = false;
                    for (auto &q : ns)
                        if (std::abs(n.dot(q)) > std::cos(pi / 180))
                            same = true;
                    if (!same)
                        ns.push_back(n);
                    if (ns.size() > 8)
                        break;
                }
                bool simple = !ns.empty() && ns.size() <= 4;
                for (auto &e : rim) {
                    bool on = false;
                    for (auto &n : ns)
                        if (std::abs((m.points[e.a] - model.a).dot(n)) <= fs.tolerance &&
                            std::abs((m.points[e.b] - model.a).dot(n)) <= fs.tolerance)
                            on = true;
                    if (!on)
                        simple = false;
                }
                for (auto &n : ns) {
                    double lo = 1e100, hi = -1e100;
                    for (auto &p : P) {
                        double v = (p - model.a).dot(n);
                        lo = std::min(lo, v);
                        hi = std::max(hi, v);
                    }
                    if (lo < -fs.tolerance && hi > fs.tolerance)
                        simple = false;
                }
                out.trimmed = !simple;
            }
        } else {
            auto arc = [](std::vector<double> v) {
                for (auto &x : v) {
                    x = std::fmod(x, 2 * pi);
                    if (x < 0)
                        x += 2 * pi;
                }
                std::sort(v.begin(), v.end());
                std::vector<double> kept;
                for (unsigned i = 0; i < v.size(); ++i)
                    if (i == 0 || v[i] - v[i - 1] > 1e-6)
                        kept.push_back(v[i]);
                v = std::move(kept);
                double gap = -1;
                unsigned at = 0;
                for (unsigned i = 0; i < v.size(); ++i) {
                    double g = (i + 1 < v.size() ? v[i + 1] : v[0] + 2 * pi) - v[i];
                    if (g > gap) {
                        gap = g;
                        at = i;
                    }
                }
                return v.size() >= 6 && gap <= pi / 4
                           ? Eigen::Vector2d(0, 2 * pi)
                           : Eigen::Vector2d(v[(at + 1) % v.size()], 2 * pi - gap);
            };
            auto local = [&](const Vec &p) {
                Vec q = p - model.a;
                double z = q.dot(model.d);
                Vec w = q - z * model.d;
                return Vec(w.norm(), z, std::atan2(w.dot(model.e2), w.dot(model.e1)));
            };
            std::vector<double> angles, profiles;
            double lo = 1e100, hi = -1e100;
            for (auto &p : P) {
                auto q = local(p);
                if (q.x() > fs.tolerance)
                    angles.push_back(q.z());
                lo = std::min(lo, q.y());
                hi = std::max(hi, q.y());
                if (model.kind == "torus" || model.sphere_cap)
                    profiles.push_back(
                        std::atan2(q.y() - model.zc, q.x() - (model.sphere_cap ? 0 : model.r)));
            }
            auto u = arc(angles);
            angular_span = u.y();
            double tube = model.sphere_cap ? model.r : model.tube;
            if (model.kind == "torus" || model.sphere_cap) {
                auto v = arc(profiles);
                lo = v.x();
                hi = v.x() + v.y();
            }
            profile_span = hi - lo;
            profile_lo = lo;
            profile_hi = hi;
            out.u0 = u.x();
            out.span = u.y();
            out.lo = lo;
            out.hi = hi;
            auto angle_gap = [](double a, double b) {
                return std::abs(std::remainder(a - b, 2 * pi));
            };
            auto on = [&](const Vec &p) {
                auto q = local(p);
                std::array<bool, 4> r{false, false, false, false};
                if (model.kind == "torus" || model.sphere_cap) {
                    double v =
                        std::atan2(q.y() - model.zc, q.x() - (model.sphere_cap ? 0 : model.r));
                    r[0] = angle_gap(v, lo) * tube <= fs.tolerance;
                    r[1] = angle_gap(v, hi) * tube <= fs.tolerance;
                } else {
                    r[0] = std::abs(q.y() - lo) <= fs.tolerance;
                    r[1] = std::abs(q.y() - hi) <= fs.tolerance;
                }
                if (u.y() < 2 * pi) {
                    r[2] = angle_gap(q.z(), u.x()) * q.x() <= fs.tolerance || q.x() <= fs.tolerance;
                    r[3] = angle_gap(q.z(), u.x() + u.y()) * q.x() <= fs.tolerance ||
                           q.x() <= fs.tolerance;
                }
                return r;
            };
            free_edges.clear();
            for (auto &e : rim) {
                auto a = on(m.points[e.a]), b = on(m.points[e.b]);
                bool along = false;
                for (unsigned k = 0; k < 4; ++k)
                    if (a[k] && b[k])
                        along = true;
                if (!along) {
                    out.trimmed = true;
                    free_edges.push_back(e);
                }
            }
        }
        if (strict) {
            if (model.kind == "sphere" && !model.sphere_cap) {
                if (!rim.empty()) {
                    if (P.size() < 8)
                        return {};
                    Vec mean = Vec::Zero();
                    for (auto &p : P)
                        mean += (p - model.a) / model.r;
                    mean /= double(P.size());
                    double least = 1;
                    for (auto &p : P)
                        least = std::min(least, ((p - model.a) / model.r).dot(mean) / mean.norm());
                    if (2 * std::acos(std::clamp(least, -1., 1.)) < 15 * pi / 180)
                        return {};
                }
            } else {
                unsigned minimum = model.kind == "cylinder" ? 8
                                   : model.kind == "cone"   ? 10
                                   : model.sphere_cap       ? 8
                                                            : 12;
                if (P.size() < minimum || angular_span < 15 * pi / 180 ||
                    ((model.kind == "torus" || model.sphere_cap) && profile_span < 15 * pi / 180))
                    return {};
                if (model.kind == "cone") {
                    bool inside = false;
                    for (auto &p : P) {
                        double z = (p - model.a).dot(model.d);
                        if (std::abs(z - profile_lo) > fs.tolerance &&
                            std::abs(z - profile_hi) > fs.tolerance)
                            inside = true;
                    }
                    if (!inside) {
                        std::set<unsigned> members(ids.begin(), ids.end());
                        std::array<bool, 2> soft{false, false};
                        for (auto a : ids)
                            for (auto b : fs.neighbors[a])
                                if (!members.count(b) &&
                                    fs.normals[a].dot(fs.normals[b]) > std::cos(20 * pi / 180)) {
                                    auto entry = fs.shared.find(std::minmax(a, b));
                                    if (entry == fs.shared.end())
                                        continue;
                                    for (auto &e : entry->second)
                                        for (unsigned end = 0; end < 2; ++end) {
                                            double level = end ? profile_hi : profile_lo;
                                            if (std::abs((m.points[e.a] - model.a).dot(model.d) -
                                                         level) <= fs.tolerance &&
                                                std::abs((m.points[e.b] - model.a).dot(model.d) -
                                                         level) <= fs.tolerance)
                                                soft[end] = true;
                                        }
                                }
                        bool all = true;
                        for (unsigned end = 0; end < 2; ++end)
                            if (model.r + model.k * (end ? profile_hi : profile_lo) >
                                    fs.tolerance &&
                                !soft[end])
                                all = false;
                        if (all)
                            return {};
                    }
                }
            }
            if (out.trimmed) {
                if (P.size() < 12 || ids.size() < 6)
                    return {};
                double least = 1;
                for (auto a : ids)
                    for (auto b : ids)
                        least = std::min(least, fs.normals[a].dot(fs.normals[b]));
                if (least > std::cos(15 * pi / 180))
                    return {};
                std::set<unsigned> members(ids.begin(), ids.end());
                double total = 0, smooth = 0;
                for (auto &e : free_edges) {
                    double length = (m.points[e.a] - m.points[e.b]).norm();
                    total += length;
                    int mine = -1, other = -1;
                    for (auto t : top.edges.at(e)) {
                        unsigned f = fs.facet_of[t];
                        if (members.count(f))
                            mine = int(f);
                        else
                            other = int(f);
                    }
                    if (mine >= 0 && other >= 0 &&
                        fs.normals[mine].dot(fs.normals[other]) > std::cos(20 * pi / 180))
                        smooth += length;
                }
                if (smooth > .25 * total)
                    return {};
            }
        }
        unsigned distinct_needed = model.kind == "cylinder" ? 3 : model.kind == "cone" ? 6 : 8;
        std::set<std::array<long long, 3>> distinct;
        for (auto f : ids) {
            auto n = fs.normals[f];
            distinct.insert({std::llround(n.x() * 1000), std::llround(n.y() * 1000),
                             std::llround(n.z() * 1000)});
        }
        if (loose && !anchored && !strict && distinct.size() < distinct_needed)
            return out;
        double least = 1;
        for (auto f : ids)
            for (auto g : ids)
                least = std::min(least, fs.normals[f].dot(fs.normals[g]));
        double turn = model.kind == "cylinder" || model.kind == "cone" ? 6 : 20;
        if (loose && !anchored && !strict && least > std::cos(turn * pi / 180))
            return out;
        out.model = model;
        out.facets = ids;
        out.convex = fs.normals[ids.front()].dot(model.normal(fs.centers[ids.front()])) > 0;
        for (auto f : ids) {
            out.area += fs.areas[f];
            for (auto t : fs.triangles[f]) {
                auto tri = m.triangles[t];
                for (unsigned j = 0; j < 3; ++j)
                    out.worst = std::max(
                        out.worst, std::abs(model.signed_at(
                                       (m.points[tri[j]] + m.points[tri[(j + 1) % 3]]) * .5)));
            }
        }
        if (out.worst < 1e-5 || (loose && !anchored && !strict && out.worst > .08))
            return {};
        return out;
    }
    Found candidate(Model model, const std::vector<unsigned> &seed, bool fixed = false) const {
        auto ids = grow(model, seed);
        if (ids.empty())
            return {};
        for (unsigned n = 0; n < 8; ++n) {
            auto better = refit(model, ids, fixed);
            if (!better || (!fixed && better.kind != model.kind))
                break;
            auto idx = grow(better, seed);
            if (idx.size() < ids.size())
                break;
            bool grew = idx.size() > ids.size();
            model = better;
            ids = std::move(idx);
            if (!grew)
                break;
        }
        std::vector<std::pair<Model, std::vector<unsigned>>> options;
        if (model.kind == "cylinder" || model.kind == "cone" || model.kind == "torus") {
            auto cut = ids;
            for (unsigned n = 0; n < 3; ++n) {
                auto p = points(m, fs, cut);
                double scale = model.kind == "torus" ? model.tube : 1.;
                auto param = [&](const Vec &q) {
                    Vec d = q - model.a;
                    double z = d.dot(model.d);
                    return model.kind == "torus"
                               ? std::atan2(z - model.zc, (d - z * model.d).norm() - model.r)
                               : z;
                };
                std::map<long long, unsigned> levels;
                for (auto &q : p)
                    ++levels[std::nearbyint(param(q) * scale / fs.tolerance)];
                unsigned most = 0;
                for (auto &entry : levels)
                    most = std::max(most, entry.second);
                double lo = 1e100, hi = -1e100;
                for (auto &entry : levels)
                    if (entry.second >= std::max(4., .2 * most)) {
                        double t = entry.first * fs.tolerance / scale;
                        lo = std::min(lo, t);
                        hi = std::max(hi, t);
                    }
                if (lo > hi)
                    break;
                lo -= 2 * fs.tolerance / scale;
                hi += 2 * fs.tolerance / scale;
                std::vector<unsigned> keep;
                for (auto f : cut) {
                    bool inside = true;
                    for (auto v : fs.verts[f]) {
                        double t = param(m.points[v]);
                        if (t < lo || t > hi)
                            inside = false;
                    }
                    if (inside)
                        keep.push_back(f);
                }
                if (keep.size() == cut.size())
                    break;
                cut = std::move(keep);
                if (cut.size() < 3)
                    break;
            }
            if (cut.size() >= 3 && cut.size() < ids.size()) {
                auto better = refit(model, cut, fixed);
                bool ok = bool(better) && (fixed || better.kind == model.kind);
                if (ok)
                    for (auto f : cut)
                        if (!test(better, f))
                            ok = false;
                options.emplace_back(ok ? better : model, cut);
            }
        }
        options.emplace_back(model, ids);
        Found result;
        for (auto &option : options) {
            auto f = describe(option.first, option.second);
            if (!f)
                continue;
            if (!f.trimmed)
                return f;
            if (!result)
                result = f;
        }
        if (!result && ids.size() >= 6)
            for (auto f : ids)
                explored[f] = true;
        return result;
    }
    static Found best(const std::vector<Found> &found) {
        Found b;
        for (auto &f : found)
            if (f && (!b || f.area > b.area))
                b = f;
        if (!b)
            return b;
        Found tight = b;
        for (auto &f : found)
            if (f && f.area >= .99 * b.area && f.worst < tight.worst)
                tight = f;
        return tight.worst < .5 * b.worst ? tight : b;
    }
    Found seed(unsigned i, unsigned j) const {
        auto p = points(m, fs, {i, j});
        std::vector<unsigned> near{i, j};
        for (auto f : fs.neighbors[i])
            if (free[f])
                near.push_back(f);
        for (auto f : fs.neighbors[j])
            if (free[f])
                near.push_back(f);
        auto nearp = points(m, fs, near);
        std::vector<Found> found;
        std::vector<Vec> axes;
        auto shared = fs.shared.find(std::minmax(i, j));
        if (shared != fs.shared.end()) {
            Vec edge = Vec::Zero();
            for (auto &e : shared->second) {
                Vec q = m.points[e.b] - m.points[e.a];
                if (q.norm() > edge.norm())
                    edge = q;
            }
            if (edge.norm() > 0)
                axes.push_back(normalized(edge));
        }
        Vec cross = fs.normals[i].cross(fs.normals[j]);
        if (cross.norm() > std::sin(pi / 180)) {
            cross = normalized(cross);
            if (axes.empty() || std::abs(cross.dot(axes.front())) < std::cos(.5 * pi / 180))
                axes.push_back(cross);
        }
        for (auto &axis : axes) {
            auto cyl = cylinder(axis, p);
            if (fits(cyl, p, fs.tolerance))
                found.push_back(candidate(cyl, {i, j}));
        }
        for (auto k : fs.neighbors[j])
            if (k != i && free[k]) {
                std::vector<Vec> ns{fs.normals[i], fs.normals[j], fs.normals[k]};
                std::set<std::array<double, 3>> distinct;
                for (auto &n : ns)
                    distinct.insert({std::nearbyint(n.x() * 1e4), std::nearbyint(n.y() * 1e4),
                                     std::nearbyint(n.z() * 1e4)});
                if (distinct.size() < 3)
                    break;
                Vec mean = (ns[0] + ns[1] + ns[2]) / 3;
                Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
                for (auto &n : ns) {
                    Vec d = n - mean;
                    cov += d * d.transpose();
                }
                auto pts = points(m, fs, {i, j, k});
                auto cone = line_on_direction(
                    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(cov).eigenvectors().col(0), pts,
                    fs.tolerance);
                if (fits(cone, pts, fs.tolerance))
                    found.push_back(candidate(cone, {i, j, k}));
                break;
            }
        bool any = false;
        for (auto &f : found)
            if (f)
                any = true;
        if (!any && near.size() >= 3) {
            Eigen::MatrixXd N(near.size(), 3);
            for (unsigned k = 0; k < near.size(); ++k)
                N.row(k) = fs.normals[near[k]];
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(N, Eigen::ComputeThinV);
            if (svd.singularValues()[2] <= .1 * svd.singularValues()[0]) {
                auto model = tube_fit(nearp, cylinder(svd.matrixV().col(2), nearp));
                if (fits(model, nearp, fs.tolerance) && model.r < 1e3)
                    found.push_back(candidate(model, {i, j}));
            }
        }
        auto before = found;
        for (auto &f : before)
            if (f) {
                auto ring = on_axis(f.model.a, f.model.d, nearp, fs.tolerance);
                if ((ring.kind == "torus" || ring.kind == "sphere") &&
                    fits(ring, nearp, fs.tolerance))
                    found.push_back(candidate(ring, {i, j}, true));
            }
        for (auto pts : {nearp, p}) {
            auto ball = sphere(pts);
            if (fits(ball, pts, fs.tolerance)) {
                found.push_back(candidate(ball, {i, j}));
                break;
            }
        }
        auto pair_best = best(found);
        found.clear();
        found.push_back(pair_best);
        Vec middle = Vec::Zero();
        for (auto &q : nearp)
            middle += q;
        middle /= double(nearp.size());
        std::vector<std::pair<double, unsigned>> close;
        for (unsigned a = 0; a < this->axes.size(); ++a) {
            double dist = (middle - this->axes[a].first.a).cross(this->axes[a].first.d).norm();
            if (dist <= this->axes[a].second)
                close.emplace_back(dist, a);
        }
        std::stable_sort(close.begin(), close.end());
        for (unsigned a = 0; a < std::min<std::size_t>(6, close.size()); ++a) {
            const auto &axis = this->axes[close[a].second].first;
            auto model = on_axis(axis.a, axis.d, nearp, fs.tolerance);
            if (fits(model, nearp, fs.tolerance))
                found.push_back(candidate(model, {i, j}, true));
        }
        return best(found);
    }
    std::vector<Found> run() {
        std::vector<Found> found;
        for (unsigned i = 0; i < free.size(); ++i) {
            if (!free[i] || explored[i])
                continue;
            for (unsigned n = 0; n < std::min<std::size_t>(3, fs.neighbors[i].size()); ++n) {
                if (explored[i])
                    break;
                auto j = fs.neighbors[i][n];
                if (!free[j])
                    continue;
                auto f = seed(i, j);
                if (f) {
                    for (auto a : f.facets)
                        free[a] = false;
                    found.push_back(std::move(f));
                    break;
                }
            }
        }
        return found;
    }
};
} // namespace stl2curves::detail::reference_features

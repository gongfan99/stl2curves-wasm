#include "internal.hpp"
#include "projection.hpp"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <Eigen/QR>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomFill_Pipe.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Standard_Failure.hxx>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>

namespace stl2curves::detail {
namespace {
constexpr double pi = 3.14159265358979323846;
gp_Pnt point(const Vec &v) {
    return {v.x(), v.y(), v.z()};
}
gp_Dir direction(const Vec &v) {
    return {v.x(), v.y(), v.z()};
}
struct Basis {
    unsigned count;
    std::vector<double> knots;
    explicit Basis(unsigned n) : count(n) {
        knots.insert(knots.end(), 4, 0);
        for (unsigned i = 1; i < n - 3; ++i)
            knots.push_back(double(i) / (n - 3));
        knots.insert(knots.end(), 4, 1);
    }
    std::vector<double> at(double t) const {
        t = std::clamp(t, 0.0, 1.0);
        std::vector<double> a(count + 3, 0), b(count + 3, 0);
        if (t == 1) {
            a[count - 1] = 1;
            return std::vector<double>(a.begin(), a.begin() + count);
        }
        for (unsigned i = 0; i < count + 3; ++i)
            if (t >= knots[i] && t < knots[i + 1])
                a[i] = 1;
        for (unsigned d = 1; d <= 3; ++d) {
            std::fill(b.begin(), b.end(), 0);
            for (unsigned i = 0; i < count + 3 - d; ++i) {
                double l = knots[i + d] - knots[i], r = knots[i + d + 1] - knots[i + 1];
                if (l > 0)
                    b[i] += (t - knots[i]) * a[i] / l;
                if (r > 0)
                    b[i] += (knots[i + d + 1] - t) * a[i + 1] / r;
            }
            a.swap(b);
        }
        a.resize(count);
        return a;
    }
    double greville(unsigned i) const { return (knots[i + 1] + knots[i + 2] + knots[i + 3]) / 3; }
};
Handle(Geom_BSplineSurface) periodic_grid(const PointGrid &grid) {
    const unsigned nu = unsigned(grid.ColLength()), nv = unsigned(grid.RowLength());
    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(nu, nu);
    for (unsigned i = 0; i < nu; ++i) {
        U(i, i) = 4. / 6;
        U(i, (i + 1) % nu) = 1. / 6;
        U(i, (i + nu - 1) % nu) = 1. / 6;
    }
    Eigen::LDLT<Eigen::MatrixXd> solve_u(U);
    Eigen::MatrixXd rhs(nv, nu * 3);
    for (unsigned j = 0; j < nv; ++j) {
        Eigen::MatrixXd points(nu, 3);
        for (unsigned i = 0; i < nu; ++i) {
            const auto &p = grid.Value(i + 1, j + 1);
            points.row(i) << p.X(), p.Y(), p.Z();
        }
        Eigen::MatrixXd c = solve_u.solve(points);
        for (unsigned i = 0; i < nu; ++i)
            for (unsigned d = 0; d < 3; ++d)
                rhs(j, i * 3 + d) = c(i, d);
    }
    // Averaged interpolation knots keep the endpoint systems well conditioned.
    Basis basis(nv);
    basis.knots.assign(4, 0);
    for (unsigned j = 1; j <= nv - 4; ++j)
        basis.knots.push_back(double(j + 1) / (nv - 1));
    basis.knots.insert(basis.knots.end(), 4, 1);
    std::vector<Eigen::Triplet<double>> entries;
    for (unsigned j = 0; j < nv; ++j) {
        auto b = basis.at(double(j) / (nv - 1));
        for (unsigned i = 0; i < nv; ++i)
            if (b[i] != 0)
                entries.emplace_back(j, i, b[i]);
    }
    Eigen::SparseMatrix<double> V(nv, nv);
    V.setFromTriplets(entries.begin(), entries.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solve_v;
    solve_v.compute(V);
    if (solve_v.info() != Eigen::Success)
        throw std::runtime_error("Periodic surface interpolation failed");
    Eigen::MatrixXd coefficients = solve_v.solve(rhs);
    if (!coefficients.allFinite())
        throw std::runtime_error("Nonfinite periodic surface coefficients");
    PointGrid poles(1, nu, 1, nv);
    for (unsigned i = 0; i < nu; ++i)
        for (unsigned j = 0; j < nv; ++j)
            poles.SetValue(i + 1, j + 1,
                           gp_Pnt(coefficients(j, i * 3), coefficients(j, i * 3 + 1),
                                  coefficients(j, i * 3 + 2)));
    RealArray uk(1, nu + 1), vk(1, nv - 2);
    IntArray um(1, nu + 1), vm(1, nv - 2);
    for (unsigned i = 0; i <= nu; ++i) {
        uk.SetValue(i + 1, double(i) / nu);
        um.SetValue(i + 1, 1);
    }
    for (unsigned i = 0; i < nv - 2; ++i) {
        vk.SetValue(i + 1, i == 0 ? 0 : (i == nv - 3 ? 1 : double(i + 1) / (nv - 1)));
        vm.SetValue(i + 1, (i == 0 || i == nv - 3) ? 4 : 1);
    }
    return new Geom_BSplineSurface(poles, uk, vk, um, vm, 3, 3, true, false);
}
bool validate(Patch &patch, const Mesh &mesh, const Topology &top, double tolerance,
              bool inner = true) {
    if (patch.surface.IsNull())
        return false;
    patch.deviation = 0;
    PointProjector projector(patch.surface);
    for (auto v : vertices(mesh, patch.triangles)) {
        const auto &p = projector.perform(mesh.points[v]);
        if (!p.NbPoints())
            return false;
        patch.deviation = std::max(patch.deviation, p.LowerDistance());
        if (patch.deviation > tolerance)
            return false;
    }
    if (inner)
        for (auto t : patch.triangles) {
            const auto &p = projector.perform(top.centers[t]);
            if (!p.NbPoints() || p.LowerDistance() > std::max(tolerance, 0.02))
                return false;
        }
    return true;
}
// P-spline height-field fitting follows upstream freeform.py: weighted rim
// corners, weak facet samples, and a full bending/twist penalty.
Patch height_field(const Mesh &mesh, const Topology &top, const std::vector<unsigned> &triangles,
                   double tolerance, Context &ctx) {
    Patch result;
    result.kind = "freeform";
    result.triangles = triangles;
    auto ids = vertices(mesh, triangles);
    Vec center = Vec::Zero(), normal = Vec::Zero();
    for (auto v : ids)
        center += mesh.points[v];
    center /= double(ids.size());
    for (auto t : triangles)
        normal += top.areas[t] * top.normals[t];
    if (normal.norm() < 1e-12)
        return {};
    normal.normalize();
    for (auto t : triangles)
        if (top.normals[t].dot(normal) < .25881904510252074)
            return {};
    Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
    for (auto v : ids) {
        Vec q = mesh.points[v] - center;
        q -= q.dot(normal) * normal;
        covariance += q * q.transpose();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(covariance);
    Vec e = eigen.eigenvectors().col(2);
    e -= e.dot(normal) * normal;
    e.normalize();
    Vec f = normal.cross(e);
    auto local = [&](const Vec &point) -> Vec {
        Vec q = point - center;
        return Vec(q.dot(e), q.dot(f), q.dot(normal));
    };
    std::map<Edge, unsigned> edges;
    std::unordered_map<unsigned, unsigned> index;
    for (unsigned i = 0; i < ids.size(); ++i)
        index.emplace(ids[i], i);
    std::vector<Vec> vertex_normals(ids.size(), Vec::Zero());
    for (auto t : triangles) {
        auto tri = mesh.triangles[t];
        for (unsigned k = 0; k < 3; ++k) {
            ++edges[Edge(tri[k], tri[(k + 1) % 3])];
            vertex_normals[index.at(tri[k])] += top.areas[t] * top.normals[t];
        }
    }
    std::vector<bool> rim(ids.size(), false);
    for (const auto &edge : edges)
        if (edge.second == 1) {
            rim[index.at(edge.first.a)] = true;
            rim[index.at(edge.first.b)] = true;
        }
    for (auto &n : vertex_normals)
        if (n.norm() > 1e-12)
            n.normalize();
    std::vector<Vec> corners;
    double ulo = 1e100, uhi = -1e100, vlo = 1e100, vhi = -1e100;
    for (auto v : ids) {
        auto s = local(mesh.points[v]);
        corners.push_back(s);
        ulo = std::min(ulo, s.x());
        uhi = std::max(uhi, s.x());
        vlo = std::min(vlo, s.y());
        vhi = std::max(vhi, s.y());
    }
    if (uhi - ulo < 1e-8 || vhi - vlo < 1e-8)
        return {};
    double umargin = .05 * (uhi - ulo) + .5, vmargin = .05 * (vhi - vlo) + .5;
    ulo -= umargin;
    uhi += umargin;
    vlo -= vmargin;
    vhi += vmargin;
    double uspan = uhi - ulo, vspan = vhi - vlo, step = std::max(uspan, vspan) / 4;
    double spacing = std::max(.25, std::sqrt(uspan * vspan / ids.size()));
    unsigned capu = std::min(128u, std::max(1u, unsigned(std::ceil(uspan / spacing))));
    unsigned capv = std::min(128u, std::max(1u, unsigned(std::ceil(vspan / spacing))));
    for (;; step /= 2) {
        if (ctx.expired())
            return {};
        unsigned intervals_u = std::min(capu, std::max(1u, unsigned(std::ceil(uspan / step))));
        unsigned intervals_v = std::min(capv, std::max(1u, unsigned(std::ceil(vspan / step))));
        unsigned mu = intervals_u + 3, mv = intervals_v + 3;
        Basis bu(mu), bv(mv);
        std::vector<Vec> samples = corners;
        std::vector<double> weights;
        for (auto is_rim : rim)
            weights.push_back(is_rim ? 10. : 1.);
        for (auto ti : triangles) {
            auto t = mesh.triangles[ti];
            Vec a = mesh.points[t[0]], b = mesh.points[t[1]], c = mesh.points[t[2]];
            double length = std::max({(a - b).norm(), (b - c).norm(), (c - a).norm()});
            unsigned n = std::clamp(unsigned(std::ceil(2 * length / step)), 2u, 16u);
            for (unsigned i = 0; i <= n; ++i)
                for (unsigned j = 0; j + i <= n; ++j)
                    if (i < n && j < n && i + j > 0) {
                        samples.push_back(local(
                            (double(i) * a + double(j) * b + double(n - i - j) * c) / double(n)));
                        weights.push_back(.001);
                    }
        }
        std::vector<Eigen::Triplet<double>> entries;
        Eigen::VectorXd heights(samples.size());
        for (unsigned row = 0; row < samples.size(); ++row) {
            auto u = bu.at((samples[row].x() - ulo) / uspan),
                 v = bv.at((samples[row].y() - vlo) / vspan);
            double w = std::sqrt(weights[row]);
            for (unsigned i = 0; i < mu; ++i)
                if (u[i])
                    for (unsigned j = 0; j < mv; ++j)
                        if (v[j])
                            entries.emplace_back(row, i * mv + j, w * u[i] * v[j]);
            heights[row] = w * samples[row].z();
        }
        Eigen::SparseMatrix<double> A(samples.size(), mu * mv);
        A.setFromTriplets(entries.begin(), entries.end());
        entries.clear();
        unsigned row = 0;
        double smoothing = 1e-7 * std::max(1., double(samples.size()) / (mu * mv)),
               w = std::sqrt(smoothing);
        for (unsigned i = 0; i < mu; ++i)
            for (unsigned j = 0; j + 2 < mv; ++j) {
                for (unsigned k = 0; k < 3; ++k)
                    entries.emplace_back(row, i * mv + j + k, (k == 1 ? -2 : 1) * w);
                ++row;
            }
        for (unsigned i = 0; i + 2 < mu; ++i)
            for (unsigned j = 0; j < mv; ++j) {
                for (unsigned k = 0; k < 3; ++k)
                    entries.emplace_back(row, (i + k) * mv + j, (k == 1 ? -2 : 1) * w);
                ++row;
            }
        for (unsigned i = 0; i + 1 < mu; ++i)
            for (unsigned j = 0; j + 1 < mv; ++j) {
                double twist = std::sqrt(2 * smoothing);
                entries.emplace_back(row, i * mv + j, twist);
                entries.emplace_back(row, (i + 1) * mv + j, -twist);
                entries.emplace_back(row, i * mv + j + 1, -twist);
                entries.emplace_back(row, (i + 1) * mv + j + 1, twist);
                ++row;
            }
        Eigen::SparseMatrix<double> penalty(row, mu * mv);
        penalty.setFromTriplets(entries.begin(), entries.end());
        Eigen::SparseMatrix<double> system = A.transpose() * A + penalty.transpose() * penalty;
        for (unsigned i = 0; i < mu * mv; ++i)
            system.coeffRef(i, i) += 1e-12;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solve;
        solve.compute(system);
        if (solve.info() == Eigen::Success) {
            Eigen::VectorXd coefficients = solve.solve(A.transpose() * heights);
            if (coefficients.allFinite()) {
                auto height = [&](const Vec &p) {
                    auto u = bu.at((p.x() - ulo) / uspan), v = bv.at((p.y() - vlo) / vspan);
                    double h = 0;
                    for (unsigned i = 0; i < mu; ++i)
                        if (u[i])
                            for (unsigned j = 0; j < mv; ++j)
                                if (v[j])
                                    h += u[i] * v[j] * coefficients[i * mv + j];
                    return h;
                };
                double deviation = 0;
                for (const auto &p : corners)
                    deviation = std::max(deviation, std::abs(height(p) - p.z()));
                bool hugs = deviation <= tolerance;
                if (hugs)
                    for (auto ti : triangles) {
                        auto t = mesh.triangles[ti];
                        double theta = 0, length = 0;
                        for (unsigned k = 0; k < 3; ++k) {
                            theta = std::max(
                                theta,
                                std::acos(std::clamp(
                                    vertex_normals[index.at(t[k])].dot(top.normals[ti]), -1., 1.)));
                            length = std::max(
                                length, (mesh.points[t[k]] - mesh.points[t[(k + 1) % 3]]).norm());
                        }
                        double allowed = length * theta / 4 + 3.5 * tolerance;
                        for (unsigned k = 0; k < 4; ++k) {
                            Vec point =
                                k == 3 ? Vec((mesh.points[t[0]] + mesh.points[t[1]] +
                                              mesh.points[t[2]]) /
                                             3)
                                       : Vec((mesh.points[t[k]] + mesh.points[t[(k + 1) % 3]]) / 2);
                            auto p = local(point);
                            if (std::abs(height(p) - p.z()) > allowed)
                                hugs = false;
                        }
                        if (!hugs)
                            break;
                    }
                if (hugs) {
                    PointGrid poles(1, mu, 1, mv);
                    for (unsigned i = 0; i < mu; ++i)
                        for (unsigned j = 0; j < mv; ++j)
                            poles.SetValue(i + 1, j + 1,
                                           point(center + (ulo + bu.greville(i) * uspan) * e +
                                                 (vlo + bv.greville(j) * vspan) * f +
                                                 coefficients[i * mv + j] * normal));
                    RealArray uk(1, intervals_u + 1), vk(1, intervals_v + 1);
                    IntArray um(1, intervals_u + 1), vm(1, intervals_v + 1);
                    for (unsigned i = 0; i <= intervals_u; ++i) {
                        uk.SetValue(i + 1, double(i) / intervals_u);
                        um.SetValue(i + 1, i == 0 || i == intervals_u ? 4 : 1);
                    }
                    for (unsigned j = 0; j <= intervals_v; ++j) {
                        vk.SetValue(j + 1, double(j) / intervals_v);
                        vm.SetValue(j + 1, j == 0 || j == intervals_v ? 4 : 1);
                    }
                    result.surface = new Geom_BSplineSurface(poles, uk, vk, um, vm, 3, 3);
                    result.origin = center;
                    result.axis = normal;
                    result.deviation = deviation;
                    return result;
                }
            }
        }
        if (intervals_u >= capu && intervals_v >= capv)
            return {};
    }
}
Patch extrusion(const Mesh &mesh, const Topology &top, const std::vector<unsigned> &triangles,
                double tolerance) {
    Patch p;
    p.kind = "extrusion";
    p.triangles = triangles;
    Eigen::Matrix3d NN = Eigen::Matrix3d::Zero();
    for (auto t : triangles)
        NN += top.normals[t] * top.normals[t].transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(NN);
    Vec axis = eigen.eigenvectors().col(0);
    for (auto t : triangles)
        if (std::abs(axis.dot(top.normals[t])) > 0.02)
            return {};
    auto ids = vertices(mesh, triangles);
    Vec center = Vec::Zero();
    for (auto i : ids)
        center += mesh.points[i];
    center /= double(ids.size());
    Vec e = axis.unitOrthogonal(), f = axis.cross(e);
    std::vector<Vec> profile;
    std::map<std::pair<long long, long long>, unsigned> seen;
    std::unordered_map<unsigned, unsigned> vertex_profile;
    for (auto i : ids) {
        Vec v = mesh.points[i] - center;
        auto key = std::make_pair(std::llround(v.dot(e) * 1e4), std::llround(v.dot(f) * 1e4));
        auto inserted = seen.emplace(key, unsigned(profile.size()));
        if (inserted.second)
            profile.push_back(v - v.dot(axis) * axis);
        vertex_profile.emplace(i, inserted.first->second);
    }
    if (profile.size() < 4)
        return {};
    // Preserve the profile's mesh connectivity. A nearest-neighbor walk can
    // connect across a concavity or leave a closed profile with an open seam.
    std::vector<std::set<unsigned>> graph(profile.size());
    for (auto ti : triangles) {
        auto t = mesh.triangles[ti];
        for (unsigned k = 0; k < 3; ++k) {
            auto a = vertex_profile.at(t[k]), b = vertex_profile.at(t[(k + 1) % 3]);
            if (a != b) {
                graph[a].insert(b);
                graph[b].insert(a);
            }
        }
    }
    unsigned start = 0, endpoints = 0;
    for (unsigned i = 0; i < graph.size(); ++i) {
        if (graph[i].empty() || graph[i].size() > 2)
            return {};
        if (graph[i].size() == 1) {
            start = i;
            ++endpoints;
        }
    }
    if (endpoints != 0 && endpoints != 2)
        return {};
    bool periodic = endpoints == 0;
    std::vector<unsigned> ordered{start};
    std::vector<bool> visited(profile.size(), false);
    visited[start] = true;
    while (ordered.size() < profile.size()) {
        unsigned next = ~0u;
        for (auto neighbor : graph[ordered.back()])
            if (!visited[neighbor]) {
                next = neighbor;
                break;
            }
        if (next == ~0u)
            return {};
        visited[next] = true;
        ordered.push_back(next);
    }
    if (periodic && !graph[ordered.back()].count(start))
        return {};
    Handle(PointHandleArray) points = new PointHandleArray(1, ordered.size());
    for (unsigned i = 0; i < ordered.size(); ++i)
        points->SetValue(i + 1, point(center + profile[ordered[i]]));
    GeomAPI_Interpolate fit(points, periodic, 1e-8);
    fit.Perform();
    if (!fit.IsDone())
        return {};
    p.axis = axis;
    p.origin = center;
    p.surface = new Geom_SurfaceOfLinearExtrusion(fit.Curve(), direction(axis));
    return validate(p, mesh, top, tolerance) ? p : Patch{};
}
Patch pipe(const Mesh &mesh, const Topology &top, const std::vector<unsigned> &triangles,
           const std::vector<double> &radii, double tolerance) {
    for (double radius : radii)
        for (double sign : {-1.0, 1.0}) {
            // Offset facet samples toward the rolling-ball center. Cross-section
            // centers cluster into a narrow curve when a constant-radius tube fits.
            std::vector<Vec> centers;
            for (auto t : triangles)
                centers.push_back(top.centers[t] + sign * radius * top.normals[t]);
            Vec mean = Vec::Zero();
            for (const auto &c : centers)
                mean += c;
            mean /= double(centers.size());
            Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
            for (const auto &c : centers) {
                Vec q = c - mean;
                covariance += q * q.transpose();
            }
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(covariance);
            Vec axis = eigen.eigenvectors().col(2);
            std::sort(centers.begin(), centers.end(),
                      [&](const Vec &a, const Vec &b) { return axis.dot(a) < axis.dot(b); });
            unsigned bins =
                std::min<unsigned>(32, std::max<unsigned>(4, unsigned(centers.size() / 8)));
            if (centers.size() < bins)
                return {};
            std::vector<Vec> spine;
            for (unsigned i = 0; i < bins; ++i) {
                unsigned begin = i * centers.size() / bins, end = (i + 1) * centers.size() / bins;
                Vec c = Vec::Zero();
                for (unsigned j = begin; j < end; ++j)
                    c += centers[j];
                spine.push_back(c / double(end - begin));
            }
            PointArray points(1, bins + 2);
            points.SetValue(1, point(spine[0] - 1.5 * (spine[1] - spine[0])));
            for (unsigned i = 0; i < bins; ++i)
                points.SetValue(i + 2, point(spine[i]));
            points.SetValue(bins + 2,
                            point(spine.back() + 1.5 * (spine.back() - spine[spine.size() - 2])));
            GeomAPI_PointsToBSpline fit(points, 3, 8, GeomAbs_C2, tolerance);
            if (!fit.IsDone())
                continue;
            GeomFill_Pipe make(fit.Curve(), radius);
            make.Perform(tolerance * 0.2, false, GeomAbs_C1, 8, 30);
            Patch p;
            p.kind = "pipe";
            p.triangles = triangles;
            p.radius = radius;
            p.surface = make.Surface();
            if (validate(p, mesh, top, tolerance))
                return p;
        }
    return {};
}
Patch thread(const Mesh &mesh, const Topology &top, const std::vector<unsigned> &triangles,
             const Vec &axis, double tolerance) {
    auto ids = vertices(mesh, triangles);
    if (ids.size() < 64)
        return {};
    Vec center = Vec::Zero();
    for (auto i : ids)
        center += mesh.points[i];
    center /= double(ids.size());
    Vec e = axis.unitOrthogonal(), f = axis.cross(e);
    struct Sample {
        double u, z, r;
    };
    std::vector<Sample> samples;
    double lo = 1e100, hi = -1e100, rlo = 1e100, rhi = 0;
    for (auto i : ids) {
        Vec q = mesh.points[i] - center;
        double z = q.dot(axis), x = q.dot(e), y = q.dot(f), r = std::hypot(x, y);
        samples.push_back({std::atan2(y, x), z, r});
        lo = std::min(lo, z);
        hi = std::max(hi, z);
        rlo = std::min(rlo, r);
        rhi = std::max(rhi, r);
    }
    if (hi - lo < 1e-4 || rhi - rlo < 0.02 || rlo < 1e-4)
        return {};
    struct Ring {
        unsigned count = 0;
        double lo = 1e100, hi = -1e100, sum = 0;
    };
    std::array<Ring, 8> rings;
    for (auto s : samples) {
        auto &ring = rings[std::min(7u, unsigned((s.z - lo) / (hi - lo) * 8))];
        ++ring.count;
        ring.lo = std::min(ring.lo, s.r);
        ring.hi = std::max(ring.hi, s.r);
        ring.sum += s.r;
    }
    unsigned populated = 0;
    double span = 0, meanlo = 1e100, meanhi = -1e100;
    for (const auto &ring : rings)
        if (ring.count >= 8) {
            ++populated;
            span += ring.hi - ring.lo;
            double average = ring.sum / ring.count;
            meanlo = std::min(meanlo, average);
            meanhi = std::max(meanhi, average);
        }
    if (populated < 4 || span / populated < (rhi - rlo) * .5 || meanhi - meanlo > (rhi - rlo) * .35)
        return {};
    std::vector<double> axial_edges;
    for (auto t : triangles) {
        auto tri = mesh.triangles[t];
        for (unsigned k = 0; k < 3; ++k) {
            double dz = std::abs((mesh.points[tri[k]] - mesh.points[tri[(k + 1) % 3]]).dot(axis));
            if (dz > (hi - lo) * 1e-6)
                axial_edges.push_back(dz);
        }
    }
    if (axial_edges.empty())
        return {};
    auto middle = axial_edges.begin() + axial_edges.size() / 2;
    std::nth_element(axial_edges.begin(), middle, axial_edges.end());
    double minimum_pitch = 2.0 * (*middle);
    // Search pitch/handedness against a periodic radial profile. A fit is only
    // accepted if a single helical phase explains every mesh corner.
    double best = 1e100, bestpitch = 0;
    Eigen::VectorXd bestprofile;
    constexpr unsigned bins = 32;
    for (int hand : {-1, 1})
        for (unsigned k = 2; k <= 512; ++k) {
            double pitch = hand * (hi - lo) * 4.0 / k;
            if (std::abs(pitch) < minimum_pitch * .999)
                continue;
            Eigen::MatrixXd A = Eigen::MatrixXd::Zero(samples.size(), bins);
            Eigen::VectorXd b(samples.size());
            for (unsigned i = 0; i < samples.size(); ++i) {
                auto s = samples[i];
                double phase = s.z / pitch - s.u / (2 * pi);
                phase -= std::floor(phase);
                double x = phase * bins;
                unsigned j = unsigned(x) % bins;
                double w = x - std::floor(x);
                A(i, j) = 1 - w;
                A(i, (j + 1) % bins) = w;
                b[i] = s.r;
            }
            Eigen::VectorXd profile = A.colPivHouseholderQr().solve(b);
            double error = (A * profile - b).cwiseAbs().maxCoeff();
            if (error < best) {
                best = error;
                bestpitch = pitch;
                bestprofile = profile;
            }
        }
    if (bestpitch == 0)
        return {};
    trace(("thread coarse pitch=" + std::to_string(bestpitch) + " residual=" + std::to_string(best))
              .c_str());
    // Refine the continuous pitch, transverse center, and axis. The periodic
    // profile is solved linearly inside each residual evaluation; this removes
    // the coarse pitch-grid quantization and bias from incomplete axial turns.
    const Vec base = center, base_e = e, base_f = f;
    const double hand = bestpitch > 0 ? 1 : -1;
    Eigen::VectorXd initial(5);
    initial << 0, 0, 0, 0, std::log(std::abs(bestpitch));
    auto profiled = [&](const Eigen::VectorXd &x,
                        Eigen::VectorXd *coefficients) -> Eigen::VectorXd {
        Vec a = (axis + x[2] * base_e + x[3] * base_f).normalized(),
            origin = base + x[0] * base_e + x[1] * base_f;
        Vec u = (base_e - base_e.dot(a) * a).normalized(), v = a.cross(u);
        double pitch = hand * std::exp(x[4]);
        if (!a.allFinite() || !origin.allFinite() || !u.allFinite() || !std::isfinite(pitch) ||
            std::abs(pitch) < std::max(1e-10, minimum_pitch * .999) ||
            std::abs(pitch) > 2 * (hi - lo))
            return Eigen::VectorXd::Constant(ids.size(), std::numeric_limits<double>::infinity());
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(ids.size(), bins);
        Eigen::VectorXd b(ids.size());
        for (unsigned i = 0; i < ids.size(); ++i) {
            Vec q = mesh.points[ids[i]] - origin;
            double z = q.dot(a), angle = std::atan2(q.dot(v), q.dot(u));
            double phase = z / pitch - angle / (2 * pi);
            if (!std::isfinite(phase))
                return Eigen::VectorXd::Constant(ids.size(),
                                                 std::numeric_limits<double>::infinity());
            phase -= std::floor(phase);
            double t = phase * bins;
            unsigned j = unsigned(t) % bins;
            double w = t - std::floor(t);
            A(i, j) = 1 - w;
            A(i, (j + 1) % bins) = w;
            b[i] = (q - z * a).norm();
        }
        Eigen::VectorXd profile = A.colPivHouseholderQr().solve(b);
        if (coefficients)
            *coefficients = profile;
        return A * profile - b;
    };
    auto optimized = least_squares(
        [&](const Eigen::VectorXd &x) -> Eigen::VectorXd { return profiled(x, nullptr); }, initial,
        40, tolerance);
    Vec fitted_axis = axis;
    if (optimized.x.allFinite()) {
        Eigen::VectorXd profile;
        double error = profiled(optimized.x, &profile).cwiseAbs().maxCoeff();
        if (error < best) {
            best = error;
            bestpitch = hand * std::exp(optimized.x[4]);
            bestprofile = std::move(profile);
            center = base + optimized.x[0] * base_e + optimized.x[1] * base_f;
            fitted_axis = (axis + optimized.x[2] * base_e + optimized.x[3] * base_f).normalized();
            e = (base_e - base_e.dot(fitted_axis) * fitted_axis).normalized();
            f = fitted_axis.cross(e);
            lo = 1e100;
            hi = -1e100;
            rlo = 1e100;
            rhi = 0;
            for (auto id : ids) {
                Vec q = mesh.points[id] - center;
                double z = q.dot(fitted_axis), r = (q - z * fitted_axis).norm();
                lo = std::min(lo, z);
                hi = std::max(hi, z);
                rlo = std::min(rlo, r);
                rhi = std::max(rhi, r);
            }
        }
    }
    trace(
        ("thread refined pitch=" + std::to_string(bestpitch) + " residual=" + std::to_string(best))
            .c_str());
    if (best > tolerance)
        return {};
    if (bestprofile.minCoeff() <= 0 || bestprofile.maxCoeff() > rhi + std::max(.1, rhi - rlo))
        return {};
    double margin = std::max(4 * tolerance, std::abs(bestpitch) / 32);
    lo -= margin;
    hi += margin;
    unsigned nu = 48,
             nv = std::min<unsigned>(
                 4097, std::max<unsigned>(33, unsigned((hi - lo) / std::abs(bestpitch) * 24) + 1));
    PointGrid grid(1, nu, 1, nv);
    for (unsigned i = 0; i < nu; ++i)
        for (unsigned j = 0; j < nv; ++j) {
            double u = 2 * pi * i / nu, z = lo + (hi - lo) * j / (nv - 1);
            double phase = z / bestpitch - u / (2 * pi);
            phase -= std::floor(phase);
            double x = phase * bins;
            unsigned k = unsigned(x) % bins;
            double w = x - std::floor(x),
                   r = (1 - w) * bestprofile[k] + w * bestprofile[(k + 1) % bins];
            grid.SetValue(
                i + 1, j + 1,
                point(center + z * fitted_axis + r * (std::cos(u) * e + std::sin(u) * f)));
        }
    Patch p;
    p.kind = bestpitch > 0 ? "thread_right" : "thread_left";
    p.origin = center;
    p.axis = fitted_axis;
    p.radius = rhi;
    p.secondary = bestpitch;
    p.triangles = triangles;
    p.surface = periodic_grid(grid);
    std::size_t stride = std::max<std::size_t>(1, triangles.size() / 512);
    PointProjector projector(p.surface);
    for (std::size_t i = 0; i < triangles.size(); i += stride) {
        auto t = triangles[i];
        const auto &projection = projector.perform(top.centers[t]);
        if (!projection.NbPoints())
            return {};
        double u, v;
        projection.LowerDistanceParameters(u, v);
        gp_Pnt q;
        gp_Vec du, dv;
        p.surface->D1(u, v, q, du, dv);
        gp_Vec cross = du.Crossed(dv);
        Vec n(cross.X(), cross.Y(), cross.Z());
        if (n.norm() < 1e-10 || std::abs(n.normalized().dot(top.normals[t])) < .90) {
            return {};
        }
    }
    if (validate(p, mesh, top, tolerance, false))
        return p;
    trace(("thread surface rejected deviation=" + std::to_string(p.deviation)).c_str());
    return {};
}
} // namespace
Patch thread_patch(const Mesh &mesh, const Topology &top, const std::vector<unsigned> &triangles,
                   const Vec &axis, double tolerance) {
    try {
        return thread(mesh, top, triangles, axis, tolerance);
    } catch (const Standard_Failure &) {
        return {};
    }
}
std::vector<Patch> advanced_patches(const Mesh &mesh, const Topology &top,
                                    const std::vector<int> &owner,
                                    const std::vector<Patch> &recognized, Context &ctx) {
    if (!ctx.options.curves)
        return {};
    const double tolerance = std::max(1e-5, std::min(0.002, ctx.options.accuracy * .02));
    std::vector<double> radii;
    std::vector<Vec> axes{Vec::UnitX(), Vec::UnitY(), Vec::UnitZ()};
    for (const auto &p : recognized) {
        if (p.kind == "cylinder" || p.kind == "torus") {
            radii.push_back(p.kind == "torus" ? p.secondary : p.radius);
            axes.push_back(p.axis);
        }
    }
    std::vector<Vec> unique_axes;
    for (const auto &axis : axes) {
        if (axis.norm() < 1e-12)
            continue;
        Vec normalized = axis.normalized();
        if (std::none_of(unique_axes.begin(), unique_axes.end(), [&](const Vec &known) {
                return std::abs(known.dot(normalized)) > 1 - 1e-12;
            }))
            unique_axes.push_back(normalized);
    }
    axes = std::move(unique_axes);
    std::sort(radii.begin(), radii.end());
    radii.erase(std::unique(radii.begin(), radii.end(),
                            [](double a, double b) { return std::abs(a - b) < 1e-4; }),
                radii.end());
    if (radii.size() > 8)
        radii.resize(8);
    std::vector<std::vector<unsigned>> graph(mesh.triangles.size());
    for (unsigned i = 0; i < owner.size(); ++i)
        if (owner[i] < 0)
            for (auto j : top.neighbors[i])
                if (owner[j] < 0 && top.normals[i].dot(top.normals[j]) > 0.72)
                    graph[i].push_back(j);
    std::vector<Patch> out;
    std::function<void(const std::vector<unsigned> &, unsigned)> fit_region;
    fit_region = [&](const std::vector<unsigned> &ids, unsigned depth) {
        if (ids.size() < 4 || ctx.expired())
            return;
        try {
            Patch p;
            // Narrow coplanar regions can be deferred by the analytic pass.
            // Keep their exact plane instead of fitting a degenerate extrusion.
            Vec origin = mesh.points[mesh.triangles[ids.front()][0]],
                axis = top.normals[ids.front()];
            bool planar = true;
            for (auto i : ids) {
                if (axis.dot(top.normals[i]) < .9999999)
                    planar = false;
                for (auto v : mesh.triangles[i])
                    if (std::abs((mesh.points[v] - origin).dot(axis)) > 1e-5)
                        planar = false;
            }
            if (planar) {
                p.kind = "plane";
                p.origin = origin;
                p.axis = axis;
                p.triangles = ids;
                p.surface = new Geom_Plane(gp_Ax3(point(origin), direction(axis)));
                out.push_back(std::move(p));
                return;
            }
            if (ids.size() >= 64 && depth == 0)
                for (const auto &axis : axes) {
                    p = thread(mesh, top, ids, axis, tolerance);
                    if (!p.surface.IsNull()) {
                        out.push_back(std::move(p));
                        return;
                    }
                    if (ctx.expired())
                        return;
                }
            if (!ctx.options.blends)
                return;
            p = extrusion(mesh, top, ids, tolerance);
            if (!p.surface.IsNull()) {
                out.push_back(std::move(p));
                return;
            }
            if (!radii.empty()) {
                p = pipe(mesh, top, ids, radii, tolerance);
                if (!p.surface.IsNull()) {
                    out.push_back(std::move(p));
                    return;
                }
            }
            p = height_field(mesh, top, ids, std::max(1e-5, tolerance * .5), ctx);
            if (p.surface.IsNull() && !ctx.expired())
                p = height_field(mesh, top, ids, tolerance, ctx);
            if (!p.surface.IsNull()) {
                if (depth > 0)
                    p.kind = "blend";
                out.push_back(std::move(p));
                return;
            }
        } catch (const Standard_Failure &) { /* Reject this candidate and retain mesh geometry. */
        }
        if (depth >= 4 || ids.size() < 16)
            return;
        // Split by facet direction, retaining connected subregions only.
        Vec average = Vec::Zero();
        for (auto i : ids)
            average += top.normals[i];
        if (average.norm() < 1e-10)
            average = Vec::UnitZ();
        else
            average.normalize();
        Vec e = average.unitOrthogonal();
        std::vector<unsigned> sorted = ids;
        std::sort(sorted.begin(), sorted.end(), [&](unsigned a, unsigned b) {
            return top.normals[a].dot(e) < top.normals[b].dot(e);
        });
        auto middle = sorted.begin() + sorted.size() / 2;
        std::vector<unsigned> a(sorted.begin(), middle), b(middle, sorted.end());
        fit_region(a, depth + 1);
        fit_region(b, depth + 1);
    };
    for (const auto &group : components(graph))
        if (group.size() >= 4 && owner[group.front()] < 0)
            fit_region(group, 0);
    return out;
}
} // namespace stl2curves::detail

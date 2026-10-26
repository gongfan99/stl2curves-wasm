// Surface recognition and robust fitting derived from the stl2curves pipeline.
#include "internal.hpp"
#include "projection.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <Eigen/QR>
#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <gp_Ax3.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>

namespace stl2curves::detail {
namespace {
gp_Pnt point(const Vec &v) {
    return gp_Pnt(v.x(), v.y(), v.z());
}
gp_Dir direction(const Vec &v) {
    return gp_Dir(v.x(), v.y(), v.z());
}
Vec mean(const Mesh &m, const std::vector<unsigned> &ids) {
    Vec c = Vec::Zero();
    for (auto i : ids)
        c += m.points[i];
    return c / double(ids.size());
}
double deviation(const Patch &p, const Vec &v) {
    Vec q = v - p.origin;
    if (p.kind == "plane")
        return std::abs(q.dot(p.axis));
    if (p.kind == "sphere")
        return std::abs(q.norm() - p.radius);
    double z = q.dot(p.axis), rho = (q - z * p.axis).norm();
    if (p.kind == "cylinder")
        return std::abs(rho - p.radius);
    if (p.kind == "cone")
        return std::abs(rho - p.radius - z * p.secondary) /
               std::sqrt(1 + p.secondary * p.secondary);
    return std::abs(std::hypot(rho - p.radius, z) - p.secondary);
}
Vec normal(const Patch &p, const Vec &v) {
    Vec q = v - p.origin;
    if (p.kind == "plane")
        return p.axis;
    if (p.kind == "sphere")
        return q.normalized();
    double z = q.dot(p.axis);
    Vec radial = q - z * p.axis;
    double rho = radial.norm();
    if (rho < 1e-12)
        return p.axis;
    if (p.kind == "cylinder")
        return radial / rho;
    if (p.kind == "cone")
        return (radial / rho - p.secondary * p.axis).normalized();
    return ((rho - p.radius) * radial / rho + z * p.axis).normalized();
}
void surface(Patch &p) {
    gp_Ax3 frame(point(p.origin), direction(p.axis));
    if (p.kind == "plane")
        p.surface = new Geom_Plane(frame);
    else if (p.kind == "cylinder")
        p.surface = new Geom_CylindricalSurface(frame, p.radius);
    else if (p.kind == "sphere")
        p.surface = new Geom_SphericalSurface(frame, p.radius);
    else if (p.kind == "cone") {
        if (std::abs(p.secondary) < 1e-8) {
            p.kind = "cylinder";
            p.surface = new Geom_CylindricalSurface(frame, p.radius);
        } else
            p.surface = new Geom_ConicalSurface(frame, std::atan(p.secondary), p.radius);
    } else if (p.kind == "torus")
        p.surface = new Geom_ToroidalSurface(frame, p.radius, p.secondary);
}
bool normal_agrees(const Patch &p, const Topology &top, unsigned t) {
    auto n = normal(p, top.centers[t]);
    return n.allFinite() && std::abs(n.dot(top.normals[t])) >= .96;
}
bool holds(Patch &p, const Mesh &m, const Topology &top, const std::vector<unsigned> &tris,
           double tolerance, bool require_turn = true) {
    if (!p.origin.allFinite() || !p.axis.allFinite() || !std::isfinite(p.radius) ||
        !std::isfinite(p.secondary))
        return false;
    p.deviation = 0;
    for (auto i : vertices(m, tris))
        p.deviation = std::max(p.deviation, deviation(p, m.points[i]));
    if (!std::isfinite(p.deviation) || p.deviation > tolerance)
        return false;
    // Normals at facet centroids prevent fitting a broad sphere to a planar strip.
    for (auto t : tris) {
        if (!normal_agrees(p, top, t))
            return false;
    }
    if (require_turn && p.kind != "plane" && tris.size() > 1) {
        Vec first = normal(p, top.centers[tris.front()]);
        double spread = 0;
        for (auto t : tris)
            spread = std::max(spread, (normal(p, top.centers[t]) - first).norm());
        if (spread < 0.1)
            return false;
    }
    return true;
}
Patch fit_plane(const Mesh &m, const std::vector<unsigned> &triangles) {
    Patch p;
    p.kind = "plane";
    p.triangles = triangles;
    auto ids = vertices(m, triangles);
    p.origin = mean(m, ids);
    Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
    for (auto i : ids) {
        Vec q = m.points[i] - p.origin;
        covariance += q * q.transpose();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(covariance);
    p.axis = eig.eigenvectors().col(0);
    return p;
}
std::vector<Patch> fits(const Mesh &m, const Topology &top, const std::vector<unsigned> &triangles,
                        bool tori) {
    auto ids = vertices(m, triangles);
    if (ids.size() < 4)
        return {};
    Vec center = mean(m, ids);
    std::vector<Patch> candidates;
    candidates.push_back(fit_plane(m, triangles));
    Eigen::MatrixXd A(ids.size(), 4);
    Eigen::VectorXd b(ids.size());
    for (unsigned i = 0; i < ids.size(); ++i) {
        Vec v = m.points[ids[i]] - center;
        A.row(i) << 2 * v.x(), 2 * v.y(), 2 * v.z(), 1;
        b[i] = v.squaredNorm();
    }
    Eigen::Vector4d s = A.colPivHouseholderQr().solve(b);
    if (s.allFinite() && s[3] + s.head<3>().squaredNorm() > 0) {
        Patch p;
        p.kind = "sphere";
        p.origin = center + s.head<3>();
        p.radius = std::sqrt(s[3] + s.head<3>().squaredNorm());
        p.triangles = triangles;
        candidates.push_back(p);
    }
    // The planar facets of a tessellated cone meet at its apex. Solving for that
    // apex first avoids the biased circle initialization of a short countersink.
    Eigen::MatrixXd planes(triangles.size(), 3);
    Eigen::VectorXd offsets(triangles.size());
    for (unsigned i = 0; i < triangles.size(); ++i) {
        planes.row(i) = top.normals[triangles[i]].transpose();
        offsets[i] = top.normals[triangles[i]].dot(top.centers[triangles[i]] - center);
    }
    auto apex_solver = planes.colPivHouseholderQr();
    if (apex_solver.rank() == 3) {
        Vec relative = apex_solver.solve(offsets);
        double extent = 0;
        for (auto i : ids)
            extent = std::max(extent, (m.points[i] - center).norm());
        if (relative.allFinite() && relative.norm() < std::max(1.0, extent) * 1e5 &&
            (planes * relative - offsets).cwiseAbs().maxCoeff() < 1e-4) {
            Vec apex = center + relative, average = Vec::Zero();
            std::vector<Vec> rays;
            for (auto i : ids) {
                Vec ray = m.points[i] - apex;
                if (ray.norm() > 1e-8 * std::max(1.0, extent)) {
                    ray.normalize();
                    rays.push_back(ray);
                    average += ray;
                }
            }
            if (rays.size() >= 3) {
                average /= double(rays.size());
                Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
                for (const auto &ray : rays) {
                    Vec q = ray - average;
                    covariance += q * q.transpose();
                }
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(covariance);
                Vec axis = eigen.eigenvectors().col(0);
                if (axis.dot(average) < 0)
                    axis = -axis;
                double cosine = axis.dot(average);
                if (eigen.eigenvalues()[1] > 1e-12 * eigen.eigenvalues()[2] && cosine > 0.01 &&
                    cosine < 1 - 1e-10) {
                    double slope = std::sqrt(1 - cosine * cosine) / cosine,
                           z = (center - apex).dot(axis);
                    Patch p;
                    p.kind = "cone";
                    p.origin = apex + z * axis;
                    p.axis = axis;
                    p.radius = z * slope;
                    p.secondary = slope;
                    p.triangles = triangles;
                    if (p.radius > 1e-6)
                        candidates.push_back(p);
                }
            }
        }
    }
    Eigen::Matrix3d NN = Eigen::Matrix3d::Zero();
    Vec avg = Vec::Zero();
    for (auto t : triangles) {
        NN += top.normals[t] * top.normals[t].transpose();
        avg += top.normals[t];
    }
    avg /= double(triangles.size());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(NN);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> coneEig(NN - double(triangles.size()) * avg *
                                                                    avg.transpose());
    std::vector<Vec> axes{eig.eigenvectors().col(0), coneEig.eigenvectors().col(0)};
    if (tori || triangles.size() <= 8) {
        axes.push_back(Vec::UnitX());
        axes.push_back(Vec::UnitY());
        axes.push_back(Vec::UnitZ());
    }
    for (auto axis : axes) {
        Vec e = axis.unitOrthogonal(), f = axis.cross(e);
        Eigen::MatrixXd C(ids.size(), 3);
        Eigen::VectorXd rr(ids.size());
        for (unsigned i = 0; i < ids.size(); ++i) {
            Vec q = m.points[ids[i]] - center;
            double x = q.dot(e), y = q.dot(f);
            C.row(i) << 2 * x, 2 * y, 1;
            rr[i] = x * x + y * y;
        }
        Eigen::Vector3d circle = C.colPivHouseholderQr().solve(rr);
        double rad2 = circle[2] + circle.head<2>().squaredNorm();
        if (!circle.allFinite() || rad2 <= 0)
            continue;
        Patch p;
        p.kind = "cylinder";
        p.axis = axis;
        p.origin = center + circle[0] * e + circle[1] * f;
        p.radius = std::sqrt(rad2);
        p.triangles = triangles;
        // Five identifiable degrees of freedom: two transverse offsets, two
        // axis tilts, and radius. Removing the axial/normalization gauges avoids
        // drifting an otherwise correct center thousands of millimeters away.
        Eigen::VectorXd initial(5);
        initial << circle[0], circle[1], 0, 0, p.radius;
        auto residual = [&](const Eigen::VectorXd &x) -> Eigen::VectorXd {
            Vec a = (axis + x[2] * e + x[3] * f).normalized(),
                origin = center + x[0] * e + x[1] * f;
            Eigen::VectorXd r(ids.size());
            for (unsigned i = 0; i < ids.size(); ++i) {
                Vec q = m.points[ids[i]] - origin;
                r[i] = (q - q.dot(a) * a).norm() - x[4];
            }
            return r;
        };
        auto solution = least_squares(residual, initial, 40);
        if (solution.x.allFinite() && solution.x[4] > 1e-6) {
            p.origin = center + solution.x[0] * e + solution.x[1] * f;
            p.axis = (axis + solution.x[2] * e + solution.x[3] * f).normalized();
            p.origin -= (p.origin - center).dot(p.axis) * p.axis;
            p.radius = solution.x[4];
            candidates.push_back(p);
        }
        // Cone fit uses the normal-space plane for its axis and a robust radial profile.
        p.kind = "cone";
        p.axis = axis;
        p.origin = center + circle[0] * e + circle[1] * f;
        Eigen::MatrixXd Z(ids.size(), 2);
        Eigen::VectorXd R(ids.size());
        for (unsigned i = 0; i < ids.size(); ++i) {
            Vec q = m.points[ids[i]] - p.origin;
            double z = q.dot(axis);
            Z.row(i) << 1, z;
            R[i] = (q - z * axis).norm();
        }
        Eigen::Vector2d profile = Z.colPivHouseholderQr().solve(R);
        if (profile[0] > 0 && std::abs(profile[1]) > 1e-5 && std::abs(profile[1]) < 100) {
            Eigen::VectorXd x(6);
            x << circle[0], circle[1], 0, 0, profile[0], profile[1];
            auto coneResidual = [&](const Eigen::VectorXd &v) -> Eigen::VectorXd {
                Vec a = (axis + v[2] * e + v[3] * f).normalized(),
                    origin = center + v[0] * e + v[1] * f;
                Eigen::VectorXd r(ids.size());
                for (unsigned i = 0; i < ids.size(); ++i) {
                    Vec q = m.points[ids[i]] - origin;
                    double z = q.dot(a);
                    r[i] = ((q - z * a).norm() - v[4] - z * v[5]) / std::sqrt(1 + v[5] * v[5]);
                }
                return r;
            };
            auto fit = least_squares(coneResidual, x, 50);
            if (fit.x.allFinite() && fit.x[4] > 0) {
                p.origin = center + fit.x[0] * e + fit.x[1] * f;
                p.axis = (axis + fit.x[2] * e + fit.x[3] * f).normalized();
                p.radius = fit.x[4];
                p.secondary = fit.x[5];
                candidates.push_back(p);
            }
        }
        if (tori && ids.size() >= 12) {
            // In a meridian plane a torus is a circle in (rho,z).
            Eigen::MatrixXd T(ids.size(), 3);
            Eigen::VectorXd V(ids.size());
            for (unsigned i = 0; i < ids.size(); ++i) {
                Vec q = m.points[ids[i]] - center;
                double z = q.dot(axis), rho = (q - z * axis).norm();
                T.row(i) << 2 * rho, 2 * z, 1;
                V[i] = rho * rho + z * z;
            }
            Eigen::Vector3d ring = T.colPivHouseholderQr().solve(V);
            double r2 = ring[2] + ring.head<2>().squaredNorm();
            if (ring[0] <= 0 || r2 <= 0)
                continue;
            Eigen::VectorXd x(7);
            x << center + ring[1] * axis, 0, 0, ring[0], std::sqrt(r2);
            auto torusResidual = [&](const Eigen::VectorXd &v) -> Eigen::VectorXd {
                Vec a = (axis + v[3] * e + v[4] * f).normalized();
                Eigen::VectorXd r(ids.size());
                for (unsigned i = 0; i < ids.size(); ++i) {
                    Vec q = m.points[ids[i]] - v.head<3>();
                    double z = q.dot(a), rho = (q - z * a).norm();
                    r[i] = std::hypot(rho - v[5], z) - v[6];
                }
                return r;
            };
            auto fit = least_squares(torusResidual, x, 50);
            if (fit.x.allFinite() && fit.x[5] > 1e-5 && fit.x[6] > 1e-5 && fit.x[5] >= fit.x[6]) {
                p.kind = "torus";
                p.origin = fit.x.head<3>();
                p.axis = (axis + fit.x[3] * e + fit.x[4] * f).normalized();
                p.radius = fit.x[5];
                p.secondary = fit.x[6];
                candidates.push_back(p);
            }
        }
    }
    auto fitted = candidates;
    for (auto p : fitted)
        if (p.kind != "plane") {
            Eigen::Index component;
            p.axis.cwiseAbs().maxCoeff(&component);
            if (std::abs(p.axis[component]) > .999999) {
                Vec axis = Vec::Zero();
                axis[component] = p.axis[component] > 0 ? 1 : -1;
                p.axis = axis;
            }
            Vec shift;
            for (unsigned k = 0; k < 3; ++k)
                shift[k] = std::round(p.origin[k] * 100) / 100 - p.origin[k];
            if (p.kind == "cylinder")
                shift -= shift.dot(p.axis) * p.axis;
            if (shift.norm() > .005)
                continue;
            p.origin += shift;
            if (p.kind != "cone") {
                double r = std::round(p.radius * 100) / 100;
                if (r > 0 && std::abs(r - p.radius) < .005)
                    p.radius = r;
            }
            if (p.kind == "torus") {
                double r = std::round(p.secondary * 100) / 100;
                if (r > 0 && std::abs(r - p.secondary) < .005)
                    p.secondary = r;
            }
            if (holds(p, m, top, triangles, 1e-4, false))
                candidates.push_back(std::move(p));
        }
    // Two axial rings of a cylinder are also cospherical. Prefer the lower-order
    // developable model before accepting that ambiguous sphere fit.
    auto rank = [](const std::string &kind) {
        if (kind == "plane")
            return 0;
        if (kind == "cylinder")
            return 1;
        if (kind == "cone")
            return 2;
        if (kind == "sphere")
            return 3;
        return 4;
    };
    std::stable_sort(candidates.begin(), candidates.end(),
                     [&](const Patch &a, const Patch &b) { return rank(a.kind) < rank(b.kind); });
    return candidates;
}
std::vector<Patch> rolling_fillets(const Mesh &mesh, const Topology &top,
                                   const std::vector<unsigned> &triangles,
                                   const std::vector<unsigned> &region,
                                   const std::vector<int> &owner, const std::vector<Patch> &known) {
    std::set<unsigned> neighbors;
    for (auto i : region) {
        if (owner[i] >= 0)
            neighbors.insert(unsigned(owner[i]));
        for (auto j : top.neighbors[i])
            if (owner[j] >= 0)
                neighbors.insert(unsigned(owner[j]));
    }
    std::vector<const Patch *> planes, cylinders;
    for (auto i : neighbors) {
        if (known[i].kind == "plane")
            planes.push_back(&known[i]);
        if (known[i].kind == "cylinder")
            cylinders.push_back(&known[i]);
    }
    auto ids = vertices(mesh, triangles);
    std::vector<Patch> out;
    Vec center = mean(mesh, ids);
    for (unsigned a = 0; a < planes.size(); ++a)
        for (unsigned b = 0; b < a; ++b) {
            Vec axis = planes[a]->axis.cross(planes[b]->axis);
            if (axis.norm() < 1e-4)
                continue;
            axis.normalize();
            Eigen::Matrix3d A;
            A.row(0) = planes[a]->axis.transpose();
            A.row(1) = planes[b]->axis.transpose();
            A.row(2) = axis.transpose();
            Eigen::Vector3d rhs(planes[a]->origin.dot(planes[a]->axis),
                                planes[b]->origin.dot(planes[b]->axis), center.dot(axis));
            Vec origin = A.colPivHouseholderQr().solve(rhs);
            for (double s1 : {-1., 1.})
                for (double s2 : {-1., 1.}) {
                    Vec shift = A.colPivHouseholderQr().solve(Vec(s1, s2, 0));
                    std::vector<double> radii;
                    for (auto id : ids) {
                        Vec q = mesh.points[id] - origin;
                        q -= q.dot(axis) * axis;
                        double aa = shift.squaredNorm() - 1, bb = -2 * q.dot(shift),
                               cc = q.squaredNorm(), disc = bb * bb - 4 * aa * cc;
                        if (disc < 0 || std::abs(aa) < 1e-10)
                            continue;
                        for (double sign : {-1., 1.}) {
                            double r = (-bb + sign * std::sqrt(disc)) / (2 * aa);
                            if (r > 1e-5 && std::isfinite(r))
                                radii.push_back(r);
                        }
                    }
                    if (radii.empty())
                        continue;
                    std::sort(radii.begin(), radii.end());
                    Eigen::VectorXd initial(1);
                    initial[0] = radii[radii.size() / 2];
                    auto residual = [&](const Eigen::VectorXd &x) -> Eigen::VectorXd {
                        Eigen::VectorXd r(ids.size());
                        for (unsigned i = 0; i < ids.size(); ++i) {
                            Vec q = mesh.points[ids[i]] - origin - x[0] * shift;
                            q -= q.dot(axis) * axis;
                            r[i] = q.norm() - x[0];
                        }
                        return r;
                    };
                    auto fit = least_squares(residual, initial, 40);
                    if (fit.x[0] > 1e-5) {
                        Patch p;
                        p.kind = "cylinder";
                        p.origin = origin + fit.x[0] * shift;
                        p.axis = axis;
                        p.radius = fit.x[0];
                        p.triangles = triangles;
                        out.push_back(p);
                    }
                }
        }
    for (const auto *cylinder : cylinders)
        for (const auto *plane : planes)
            if (std::abs(cylinder->axis.dot(plane->axis)) > .9999) {
                Vec axis = cylinder->axis,
                    origin =
                        cylinder->origin + ((plane->origin - cylinder->origin).dot(axis)) * axis;
                for (double sr : {-1., 1.})
                    for (double sz : {-1., 1.}) {
                        std::vector<double> radii;
                        for (auto id : ids) {
                            Vec q = mesh.points[id] - origin;
                            double z = q.dot(axis), rho = (q - z * axis).norm() - cylinder->radius,
                                   a = sr * rho + sz * z, disc = a * a - rho * rho - z * z;
                            if (disc < 0)
                                continue;
                            for (double sign : {-1., 1.}) {
                                double r = a + sign * std::sqrt(disc);
                                if (r > 1e-5 && cylinder->radius + sr * r > r)
                                    radii.push_back(r);
                            }
                        }
                        if (radii.empty())
                            continue;
                        std::sort(radii.begin(), radii.end());
                        Eigen::VectorXd initial(1);
                        initial[0] = radii[radii.size() / 2];
                        auto residual = [&](const Eigen::VectorXd &x) -> Eigen::VectorXd {
                            Eigen::VectorXd r(ids.size());
                            for (unsigned i = 0; i < ids.size(); ++i) {
                                Vec q = mesh.points[ids[i]] - origin;
                                double z = q.dot(axis), rho = (q - z * axis).norm();
                                r[i] =
                                    std::hypot(rho - cylinder->radius - sr * x[0], z - sz * x[0]) -
                                    x[0];
                            }
                            return r;
                        };
                        auto fit = least_squares(residual, initial, 40);
                        double r = fit.x[0];
                        if (r > 1e-5 && cylinder->radius + sr * r > r) {
                            Patch p;
                            p.kind = "torus";
                            p.axis = axis;
                            p.origin = origin + sz * r * axis;
                            p.radius = cylinder->radius + sr * r;
                            p.secondary = r;
                            p.triangles = triangles;
                            out.push_back(p);
                        }
                    }
            }
    return out;
}
} // namespace
std::vector<Patch> detect(const Mesh &m, Context &ctx) {
    trace("analytic detection");
    const auto top = topology(m);
    const unsigned n = unsigned(m.triangles.size());
    std::vector<int> owner(n, -1);
    std::vector<unsigned> bodies(n);
    unsigned body = 0;
    for (const auto &group : components(top.neighbors)) {
        for (auto i : group)
            bodies[i] = body;
        ++body;
    }
    std::vector<Patch> patches;
    auto supported_region = [&](const Patch &p) {
        // Local six-parameter fits can interpolate a tiny piece of a general
        // blend. Reserve such pieces for the freeform pass instead of claiming
        // them as independent analytic faces.
        double area = 0;
        for (auto t : p.triangles)
            area += top.areas[t];
        if (p.kind == "sphere")
            return p.triangles.size() >= 48 && area >= p.radius * p.radius;
        if (p.kind == "cone")
            return p.triangles.size() >= 48 && area >= .2 * p.radius * p.radius;
        if (p.kind == "cylinder")
            return p.triangles.size() >= 24 && area >= p.radius * p.radius;
        if (p.kind == "torus")
            return area >= p.secondary * p.secondary;
        return true;
    };
    auto prune_weak_cylinder = [&](Patch &p) {
        if (p.kind != "cylinder" || p.radius <= 0 || p.triangles.empty())
            return;
        std::set<unsigned> weak;
        double total = 0, weak_area = 0;
        double resolvable = std::sqrt(8e-4 / p.radius);
        for (auto t : p.triangles) {
            total += top.areas[t];
            Vec radial[3];
            double z[3];
            auto tri = m.triangles[t];
            bool valid = true;
            for (unsigned k = 0; k < 3; ++k) {
                Vec q = m.points[tri[k]] - p.origin;
                z[k] = q.dot(p.axis);
                radial[k] = q - z[k] * p.axis;
                if (radial[k].norm() < 1e-12)
                    valid = false;
            }
            if (!valid)
                continue;
            double turn = 0, axial_edge = 0;
            for (unsigned k = 0; k < 3; ++k) {
                unsigned j = (k + 1) % 3;
                turn = std::max(turn, (radial[k].normalized() - radial[j].normalized()).norm());
                if ((radial[k] - radial[j]).norm() < 1e-5)
                    axial_edge = std::max(axial_edge, std::abs(z[k] - z[j]));
            }
            if (turn < resolvable && axial_edge < p.radius * .25) {
                weak.insert(t);
                weak_area += top.areas[t];
            }
        }
        // A well-supported curved strip must not claim nearly planar needles
        // from an adjoining blend merely because they fall inside its fit band.
        if (weak_area > 0 && weak_area < .1 * total)
            p.triangles.erase(std::remove_if(p.triangles.begin(), p.triangles.end(),
                                             [&](unsigned t) { return weak.count(t) > 0; }),
                              p.triangles.end());
    };
    auto add = [&](Patch p) {
        surface(p);
        p.body = bodies[p.triangles.front()];
        for (auto t : p.triangles)
            owner[t] = int(patches.size());
        patches.push_back(std::move(p));
    };
    // Flat connected components first. Narrow quads in curved bands stay available
    // to the analytic recognizer rather than becoming one plane per mesh strip.
    std::vector<std::vector<unsigned>> flat(n);
    for (unsigned i = 0; i < n; ++i)
        for (auto j : top.neighbors[i])
            if (top.normals[i].dot(top.normals[j]) > 1 - 1e-7)
                flat[i].push_back(j);
    for (const auto &group : components(flat)) {
        auto p = fit_plane(m, group);
        if (!holds(p, m, top, group, 1e-4))
            continue;
        bool curved_border = false;
        std::set<unsigned> surround;
        double area = 0;
        for (auto i : group) {
            area += top.areas[i];
            for (auto j : top.neighbors[i]) {
                double d = top.normals[i].dot(top.normals[j]);
                if (d < 1 - 1e-7 && d > 0.75) {
                    curved_border = true;
                    surround.insert(j);
                }
            }
        }
        double border_area = 0;
        for (auto i : surround)
            border_area += top.areas[i];
        auto ids = vertices(m, group);
        Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
        for (auto i : ids) {
            Vec q = m.points[i] - p.origin;
            cov += q * q.transpose();
        }
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(cov);
        bool broad = eigen.eigenvalues()[1] > 0.04 * eigen.eigenvalues()[2];
        if (!ctx.options.curves || !curved_border || group.size() > 32 ||
            (broad && group.size() > 4) || area > 4 * border_area)
            add(p);
    }
    if (ctx.options.curves) {
        std::vector<std::vector<unsigned>> smooth(n);
        for (unsigned i = 0; i < n; ++i)
            if (owner[i] < 0)
                for (auto j : top.neighbors[i])
                    if (owner[j] < 0 && top.normals[i].dot(top.normals[j]) > 0.72)
                        smooth[i].push_back(j);
        auto seed_neighborhood = [&](unsigned seed, unsigned limit) {
            // Cross long shared edges first. On an extruded curved strip these
            // stay on the same support; short end edges lead into its blends.
            std::priority_queue<std::pair<double, unsigned>> queue;
            queue.push({0, seed});
            std::set<unsigned> visited;
            std::vector<unsigned> local;
            while (!queue.empty() && local.size() < limit) {
                auto i = queue.top().second;
                queue.pop();
                if (owner[i] >= 0 || !visited.insert(i).second)
                    continue;
                local.push_back(i);
                for (auto j : smooth[i])
                    if (owner[j] < 0 && !visited.count(j)) {
                        unsigned common[2], count = 0;
                        for (auto a : m.triangles[i])
                            for (auto b : m.triangles[j])
                                if (a == b && count < 2)
                                    common[count++] = a;
                        if (count == 2)
                            queue.push(
                                {(m.points[common[0]] - m.points[common[1]]).squaredNorm(), j});
                    }
            }
            return local;
        };
        auto groups = components(smooth);
        std::vector<std::vector<Patch>> proposals(groups.size());
        ctx.diagnostics.parallel_tasks += parallel_for(groups.size(), [&](std::size_t k) {
            auto &g = groups[k];
            if (g.size() < 4 || owner[g.front()] >= 0)
                return;
            proposals[k] = fits(m, top, g, true);
        });
        const double strict = std::max(1e-4, std::min(0.002, ctx.options.accuracy * 0.02));
        const double seed_tolerance = 1e-4;
        for (unsigned k = 0; k < groups.size(); ++k) {
            const auto &g = groups[k];
            if (g.size() < 4 || owner[g.front()] >= 0)
                continue;
            bool accepted = false;
            for (auto p : proposals[k])
                if (p.kind != "plane" && holds(p, m, top, g, strict)) {
                    add(p);
                    accepted = true;
                    break;
                }
            if (accepted)
                continue;
            // Preserve helical regions before local analytic seeds fragment them
            // into short strips. A cheap ring test rejects ordinary fillets.
            std::vector<Vec> thread_axes;
            for (const auto &p : proposals[k])
                if (p.kind == "cylinder") {
                    bool duplicate = false;
                    for (const auto &a : thread_axes)
                        if (std::abs(a.dot(p.axis)) > .9999)
                            duplicate = true;
                    if (!duplicate)
                        thread_axes.push_back(p.axis);
                }
            thread_axes.push_back(Vec::UnitX());
            thread_axes.push_back(Vec::UnitY());
            thread_axes.push_back(Vec::UnitZ());
            if (g.size() >= 128)
                for (const Vec &axis : thread_axes) {
                    auto helical = thread_patch(
                        m, top, g, axis, std::max(1e-5, std::min(.002, ctx.options.accuracy)));
                    if (!helical.surface.IsNull()) {
                        add(helical);
                        accepted = true;
                        break;
                    }
                }
            if (accepted)
                continue;
            // Compare support across the whole smooth island before committing
            // local seeds. A seed on a blend junction otherwise fragments the
            // adjacent cylinder into many convincing but short approximations.
            unsigned stride = std::max<unsigned>(1, unsigned(g.size() / 256));
            std::vector<unsigned> seeds;
            for (unsigned i = 0; i < g.size(); i += stride)
                if (owner[g[i]] < 0)
                    seeds.push_back(g[i]);
            std::vector<std::vector<Patch>> island_proposals(seeds.size());
            ctx.diagnostics.parallel_tasks += parallel_for(seeds.size(), [&](std::size_t index) {
                unsigned seed = seeds[index];
                auto local = seed_neighborhood(seed, 24);
                if (local.size() < 6)
                    return;
                auto candidates = fits(m, top, local, true);
                for (unsigned size : {4u, 8u})
                    if (local.size() > size) {
                        auto small =
                            fits(m, top, std::vector<unsigned>(local.begin(), local.begin() + size),
                                 false);
                        candidates.insert(candidates.end(), small.begin(), small.end());
                    }
                auto rolling = rolling_fillets(m, top, local, g, owner, patches);
                candidates.insert(candidates.end(), rolling.begin(), rolling.end());
                for (auto p : candidates) {
                    if (p.kind == "plane")
                        continue;
                    auto on = [&](unsigned t) {
                        for (auto v : m.triangles[t])
                            if (deviation(p, m.points[v]) > seed_tolerance)
                                return false;
                        return normal_agrees(p, top, t);
                    };
                    if (!on(seed))
                        continue;
                    std::vector<unsigned> grown{seed};
                    std::set<unsigned> visited{seed};
                    for (unsigned i = 0; i < grown.size(); ++i)
                        for (auto j : smooth[grown[i]])
                            if (owner[j] < 0 && visited.insert(j).second && on(j))
                                grown.push_back(j);
                    if (grown.size() < 24)
                        continue;
                    p.triangles = std::move(grown);
                    prune_weak_cylinder(p);
                    if (p.triangles.size() >= 24 && holds(p, m, top, p.triangles, seed_tolerance) &&
                        supported_region(p))
                        island_proposals[index].push_back(std::move(p));
                }
            });
            std::vector<Patch> supported;
            for (auto &proposals : island_proposals)
                for (auto &p : proposals)
                    supported.push_back(std::move(p));
            std::stable_sort(supported.begin(), supported.end(),
                             [](const Patch &a, const Patch &b) {
                                 if ((a.deviation < 1e-4) != (b.deviation < 1e-4))
                                     return a.deviation < 1e-4;
                                 return a.triangles.size() != b.triangles.size()
                                            ? a.triangles.size() > b.triangles.size()
                                            : a.deviation < b.deviation;
                             });
            for (auto p : supported) {
                auto old_size = p.triangles.size();
                p.triangles.erase(std::remove_if(p.triangles.begin(), p.triangles.end(),
                                                 [&](unsigned t) { return owner[t] >= 0; }),
                                  p.triangles.end());
                if (p.triangles.size() < 24 || p.triangles.size() * 5 < old_size * 4)
                    continue;
                if (holds(p, m, top, p.triangles, strict))
                    add(std::move(p));
            }
            // Local seed/grow passes for smooth islands spanning several surfaces.
            for (unsigned seed : g) {
                if (owner[seed] >= 0)
                    continue;
                auto neighborhood = seed_neighborhood(seed, 24);
                if (neighborhood.size() < 4)
                    continue;
                double turning = 0;
                for (auto i : neighborhood)
                    turning = std::max(turning, (top.normals[i] - top.normals[seed]).norm());
                if (turning < .03)
                    continue;
                Patch best;
                std::size_t best_size = 0;
                auto candidates = fits(m, top, neighborhood, false);
                for (unsigned size : {4u, 8u})
                    if (neighborhood.size() > size) {
                        auto small = fits(m, top,
                                          std::vector<unsigned>(neighborhood.begin(),
                                                                neighborhood.begin() + size),
                                          false);
                        candidates.insert(candidates.end(), small.begin(), small.end());
                    }
                auto fillets = rolling_fillets(m, top, neighborhood, g, owner, patches);
                candidates.insert(candidates.end(), fillets.begin(), fillets.end());
                for (auto p : candidates) {
                    if (p.kind == "plane")
                        continue;
                    std::vector<unsigned> grown{seed};
                    std::set<unsigned> visited{seed};
                    auto on = [&](unsigned t) {
                        for (auto v : m.triangles[t])
                            if (deviation(p, m.points[v]) > seed_tolerance)
                                return false;
                        return normal_agrees(p, top, t);
                    };
                    if (!on(seed))
                        continue;
                    for (std::size_t q = 0; q < grown.size(); ++q)
                        for (auto j : smooth[grown[q]])
                            if (owner[j] < 0 && visited.insert(j).second && on(j))
                                grown.push_back(j);
                    if (grown.size() >= 6 && grown.size() > best_size) {
                        p.triangles = std::move(grown);
                        prune_weak_cylinder(p);
                        if (p.triangles.size() > best_size &&
                            holds(p, m, top, p.triangles, seed_tolerance) && supported_region(p) &&
                            (p.triangles.size() >= 64 || p.deviation < 1e-4)) {
                            best_size = p.triangles.size();
                            best = std::move(p);
                        }
                    }
                }
                if (best_size) {
                    for (unsigned iteration = 0; iteration < 2 && best.kind != "torus";
                         ++iteration) {
                        auto refined = fits(m, top, best.triangles, false);
                        Patch chosen = best;
                        for (auto q : refined)
                            if (q.kind == best.kind &&
                                holds(q, m, top, best.triangles, seed_tolerance)) {
                                std::vector<unsigned> grown{seed};
                                std::set<unsigned> visited{seed};
                                auto on = [&](unsigned t) {
                                    for (auto v : m.triangles[t])
                                        if (deviation(q, m.points[v]) > seed_tolerance)
                                            return false;
                                    return normal_agrees(q, top, t);
                                };
                                for (std::size_t j = 0; j < grown.size(); ++j)
                                    for (auto next : smooth[grown[j]])
                                        if (owner[next] < 0 && visited.insert(next).second &&
                                            on(next))
                                            grown.push_back(next);
                                q.triangles = std::move(grown);
                                prune_weak_cylinder(q);
                                if (q.triangles.size() >= chosen.triangles.size() &&
                                    holds(q, m, top, q.triangles, seed_tolerance) &&
                                    (q.triangles.size() > chosen.triangles.size() ||
                                     q.deviation < chosen.deviation))
                                    chosen = std::move(q);
                            }
                        best = std::move(chosen);
                    }
                    add(best);
                }
            }
        }
    }
    // A verified mesh symmetry can propose an independently checked patch on the
    // opposite side. It never bypasses the point/normal acceptance tests.
    if (ctx.options.curves) {
        auto mirrors = symmetry_planes(m);
        auto known = patches;
        for (const auto &mirror : mirrors)
            for (auto p : known)
                if (p.kind != "plane") {
                    const auto &n = mirror.first;
                    p.origin -= 2 * (n.dot(p.origin) - mirror.second) * n;
                    p.axis -= 2 * n.dot(p.axis) * n;
                    std::vector<bool> on(owner.size(), false);
                    for (unsigned i = 0; i < owner.size(); ++i)
                        if (owner[i] < 0) {
                            Patch one = p;
                            if (holds(one, m, top, {i},
                                      std::max(1e-4, std::min({0.002, ctx.options.accuracy * 0.02,
                                                               2 * p.deviation + 1e-5}))))
                                on[i] = true;
                        }
                    std::vector<std::vector<unsigned>> graph(owner.size());
                    for (unsigned i = 0; i < on.size(); ++i)
                        if (on[i])
                            for (auto j : top.neighbors[i])
                                if (on[j])
                                    graph[i].push_back(j);
                    for (auto ids : components(graph))
                        if (ids.size() >= 6 && on[ids.front()]) {
                            auto candidate = p;
                            candidate.triangles = std::move(ids);
                            add(candidate);
                        }
                }
    }
    // Remove only sub-weld noise before these supports constrain blend rims.
    for (auto &p : patches)
        if (p.kind == "cylinder" || p.kind == "sphere" || p.kind == "torus") {
            auto candidate = p;
            Eigen::Index k;
            candidate.axis.cwiseAbs().maxCoeff(&k);
            if (std::abs(candidate.axis[k]) > .999999) {
                Vec axis = Vec::Zero();
                axis[k] = candidate.axis[k] > 0 ? 1 : -1;
                candidate.axis = axis;
            }
            Vec shift;
            for (unsigned j = 0; j < 3; ++j)
                shift[j] = std::round(p.origin[j] * 100) / 100 - p.origin[j];
            if (p.kind == "cylinder")
                shift -= shift.dot(candidate.axis) * candidate.axis;
            double radius = std::round(p.radius * 100) / 100;
            if (shift.norm() > 1e-4 || radius <= 0 || std::abs(radius - p.radius) > 1e-4)
                continue;
            candidate.origin += shift;
            candidate.radius = radius;
            if (p.kind == "torus") {
                double r = std::round(p.secondary * 100) / 100;
                if (r > 0 && std::abs(r - p.secondary) < 1e-4)
                    candidate.secondary = r;
            }
            if (holds(candidate, m, top, candidate.triangles, 2e-4, false)) {
                surface(candidate);
                p = std::move(candidate);
            }
        }
    // Recover true trimmed cylinder facets after the needle test, using a much
    // tighter point test than the initial search band.
    for (unsigned pi = 0; pi < patches.size(); ++pi) {
        auto &p = patches[pi];
        if (p.kind != "cylinder")
            continue;
        for (unsigned i = 0; i < p.triangles.size(); ++i)
            for (auto t : top.neighbors[p.triangles[i]])
                if (owner[t] < 0) {
                    bool on = normal_agrees(p, top, t);
                    for (auto v : m.triangles[t])
                        if (deviation(p, m.points[v]) > 1e-5)
                            on = false;
                    if (on) {
                        owner[t] = pi;
                        p.triangles.push_back(t);
                    }
                }
    }
    // Optional reconstruction starts after complete analytic analysis.
    trace("advanced fitting");
    ctx.deadline = ctx.options.time_limit
                       ? Clock::now() + std::chrono::duration_cast<Clock::duration>(
                                            std::chrono::duration<double>(ctx.options.time_limit))
                       : Clock::time_point::max();
    auto virtual_owner = owner;
    std::vector<bool> eligible(n, false);
    for (unsigned i = 0; i < n; ++i)
        eligible[i] = owner[i] < 0;
    for (const auto &p : patches) {
        double area = 0;
        for (auto i : p.triangles)
            area += top.areas[i];
        bool small = false;
        if (p.kind == "plane")
            small =
                p.triangles.size() <= 2 && area < ctx.options.accuracy * ctx.options.accuracy * 10;
        else if (p.kind == "sphere")
            small = p.triangles.size() <= 128 && area < 4 * p.radius * p.radius;
        else if (p.kind == "cylinder" || p.kind == "cone")
            small = p.triangles.size() <= 64 && area < 3 * p.radius * p.radius;
        else if (p.kind == "torus")
            small = p.triangles.size() <= 64 && area < 4 * p.secondary * p.secondary;
        if (small)
            for (auto i : p.triangles)
                eligible[i] = true;
    }
    std::vector<std::vector<unsigned>> transition_graph(n);
    for (unsigned i = 0; i < n; ++i)
        if (eligible[i])
            for (auto j : top.neighbors[i])
                if (eligible[j] && top.normals[i].dot(top.normals[j]) > .72)
                    transition_graph[i].push_back(j);
    std::vector<Patch> rescued;
    for (const auto &group : components(transition_graph)) {
        if (group.size() < 6 || !eligible[group.front()])
            continue;
        std::vector<double> adjoining;
        for (auto i : group)
            for (auto j : top.neighbors[i])
                if (!eligible[j] && owner[j] >= 0 && top.normals[i].dot(top.normals[j]) > .95) {
                    const auto &p = patches[owner[j]];
                    if (p.kind == "cylinder" || p.kind == "sphere")
                        adjoining.push_back(p.radius);
                    else if (p.kind == "torus")
                        adjoining.push_back(p.secondary);
                }
        if (adjoining.size() < 2)
            continue;
        auto limits = std::minmax_element(adjoining.begin(), adjoining.end());
        if (*limits.second - *limits.first < std::max(.01, ctx.options.accuracy * .1))
            continue;
        bool exact = false;
        for (auto p : fits(m, top, group, true))
            if (p.kind != "plane" &&
                holds(p, m, top, group,
                      std::max(1e-4, std::min(.002, ctx.options.accuracy * .02)))) {
                rescued.push_back(p);
                exact = true;
                break;
            }
        for (auto i : group)
            virtual_owner[i] = exact ? int(patches.size()) : -1;
    }
    auto advanced = advanced_patches(m, top, virtual_owner, patches, ctx);
    advanced.insert(advanced.end(), rescued.begin(), rescued.end());
    if (!advanced.empty()) {
        std::vector<bool> replaced(n, false);
        for (const auto &p : advanced)
            for (auto i : p.triangles)
                replaced[i] = true;
        for (auto &p : patches)
            p.triangles.erase(std::remove_if(p.triangles.begin(), p.triangles.end(),
                                             [&](unsigned i) { return replaced[i]; }),
                              p.triangles.end());
        patches.erase(std::remove_if(patches.begin(), patches.end(),
                                     [](const Patch &p) { return p.triangles.empty(); }),
                      patches.end());
        std::fill(owner.begin(), owner.end(), -1);
        for (unsigned i = 0; i < patches.size(); ++i)
            for (auto j : patches[i].triangles)
                owner[j] = int(i);
        for (auto p : advanced)
            add(std::move(p));
    }
    // Unrecognized facets remain exact planar geometry. They are never discarded.
    for (const auto &group : components(flat)) {
        std::vector<unsigned> left;
        for (auto i : group)
            if (owner[i] < 0)
                left.push_back(i);
        if (left.empty())
            continue;
        auto p = fit_plane(m, left);
        if (holds(p, m, top, left, 1e-5))
            add(p);
        else
            for (auto i : left)
                add(fit_plane(m, {i}));
    }
    // A very slender tessellation facet can have an unreliable normal even
    // though its corners lie on an adjoining curved surface. Do not retain it
    // as an independent planar face if a verified neighbor explains its points
    // and preserves the orientation of the projected triangle.
    unsigned absorbed = 0;
    for (unsigned pi = 0; pi < patches.size(); ++pi) {
        auto &p = patches[pi];
        if (p.kind != "plane" || p.triangles.empty() || p.triangles.size() > 2)
            continue;
        double area = 0, longest = 0;
        bool slender = true;
        for (auto t : p.triangles) {
            auto tri = m.triangles[t];
            double length2 = 0;
            for (unsigned k = 0; k < 3; ++k)
                length2 = std::max(length2,
                                   (m.points[tri[k]] - m.points[tri[(k + 1) % 3]]).squaredNorm());
            area += top.areas[t];
            longest = std::max(longest, std::sqrt(length2));
            if (length2 <= 0 || 2 * top.areas[t] > .05 * length2)
                slender = false;
        }
        if (!slender || area > 10 * ctx.options.accuracy * ctx.options.accuracy ||
            longest > std::max(4 * ctx.options.accuracy, 10 * ctx.options.tolerance))
            continue;
        std::set<unsigned> neighbors;
        for (auto t : p.triangles)
            for (auto j : top.neighbors[t])
                if (owner[j] >= 0 && unsigned(owner[j]) != pi && patches[owner[j]].kind != "plane")
                    neighbors.insert(unsigned(owner[j]));
        unsigned best = unsigned(patches.size());
        double best_error = 1e100;
        const double limit = std::max(1e-5, std::min(.002, ctx.options.accuracy * .02));
        for (auto candidate : neighbors) {
            const auto &support = patches[candidate];
            PointProjector project(support.surface);
            std::map<unsigned, Vec> projected;
            double error = 0;
            bool fits = true;
            for (auto v : vertices(m, p.triangles)) {
                const auto &q = project.perform(m.points[v]);
                if (!q.NbPoints()) {
                    fits = false;
                    break;
                }
                error = std::max(error, q.LowerDistance());
                auto at = q.NearestPoint();
                projected.emplace(v, Vec(at.X(), at.Y(), at.Z()));
            }
            if (!fits || error > limit || error >= best_error)
                continue;
            for (auto t : p.triangles) {
                auto tri = m.triangles[t];
                Vec cross = (projected.at(tri[1]) - projected.at(tri[0]))
                                .cross(projected.at(tri[2]) - projected.at(tri[0]));
                if (cross.dot(top.normals[t]) <= std::max(1e-14, 2 * top.areas[t] * 1e-4)) {
                    fits = false;
                    break;
                }
                const auto &q = project.perform(top.centers[t]);
                if (!q.NbPoints()) {
                    fits = false;
                    break;
                }
                double u, v;
                q.LowerDistanceParameters(u, v);
                gp_Pnt at;
                gp_Vec du, dv;
                support.surface->D1(u, v, at, du, dv);
                auto normal = du.Crossed(dv);
                Vec n(normal.X(), normal.Y(), normal.Z());
                if (n.norm() < 1e-12 || std::abs(n.normalized().dot(top.normals[t])) > .72) {
                    fits = false;
                    break;
                }
            }
            if (fits) {
                best = candidate;
                best_error = error;
            }
        }
        if (best < patches.size()) {
            auto &target = patches[best];
            target.triangles.insert(target.triangles.end(), p.triangles.begin(), p.triangles.end());
            target.deviation = std::max(target.deviation, best_error);
            for (auto t : p.triangles)
                owner[t] = int(best);
            absorbed += unsigned(p.triangles.size());
            p.triangles.clear();
        }
    }
    if (absorbed) {
        patches.erase(std::remove_if(patches.begin(), patches.end(),
                                     [](const Patch &p) { return p.triangles.empty(); }),
                      patches.end());
        std::fill(owner.begin(), owner.end(), -1);
        for (unsigned i = 0; i < patches.size(); ++i)
            for (auto t : patches[i].triangles)
                owner[t] = int(i);
        ctx.diagnostics.repairs.push_back("Reconstructed " + std::to_string(absorbed) +
                                          " slender transition facets on adjoining surfaces");
    }
    // Snap numerical noise in repeated radii and almost axial frames, retaining
    // each change only if the original mesh points still satisfy the model.
    std::vector<double> radii;
    for (const auto &p : patches)
        if (p.radius > 0 && (p.kind == "cylinder" || p.kind == "sphere" || p.kind == "cone"))
            radii.push_back(p.radius);
    std::sort(radii.begin(), radii.end());
    for (auto &p : patches)
        if (p.kind == "cylinder" || p.kind == "cone" || p.kind == "sphere" || p.kind == "torus") {
            auto candidate = p;
            double closest = p.radius;
            std::vector<double> equal;
            for (auto r : radii)
                if (std::abs(r - p.radius) < 1e-4)
                    equal.push_back(r);
            if (!equal.empty())
                closest = equal[equal.size() / 2];
            double nominal = std::round(closest * 10) / 10;
            if (nominal > 0 && std::abs(nominal - closest) < 1e-5)
                closest = nominal;
            candidate.radius = closest;
            if (p.kind != "sphere") {
                Eigen::Index axis;
                candidate.axis.cwiseAbs().maxCoeff(&axis);
                if (std::abs(candidate.axis[axis]) > 0.999999) {
                    Vec axial = Vec::Zero();
                    axial[axis] = candidate.axis[axis] > 0 ? 1 : -1;
                    candidate.axis = axial;
                }
            }
            if (holds(candidate, m, top, candidate.triangles,
                      std::max(1e-4, std::min(.002, ctx.options.accuracy * .02)))) {
                surface(candidate);
                p = std::move(candidate);
            }
        }
    // Remove sub-quantization noise from support positions. Small transverse
    // center errors become large parameter-angle errors near a sphere's poles.
    for (auto &p : patches)
        if (p.kind == "plane" || p.kind == "sphere" || p.kind == "cylinder" || p.kind == "torus") {
            Vec shift;
            for (unsigned k = 0; k < 3; ++k)
                shift[k] = std::round(p.origin[k] * 100) / 100 - p.origin[k];
            if (p.kind == "cylinder")
                shift -= shift.dot(p.axis) * p.axis;
            if (p.kind == "plane")
                shift = shift.dot(p.axis) * p.axis;
            if (shift.norm() > .002)
                continue;
            auto candidate = p;
            candidate.origin += shift;
            if (p.kind != "plane") {
                double radius = std::round(p.radius * 100) / 100;
                if (radius > 0 && std::abs(radius - p.radius) < .002)
                    candidate.radius = radius;
            }
            if (p.kind == "torus") {
                double radius = std::round(p.secondary * 100) / 100;
                if (radius > 0 && std::abs(radius - p.secondary) < .002)
                    candidate.secondary = radius;
            }
            if (holds(candidate, m, top, candidate.triangles, 1e-4, false)) {
                surface(candidate);
                p = std::move(candidate);
            }
        }
    for (auto &p : patches)
        if (p.kind == "torus") {
            auto candidate = p;
            double radius = std::round(p.secondary * 100) / 100;
            if (radius > 0 && std::abs(radius - p.secondary) < 1e-4) {
                candidate.secondary = radius;
                if (holds(candidate, m, top, candidate.triangles, 1e-4)) {
                    surface(candidate);
                    p = std::move(candidate);
                }
            }
        }
    // The spine circle of a torus constrains adjoining cylindrical fillets.
    for (auto &p : patches)
        if (p.kind == "cylinder")
            for (const auto &q : patches)
                if (q.kind == "torus" && q.deviation < 1e-4 &&
                    std::abs(p.axis.dot(q.axis)) < 1e-4 && std::abs(p.radius - q.secondary) < .01) {
                    Vec shift = p.origin - q.origin;
                    Vec radial = shift - shift.dot(p.axis) * p.axis;
                    radial -= radial.dot(q.axis) * q.axis;
                    if (radial.norm() < 1e-6)
                        continue;
                    radial.normalize();
                    auto candidate = p;
                    candidate.origin = q.origin + q.radius * radial + shift.dot(p.axis) * p.axis;
                    candidate.radius = q.secondary;
                    if ((candidate.origin - p.origin).norm() < .01 &&
                        holds(candidate, m, top, candidate.triangles,
                              std::max(1e-4, std::min(.002, ctx.options.accuracy * .02))) &&
                        candidate.deviation < p.deviation) {
                        surface(candidate);
                        p = std::move(candidate);
                    }
                }
    // A spherical corner gives an exact center on each adjoining fillet axis.
    // Retain the constraint only when every sampled vertex still fits.
    for (auto &p : patches)
        if (p.kind == "cylinder")
            for (const auto &q : patches)
                if (q.kind == "sphere" && q.deviation < 1e-4 &&
                    std::abs(q.radius - p.radius) < .01) {
                    Vec shift = q.origin - p.origin;
                    shift -= shift.dot(p.axis) * p.axis;
                    if (shift.norm() > .01)
                        continue;
                    auto candidate = p;
                    candidate.origin += shift;
                    candidate.radius = q.radius;
                    if (holds(candidate, m, top, candidate.triangles,
                              std::max(1e-4, std::min(.002, ctx.options.accuracy * .02)))) {
                        surface(candidate);
                        p = std::move(candidate);
                        break;
                    }
                }
    std::vector<double> clean_radii;
    for (const auto &p : patches)
        if ((p.kind == "cylinder" || p.kind == "sphere") && p.deviation < 1e-4)
            clean_radii.push_back(p.radius);
    for (auto &p : patches)
        if (p.kind == "torus") {
            Patch best = p;
            for (double radius : clean_radii)
                if (std::abs(radius - p.secondary) < std::min(.01, ctx.options.accuracy * .05))
                    for (double sr : {-1., 1.})
                        for (double sz : {-1., 1.}) {
                            Patch candidate = p;
                            double delta = radius - p.secondary;
                            candidate.secondary = radius;
                            candidate.radius += sr * delta;
                            candidate.origin += sz * delta * p.axis;
                            if (candidate.radius > candidate.secondary &&
                                holds(candidate, m, top, candidate.triangles,
                                      std::max(1e-4, std::min(.002, ctx.options.accuracy * .02))) &&
                                candidate.deviation < best.deviation)
                                best = std::move(candidate);
                        }
            surface(best);
            p = std::move(best);
        }
    std::vector<Patch> merged;
    for (auto &p : patches) {
        bool joined = false;
        if (p.kind == "plane" || p.kind == "cylinder" || p.kind == "cone" || p.kind == "sphere" ||
            p.kind == "torus")
            for (auto &q : merged)
                if (q.kind == p.kind && q.body == p.body) {
                    bool same = false;
                    Vec shift = p.origin - q.origin;
                    if (p.kind == "plane")
                        same = std::abs(p.axis.dot(q.axis)) > 0.9999999 &&
                               std::abs(shift.dot(q.axis)) < 1e-4;
                    else if (p.kind == "sphere")
                        same = shift.norm() < .02 && std::abs(p.radius - q.radius) < .02;
                    else if (p.kind == "cylinder")
                        same = std::abs(p.axis.dot(q.axis)) > 0.9999 &&
                               (shift - shift.dot(q.axis) * q.axis).norm() < .02 &&
                               std::abs(p.radius - q.radius) < .02;
                    else if (p.kind == "torus")
                        same = std::abs(p.axis.dot(q.axis)) > 0.9999 && shift.norm() < .02 &&
                               std::abs(p.radius - q.radius) < .02 &&
                               std::abs(p.secondary - q.secondary) < .02;
                    else
                        same =
                            p.axis.dot(q.axis) > 0.9999 &&
                            (shift - shift.dot(q.axis) * q.axis).norm() < .02 &&
                            std::abs(p.secondary - q.secondary) < .001 &&
                            std::abs(p.radius - q.radius - q.secondary * shift.dot(q.axis)) < .02;
                    if (same) {
                        auto candidate = q;
                        candidate.triangles.insert(candidate.triangles.end(), p.triangles.begin(),
                                                   p.triangles.end());
                        double limit =
                            p.kind == "plane"
                                ? 1e-4
                                : std::max(1e-4, std::min(.002, ctx.options.accuracy * .02));
                        if (!holds(candidate, m, top, candidate.triangles, limit)) {
                            auto refits =
                                p.kind == "plane"
                                    ? std::vector<Patch>{fit_plane(m, candidate.triangles)}
                                    : fits(m, top, candidate.triangles, p.kind == "torus");
                            for (auto refit : refits)
                                if (refit.kind == p.kind &&
                                    holds(refit, m, top, refit.triangles, limit)) {
                                    candidate = std::move(refit);
                                    break;
                                }
                        }
                        if (holds(candidate, m, top, candidate.triangles, limit)) {
                            candidate.body = q.body;
                            surface(candidate);
                            q = std::move(candidate);
                            joined = true;
                            break;
                        }
                    }
                }
        if (!joined)
            merged.push_back(std::move(p));
    }
    // Reclaim adjacent loose facets and short stand-in surfaces once a larger
    // exact support is known. Signed normal agreement preserves opposing caps.
    // A short conical band can have only two vertex rings and is then also
    // cospherical. Refit the complete band after merging the local seeds.
    for (auto &p : merged)
        if (p.kind == "sphere") {
            double sphere_error = 0;
            for (auto t : p.triangles)
                sphere_error += 1 - std::abs(normal(p, top.centers[t]).dot(top.normals[t]));
            auto candidates = fits(m, top, p.triangles, true);
            auto ids = vertices(m, p.triangles);
            // For two cospherical rings their sphere center lies on the cone axis.
            // Use that center to avoid the biased circle fit of unequal-radius rings.
            Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
            Vec average = Vec::Zero();
            for (auto t : p.triangles)
                average += top.normals[t];
            average /= double(p.triangles.size());
            for (auto t : p.triangles) {
                Vec d = top.normals[t] - average;
                covariance += d * d.transpose();
            }
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(covariance);
            for (const Vec &axis : std::vector<Vec>{eigen.eigenvectors().col(0), Vec::UnitX(),
                                                    Vec::UnitY(), Vec::UnitZ()}) {
                Eigen::MatrixXd z(ids.size(), 2);
                Eigen::VectorXd r(ids.size());
                for (unsigned i = 0; i < ids.size(); ++i) {
                    Vec d = m.points[ids[i]] - p.origin;
                    double h = d.dot(axis);
                    z.row(i) << 1, h;
                    r[i] = (d - h * axis).norm();
                }
                Eigen::Vector2d profile = z.colPivHouseholderQr().solve(r);
                if (profile[0] > 0 && std::abs(profile[1]) > 1e-5 && std::abs(profile[1]) < 100) {
                    Patch candidate = p;
                    candidate.kind = "cone";
                    candidate.axis = axis;
                    candidate.radius = profile[0];
                    candidate.secondary = profile[1];
                    candidates.push_back(candidate);
                }
            }
            for (auto candidate : candidates)
                if (candidate.kind == "cone" || candidate.kind == "cylinder") {
                    bool accepted = holds(candidate, m, top, p.triangles, 1e-4);
                    if (!accepted)
                        continue;
                    double error = 0;
                    for (auto t : p.triangles)
                        error +=
                            1 - std::abs(normal(candidate, top.centers[t]).dot(top.normals[t]));
                    if (error < sphere_error * .8) {
                        candidate.body = p.body;
                        candidate.part = p.part;
                        surface(candidate);
                        p = std::move(candidate);
                        break;
                    }
                }
        }
    std::vector<unsigned> assigned(n);
    for (unsigned i = 0; i < merged.size(); ++i)
        for (auto t : merged[i].triangles)
            assigned[t] = i;
    std::vector<bool> exact;
    for (const auto &p : merged)
        exact.push_back(p.deviation < 1e-4);
    for (unsigned pass = 0; pass < 2; ++pass)
        for (const std::string kind : {"torus", "sphere", "cylinder", "cone"})
            for (unsigned qi = 0; qi < merged.size(); ++qi) {
                auto &q = merged[qi];
                if (q.kind != kind || q.triangles.empty())
                    continue;
                double sign = normal(q, top.centers[q.triangles.front()])
                                          .dot(top.normals[q.triangles.front()]) >= 0
                                  ? 1
                                  : -1;
                for (unsigned pi = 0; pi < merged.size(); ++pi) {
                    auto &p = merged[pi];
                    if (pi == qi || p.body != q.body || p.triangles.empty())
                        continue;
                    if (p.kind != "plane" && p.kind != "cylinder" && p.kind != "cone" &&
                        p.kind != "sphere" && p.kind != "torus" && p.kind != "freeform" &&
                        p.kind != "blend" && p.kind != "extrusion")
                        continue;
                    bool carve = exact[qi] && (!exact[pi] || p.kind == "freeform" ||
                                               p.kind == "blend" || p.kind == "extrusion");
                    if (!carve && p.triangles.size() > q.triangles.size())
                        continue;
                    bool adjacent = false;
                    for (auto t : p.triangles)
                        for (auto j : top.neighbors[t])
                            if (assigned[j] == qi)
                                adjacent = true;
                    if (!adjacent)
                        continue;
                    if (carve) {
                        double limit = std::max(1e-4, std::min(.002, ctx.options.accuracy * .02));
                        // Preserve the precision of an exact cylindrical seed.
                        // A looser growth band steals almost tangent blend
                        // facets and creates oscillating shared outlines.
                        if (q.kind == "cylinder")
                            limit = std::min(limit, 1e-4);
                        bool changed = true;
                        std::vector<unsigned> remaining = p.triangles;
                        while (changed) {
                            changed = false;
                            std::vector<unsigned> keep;
                            for (auto t : remaining) {
                                bool touch = false;
                                for (auto j : top.neighbors[t])
                                    if (assigned[j] == qi)
                                        touch = true;
                                bool fits =
                                    touch &&
                                    sign * normal(q, top.centers[t]).dot(top.normals[t]) > .96;
                                auto tri = m.triangles[t];
                                for (unsigned k = 0; k < 3 && fits; ++k)
                                    fits = deviation(q, m.points[tri[k]]) <= limit &&
                                           deviation(
                                               q, (m.points[tri[k]] + m.points[tri[(k + 1) % 3]]) /
                                                      2) <= std::min(.02, ctx.options.accuracy);
                                if (fits) {
                                    assigned[t] = qi;
                                    q.triangles.push_back(t);
                                    changed = true;
                                } else
                                    keep.push_back(t);
                            }
                            remaining = std::move(keep);
                        }
                        p.triangles = std::move(remaining);
                        holds(q, m, top, q.triangles, limit);
                        if (!p.triangles.empty()) {
                            auto plane = fit_plane(m, p.triangles);
                            if (holds(plane, m, top, plane.triangles, 1e-4, false)) {
                                plane.body = p.body;
                                plane.part = p.part;
                                surface(plane);
                                p = std::move(plane);
                            }
                        }
                        continue;
                    }
                    auto candidate = q;
                    double limit = std::max(1e-4, std::min(.002, ctx.options.accuracy * .02));
                    if (!holds(candidate, m, top, p.triangles, limit, false))
                        continue;
                    bool aligned = true;
                    for (auto t : p.triangles) {
                        if (sign * normal(q, top.centers[t]).dot(top.normals[t]) < .96)
                            aligned = false;
                        if (p.kind == "plane") {
                            auto tri = m.triangles[t];
                            for (unsigned k = 0; k < 3; ++k)
                                if (deviation(q, (m.points[tri[k]] + m.points[tri[(k + 1) % 3]]) /
                                                     2) > std::min(.02, ctx.options.accuracy))
                                    aligned = false;
                        }
                    }
                    if (!aligned)
                        continue;
                    q.triangles.insert(q.triangles.end(), p.triangles.begin(), p.triangles.end());
                    for (auto t : p.triangles)
                        assigned[t] = qi;
                    p.triangles.clear();
                    holds(q, m, top, q.triangles, limit);
                }
            }
    // Tangent surfaces can both explain a border facet within the fitting
    // tolerance. Prefer the support with decisively smaller vertex residuals.
    // This removes one-facet protrusions without shifting any support surface.
    auto analytic = [](const Patch &p) {
        return p.kind == "cylinder" || p.kind == "cone" || p.kind == "sphere" || p.kind == "torus";
    };
    for (unsigned pass = 0; pass < 2; ++pass) {
        std::vector<std::pair<unsigned, unsigned>> moves;
        for (unsigned t = 0; t < n; ++t) {
            unsigned pi = assigned[t];
            if (!analytic(merged[pi]))
                continue;
            double current = 0;
            for (auto v : m.triangles[t])
                current = std::max(current, deviation(merged[pi], m.points[v]));
            if (current <= 1e-5)
                continue;
            double best = current;
            unsigned target = pi;
            for (auto j : top.neighbors[t]) {
                unsigned qi = assigned[j];
                if (qi == pi || !analytic(merged[qi]))
                    continue;
                double error = 0;
                for (auto v : m.triangles[t])
                    error = std::max(error, deviation(merged[qi], m.points[v]));
                if (error < best && error < current * .25 && error < 1e-4 &&
                    std::abs(normal(merged[qi], top.centers[t]).dot(top.normals[t])) > .96) {
                    target = qi;
                    best = error;
                }
            }
            if (target != pi)
                moves.emplace_back(t, target);
        }
        if (moves.empty())
            break;
        for (auto move : moves)
            assigned[move.first] = move.second;
        for (unsigned i = 0; i < merged.size(); ++i) {
            auto &p = merged[i];
            p.triangles.erase(std::remove_if(p.triangles.begin(), p.triangles.end(),
                                             [&](unsigned t) { return assigned[t] != i; }),
                              p.triangles.end());
        }
        for (auto move : moves)
            merged[move.second].triangles.push_back(move.first);
    }
    merged.erase(std::remove_if(merged.begin(), merged.end(),
                                [](const Patch &p) { return p.triangles.empty(); }),
                 merged.end());
    for (unsigned i = 0; i < merged.size(); ++i)
        if (merged[i].kind == "plane" || merged[i].kind == "cone")
            for (unsigned j = i + 1; j < merged.size();) {
                auto &p = merged[i];
                const auto &q = merged[j];
                Vec shift = p.origin - q.origin;
                bool same = false;
                if (p.kind == q.kind && p.body == q.body) {
                    if (p.kind == "plane")
                        same = std::abs(p.axis.dot(q.axis)) > .9999999 &&
                               std::abs(shift.dot(q.axis)) < 1e-4;
                    else
                        same =
                            p.axis.dot(q.axis) > .9999999 &&
                            (shift - shift.dot(q.axis) * q.axis).norm() < 1e-4 &&
                            std::abs(p.secondary - q.secondary) < 1e-4 &&
                            std::abs(p.radius - q.radius - q.secondary * shift.dot(q.axis)) < 1e-4;
                }
                if (same) {
                    auto ids = p.triangles;
                    ids.insert(ids.end(), q.triangles.begin(), q.triangles.end());
                    auto candidate = p.kind == "plane" ? fit_plane(m, ids) : p;
                    candidate.triangles = ids;
                    if (holds(candidate, m, top, ids, 1e-4, false)) {
                        candidate.body = p.body;
                        candidate.part = p.part;
                        surface(candidate);
                        p = std::move(candidate);
                        merged.erase(merged.begin() + j);
                        continue;
                    }
                }
                ++j;
            }
    trace("detection complete");
    return merged;
}
void design_size(Mesh &mesh, std::vector<Patch> &patches, Context &ctx) {
    // Port of sizing.py: weighted distinct measurements, explicit slicer/unit
    // hypotheses, and a clear-margin requirement against the unscaled baseline.
    struct Measurement {
        double value, weight;
    };
    std::vector<Measurement> values;
    std::vector<double> areas(patches.size());
    std::vector<Vec> normals(patches.size());
    std::map<std::pair<unsigned, unsigned>, std::vector<unsigned>> bodies;
    auto full_turn = [&](const Patch &p) {
        struct Rim {
            unsigned count = 0, a = 0, b = 0;
        };
        std::map<Edge, Rim> edges;
        for (auto i : p.triangles) {
            auto t = mesh.triangles[i];
            for (unsigned k = 0; k < 3; ++k) {
                auto &edge = edges[Edge(t[k], t[(k + 1) % 3])];
                ++edge.count;
                edge.a = t[k];
                edge.b = t[(k + 1) % 3];
            }
        }
        std::map<unsigned, unsigned> next;
        for (const auto &item : edges)
            if (item.second.count == 1)
                if (!next.emplace(item.second.a, item.second.b).second)
                    return false;
        if (next.empty())
            return p.kind != "sphere";
        Vec e = p.axis.unitOrthogonal(), f = p.axis.cross(e), first_normal = Vec::Zero();
        bool planar = true, winding = false;
        while (!next.empty()) {
            unsigned start = next.begin()->first, current = start;
            std::vector<unsigned> loop;
            double turn = 0;
            do {
                auto found = next.find(current);
                if (found == next.end())
                    return false;
                unsigned to = found->second;
                next.erase(found);
                loop.push_back(current);
                Vec a = mesh.points[current] - p.origin, b = mesh.points[to] - p.origin;
                double delta = std::atan2(b.dot(f), b.dot(e)) - std::atan2(a.dot(f), a.dot(e));
                turn += std::atan2(std::sin(delta), std::cos(delta));
                current = to;
            } while (current != start);
            winding = winding || std::abs(turn) > 3.141592653589793;
            if (p.kind == "sphere") {
                Vec c = Vec::Zero();
                for (auto i : loop)
                    c += mesh.points[i];
                c /= double(loop.size());
                Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
                for (auto i : loop) {
                    Vec q = mesh.points[i] - c;
                    cov += q * q.transpose();
                }
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(cov);
                Vec n = eigen.eigenvectors().col(0);
                if (eigen.eigenvalues()[0] > 1e-8 * eigen.eigenvalues()[2] + 1e-10)
                    planar = false;
                if (first_normal.isZero())
                    first_normal = n;
                else if (std::abs(first_normal.dot(n)) < .99999)
                    planar = false;
            }
        }
        return p.kind == "sphere" ? planar : winding;
    };
    for (unsigned i = 0; i < patches.size(); ++i) {
        const auto &p = patches[i];
        for (auto t : p.triangles) {
            auto tri = mesh.triangles[t];
            areas[i] += (mesh.points[tri[1]] - mesh.points[tri[0]])
                            .cross(mesh.points[tri[2]] - mesh.points[tri[0]])
                            .norm() /
                        2;
        }
        auto tri = mesh.triangles[p.triangles.front()];
        normals[i] = (mesh.points[tri[1]] - mesh.points[tri[0]])
                         .cross(mesh.points[tri[2]] - mesh.points[tri[0]])
                         .normalized();
        bodies[{p.part, p.body}].push_back(i);
        if (p.kind == "cylinder" || p.kind == "sphere" || p.kind == "torus" || p.kind == "pipe") {
            double radius = p.kind == "torus" ? p.secondary : p.radius;
            if (radius > 0)
                values.push_back(
                    {radius * ((p.kind != "pipe" && full_turn(p)) ? 2 : 1), std::sqrt(areas[i])});
        }
    }
    for (const auto &body : bodies) {
        double total = 0;
        std::vector<unsigned> planes;
        for (auto i : body.second) {
            total += areas[i];
            if (patches[i].kind == "plane")
                planes.push_back(i);
        }
        std::sort(planes.begin(), planes.end(),
                  [&](unsigned a, unsigned b) { return areas[a] > areas[b]; });
        if (planes.size() > 40)
            planes.resize(40);
        for (unsigned a = 0; a < planes.size(); ++a)
            if (areas[planes[a]] > 1e-3 * total)
                for (unsigned b = 0; b < a; ++b)
                    if (normals[planes[a]].dot(normals[planes[b]]) < -.99999) {
                        double gap =
                            std::abs((patches[planes[a]].origin - patches[planes[b]].origin)
                                         .dot(normals[planes[a]]));
                        if (gap > 0)
                            values.push_back(
                                {gap, std::sqrt(std::min(areas[planes[a]], areas[planes[b]]))});
                    }
    }
    std::sort(values.begin(), values.end(),
              [](const Measurement &a, const Measurement &b) { return a.value < b.value; });
    std::vector<Measurement> distinct;
    for (auto v : values)
        if (std::isfinite(v.value) && v.value > 0 && v.weight > 0) {
            if (!distinct.empty() &&
                std::abs(v.value - distinct.back().value) <= 1e-3 * std::max(v.value, 1.0))
                distinct.back().weight = std::max(distinct.back().weight, v.weight);
            else
                distinct.push_back(v);
        }
    if (distinct.size() < 3)
        return;
    double weight = 0;
    for (auto v : distinct)
        weight += v.weight;
    auto roundness = [](double v, bool inch) {
        double x = v / (inch ? 25.4 : 1), relative = std::max(.0015 * x, .002 / (inch ? 25.4 : 1));
        std::vector<std::pair<double, double>> steps =
            inch ? std::vector<std::pair<double, double>>{{1. / 16, 1},
                                                          {.05, .5},
                                                          {1. / 64, .5},
                                                          {.01, .5}}
                 : std::vector<std::pair<double, double>>{{.5, 1}, {.1, .5}};
        for (auto step : steps) {
            double n = std::nearbyint(x / step.first);
            if (n != 0 && std::abs(x - step.first * n) <= std::min(relative, .05 * step.first))
                return step.second;
        }
        return 0.;
    };
    auto score = [&](double factor, bool inch) {
        double sum = 0;
        for (auto v : distinct)
            sum += v.weight * roundness(v.value * factor, inch);
        return sum / weight;
    };
    Vec lo = mesh.points[mesh.triangles.front()[0]], hi = lo;
    for (auto t : mesh.triangles)
        for (auto v : t) {
            lo = lo.cwiseMin(mesh.points[v]);
            hi = hi.cwiseMax(mesh.points[v]);
        }
    double size = (hi - lo).maxCoeff(), baseline = score(1, false);
    bool unusual = baseline < .35 || size < 5 || size > 1000;
    struct Candidate {
        std::string label;
        double factor;
        bool inch = false, file = false;
    };
    std::vector<Candidate> candidates{{"as is", 1}};
    for (double percent : {90., 95., 97., 97.5, 98., 98.5, 99., 99.5, 100.5, 101., 101.5, 102.,
                           102.5, 103., 105., 110.}) {
        std::ostringstream text;
        text << "scaled to " << percent << '%';
        candidates.push_back({text.str(), 100 / percent});
    }
    candidates.push_back({"drawn in inches", 1, true});
    candidates.push_back({"file in inches", 25.4, false, true});
    candidates.push_back({"file in centimetres", 10, false, true});
    candidates.push_back({"file in metres", 1000, false, true});
    candidates.push_back({"file in tenths of a mm", .1, false, true});
    Candidate best{"as is", 1};
    double best_score = -1;
    for (const auto &candidate : candidates) {
        if (size * candidate.factor < 1 || size * candidate.factor > 2000 ||
            (candidate.file && !unusual))
            continue;
        double s = score(candidate.factor, candidate.inch);
        if (s > best_score || (s == best_score && candidate.label > best.label)) {
            best = candidate;
            best_score = s;
        }
    }
    if (best_score < 0)
        return;
    if (best.label != "as is" &&
        !(best_score >= std::max(.5, 1.5 * baseline) && best_score - baseline >= .2))
        best = {"as is", 1};
    ctx.diagnostics.scale = best.factor;
    ctx.diagnostics.design_unit = best.inch ? "inch" : "mm";
    ctx.diagnostics.design_label = best.label;
    if (!ctx.options.true_size)
        return;
    for (auto &v : mesh.points)
        v *= best.factor;
    for (auto &p : patches) {
        p.origin *= best.factor;
        p.radius *= best.factor;
        if (p.kind == "torus" || p.kind.rfind("thread_", 0) == 0)
            p.secondary *= best.factor;
        p.deviation *= best.factor;
        if (p.kind == "plane" || p.kind == "cylinder" || p.kind == "cone" || p.kind == "sphere" ||
            p.kind == "torus")
            surface(p);
        else {
            gp_Trsf transform;
            transform.SetScale(gp_Pnt(0, 0, 0), best.factor);
            p.surface = Handle(Geom_Surface)::DownCast(p.surface->Transformed(transform));
        }
    }
    for (auto &p : patches)
        if (p.kind == "cylinder" || p.kind == "sphere") {
            double unit = best.inch ? 25.4 : 1, x = p.radius / unit;
            std::vector<double> steps = best.inch ? std::vector<double>{1. / 16, .05, 1. / 64, .01}
                                                  : std::vector<double>{.5, .1};
            for (double step : steps) {
                double target = std::nearbyint(x / step) * step * unit;
                if (target <= 0 ||
                    std::abs(target - p.radius) >
                        std::min(std::max(.0015 * p.radius, .002), .05 * step * unit))
                    continue;
                auto candidate = p;
                candidate.radius = target;
                bool fits = true;
                for (auto v : vertices(mesh, p.triangles))
                    if (deviation(candidate, mesh.points[v]) >
                        std::max(1e-4, std::min(.002, ctx.options.accuracy * .02)))
                        fits = false;
                if (fits) {
                    surface(candidate);
                    p = std::move(candidate);
                }
                break;
            }
        }
    ctx.diagnostics.repairs.push_back("Applied inferred design scale " +
                                      std::to_string(best.factor));
}
} // namespace stl2curves::detail

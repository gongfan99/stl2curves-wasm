#include "internal.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>

// Exercise the same private operations as repair_mesh, independently, to avoid
// a later repair hiding an earlier mismatch (e.g. refilling a canceled wall).
namespace stl2curves::detail {
#include "mesh_bodies.hpp"
#include "mesh_repair.hpp"
std::vector<Tri> oriented_triangles(std::vector<Tri> triangles) {
    for (auto &t : triangles)
        std::rotate(t.begin(), std::min_element(t.begin(), t.end()), t.end());
    std::sort(triangles.begin(), triangles.end());
    return triangles;
}
// A planar polygon may have several equally valid triangulations: Eigen and
// NumPy can choose opposite signs for their SVD basis. Compare the oriented
// boundary on each plane, so internal diagonals cancel but missing material,
// changed facing, folds, or holes still fail.
auto planar_boundaries(const Mesh &m) {
    using Plane = std::array<std::int64_t, 4>;
    std::map<Plane, std::map<Edge, int>> result;
    for (auto t : m.triangles) {
        const Vec &a = m.points[t[0]], &b = m.points[t[1]], &c = m.points[t[2]];
        Vec n = (b - a).cross(c - a).normalized();
        Plane plane{std::llround(n.x() * 1e9), std::llround(n.y() * 1e9), std::llround(n.z() * 1e9),
                    std::llround(n.dot(a) * 1e9)};
        auto &edges = result[plane];
        for (unsigned k = 0; k < 3; ++k) {
            auto i = t[k], j = t[(k + 1) % 3];
            Edge e(i, j);
            edges[e] += i < j ? 1 : -1;
            if (!edges[e])
                edges.erase(e);
        }
    }
    return result;
}
double point_triangle_gap(const Vec &p, const Vec &a, const Vec &b, const Vec &c) {
    Vec n = (b - a).cross(c - a);
    double gap = std::numeric_limits<double>::infinity();
    if (n.squaredNorm() > 1e-30) {
        Vec q = p - n * ((p - a).dot(n) / n.squaredNorm());
        if ((b - a).cross(q - a).dot(n) >= 0 && (c - b).cross(q - b).dot(n) >= 0 &&
            (a - c).cross(q - c).dot(n) >= 0)
            gap = (p - q).norm();
    }
    const Vec corners[]{a, b, c};
    for (unsigned k = 0; k < 3; ++k) {
        Vec edge = corners[(k + 1) % 3] - corners[k];
        double t =
            std::clamp((p - corners[k]).dot(edge) / std::max(edge.squaredNorm(), 1e-30), 0., 1.);
        gap = std::min(gap, (p - corners[k] - t * edge).norm());
    }
    return gap;
}
int mesh_parity(const char *file) {
    std::ifstream in(file);
    if (!in)
        throw std::runtime_error("Cannot open Python mesh fixtures");
    unsigned count, failures = 0;
    in >> count;
    for (unsigned i = 0; i < count; ++i) {
        std::string name, stage;
        unsigned np, nt, nq, nu, ng;
        in >> name >> stage >> np >> nt >> nq >> nu >> ng;
        Mesh m, expected;
        auto read = [&](Mesh &mesh, unsigned points, unsigned triangles) {
            mesh.points.resize(points);
            mesh.triangles.resize(triangles);
            for (auto &p : mesh.points)
                in >> p.x() >> p.y() >> p.z();
            for (auto &t : mesh.triangles)
                in >> t[0] >> t[1] >> t[2];
        };
        read(m, np, nt);
        read(expected, nq, nu);
        std::vector<unsigned> sizes(ng);
        for (auto &n : sizes)
            in >> n;
        if (!in)
            throw std::runtime_error("Truncated Python mesh fixture " + name);
        bool body_match = true;
        bool simplification_bound = true;
        if (stage == "bodies") {
            std::vector<std::vector<Tri>> actual, reference;
            for (const auto &group : mesh_body_groups(m)) {
                std::vector<Tri> tris;
                for (auto i : group)
                    tris.push_back(m.triangles[i]);
                actual.push_back(oriented_triangles(std::move(tris)));
            }
            unsigned offset = 0;
            for (auto n : sizes) {
                reference.push_back(oriented_triangles(std::vector<Tri>(
                    expected.triangles.begin() + offset, expected.triangles.begin() + offset + n)));
                offset += n;
            }
            std::sort(actual.begin(), actual.end());
            std::sort(reference.begin(), reference.end());
            body_match = actual == reference;
            m.triangles = expected.triangles;
        } else if (stage == "duplicates")
            remove_duplicates(m);
        else if (stage == "orient")
            orient_mesh(m);
        else if (stage == "slivers")
            split_slivers(m);
        else if (stage == "flaps")
            remove_flaps(m);
        else if (stage == "knots")
            untie_knots(m);
        else if (stage == "simplify") {
            Context ctx;
            simplify_mesh(m, .01, ctx);
            // Check against the original corners globally, independently of
            // the simplifier's per-collapse attachment bookkeeping.
            for (const auto &p : m.points) {
                double nearest = std::numeric_limits<double>::infinity();
                for (auto t : m.triangles)
                    nearest = std::min(nearest, point_triangle_gap(p, m.points[t[0]],
                                                                   m.points[t[1]], m.points[t[2]]));
                simplification_bound = simplification_bound && nearest <= .01 + 1e-10;
            }
        } else if (stage == "cracks")
            close_cracks(m, .001);
        else if (stage == "holes")
            fill_holes(m);
        else if (stage == "repair") {
            Context ctx;
            repair_mesh(m, ctx);
        } else
            throw std::runtime_error("Unknown parity stage " + stage);
        bool same_triangles =
            oriented_triangles(m.triangles) == oriented_triangles(expected.triangles);
        if (!same_triangles && name.rfind("large_hole_", 0) == 0 &&
            m.triangles.size() == expected.triangles.size())
            same_triangles = planar_boundaries(m) == planar_boundaries(expected);
        bool ok = m.points.size() == expected.points.size() && same_triangles && body_match &&
                  simplification_bound;
        for (unsigned j = 0; ok && j < m.points.size(); ++j)
            ok = (m.points[j] - expected.points[j]).norm() <= 1e-12;
        std::cout << (ok ? "PASS " : "FAIL ") << stage << ' ' << name << " triangles "
                  << m.triangles.size() << '/' << expected.triangles.size() << '\n';
        failures += !ok;
    }
    std::cout << count - failures << '/' << count << " Python mesh parity cases passed\n";
    return failures ? 1 : 0;
}
} // namespace stl2curves::detail
int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Usage: test_python_mesh_parity mesh_cases.txt\n";
        return 2;
    }
    try {
        return stl2curves::detail::mesh_parity(argv[1]);
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

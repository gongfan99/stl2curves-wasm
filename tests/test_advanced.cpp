#include "feature_islands.hpp"
#include "internal.hpp"
#include "projection.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <Eigen/Geometry>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Curve.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
using namespace stl2curves::detail;
namespace {
constexpr double pi = 3.14159265358979323846;
void check(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
void connect(Mesh &mesh, unsigned rings, unsigned around) {
    for (unsigned i = 0; i + 1 < rings; ++i)
        for (unsigned j = 0; j < around; ++j) {
            unsigned k = (j + 1) % around, a = i * around + j, b = i * around + k,
                     c = (i + 1) * around + j, d = (i + 1) * around + k;
            mesh.triangles.push_back({a, b, d});
            mesh.triangles.push_back({a, d, c});
        }
}
void thread_test(double handedness, bool rotated) {
    Mesh mesh;
    constexpr unsigned around = 24, rings = 65;
    for (unsigned i = 0; i < rings; ++i)
        for (unsigned j = 0; j < around; ++j) {
            double z = 4. * i / (rings - 1), u = 2 * pi * j / around,
                   r = 3 + .08 * std::cos(2 * pi * z / (handedness * 1.07) - u);
            mesh.points.emplace_back(r * std::cos(u), r * std::sin(u), z);
        }
    connect(mesh, rings, around);
    Vec axis = Vec::UnitZ();
    if (rotated) {
        Eigen::AngleAxisd rotation(.73, Vec(1, 2, 3).normalized());
        axis = rotation * axis;
        for (auto &p : mesh.points)
            p = rotation * p + Vec(4, -7, 11);
    }
    auto top = topology(mesh);
    std::vector<unsigned> ids(mesh.triangles.size());
    std::iota(ids.begin(), ids.end(), 0);
    auto feature = thread_patch(mesh, top, ids, axis, .002);
    check(!feature.surface.IsNull(), "helical surface fit");
    check(feature.kind == (handedness > 0 ? "thread_right" : "thread_left") &&
              std::abs(feature.secondary - handedness * 1.07) < .005,
          "thread pitch and handedness");
    check(feature.deviation <= .002, "thread deviation");
}
void extrusion_test() {
    Mesh mesh;
    constexpr unsigned around = 64, rings = 3;
    Vec axis = Vec(1, 2, 3).normalized(), e = axis.unitOrthogonal(), f = axis.cross(e);
    for (unsigned i = 0; i < rings; ++i)
        for (unsigned j = 0; j < around; ++j) {
            double u = 2 * pi * j / around, r = 1 + .1 * std::sin(3 * u) + .07 * std::cos(5 * u);
            mesh.points.push_back(Vec(3, -2, 7) + double(i) * axis +
                                  r * (std::cos(u) * e + std::sin(u) * f));
        }
    connect(mesh, rings, around);
    auto top = topology(mesh);
    std::vector<int> owner(mesh.triangles.size(), -1);
    Context ctx;
    auto features = advanced_patches(mesh, top, owner, {}, ctx);
    bool fitted = false;
    for (const auto &p : features)
        if (p.kind == "extrusion" && p.deviation <= .002) {
            fitted = true;
            check(p.surface->IsUPeriodic(), "closed extrusion profile must retain its seam");
        }
    check(fitted, "noncircular rotated extrusion recognition");
}
void pipe_test() {
    Mesh mesh;
    constexpr unsigned around = 32, rings = 49;
    for (unsigned i = 0; i < rings; ++i)
        for (unsigned j = 0; j < around; ++j) {
            double x = 6. * i / (rings - 1), u = 2 * pi * j / around;
            Vec tangent(1, .05 * std::cos(x), 0);
            tangent.normalize();
            Vec e = Vec::UnitZ(), f = tangent.cross(e);
            mesh.points.push_back(Vec(x, .05 * std::sin(x), 0) +
                                  .2 * (std::cos(u) * e + std::sin(u) * f));
        }
    connect(mesh, rings, around);
    auto top = topology(mesh);
    std::vector<int> owner(mesh.triangles.size(), -1);
    Patch known;
    known.kind = "cylinder";
    known.radius = .2;
    Context context;
    auto features = advanced_patches(mesh, top, owner, {known}, context);
    bool fitted = false;
    for (const auto &feature : features)
        if (feature.kind == "pipe" && feature.deviation <= .002)
            fitted = true;
    check(fitted, "constant-radius swept pipe fit");
}
void steep_freeform_test(bool rotated) {
    // A known quadratic roof, sampled coarsely enough that vertical error
    // overstates its normal distance near the steep sides. The true smooth
    // surface should remain a light height field rather than follow each facet.
    constexpr unsigned n = 17;
    const Eigen::AngleAxisd rotation(rotated ? .73 : 0., Vec(1, 2, 3).normalized());
    const Vec offset = rotated ? Vec(4, -7, 11) : Vec::Zero();
    auto point = [&](double x, double y) -> Vec {
        return rotation * Vec(x, y, .6 * x * x + .1 * y * y + .05 * x * y) + offset;
    };
    Mesh mesh;
    for (unsigned i = 0; i < n; ++i)
        for (unsigned j = 0; j < n; ++j)
            mesh.points.push_back(point(-2. + 4. * i / (n - 1), -1. + 2. * j / (n - 1)));
    for (unsigned i = 0; i + 1 < n; ++i)
        for (unsigned j = 0; j + 1 < n; ++j) {
            unsigned a = i * n + j, b = (i + 1) * n + j, c = (i + 1) * n + j + 1, d = i * n + j + 1;
            mesh.triangles.push_back({a, b, c});
            mesh.triangles.push_back({a, c, d});
        }
    auto top = topology(mesh);
    Context context;
    std::vector<int> owner(mesh.triangles.size(), -1);
    auto features = advanced_patches(mesh, top, owner, {}, context);
    check(features.size() == 1 && features[0].triangles.size() == mesh.triangles.size(),
          "steep quadratic roof must remain one complete smooth patch");
    auto spline = Handle(Geom_BSplineSurface)::DownCast(features[0].surface);
    check(!spline.IsNull() && spline->NbUPoles() * spline->NbVPoles() <= 64,
          "steep quadratic roof must not acquire facet-scale spline oscillations");
    PointProjector projector(spline);
    for (unsigned i = 0; i < 51; ++i)
        for (unsigned j = 0; j < 31; ++j) {
            const auto &projection = projector.perform(point(-2. + 4. * i / 50, -1. + 2. * j / 30));
            check(projection.NbPoints() && projection.LowerDistance() < .001,
                  "steep height field must follow the known smooth surface between mesh corners");
        }
}
void compact_freeform_test() {
    // Sparse samples can force precise fits to split one smooth quartic roof.
    // At the requested accuracy a lighter whole surface should replace those
    // artificial boundaries without losing any input triangles.
    constexpr unsigned n = 13;
    auto point = [](double x, double y) -> Vec {
        return Vec(x, y, .005 * x * x * x * x + .1 * y * y + .03 * x * y);
    };
    Mesh mesh;
    for (unsigned i = 0; i < n; ++i)
        for (unsigned j = 0; j < n; ++j)
            mesh.points.push_back(point(-4. + 8. * i / (n - 1), -1. + 2. * j / (n - 1)));
    for (unsigned i = 0; i + 1 < n; ++i)
        for (unsigned j = 0; j + 1 < n; ++j) {
            unsigned a = i * n + j, b = (i + 1) * n + j, c = (i + 1) * n + j + 1, d = i * n + j + 1;
            mesh.triangles.push_back({a, b, c});
            mesh.triangles.push_back({a, c, d});
        }
    Context context;
    auto top = topology(mesh);
    std::vector<int> owner(mesh.triangles.size(), -1);
    auto features = advanced_patches(mesh, top, owner, {}, context);
    check(features.size() == 1 && features[0].triangles.size() == mesh.triangles.size(),
          "smooth quartic roof must not retain artificial freeform partition boundaries");
    PointProjector projector(features[0].surface);
    for (unsigned i = 0; i < 51; ++i)
        for (unsigned j = 0; j < 31; ++j) {
            const auto &projection = projector.perform(point(-4. + 8. * i / 50, -1. + 2. * j / 30));
            check(projection.NbPoints() && projection.LowerDistance() < .01,
                  "compact quartic fit must preserve the known smooth shape");
        }
}
void python_feature_islands_test() {
    auto file = std::filesystem::path(__FILE__).parent_path() / "fixtures/python_feature_cases.txt";
    std::ifstream input(file);
    check(bool(input), "Python feature oracle fixture must be readable");
    unsigned cases;
    input >> cases;
    auto indices = [&] {
        unsigned count;
        input >> count;
        std::vector<unsigned> result(count);
        for (auto &i : result)
            input >> i;
        return result;
    };
    for (unsigned c = 0; c < cases; ++c) {
        std::string name, kind;
        double size;
        unsigned nv, nt, np, nf;
        input >> name >> kind >> size >> nv >> nt >> np >> nf;
        Mesh mesh;
        for (unsigned i = 0; i < nv; ++i) {
            double x, y, z;
            input >> x >> y >> z;
            mesh.points.emplace_back(x, y, z);
        }
        mesh.triangles.resize(nt);
        for (auto &t : mesh.triangles)
            input >> t[0] >> t[1] >> t[2];
        std::vector<std::vector<unsigned>> expected_facets;
        for (unsigned f = 0; f < nf; ++f)
            expected_facets.push_back(indices());
        std::vector<std::vector<unsigned>> expected_neighbors;
        for (unsigned f = 0; f < nf; ++f)
            expected_neighbors.push_back(indices());
        std::vector<Patch> patches;
        std::vector<std::vector<unsigned>> expected;
        for (unsigned p = 0; p < np; ++p) {
            Patch patch;
            patch.kind = kind;
            patch.triangles = indices();
            expected.push_back(indices());
            patch.radius = size;
            if (kind == "freeform") {
                PointGrid poles(1, 4, 1, 4);
                const double q[] = {0, 0, size * size / 3, size * size};
                for (unsigned i = 0; i < 4; ++i)
                    for (unsigned j = 0; j < 4; ++j)
                        poles.SetValue(
                            i + 1, j + 1,
                            gp_Pnt(size * i / 3, size * j / 3, .005 * q[i] + .007 * q[j]));
                RealArray knots(1, 2);
                knots.SetValue(1, 0);
                knots.SetValue(2, size);
                IntArray mult(1, 2);
                mult.SetValue(1, 4);
                mult.SetValue(2, 4);
                patch.surface = new Geom_BSplineSurface(poles, knots, knots, mult, mult, 3, 3);
            }
            patches.push_back(std::move(patch));
        }
        check(bool(input), "Python feature oracle fixture must be complete");
        const auto top = topology(mesh);
        FeatureFacets facets(mesh, top);
        if (facets.triangles != expected_facets)
            throw std::runtime_error(name + ": facet grouping differs from Python");
        if (facets.neighbors != expected_neighbors)
            throw std::runtime_error(name + ": facet neighbor order differs from Python");
        absorb_feature_islands(mesh, top, patches);
        if (patches.size() != expected.size())
            throw std::runtime_error(name + ": island ownership changed the patch count");
        for (unsigned p = 0; p < patches.size(); ++p) {
            auto actual = patches[p].triangles;
            std::sort(actual.begin(), actual.end());
            std::sort(expected[p].begin(), expected[p].end());
            if (actual != expected[p])
                throw std::runtime_error(name + ": island triangle ownership differs from Python");
        }
    }
    unsigned order_cases;
    input >> order_cases;
    check(bool(input), "Python neighbor-order oracle must be readable");
    for (unsigned c = 0; c < order_cases; ++c) {
        auto inserted = indices();
        auto expected = indices();
        FeatureNeighborOrder order;
        for (auto neighbor : inserted)
            order.add(neighbor);
        check(order.values() == expected, "Neighbor insertion/resize order differs from Python");
    }
    check(bool(input), "Python neighbor-order oracle must be complete");
}
void python_corner_reference_test() {
    auto path =
        std::filesystem::path(__FILE__).parent_path() / "fixtures/python_corner_reference.txt";
    std::ifstream input(path);
    check(bool(input), "Python corner reference fixture must be readable");
    unsigned nv, nt, nf;
    input >> nv >> nt >> nf;
    Mesh mesh;
    for (unsigned i = 0; i < nv; ++i) {
        double x, y, z;
        input >> x >> y >> z;
        mesh.points.emplace_back(x, y, z);
    }
    mesh.triangles.resize(nt);
    for (auto &triangle : mesh.triangles)
        input >> triangle[0] >> triangle[1] >> triangle[2];
    Context context;
    context.options.fuse = false;
    context.options.time_limit = 0;
    auto patches = detect(mesh, context);
    for (const auto &patch : patches)
        check(!patch.surface.IsNull(),
              "Corner output patch must have a surface for reconstruction");
    std::vector<bool> matched(patches.size(), false);
    double largest_gap = 0;
    for (unsigned f = 0; f < nf; ++f) {
        std::string kind;
        bool trimmed;
        double radius, secondary;
        Vec origin, axis;
        input >> kind >> trimmed >> radius >> secondary >> origin[0] >> origin[1] >> origin[2] >>
            axis[0] >> axis[1] >> axis[2];
        unsigned count;
        input >> count;
        std::vector<unsigned> expected(count);
        for (auto &triangle : expected)
            input >> triangle;
        std::sort(expected.begin(), expected.end());
        std::size_t found = patches.size();
        for (std::size_t k = 0; k < patches.size(); ++k) {
            auto owned = patches[k].triangles;
            std::sort(owned.begin(), owned.end());
            if (patches[k].kind == kind && owned == expected) {
                check(!matched[k], "Corner reference feature matched more than once");
                found = k;
                break;
            }
        }
        if (found == patches.size())
            throw std::runtime_error("Corner reference ownership differs for feature " +
                                     std::to_string(f) + " (" + kind + ")");
        matched[found] = true;
        const auto &patch = patches[found];
        check(patch.trimmed == trimmed, "Corner reference trimmed/natural classification differs");
        if (kind != "freeform") {
            check(std::abs(patch.radius - radius) < 1e-9 &&
                      std::abs(patch.secondary - secondary) < 1e-9,
                  "Corner reference analytic radii differ");
            check((patch.origin - origin).norm() < 1e-9 && (patch.axis - axis).norm() < 1e-9,
                  "Corner reference analytic frame differs");
        }
        unsigned samples;
        input >> samples;
        for (unsigned s = 0; s < samples; ++s) {
            Vec p, n;
            input >> p[0] >> p[1] >> p[2] >> n[0] >> n[1] >> n[2];
            PointProjector projector(patch.surface);
            const auto &projection = projector.perform(p);
            check(projection.NbPoints() != 0, "Corner reference surface projection failed");
            largest_gap = std::max(largest_gap, projection.LowerDistance());
            double u, v;
            projection.LowerDistanceParameters(u, v);
            gp_Pnt at;
            gp_Vec du, dv;
            patch.surface->D1(u, v, at, du, dv);
            auto cross = du.Crossed(dv);
            Vec actual(cross.X(), cross.Y(), cross.Z());
            check(std::abs(actual.normalized().dot(n)) > 1 - 1e-9,
                  "Corner reference surface normal differs");
        }
    }
    check(bool(input), "Python corner reference fixture must be complete");
    for (std::size_t k = 0; k < patches.size(); ++k)
        check(matched[k] || patches[k].kind == "plane", "Unexpected curved corner patch");
    std::cout << "Python corner reference: " << nf << " exact ownership groups, max surface gap "
              << largest_gap << " mm\n";
    check(largest_gap < 1e-7, "Corner reference B-spline surfaces differ");
    auto shape = reconstruct(mesh, patches, context);
    stl2curves::Diagnostics diagnostics;
    check(inspect(shape, diagnostics) && diagnostics.faces == 30 && diagnostics.solids == 1,
          "Corner reference reconstruction must be one valid 30-face solid");
    ShapeMap edges;
    TopExp::MapShapes(shape, TopAbs_EDGE, edges);
    struct Curve {
        int kind;
        Handle(Geom_Curve) geometry;
        double first, last;
    };
    std::vector<Curve> curves;
    std::map<int, unsigned> actual_counts, expected_counts;
    for (int i = 1; i <= edges.Extent(); ++i) {
        BRepAdaptor_Curve curve(TopoDS::Edge(edges(i)));
        auto geometry = Handle(Geom_Curve)::DownCast(curve.Curve().Curve()->Copy());
        geometry->Transform(curve.Trsf());
        curves.push_back(
            {int(curve.GetType()), geometry, curve.FirstParameter(), curve.LastParameter()});
        ++actual_counts[int(curve.GetType())];
    }
    unsigned expected_edges;
    input >> expected_edges;
    check(bool(input), "Python boundary reference must be readable");
    check(unsigned(edges.Extent()) == expected_edges, "Corner reference edge count differs");
    double edge_gap = 0;
    for (unsigned e = 0; e < expected_edges; ++e) {
        int kind;
        unsigned samples;
        input >> kind >> samples;
        ++expected_counts[kind];
        for (unsigned s = 0; s < samples; ++s) {
            double x, y, z;
            input >> x >> y >> z;
            double best = std::numeric_limits<double>::infinity();
            for (const auto &curve : curves)
                if (curve.kind == kind) {
                    GeomAPI_ProjectPointOnCurve project(gp_Pnt(x, y, z), curve.geometry,
                                                        curve.first, curve.last);
                    if (project.NbPoints())
                        best = std::min(best, project.LowerDistance());
                }
            edge_gap = std::max(edge_gap, best);
        }
    }
    check(bool(input), "Python boundary reference must be complete");
    check(actual_counts == expected_counts, "Corner reference curve types/counts differ");
    std::cout << "Python corner boundaries: " << expected_edges << " edges, max gap " << edge_gap
              << " mm\n";
    check(edge_gap < 1e-7, "Corner reference boundary geometry differs");
    double expected_volume, expected_epsilon;
    input >> expected_volume >> expected_epsilon;
    check(bool(input), "Python corner reference volume must be readable");
    std::cout << "Python corner volume difference: " << diagnostics.volume - expected_volume
              << " mm^3\n";
    check(volume_epsilon == expected_epsilon, "Corner reference integration precision differs");
    // Adaptive integration is an estimate at this relative precision;
    // surface and boundary geometry are checked independently above.
    check(std::abs(diagnostics.volume - expected_volume) < expected_epsilon * expected_volume,
          "Corner reference volume differs at the reference integration precision");
    Context tighter;
    tighter.options.fuse = false;
    tighter.options.time_limit = 0;
    tighter.options.accuracy = .0001;
    auto tight_patches = detect(mesh, tighter);
    std::vector<bool> covered(mesh.triangles.size(), false);
    for (const auto &patch : tight_patches) {
        check(!patch.surface.IsNull(), "Tight corner fit must retain reconstructable surfaces");
        for (auto triangle : patch.triangles) {
            check(!covered.at(triangle), "Tight corner fit duplicates a triangle");
            covered[triangle] = true;
        }
    }
    check(std::all_of(covered.begin(), covered.end(), [](bool owned) { return owned; }),
          "Tight corner fit loses mesh triangles");
}
} // namespace
int main() {
    try {
        thread_test(1, false);
        thread_test(-1, true);
        pipe_test();
        extrusion_test();
        steep_freeform_test(false);
        steep_freeform_test(true);
        compact_freeform_test();
        python_feature_islands_test();
        python_corner_reference_test();
        std::cout << "Advanced fitting tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

#include "internal.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepLib.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GeomConvert.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <STEPControl_Reader.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Ax3.hxx>
#include <iostream>
#include <sstream>

using namespace stl2curves;
using namespace stl2curves::detail;
namespace {
constexpr double pi = 3.14159265358979323846;
void check(bool ok, const std::string &message) {
    if (!ok)
        throw std::runtime_error(message);
}
void analytic_edges(const TopoDS_Shape &shape) {
    ShapeMap edges;
    TopExp::MapShapes(shape, TopAbs_EDGE, edges);
    unsigned circles = 0;
    for (int i = 1; i <= edges.Extent(); ++i) {
        auto type = BRepAdaptor_Curve(TopoDS::Edge(edges(i))).GetType();
        check(type == GeomAbs_Line || type == GeomAbs_Circle,
              "analytic cylinder boundaries must remain lines and circles");
        circles += type == GeomAbs_Circle;
    }
    check(circles > 0, "circular boundaries must be retained");
}
// A quarter cylinder, optionally with a step in its top. The step forces
// arbitrary trimming: its cylinder cannot use a rectangular parameter domain.
void cylinder_boundaries(double start, bool stepped, bool recovery, bool bad_plane = false,
                         bool noisy_rim = false, bool short_edge = false) {
    constexpr unsigned segments = 36, ring = segments + 1;
    constexpr double radius = 10, high = 3.5, low = 2;
    Mesh mesh;
    const std::vector<double> heights =
        stepped ? std::vector<double>{0, low, high} : std::vector<double>{0, high};
    for (double z : heights)
        for (unsigned i = 0; i <= segments; ++i) {
            double angle = start + pi / 2 * i / segments;
            if (short_edge && i == 1)
                angle = start + 1e-4;
            double axial_noise = noisy_rim ? (i % 2 ? 2e-5 : -2e-5) : 0;
            mesh.points.emplace_back(radius * std::cos(angle), radius * std::sin(angle),
                                     z + axial_noise);
        }
    const unsigned centers = unsigned(mesh.points.size());
    for (double z : heights)
        mesh.points.emplace_back(0, 0, z);
    auto vertex = [](unsigned layer, unsigned i) { return layer * ring + i; };
    std::vector<Patch> patches;
    Patch curved;
    curved.kind = "cylinder";
    curved.radius = radius;
    curved.surface = new Geom_CylindricalSurface(gp_Ax3(), radius);
    patches.push_back(curved);
    auto plane = [&](double z, Vec normal) {
        Patch p;
        p.kind = "plane";
        p.origin = Vec(0, 0, z);
        p.axis = normal;
        p.surface = new Geom_Plane(gp_Pnt(0, 0, z), gp_Dir(normal.x(), normal.y(), normal.z()));
        patches.push_back(p);
        return unsigned(patches.size() - 1);
    };
    auto bottom = plane(0, -Vec::UnitZ());
    auto top = plane(high, Vec::UnitZ());
    auto start_face = plane(0, Vec(std::sin(start), -std::cos(start), 0));
    double end = start + pi / 2;
    auto end_face = plane(0, Vec(-std::sin(end), std::cos(end), 0));
    unsigned low_top = 0, step_face = 0;
    if (stepped) {
        low_top = plane(low, Vec::UnitZ());
        double middle = start + pi / 4;
        step_face = plane(low, Vec(-std::sin(middle), std::cos(middle), 0));
    }
    auto triangle = [&](unsigned patch, unsigned a, unsigned b, unsigned c) {
        patches[patch].triangles.push_back(unsigned(mesh.triangles.size()));
        mesh.triangles.push_back({a, b, c});
    };
    auto side = [&](unsigned patch, unsigned layer, unsigned i, bool radial) {
        unsigned a = radial ? centers + layer : vertex(layer, i);
        unsigned b = vertex(layer, radial ? i : i + 1);
        unsigned c = radial ? centers + layer + 1 : vertex(layer + 1, i);
        unsigned d = vertex(layer + 1, radial ? i : i + 1);
        triangle(patch, a, b, d);
        triangle(patch, a, d, c);
    };
    for (unsigned i = 0; i < segments; ++i) {
        triangle(bottom, centers, vertex(0, i + 1), vertex(0, i));
        side(0, 0, i, false);
        unsigned layer = 1;
        unsigned cap = top;
        if (stepped) {
            if (i < segments / 2) {
                side(0, 1, i, false);
                layer = 2;
            } else
                cap = low_top;
        }
        triangle(cap, centers + layer, vertex(layer, i), vertex(layer, i + 1));
    }
    side(start_face, 0, 0, true);
    // End and step walls face towards increasing angle.
    auto end_wall = [&](unsigned patch, unsigned layer, unsigned i) {
        triangle(patch, centers + layer, centers + layer + 1, vertex(layer + 1, i));
        triangle(patch, centers + layer, vertex(layer + 1, i), vertex(layer, i));
    };
    end_wall(end_face, 0, segments);
    if (stepped) {
        side(start_face, 1, 0, true);
        end_wall(step_face, 1, segments / 2);
    }
    if (bad_plane) {
        // A defective fitted support must not prevent restoration of the
        // independently valid cylinder beside it.
        auto &patch = patches[bottom];
        patch.origin.z() -= .05;
        patch.surface = new Geom_Plane(gp_Pnt(0, 0, -.05), gp_Dir(0, 0, -1));
    }
    Context ctx;
    auto shape = recovery ? recover(mesh, patches, ctx) : reconstruct(mesh, patches, ctx);
    Diagnostics info;
    check(inspect(shape, info), "reconstructed cylinder must be a valid closed solid");
    check(info.faces == (stepped ? 7u : 5u), "cylinder face count");
    double expected_volume = radius * radius * pi / 4 * (stepped ? (high + low) / 2 : high);
    // Axial rounding on both rims can change the cylinder's bound by 40 nm
    // in height (0.0032 mm^3 here); exact inputs retain the tighter check.
    check(std::abs(info.volume - expected_volume) < (noisy_rim ? .005 : .001),
          "analytic cylinder volume: expected=" + std::to_string(expected_volume) +
              " actual=" + std::to_string(info.volume) + " noisy=" + std::to_string(noisy_rim) +
              " bad_plane=" + std::to_string(bad_plane));
    analytic_edges(shape);
    if (!stepped) {
        ShapeMap edges;
        TopExp::MapShapes(shape, TopAbs_EDGE, edges);
        check(edges.Extent() == 9, "quarter cylinder must have nine unsplit edges");
    }
    auto bytes = write_step(shape, ctx.options.tolerance);
    check(verify_step(bytes, shape, ctx.options.tolerance), "cylinder STEP round trip");
    STEPControl_Reader reader;
    std::istringstream stream(std::string(bytes.begin(), bytes.end()));
    check(reader.ReadStream("boundaries.step", stream) == IFSelect_RetDone &&
              reader.TransferRoots() > 0,
          "read cylinder STEP");
    analytic_edges(reader.OneShape());
}
void freeform_solid() {
    // A quadratic roof has curved intersections with four planar walls.
    // It should remain a six-face solid, including after STEP serialization.
    constexpr unsigned n = 17, layer = n * n;
    Mesh mesh;
    for (unsigned z = 0; z < 2; ++z)
        for (unsigned y = 0; y < n; ++y)
            for (unsigned x = 0; x < n; ++x) {
                double u = -5 + 10. * x / (n - 1), v = -5 + 10. * y / (n - 1);
                mesh.points.emplace_back(u, v, z ? 3 + .02 * u * u + .03 * v * v + .01 * u * v : 0);
            }
    auto id = [](unsigned z, unsigned x, unsigned y) { return z * layer + y * n + x; };
    auto triangle = [&](unsigned a, unsigned b, unsigned c) {
        mesh.triangles.push_back({a, b, c});
    };
    for (unsigned y = 0; y + 1 < n; ++y)
        for (unsigned x = 0; x + 1 < n; ++x) {
            auto a = id(0, x, y), b = id(0, x + 1, y), c = id(0, x, y + 1), d = id(0, x + 1, y + 1);
            triangle(a, d, b);
            triangle(a, c, d);
            triangle(a + layer, b + layer, d + layer);
            triangle(a + layer, d + layer, c + layer);
        }
    auto wall = [&](unsigned a, unsigned b, bool reverse) {
        if (reverse) {
            triangle(a, a + layer, b + layer);
            triangle(a, b + layer, b);
        } else {
            triangle(a, b, b + layer);
            triangle(a, b + layer, a + layer);
        }
    };
    for (unsigned i = 0; i + 1 < n; ++i) {
        wall(id(0, i, 0), id(0, i + 1, 0), false);
        wall(id(0, i, n - 1), id(0, i + 1, n - 1), true);
        wall(id(0, 0, i), id(0, 0, i + 1), true);
        wall(id(0, n - 1, i), id(0, n - 1, i + 1), false);
    }
    Context ctx;
    auto patches = detect(mesh, ctx);
    auto shape = reconstruct(mesh, patches, ctx);
    Diagnostics info;
    check(inspect(shape, info), "freeform roof must be a valid closed solid");
    check(info.faces == 6 && info.surface_counts["bspline"] == 1,
          "freeform roof must retain one curved face and five planar faces");
    check(std::abs(info.volume - (300 + 125. / 3)) < .05,
          "quadratic roof volume: " + std::to_string(info.volume));
    auto bytes = write_step(shape, ctx.options.tolerance);
    Diagnostics restored;
    check(verify_step(bytes, shape, ctx.options.tolerance, &restored),
          "freeform roof STEP round trip");
    check(restored.faces == 6 && restored.surface_counts["bspline"] == 1,
          "freeform roof STEP topology");

    // Force the recovery builder, rather than allowing the initial shared-run
    // builder to hide a per-chord trim regression. Both sides of each curved
    // roof/wall junction must retain one edge, as Python's _wire does.
    auto recovered = recover(mesh, patches, ctx);
    ShapeMap input_vertices;
    TopExp::MapShapes(recovered, TopAbs_VERTEX, input_vertices);
    BRep_Builder builder;
    std::vector<double> input_tolerances;
    for (int i = 1; i <= input_vertices.Extent(); ++i) {
        auto vertex = TopoDS::Vertex(input_vertices(i));
        builder.UpdateVertex(vertex, .1);
        input_tolerances.push_back(BRep_Tool::Tolerance(vertex));
    }
    recovered = join_boundary_runs(recovered, mesh, patches, ctx.options.tolerance);
    for (int i = 1; i <= input_vertices.Extent(); ++i)
        check(BRep_Tool::Tolerance(TopoDS::Vertex(input_vertices(i))) == input_tolerances[i - 1],
              "boundary joining must preserve the input's sewing tolerances");
    check(inspect(recovered, info) && info.faces == 6,
          "recovered freeform roof retains all six faces");
    ShapeMap recovered_edges;
    TopExp::MapShapes(recovered, TopAbs_EDGE, recovered_edges);
    check(recovered_edges.Extent() == 12,
          "recovered roof must use twelve whole boundary edges, not mesh-chord segments (got " +
              std::to_string(recovered_edges.Extent()) + ")");
    unsigned splines = 0;
    for (int i = 1; i <= recovered_edges.Extent(); ++i)
        splines +=
            BRepAdaptor_Curve(TopoDS::Edge(recovered_edges(i))).GetType() == GeomAbs_BSplineCurve;
    check(splines == 4, "each recovered curved roof rim must be one smooth spline");
}
void planar_strip_roof() {
    // A curved roof exported as broad, planar strips, each with six triangles.
    // Recognizing those strips as permanent planes prevents the smooth fitter
    // from recovering the roof. The large side walls establish its facet scale.
    constexpr unsigned nx = 16, ny = 4, layer = nx * ny;
    Mesh mesh;
    for (unsigned z = 0; z < 2; ++z)
        for (unsigned y = 0; y < ny; ++y)
            for (unsigned x = 0; x < nx; ++x) {
                double u = 10. * x / (nx - 1), v = .3 * y / (ny - 1);
                mesh.points.emplace_back(u, v, z ? 3 + .05 * u * u : 0);
            }
    auto id = [](unsigned x, unsigned y) { return y * nx + x; };
    auto triangle = [&](unsigned a, unsigned b, unsigned c) {
        mesh.triangles.push_back({a, b, c});
    };
    for (unsigned y = 0; y + 1 < ny; ++y)
        for (unsigned x = 0; x + 1 < nx; ++x) {
            unsigned a = id(x, y), b = id(x + 1, y), c = id(x, y + 1), d = id(x + 1, y + 1);
            triangle(a, d, b);
            triangle(a, c, d);
            triangle(a + layer, b + layer, d + layer);
            triangle(a + layer, d + layer, c + layer);
        }
    auto wall = [&](unsigned a, unsigned b) {
        triangle(a, b, b + layer);
        triangle(a, b + layer, a + layer);
    };
    for (unsigned x = 0; x + 1 < nx; ++x) {
        wall(id(x, 0), id(x + 1, 0));
        wall(id(x + 1, ny - 1), id(x, ny - 1));
    }
    for (unsigned y = 0; y + 1 < ny; ++y) {
        wall(id(0, y + 1), id(0, y));
        wall(id(nx - 1, y), id(nx - 1, y + 1));
    }
    Context ctx;
    auto patches = detect(mesh, ctx);
    auto shape = reconstruct(mesh, patches, ctx);
    Diagnostics info;
    check(inspect(shape, info), "planar-strip roof must close");
    check(info.faces == 6 && info.surface_counts["plane"] == 5,
          "planar mesh strips must recover one smooth roof and five planes");
    check(std::abs(info.volume - 14.) < .01, "planar-strip roof volume");
}
void periodic_spline_boundaries() {
    // Generic periodic B-spline trims create separate charts at their seam.
    // They must merge on their shared support instead of retaining the charts
    // as separate editable faces (even when all geometry is copied for repair).
    constexpr unsigned count = 192;
    constexpr double radius = 2, height = 3;
    Mesh mesh;
    for (unsigned layer = 0; layer < 2; ++layer)
        for (unsigned i = 0; i < count; ++i) {
            double angle = 2 * pi * i / count;
            mesh.points.emplace_back(radius * std::cos(angle), radius * std::sin(angle),
                                     layer * height);
        }
    mesh.points.emplace_back(0, 0, 0);
    mesh.points.emplace_back(0, 0, height);
    std::vector<Patch> patches(3);
    auto spline = GeomConvert::SurfaceToBSplineSurface(new Geom_RectangularTrimmedSurface(
        new Geom_CylindricalSurface(gp_Ax3(), radius), 0, 2 * pi, 0, height));
    if (!spline->IsUPeriodic())
        spline->SetUPeriodic();
    patches[0].kind = "freeform";
    patches[0].surface = spline;
    for (unsigned i = 1; i < 3; ++i) {
        double z = i == 1 ? 0 : height;
        patches[i].kind = "plane";
        patches[i].origin = Vec(0, 0, z);
        patches[i].axis = Vec(0, 0, i == 1 ? -1 : 1);
        patches[i].surface = new Geom_Plane(gp_Pnt(0, 0, z), gp_Dir(0, 0, i == 1 ? -1 : 1));
    }
    auto triangle = [&](unsigned patch, unsigned a, unsigned b, unsigned c) {
        patches[patch].triangles.push_back(unsigned(mesh.triangles.size()));
        mesh.triangles.push_back({a, b, c});
    };
    for (unsigned i = 0; i < count; ++i) {
        unsigned j = (i + 1) % count;
        triangle(0, i, j, count + j);
        triangle(0, i, count + j, count + i);
        triangle(1, 2 * count, j, i);
        triangle(2, 2 * count + 1, count + i, count + j);
    }
    Context ctx;
    auto shape = reconstruct(mesh, patches, ctx);
    Diagnostics info;
    check(inspect(shape, info), "periodic spline cylinder must close");
    check(info.faces == 3 && info.surface_counts["bspline"] == 1,
          "periodic charts must merge into one B-spline support face");
    check(std::abs(info.volume - pi * radius * radius * height) < .002,
          "periodic spline cylinder volume");
    auto bytes = write_step(shape, ctx.options.tolerance);
    check(verify_step(bytes, shape, ctx.options.tolerance), "periodic spline STEP round trip");

    // Two editable faces deliberately share one spline support. A cleanup
    // copy must retain that identity; copying each face's geometry separately
    // silently leaves the artificial partition in the final solid.
    BRepBuilderAPI_Sewing sewing(1e-5);
    for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next()) {
        auto face = TopoDS::Face(ex.Current());
        if (BRepAdaptor_Surface(face).GetType() != GeomAbs_BSplineSurface) {
            sewing.Add(face);
            continue;
        }
        auto support = BRep_Tool::Surface(face);
        double u0, u1, v0, v1;
        BRepTools::UVBounds(face, u0, u1, v0, v1);
        double middle = (v0 + v1) / 2;
        for (auto range : {std::pair<double, double>{v0, middle}, {middle, v1}}) {
            auto piece =
                BRepBuilderAPI_MakeFace(support, u0, u1, range.first, range.second, 1e-7).Face();
            piece.Orientation(face.Orientation());
            sewing.Add(piece);
        }
    }
    sewing.Perform();
    TopExp_Explorer shells(sewing.SewedShape(), TopAbs_SHELL);
    check(shells.More(), "split spline cylinder shell");
    auto split = BRepBuilderAPI_MakeSolid(TopoDS::Shell(shells.Current())).Solid();
    check(BRepLib::OrientClosedSolid(split) && inspect(split, info) && info.faces == 4,
          "split spline cylinder must start with four valid faces");
    auto tidied = tidy_shape(split, mesh, patches, ctx.options.tolerance);
    check(inspect(tidied, info) && info.faces == 3,
          "shared spline support must survive cleanup copying and merge");
}
} // namespace
int main() {
    unsigned failures = 0;
    for (bool recovery : {false, true})
        for (bool stepped : {false, true})
            for (double start : {.7, .1, pi, 5.7}) {
                try {
                    cylinder_boundaries(start, stepped, recovery);
                } catch (const std::exception &error) {
                    std::cerr << "start=" << start << " stepped=" << stepped
                              << " recovery=" << recovery << ": " << error.what() << '\n';
                    ++failures;
                }
            }
    if (failures)
        return 1;
    try {
        cylinder_boundaries(.7, false, true, true);
        cylinder_boundaries(5.7, true, true, true);
        cylinder_boundaries(.7, false, false, false, true);
        cylinder_boundaries(.7, false, false, false, false, true);
        cylinder_boundaries(5.7, true, true, false, false, true);
        freeform_solid();
        planar_strip_roof();
        periodic_spline_boundaries();
    } catch (const std::exception &error) {
        std::cerr << "boundary recovery and freeform reconstruction: " << error.what() << '\n';
        return 1;
    }
    std::cout << "Analytic boundary reconstruction and STEP tests passed\n";
    return 0;
}

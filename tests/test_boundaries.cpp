#include "internal.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Plane.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <STEPControl_Reader.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <gp_Ax3.hxx>
#include <iostream>
#include <sstream>

using namespace stl2curves;
using namespace stl2curves::detail;
namespace {
constexpr double pi = 3.14159265358979323846;
void check(bool ok, const char *message) {
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
void cylinder_boundaries(double start, bool stepped) {
    constexpr unsigned segments = 36, ring = segments + 1;
    constexpr double radius = 10, high = 3.5, low = 2;
    Mesh mesh;
    const std::vector<double> heights = stepped ? std::vector<double>{0, low, high}
                                                : std::vector<double>{0, high};
    for (double z : heights)
        for (unsigned i = 0; i <= segments; ++i) {
            double angle = start + pi / 2 * i / segments;
            mesh.points.emplace_back(radius * std::cos(angle), radius * std::sin(angle), z);
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
    Context ctx;
    auto shape = reconstruct(mesh, patches, ctx);
    Diagnostics info;
    check(inspect(shape, info), "reconstructed cylinder must be a valid closed solid");
    check(info.faces == (stepped ? 7u : 5u), "cylinder face count");
    double expected_volume = radius * radius * pi / 4 * (stepped ? (high + low) / 2 : high);
    check(std::abs(info.volume - expected_volume) < 1e-3, "analytic cylinder volume");
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
} // namespace
int main() {
    unsigned failures = 0;
    for (bool stepped : {false, true})
        for (double start : {.7, .1, pi, 5.7}) {
            try {
                cylinder_boundaries(start, stepped);
            } catch (const std::exception &error) {
                std::cerr << "start=" << start << " stepped=" << stepped << ": " << error.what()
                          << '\n';
                ++failures;
            }
        }
    if (failures)
        return 1;
    std::cout << "Analytic boundary reconstruction and STEP tests passed\n";
    return 0;
}

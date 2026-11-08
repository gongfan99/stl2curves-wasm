#include "internal.hpp"
#include <Eigen/Geometry>
#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <fstream>
#include <iostream>
#include <set>

using namespace stl2curves::detail;
namespace {
Patch read_patch(std::istream &in) {
    Patch p;
    unsigned count;
    in >> p.kind >> count >> p.origin.x() >> p.origin.y() >> p.origin.z() >> p.axis.x() >>
        p.axis.y() >> p.axis.z() >> p.radius >> p.secondary;
    p.triangles.resize(count);
    for (auto &t : p.triangles)
        in >> t;
    gp_Ax3 frame(gp_Pnt(p.origin.x(), p.origin.y(), p.origin.z()),
                 gp_Dir(p.axis.x(), p.axis.y(), p.axis.z()));
    if (p.kind == "cylinder")
        p.surface = new Geom_CylindricalSurface(frame, p.radius);
    else if (p.kind == "cone")
        p.surface = new Geom_ConicalSurface(frame, std::atan(p.secondary), p.radius);
    else if (p.kind == "sphere")
        p.surface = new Geom_SphericalSurface(frame, p.radius);
    else if (p.kind == "torus")
        p.surface = new Geom_ToroidalSurface(frame, p.radius, p.secondary);
    else
        throw std::runtime_error("Unsupported oracle surface " + p.kind);
    return p;
}
} // namespace
int main(int argc, char **argv) {
    if (argc != 2)
        return 2;
    try {
        std::ifstream in(argv[1]);
        if (!in)
            throw std::runtime_error("Cannot read sizing fixtures");
        unsigned count, failures = 0;
        in >> count;
        for (unsigned c = 0; c < count; ++c) {
            std::string name, unit, design_unit;
            unsigned np, nt, nf;
            double scale;
            in >> name >> unit >> np >> nt >> nf >> scale >> design_unit;
            Mesh m;
            m.points.resize(np);
            m.triangles.resize(nt);
            for (auto &p : m.points)
                in >> p.x() >> p.y() >> p.z();
            for (auto &t : m.triangles)
                in >> t[0] >> t[1] >> t[2];
            std::vector<Patch> patches, expected;
            for (unsigned i = 0; i < nf; ++i)
                patches.push_back(read_patch(in));
            for (unsigned i = 0; i < nf; ++i)
                expected.push_back(read_patch(in));
            if (!in)
                throw std::runtime_error("Truncated sizing fixture " + name);
            Context ctx;
            design_size(m, patches, ctx);
            bool ok = std::abs(ctx.diagnostics.scale - scale) < 1e-12 &&
                      ctx.diagnostics.design_unit == (design_unit == "-" ? "" : design_unit);
            snap_patches(m, patches, ctx, unit == "-" ? "" : unit);
            for (unsigned i = 0; i < nf; ++i) {
                const auto &a = patches[i], &b = expected[i];
                Vec offset = a.origin - b.origin;
                if (a.kind == "cylinder")
                    offset -= offset.dot(b.axis) * b.axis;
                ok = ok && a.kind == b.kind && std::abs(a.radius - b.radius) < 1e-8 &&
                     std::abs(a.secondary - b.secondary) < 1e-8 && offset.norm() < 1e-8 &&
                     std::abs(a.axis.dot(b.axis)) > 1 - 1e-10;
            }
            std::cout << (ok ? "PASS " : "FAIL ") << name << '\n';
            failures += !ok;
        }
        return failures ? 1 : 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

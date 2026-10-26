#include "internal.hpp"
#include <Eigen/Geometry>
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
} // namespace
int main() {
    try {
        thread_test(1, false);
        thread_test(-1, true);
        pipe_test();
        extrusion_test();
        std::cout << "Advanced fitting tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

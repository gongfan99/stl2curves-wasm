#include "internal.hpp"
#include <Eigen/Geometry>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
using namespace stl2curves;
using namespace stl2curves::detail;
int main(int argc, char **argv) {
    try {
        if (argc < 2 || argc > 3)
            throw std::runtime_error("Usage: inspect_detection FILE.stl [ACCURACY]");
        std::ifstream in(argv[1], std::ios::binary);
        if (!in)
            throw std::runtime_error("Cannot read STL");
        Bytes bytes((std::istreambuf_iterator<char>(in)), {});
        auto mesh = read_stl(bytes);
        Context context;
        if (argc == 3)
            context.options.accuracy = std::stod(argv[2]);
        repair_mesh(mesh, context);
        auto patches = detect(mesh, context);
        auto top = topology(mesh);
        std::vector<unsigned> uses(mesh.triangles.size());
        std::map<std::string, unsigned> counts;
        std::cout.precision(9);
        for (const auto &p : patches) {
            ++counts[p.kind];
            for (auto i : p.triangles)
                ++uses.at(i);
            if (p.kind != "plane" || std::getenv("STL2CURVES_INSPECT_PLANES")) {
                Vec lo = mesh.points[mesh.triangles[p.triangles.front()][0]], hi = lo;
                for (auto v : vertices(mesh, p.triangles)) {
                    lo = lo.cwiseMin(mesh.points[v]);
                    hi = hi.cwiseMax(mesh.points[v]);
                }
                std::cout << p.kind << " patch=" << (&p - patches.data())
                          << " triangles=" << p.triangles.size() << " r=" << p.radius
                          << " secondary=" << p.secondary << " dev=" << p.deviation
                          << " origin=" << p.origin.transpose() << " axis=" << p.axis.transpose()
                          << " bounds=" << lo.transpose() << " / " << hi.transpose() << '\n';
                if (p.kind == "cylinder") {
                    Vec e = p.axis.unitOrthogonal(), f = p.axis.cross(e);
                    double minimum = 1e100;
                    unsigned folded = 0;
                    for (auto ti : p.triangles) {
                        auto tri = mesh.triangles[ti];
                        double u[3], z[3];
                        for (unsigned k = 0; k < 3; ++k) {
                            Vec q = mesh.points[tri[k]] - p.origin;
                            u[k] = std::atan2(q.dot(f), q.dot(e));
                            z[k] = q.dot(p.axis);
                        }
                        for (unsigned k = 1; k < 3; ++k)
                            u[k] += std::round((u[0] - u[k]) / (2 * 3.141592653589793)) * 2 *
                                    3.141592653589793;
                        Vec radial = top.centers[ti] - p.origin;
                        radial -= radial.dot(p.axis) * p.axis;
                        double sign = radial.dot(top.normals[ti]) >= 0 ? 1 : -1;
                        double ratio =
                            p.radius * sign *
                            ((u[1] - u[0]) * (z[2] - z[0]) - (u[2] - u[0]) * (z[1] - z[0])) /
                            (2 * top.areas[ti]);
                        minimum = std::min(minimum, ratio);
                        if (ratio < .25)
                            ++folded;
                    }
                    std::cout << "  projected-area ratio min=" << minimum << " folded=" << folded
                              << '\n';
                }
            }
        }
        for (auto n : uses)
            if (n != 1)
                throw std::runtime_error("Detection does not partition the mesh");
        for (const auto &item : counts)
            std::cout << item.first << ": " << item.second << '\n';
        std::cout << "triangles " << mesh.triangles.size() << " patches " << patches.size() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

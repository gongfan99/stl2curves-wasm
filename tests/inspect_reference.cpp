#include "internal.hpp"
#include "projection.hpp"
#include <BRepAdaptor_Surface.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Eigen/Geometry>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <fstream>
#include <gp_Pnt2d.hxx>
#include <iostream>
#include <map>
#include <memory>
using namespace stl2curves;
using namespace stl2curves::detail;

double triangle_distance(const Vec &p, const Mesh &mesh, const Topology &top, unsigned i) {
    auto t = mesh.triangles[i];
    const auto &a = mesh.points[t[0]];
    const auto &b = mesh.points[t[1]];
    const auto &c = mesh.points[t[2]];
    const auto &n = top.normals[i];
    double d = (p - a).dot(n);
    Vec projected = p - d * n;
    if ((b - a).cross(projected - a).dot(n) >= -1e-12 &&
        (c - b).cross(projected - b).dot(n) >= -1e-12 &&
        (a - c).cross(projected - c).dot(n) >= -1e-12)
        return d * d;
    double best = 1e100;
    for (unsigned k = 0; k < 3; ++k) {
        Vec start = mesh.points[t[k]], edge = mesh.points[t[(k + 1) % 3]] - start;
        double s = std::clamp((p - start).dot(edge) / edge.squaredNorm(), 0., 1.);
        best = std::min(best, (p - start - s * edge).squaredNorm());
    }
    return best;
}
int main(int argc, char **argv) {
    try {
        if (argc != 3)
            throw std::runtime_error("Usage: inspect_reference MESH.stl REFERENCE.brep");
        std::ifstream file(argv[1], std::ios::binary);
        if (!file)
            throw std::runtime_error("Cannot read STL");
        Bytes bytes((std::istreambuf_iterator<char>(file)), {});
        auto mesh = read_stl(bytes);
        Context ctx;
        repair_mesh(mesh, ctx);
        auto patches = detect(mesh, ctx);
        auto top = topology(mesh);
        std::vector<unsigned> owner(mesh.triangles.size());
        std::vector<std::unique_ptr<PointProjector>> projections;
        for (unsigned i = 0; i < patches.size(); ++i) {
            for (auto t : patches[i].triangles)
                owner[t] = i;
            projections.push_back(std::make_unique<PointProjector>(patches[i].surface));
        }
        TopoDS_Shape shape;
        BRep_Builder builder;
        if (!BRepTools::Read(shape, argv[2], builder))
            throw std::runtime_error("Cannot read reference BREP");
        auto built = reconstruct(mesh, patches, ctx);
        for (bool original : {true, false}) {
            std::cout << (original ? "REFERENCE" : "REBUILT") << " face contributions\n";
            unsigned fi = 0;
            for (TopExp_Explorer ex(original ? shape : built, TopAbs_FACE); ex.More(); ex.Next()) {
                GProp_GProps mass, area;
                BRepGProp::VolumeProperties(ex.Current(), mass, 1e-7);
                BRepGProp::SurfaceProperties(ex.Current(), area);
                auto center = area.CentreOfMass();
                std::cout << ++fi << " type="
                          << int(BRepAdaptor_Surface(TopoDS::Face(ex.Current())).GetType())
                          << " area=" << area.Mass() << " volume=" << mass.Mass()
                          << " center=" << center.X() << "," << center.Y() << "," << center.Z()
                          << '\n';
            }
        }
        unsigned face_number = 0;
        for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next()) {
            ++face_number;
            auto face = TopoDS::Face(ex.Current());
            BRepAdaptor_Surface surface(face);
            double ulo, uhi, vlo, vhi;
            BRepTools::UVBounds(face, ulo, uhi, vlo, vhi);
            double worst = 0, squares = 0, signed_sum = 0;
            unsigned count = 0, worst_owner = 0;
            Vec worst_at = Vec::Zero();
            std::map<std::string, unsigned> owners;
            for (unsigned i = 0; i < 11; ++i)
                for (unsigned j = 0; j < 11; ++j) {
                    double u = ulo + (uhi - ulo) * (i + .5) / 11,
                           v = vlo + (vhi - vlo) * (j + .5) / 11;
                    BRepClass_FaceClassifier classifier(face, gp_Pnt2d(u, v), 1e-7);
                    if (classifier.State() != TopAbs_IN && classifier.State() != TopAbs_ON)
                        continue;
                    gp_Pnt at;
                    gp_Vec du, dv;
                    surface.D1(u, v, at, du, dv);
                    Vec p(at.X(), at.Y(), at.Z());
                    double nearest = 1e100;
                    unsigned triangle = 0;
                    for (unsigned t = 0; t < mesh.triangles.size(); ++t) {
                        double d = triangle_distance(p, mesh, top, t);
                        if (d < nearest) {
                            nearest = d;
                            triangle = t;
                        }
                    }
                    unsigned patch = owner[triangle];
                    ++owners[patches[patch].kind];
                    const auto &projection = projections[patch]->perform(p);
                    if (!projection.NbPoints())
                        continue;
                    double gap = projection.LowerDistance();
                    squares += gap * gap;
                    ++count;
                    auto q = projection.NearestPoint();
                    auto n = du.Crossed(dv).Normalized();
                    if (face.Orientation() == TopAbs_REVERSED)
                        n.Reverse();
                    signed_sum +=
                        (q.X() - p.x()) * n.X() + (q.Y() - p.y()) * n.Y() + (q.Z() - p.z()) * n.Z();
                    if (gap > worst) {
                        worst = gap;
                        worst_at = p;
                        worst_owner = patch;
                    }
                }
            std::cout << "face " << face_number << " type=" << int(surface.GetType())
                      << " samples=" << count << " rms=" << std::sqrt(squares / std::max(1u, count))
                      << " signed=" << signed_sum / std::max(1u, count) << " max=" << worst
                      << " at=" << worst_at.transpose() << " owner=" << patches[worst_owner].kind
                      << " vertex-deviation=" << patches[worst_owner].deviation << '\n';
            for (const auto &entry : owners)
                std::cout << "  " << entry.first << ": " << entry.second << '\n';
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

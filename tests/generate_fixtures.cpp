// C++ equivalents of upstream tests/make_*.py at 9f71aaa16eb689e759435d2e1fbf2ba7468cadfb.
// Copyright (c) 2026 Duvdivan. MIT license.
#include "occt_collections.hpp"
#include "occt_error.hpp"
using stl2curves::detail::ShapeMap;
#include <BRepAdaptor_Curve.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <StlAPI_Writer.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <gp_Ax1.hxx>
#include <gp_Trsf.hxx>
#include <iostream>
#include <stdexcept>
#include <vector>

using Shape = TopoDS_Shape;
using Predicate = std::function<bool(const BRepAdaptor_Curve &, const gp_Pnt &, const gp_Pnt &)>;
constexpr double pi = 3.14159265358979323846;
Shape cut(const Shape &a, const Shape &b) {
    return BRepAlgoAPI_Cut(a, b).Shape();
}
Shape fuse(const Shape &a, const Shape &b) {
    return BRepAlgoAPI_Fuse(a, b).Shape();
}
Shape box(double x, double y, double z, double a, double b, double c) {
    return BRepPrimAPI_MakeBox(gp_Pnt(x, y, z), a, b, c).Shape();
}
Shape cylinder(double x, double y, double z, double radius, double height,
               gp_Dir d = gp_Dir(0, 0, 1)) {
    return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x, y, z), d), radius, height).Shape();
}
Shape cone(double x, double y, double z, double r1, double r2, double height) {
    return BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(x, y, z), gp_Dir(0, 0, 1)), r1, r2, height).Shape();
}
Shape sphere(double x, double y, double z, double radius) {
    return BRepPrimAPI_MakeSphere(gp_Pnt(x, y, z), radius).Shape();
}
std::vector<TopoDS_Edge> edges(const Shape &s, const Predicate &predicate) {
    ShapeMap map;
    TopExp::MapShapes(s, TopAbs_EDGE, map);
    std::vector<TopoDS_Edge> out;
    for (int i = 1; i <= map.Extent(); ++i) {
        auto e = TopoDS::Edge(map(i));
        BRepAdaptor_Curve c(e);
        if (predicate(c, c.Value(c.FirstParameter()), c.Value(c.LastParameter())))
            out.push_back(e);
    }
    return out;
}
Predicate at(double z) {
    return [=](const BRepAdaptor_Curve &, const gp_Pnt &a, const gp_Pnt &b) {
        return std::abs(a.Z() - z) < 1e-6 && std::abs(b.Z() - z) < 1e-6;
    };
}
bool upright(const BRepAdaptor_Curve &c, const gp_Pnt &a, const gp_Pnt &b) {
    return c.GetType() == GeomAbs_Line && std::abs(a.Z() - b.Z()) > 1 &&
           std::abs(a.X() - b.X()) < 1e-6 && std::abs(a.Y() - b.Y()) < 1e-6;
}
Shape fillet(const Shape &s, double radius, const Predicate &predicate) {
    BRepFilletAPI_MakeFillet f(s);
    for (const auto &e : edges(s, predicate))
        f.Add(radius, e);
    f.Build();
    if (!f.IsDone())
        throw std::runtime_error("fixture fillet failed");
    return f.Shape();
}
void write(const Shape &s, const std::filesystem::path &directory, const char *name,
           double angle = .25) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(s, props);
    std::cout.precision(12);
    GProp_GProps precise;
    BRepGProp::VolumeProperties(s, precise, 1e-7);
    std::cout << name << " design volume " << props.Mass() << " adaptive " << precise.Mass()
              << '\n';
    std::ofstream reference(directory / (std::string(name) + ".brep"));
    BRepTools::Write(s, reference);
    if (!reference)
        throw std::runtime_error("reference BREP write failed");
    BRepTools::Clean(s);
    BRepMesh_IncrementalMesh mesh(s, .01, false, angle, true);
    StlAPI_Writer writer;
    writer.ASCIIMode() = false;
    if (!writer.Write(s, (directory / name).string().c_str()))
        throw std::runtime_error("fixture write failed");
}
void plate(const std::filesystem::path &d) {
    auto s = fillet(box(0, 0, 0, 60, 40, 15), 4, upright);
    s = cut(s, cylinder(10, 10, -1, 2.5, 17));
    s = cut(s, cylinder(20, 10, 7, 1.6, 9));
    s = cut(s, cylinder(40, 10, -1, 1.7, 17));
    s = cut(s, cylinder(40, 10, 12, 3, 4));
    s = cut(s, cylinder(-1, 30, 7.5, 1.5, 62, gp_Dir(1, 0, 0)));
    BRepBuilderAPI_MakePolygon polygon;
    for (int k = 0; k < 6; ++k)
        polygon.Add(gp_Pnt(50 + 3.2 * std::cos(k * pi / 3), 28 + 3.2 * std::sin(k * pi / 3), 11));
    polygon.Close();
    auto hex = BRepBuilderAPI_MakeFace(polygon.Wire()).Face();
    s = cut(s, BRepPrimAPI_MakePrism(hex, gp_Vec(0, 0, 5)).Shape());
    s = cut(s, cylinder(30, 25, -1, 2, 17));
    s = fuse(s, cylinder(15, 30, 15, 2, 6));
    s = fuse(s, cylinder(28, 33, 15, 1.2, 4));
    write(s, d, "test_fine.stl", .2);
    write(s, d, "test_coarse.stl", .5);
}
void hard(const std::filesystem::path &d) {
    auto s = box(0, 0, 0, 80, 40, 12);
    s = cut(s, cylinder(10, 10, -1, 2, 14));
    s = cut(s, cone(10, 10, -.001, 2.6, 2, .601));
    s = cut(s, cone(10, 10, 11.4, 2, 2.6, .601));
    s = cut(s, cylinder(30, -1, 6, 1.5, 42, gp_Dir(0, 1, 0)));
    s = cut(s, cylinder(25, 20, 6, 1.5, 10, gp_Dir(1, 0, 0)));
    s = cut(s, cylinder(50, 10, -1, 2, 14));
    s = cut(s, cylinder(56, 10, -1, 2, 14));
    s = cut(s, box(50, 8, -1, 6, 4, 14));
    s = fuse(s, cylinder(15, 30, 12, 2.5, 5));
    s = fuse(s, cone(15, 30, 17, 2.5, 1.9, .6));
    s = fuse(s, cylinder(35, 30, 11, 2, 7));
    s = fillet(s, .8, [](const BRepAdaptor_Curve &c, const gp_Pnt &, const gp_Pnt &) {
        return c.GetType() == GeomAbs_Circle && std::abs(c.Circle().Radius() - 2) < 1e-6 &&
               std::abs(c.Circle().Location().Z() - 12) < 1e-6 &&
               std::abs(c.Circle().Location().X() - 35) < 1e-6;
    });
    auto wedge = box(60, 20, 0, 30, 30, 30);
    gp_Trsf tr;
    tr.SetRotation(gp_Ax1(gp_Pnt(60, 20, 12), gp_Dir(0, 1, 0)), -pi / 6);
    wedge = BRepBuilderAPI_Transform(wedge, tr).Shape();
    s = cut(s, wedge);
    s = cut(s, cylinder(66, 30, -1, 1.8, 20));
    s = cut(s, cylinder(-1, 34, 5, 1.25, 20, gp_Dir(1, 0, 0)));
    write(s, d, "test_hard.stl", .3);
}
void lid(const std::filesystem::path &d) {
    auto s = fillet(box(0, 0, 0, 80, 50, 20), 6, upright);
    auto pocket = fillet(box(10, 10, 8, 60, 30, 20), 4, upright);
    s = cut(s, pocket);
    s = fillet(s, 1.5, at(8));
    s = fillet(s, 2, [](const BRepAdaptor_Curve &c, const gp_Pnt &a, const gp_Pnt &b) {
        if (std::abs(a.Z() - 20) > 1e-6 || std::abs(b.Z() - 20) > 1e-6)
            return false;
        auto m = c.Value((c.FirstParameter() + c.LastParameter()) / 2);
        return m.X() < 5 || m.X() > 75 || m.Y() < 5 || m.Y() > 45;
    });
    BRepFilletAPI_MakeChamfer ch(s);
    for (const auto &e : edges(s, at(0)))
        ch.Add(1, e);
    ch.Build();
    if (!ch.IsDone())
        throw std::runtime_error("fixture chamfer failed");
    s = ch.Shape();
    s = cut(s, cylinder(40, 25, -1, 1.7, 12));
    s = cut(s, cone(40, 25, -.001, 3.4, 1.7, 1.701));
    write(s, d, "test_lid.stl");
}
void spheres(const std::filesystem::path &d) {
    auto s = fillet(box(0, 0, 0, 70, 40, 12), 2,
                    [](const BRepAdaptor_Curve &, const gp_Pnt &, const gp_Pnt &) { return true; });
    s = fuse(s, cut(sphere(12, 12, 12, 5), box(0, 0, 0, 30, 30, 12)));
    s = cut(s, sphere(30, 12, 14.5, 4));
    s = fuse(s, cylinder(45, 12, 12, 2, 6));
    s = fuse(s, sphere(45, 12, 18, 2));
    s = cut(s, cylinder(58, 12, 5, 1.5, 8));
    s = cut(s, sphere(58, 12, 5, 1.5));
    s = cut(s, sphere(20, 30, 6, 2.5));
    TopoDS_Compound comp;
    BRep_Builder b;
    b.MakeCompound(comp);
    b.Add(comp, s);
    b.Add(comp, sphere(90, 20, 6, 6));
    write(comp, d, "test_sphere.stl");
}
void corners(const std::filesystem::path &d) {
    auto s = fillet(box(0, 0, 0, 60, 40, 14), 3, upright);
    s = fillet(s, 1, at(14));
    s = fillet(s, 4, at(0));
    write(s, d, "test_corners.stl");
    s = box(0, 0, 0, 40, 30, 12);
    BRepFilletAPI_MakeFillet f(s);
    for (const auto &e : edges(s, at(12)))
        f.Add(1, e);
    for (const auto &e : edges(s, upright))
        f.Add(2, e);
    for (const auto &e : edges(s, at(0)))
        f.Add(3, e);
    f.Build();
    if (!f.IsDone())
        throw std::runtime_error("corner fixture failed");
    write(f.Shape(), d, "test_corners2.stl");
}
int main(int argc, char **argv) {
    try {
        if (argc != 2)
            throw std::runtime_error("Usage: generate_fixtures DIRECTORY");
        std::filesystem::path d(argv[1]);
        std::filesystem::create_directories(d);
        plate(d);
        hard(d);
        lid(d);
        spheres(d);
        corners(d);
        return 0;
    } catch (const Standard_Failure &e) {
        std::cerr << stl2curves::detail::failure_message(e) << '\n';
        return 1;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

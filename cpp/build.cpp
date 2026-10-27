#include "internal.hpp"
#include "occt_error.hpp"
#include "projection.hpp"
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Eigen/Geometry>
#include <GProp_GProps.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <GeomAPI_IntSS.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomProjLib.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <ShapeFix_Edge.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <cstdlib>
#include <gp_Pnt2d.hxx>
#include <gp_Vec2d.hxx>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <tuple>

namespace stl2curves::detail {
namespace {
gp_Pnt point(const Vec &v) {
    return gp_Pnt(v.x(), v.y(), v.z());
}
TopoDS_Shape compound(const std::vector<TopoDS_Shape> &items) {
    TopoDS_Compound c;
    BRep_Builder b;
    b.MakeCompound(c);
    for (const auto &s : items)
        if (!s.IsNull())
            b.Add(c, s);
    return c;
}
double volume(const TopoDS_Shape &s) {
    // Integrate independent solids separately, preserving signed cavity shells.
    double total = 0;
    bool found = false;
    for (TopExp_Explorer ex(s, TopAbs_SOLID); ex.More(); ex.Next()) {
        GProp_GProps props;
        BRepGProp::VolumeProperties(ex.Current(), props, volume_epsilon);
        total += props.Mass();
        found = true;
    }
    if (found)
        return total;
    GProp_GProps props;
    BRepGProp::VolumeProperties(s, props, volume_epsilon);
    return props.Mass();
}
void trace_invalid_faces(const TopoDS_Shape &shape, const std::vector<Patch> *patches = nullptr) {
#ifndef __EMSCRIPTEN__
    if (!std::getenv("STL2CURVES_TRACE"))
        return;
    unsigned number = 0, reported = 0;
    for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More() && reported < 16; ex.Next()) {
        ++number;
        auto face = TopoDS::Face(ex.Current());
        BRepCheck_Analyzer check(face);
        if (check.IsValid())
            continue;
        ++reported;
        std::string message = "invalid face=" + std::to_string(number) + " surface=" +
                              std::to_string(int(BRepAdaptor_Surface(face).GetType()));
        auto support = BRep_Tool::Surface(face);
        if (patches)
            for (unsigned i = 0; i < patches->size(); ++i)
                if (support == (*patches)[i].surface)
                    message += " patch=" + std::to_string(i) + ":" + (*patches)[i].kind;
        for (auto kind : {TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE})
            for (TopExp_Explorer child(face, kind); child.More(); child.Next()) {
                auto result = check.Result(child.Current());
                if (result.IsNull())
                    continue;
                for (auto status : result->Status())
                    if (status != BRepCheck_NoError)
                        message += " status=" + std::to_string(int(status));
                for (result->InitContextIterator(); result->MoreShapeInContext();
                     result->NextShapeInContext())
                    for (auto status : result->StatusOnShape())
                        if (status != BRepCheck_NoError)
                            message += " context-status=" + std::to_string(int(status));
            }
        trace(message.c_str());
    }
#else
    (void)shape;
    (void)patches;
#endif
}
TopoDS_Shape close_shells(const TopoDS_Shape &sewn, double tolerance) {
    std::vector<TopoDS_Solid> solids;
    std::vector<double> original_volumes;
    std::vector<TopoDS_Shell> candidates;
    for (TopExp_Explorer ex(sewn, TopAbs_SHELL); ex.More(); ex.Next())
        candidates.push_back(TopoDS::Shell(ex.Current()));
    // Sewing may leave an already closed one-face sphere outside any shell.
    // Include such faces so isolated bodies and spherical cavities survive.
    for (TopExp_Explorer ex(sewn, TopAbs_FACE, TopAbs_SHELL); ex.More(); ex.Next()) {
        TopoDS_Shell shell;
        BRep_Builder builder;
        builder.MakeShell(shell);
        builder.Add(shell, ex.Current());
        candidates.push_back(shell);
    }
    for (auto shell : candidates) {
        auto solid = BRepBuilderAPI_MakeSolid(shell).Solid();
        original_volumes.push_back(volume(solid));
        if (!BRepLib::OrientClosedSolid(solid))
            throw std::runtime_error("Sewn shell is not closed");
        if (volume(solid) < 0)
            solid.Reverse();
        solids.push_back(solid);
    }
    if (solids.empty())
        throw std::runtime_error("Reconstruction produced no closed shells");
    // Nest shells geometrically. Odd nesting levels are cavities; even levels are
    // independent material regions, including islands inside a cavity.
    std::vector<int> parent(solids.size(), -1);
    std::vector<double> volumes;
    for (auto v : original_volumes)
        volumes.push_back(std::abs(v));
    for (unsigned i = 0; i < solids.size(); ++i) {
        TopExp_Explorer vertices(solids[i], TopAbs_VERTEX);
        if (!vertices.More())
            continue;
        gp_Pnt sample = BRep_Tool::Pnt(TopoDS::Vertex(vertices.Current()));
        double smallest = std::numeric_limits<double>::infinity();
        for (unsigned j = 0; j < solids.size(); ++j)
            if (i != j && volumes[j] > volumes[i] && volumes[j] < smallest &&
                original_volumes[i] * original_volumes[j] < 0) {
                BRepClass3d_SolidClassifier classifier(solids[j], sample, tolerance);
                if (classifier.State() == TopAbs_IN) {
                    bool contained = true;
                    for (TopExp_Explorer v(solids[i], TopAbs_VERTEX); v.More(); v.Next()) {
                        classifier.Perform(BRep_Tool::Pnt(TopoDS::Vertex(v.Current())), tolerance);
                        if (classifier.State() == TopAbs_OUT ||
                            classifier.State() == TopAbs_UNKNOWN) {
                            contained = false;
                            break;
                        }
                    }
                    if (contained) {
                        parent[i] = int(j);
                        smallest = volumes[j];
                    }
                }
            }
    }
    auto depth = [&](unsigned i) {
        unsigned d = 0;
        for (int p = parent[i]; p >= 0; p = parent[p])
            ++d;
        return d;
    };
    std::vector<TopoDS_Shape> out;
    for (unsigned i = 0; i < solids.size(); ++i)
        if (depth(i) % 2 == 0) {
            BRepBuilderAPI_MakeSolid make(solids[i]);
            for (unsigned j = 0; j < solids.size(); ++j)
                if (parent[j] == int(i)) {
                    TopExp_Explorer ex(solids[j], TopAbs_SHELL);
                    auto inner = TopoDS::Shell(ex.Current());
                    if (volume(BRepBuilderAPI_MakeSolid(inner).Solid()) > 0)
                        inner.Reverse();
                    make.Add(inner);
                }
            out.push_back(make.Solid());
        }
    return out.size() == 1 ? out.front() : compound(out);
}
using Boundaries = std::map<Edge, TopoDS_Edge>;
struct Rim {
    unsigned from, to, count = 0;
};
std::map<Edge, Rim> patch_rim(const Mesh &mesh, const std::vector<unsigned> &ids) {
    std::map<Edge, Rim> rim;
    for (auto i : ids) {
        auto t = mesh.triangles[i];
        for (unsigned k = 0; k < 3; ++k) {
            auto &r = rim[Edge(t[k], t[(k + 1) % 3])];
            r.from = t[k];
            r.to = t[(k + 1) % 3];
            ++r.count;
        }
    }
    return rim;
}
std::vector<TopoDS_Face> planar_patch(const Mesh &mesh, const Patch &patch,
                                      const Boundaries &boundaries) {
    auto top = topology(mesh);
    std::set<unsigned> available(patch.triangles.begin(), patch.triangles.end());
    std::vector<TopoDS_Face> faces;
    while (!available.empty()) {
        std::vector<unsigned> group{*available.begin()};
        available.erase(group[0]);
        for (std::size_t q = 0; q < group.size(); ++q)
            for (auto j : top.neighbors[group[q]])
                if (available.erase(j))
                    group.push_back(j);
        auto rim = patch_rim(mesh, group);
        std::map<unsigned, unsigned> next;
        for (const auto &kv : rim)
            if (kv.second.count == 1)
                if (!next.emplace(kv.second.from, kv.second.to).second)
                    throw std::runtime_error("Branching planar boundary");
        struct Loop {
            TopoDS_Wire wire;
            double area;
        };
        std::vector<Loop> loops;
        Vec e = patch.axis.unitOrthogonal(), f = patch.axis.cross(e);
        while (!next.empty()) {
            unsigned start = next.begin()->first, current = start;
            BRepBuilderAPI_MakeWire wire;
            double area = 0;
            do {
                auto found = next.find(current);
                if (found == next.end())
                    throw std::runtime_error("Open planar patch boundary");
                unsigned to = found->second;
                next.erase(found);
                const auto &a = mesh.points[current];
                const auto &b = mesh.points[to];
                area += a.dot(e) * b.dot(f) - b.dot(e) * a.dot(f);
                auto boundary = boundaries.find(Edge(current, to));
                if (boundary == boundaries.end())
                    wire.Add(BRepBuilderAPI_MakeEdge(point(a), point(b)).Edge());
                else {
                    auto edge = boundary->second;
                    if (current > to)
                        edge.Reverse();
                    wire.Add(edge);
                }
                current = to;
            } while (current != start);
            if (!wire.IsDone())
                throw std::runtime_error("Cannot connect planar patch boundary");
            loops.push_back({wire.Wire(), area / 2});
        }
        if (loops.empty())
            throw std::runtime_error("Empty planar patch boundary");
        auto outer = std::max_element(loops.begin(), loops.end(), [](const Loop &a, const Loop &b) {
            return std::abs(a.area) < std::abs(b.area);
        });
        std::iter_swap(loops.begin(), outer);
        if (loops[0].area < 0)
            loops[0].wire.Reverse();
        BRepBuilderAPI_MakeFace make(patch.surface, loops[0].wire, true);
        for (std::size_t i = 1; i < loops.size(); ++i) {
            if (loops[i].area > 0)
                loops[i].wire.Reverse();
            make.Add(loops[i].wire);
        }
        if (!make.IsDone())
            throw std::runtime_error("Cannot trim planar patch");
        auto face = make.Face();
        if (patch.axis.dot(top.normals[group[0]]) < 0)
            face.Reverse();
        faces.push_back(face);
    }
    return faces;
}
TopoDS_Face natural_patch(const Mesh &mesh, const Patch &patch, double tolerance) {
    if (patch.kind != "cylinder" && patch.kind != "cone" && patch.kind != "sphere" &&
        patch.kind != "torus")
        return {};
    auto ids = vertices(mesh, patch.triangles);
    std::map<unsigned, gp_Pnt2d> uv;
    std::set<unsigned> poles;
    PointProjector projector(patch.surface);
    for (auto id : ids) {
        const auto &project = projector.perform(mesh.points[id]);
        if (!project.NbPoints())
            return {};
        double u, v;
        project.LowerDistanceParameters(u, v);
        uv.emplace(id, gp_Pnt2d(u, v));
        gp_Pnt p;
        gp_Vec du, dv;
        patch.surface->D1(u, v, p, du, dv);
        if (du.Magnitude() < std::max(1e-9, tolerance * .001))
            poles.insert(id);
    }
    auto unwrap = [&](bool u, double period) {
        std::vector<double> values;
        for (const auto &item : uv) {
            if (u && poles.count(item.first))
                continue;
            double t = u ? item.second.X() : item.second.Y();
            t = std::fmod(t, period);
            if (t < 0)
                t += period;
            values.push_back(t);
        }
        if (values.empty())
            return;
        std::sort(values.begin(), values.end());
        double gap = -1, cut = values.front();
        for (unsigned i = 0; i < values.size(); ++i) {
            double next = values[(i + 1) % values.size()] + (i + 1 == values.size() ? period : 0);
            if (next - values[i] > gap) {
                gap = next - values[i];
                cut = std::fmod(next, period);
            }
        }
        for (auto &item : uv) {
            double t = u ? item.second.X() : item.second.Y();
            t -= std::floor((t - cut) / period) * period;
            if (u)
                item.second.SetX(t);
            else
                item.second.SetY(t);
        }
    };
    if (patch.surface->IsUPeriodic())
        unwrap(true, patch.surface->UPeriod());
    if (patch.surface->IsVPeriodic())
        unwrap(false, patch.surface->VPeriod());
    double ulo = 1e100, uhi = -1e100, vlo = 1e100, vhi = -1e100;
    for (const auto &item : uv) {
        if (!poles.count(item.first)) {
            ulo = std::min(ulo, item.second.X());
            uhi = std::max(uhi, item.second.X());
        }
        vlo = std::min(vlo, item.second.Y());
        vhi = std::max(vhi, item.second.Y());
    }
    if (uhi - ulo < 1e-8 || vhi - vlo < 1e-8)
        return {};
    auto rim = patch_rim(mesh, patch.triangles);
    double eps = 1e-5;
    bool all_v = true, all_u = true;
    unsigned edges = 0;
    for (const auto &kv : rim)
        if (kv.second.count == 1) {
            ++edges;
            auto a = uv.at(kv.first.a), b = uv.at(kv.first.b);
            bool v = (std::abs(a.Y() - vlo) < eps && std::abs(b.Y() - vlo) < eps) ||
                     (std::abs(a.Y() - vhi) < eps && std::abs(b.Y() - vhi) < eps);
            bool u = ((poles.count(kv.first.a) || std::abs(a.X() - ulo) < eps) &&
                      (poles.count(kv.first.b) || std::abs(b.X() - ulo) < eps)) ||
                     ((poles.count(kv.first.a) || std::abs(a.X() - uhi) < eps) &&
                      (poles.count(kv.first.b) || std::abs(b.X() - uhi) < eps));
            if (!u && !v)
                return {};
            all_v = all_v && v;
            all_u = all_u && u;
        }
    if (!edges) {
        auto face = BRepBuilderAPI_MakeFace(patch.surface, tolerance).Face();
        auto t = mesh.triangles[patch.triangles.front()];
        Vec centroid = (mesh.points[t[0]] + mesh.points[t[1]] + mesh.points[t[2]]) / 3;
        GeomAPI_ProjectPointOnSurf projection(point(centroid), patch.surface);
        double u, v;
        projection.LowerDistanceParameters(u, v);
        gp_Pnt q;
        gp_Vec du, dv;
        patch.surface->D1(u, v, q, du, dv);
        gp_Vec normal = du.Crossed(dv);
        Vec original =
            (mesh.points[t[1]] - mesh.points[t[0]]).cross(mesh.points[t[2]] - mesh.points[t[0]]);
        if (original.dot(Vec(normal.X(), normal.Y(), normal.Z())) < 0)
            face.Reverse();
        return face;
    }
    if (all_v && patch.surface->IsUPeriodic()) {
        ulo = 0;
        uhi = patch.surface->UPeriod();
    }
    if (all_u && patch.surface->IsVPeriodic()) {
        vlo = 0;
        vhi = patch.surface->VPeriod();
    }
    BRepBuilderAPI_MakeFace make(patch.surface, ulo, uhi, vlo, vhi, tolerance);
    if (!make.IsDone())
        return {};
    auto face = make.Face();
    auto t = mesh.triangles[patch.triangles.front()];
    Vec centroid = (mesh.points[t[0]] + mesh.points[t[1]] + mesh.points[t[2]]) / 3;
    GeomAPI_ProjectPointOnSurf projection(point(centroid), patch.surface);
    double u, v;
    projection.LowerDistanceParameters(u, v);
    gp_Pnt value;
    gp_Vec du, dv;
    patch.surface->D1(u, v, value, du, dv);
    gp_Vec n = du.Crossed(dv);
    Vec original =
        (mesh.points[t[1]] - mesh.points[t[0]]).cross(mesh.points[t[2]] - mesh.points[t[0]]);
    if (original.dot(Vec(n.X(), n.Y(), n.Z())) < 0)
        face.Reverse();
    return face;
}
std::vector<TopoDS_Face> trimmed_patch(const Mesh &mesh, const Patch &patch, double tolerance,
                                       const Boundaries &boundaries) {
    // Cut the triangulated parameter chart along periodic seams, then remove
    // interior edges in 2D. This retains arbitrary outlines and holes while
    // producing one trimmed support face per connected chart.
    using Key = std::tuple<unsigned, int, int>;
    struct Vertex {
        gp_Pnt2d uv;
        unsigned original;
        TopoDS_Vertex shape;
    };
    std::map<unsigned, gp_Pnt2d> projected;
    auto ids = vertices(mesh, patch.triangles);
    PointProjector projector(patch.surface);
    for (auto id : ids) {
        const auto &p = projector.perform(mesh.points[id]);
        if (!p.NbPoints())
            return {};
        double u, v;
        p.LowerDistanceParameters(u, v);
        if (!std::isfinite(u) || !std::isfinite(v) || std::abs(u) > 1e5 || std::abs(v) > 1e5)
            return {};
        projected.emplace(id, gp_Pnt2d(u, v));
    }
    double pu = patch.surface->IsUPeriodic() ? patch.surface->UPeriod() : 0,
           pv = patch.surface->IsVPeriodic() ? patch.surface->VPeriod() : 0;
    std::map<Key, unsigned> lookup;
    std::vector<Vertex> chart;
    std::vector<Tri> triangles;
    for (auto ti : patch.triangles) {
        auto t = mesh.triangles[ti];
        gp_Pnt2d uv[3];
        int ku[3] = {0, 0, 0}, kv[3] = {0, 0, 0};
        for (unsigned k = 0; k < 3; ++k)
            uv[k] = projected.at(t[k]);
        auto unwrap = [&](double period, bool u, int *shifts) {
            if (!period)
                return;
            double values[3];
            for (unsigned k = 0; k < 3; ++k) {
                double x = u ? uv[k].X() : uv[k].Y();
                shifts[k] = -int(std::floor(x / period));
                values[k] = x + shifts[k] * period;
            }
            double lo = *std::min_element(values, values + 3),
                   hi = *std::max_element(values, values + 3);
            if (hi - lo > period / 2)
                for (unsigned k = 0; k < 3; ++k)
                    if (values[k] < period / 2) {
                        values[k] += period;
                        ++shifts[k];
                    }
            for (unsigned k = 0; k < 3; ++k)
                if (u)
                    uv[k].SetX(values[k]);
                else
                    uv[k].SetY(values[k]);
        };
        unwrap(pu, true, ku);
        unwrap(pv, false, kv);
        Tri mapped;
        for (unsigned k = 0; k < 3; ++k) {
            Key key{t[k], ku[k], kv[k]};
            auto found = lookup.emplace(key, unsigned(chart.size()));
            if (found.second) {
                auto v =
                    BRepBuilderAPI_MakeVertex(patch.surface->Value(uv[k].X(), uv[k].Y())).Vertex();
                BRep_Builder().UpdateVertex(v, std::max(1e-6, tolerance * .1));
                chart.push_back({uv[k], t[k], v});
            }
            mapped[k] = found.first->second;
        }
        triangles.push_back(mapped);
    }
    std::map<Edge, Rim> edges;
    std::vector<std::vector<unsigned>> neighbors(triangles.size());
    std::map<Edge, std::vector<unsigned>> uses;
    for (unsigned i = 0; i < triangles.size(); ++i) {
        auto t = triangles[i];
        for (unsigned k = 0; k < 3; ++k)
            uses[Edge(t[k], t[(k + 1) % 3])].push_back(i);
    }
    for (const auto &item : uses)
        if (item.second.size() == 2) {
            neighbors[item.second[0]].push_back(item.second[1]);
            neighbors[item.second[1]].push_back(item.second[0]);
        }
    std::vector<TopoDS_Face> faces;
    for (const auto &group : components(neighbors)) {
        edges.clear();
        for (auto i : group) {
            auto t = triangles[i];
            for (unsigned k = 0; k < 3; ++k) {
                auto &edge = edges[Edge(t[k], t[(k + 1) % 3])];
                edge.from = t[k];
                edge.to = t[(k + 1) % 3];
                ++edge.count;
            }
        }
        std::map<unsigned, unsigned> next;
        for (const auto &item : edges)
            if (item.second.count == 1)
                if (!next.emplace(item.second.from, item.second.to).second)
                    return {};
        struct Loop {
            TopoDS_Wire wire;
            double area;
        };
        std::vector<Loop> loops;
        while (!next.empty()) {
            unsigned start = next.begin()->first, current = start;
            TopoDS_Wire wire;
            BRep_Builder builder;
            builder.MakeWire(wire);
            double area = 0;
            do {
                auto found = next.find(current);
                if (found == next.end())
                    return {};
                unsigned to = found->second;
                next.erase(found);
                const auto &a = chart[current];
                const auto &b = chart[to];
                area += a.uv.X() * b.uv.Y() - b.uv.X() * a.uv.Y();
                gp_Vec2d delta(a.uv, b.uv);
                if (delta.Magnitude() < 1e-12)
                    return {};
                TopoDS_Edge edge;
                auto exact = boundaries.find(Edge(a.original, b.original));
                if (exact != boundaries.end()) {
                    double first, last;
                    auto curve = BRep_Tool::Curve(exact->second, first, last);
                    double error = tolerance * .01;
                    auto pc = GeomProjLib::Curve2d(curve, first, last, patch.surface, error);
                    if (!pc.IsNull()) {
                        auto middle = pc->Value((first + last) / 2);
                        double du = pu ? std::round(((a.uv.X() + b.uv.X()) / 2 - middle.X()) / pu) *
                                             pu
                                       : 0,
                               dv = pv ? std::round(((a.uv.Y() + b.uv.Y()) / 2 - middle.Y()) / pv) *
                                             pv
                                       : 0;
                        pc->Translate(gp_Vec2d(du, dv));
                        bool reverse = pc->Value(last).SquareDistance(a.uv) <
                                       pc->Value(first).SquareDistance(a.uv);
                        BRepBuilderAPI_MakeEdge make(pc, patch.surface, reverse ? b.shape : a.shape,
                                                     reverse ? a.shape : b.shape, first, last);
                        if (make.IsDone()) {
                            edge = make.Edge();
                            if (reverse)
                                edge.Reverse();
                        }
                    }
                }
                if (edge.IsNull()) {
                    Handle(Geom2d_Line) line = new Geom2d_Line(a.uv, gp_Dir2d(delta));
                    BRepBuilderAPI_MakeEdge make(line, patch.surface, a.shape, b.shape, 0,
                                                 delta.Magnitude());
                    if (!make.IsDone())
                        return {};
                    edge = make.Edge();
                }
                BRepLib::BuildCurve3d(edge, tolerance * .01);
                builder.Add(wire, edge);
                current = to;
            } while (current != start);
            if (std::abs(area) < 1e-16)
                return {};
            loops.push_back({wire, area / 2});
        }
        if (loops.empty())
            return {};
        auto outer = std::max_element(loops.begin(), loops.end(), [](const Loop &a, const Loop &b) {
            return std::abs(a.area) < std::abs(b.area);
        });
        std::iter_swap(loops.begin(), outer);
        bool reversed = loops[0].area < 0;
        if (reversed)
            loops[0].wire.Reverse();
        BRepBuilderAPI_MakeFace make(patch.surface, loops[0].wire, true);
        for (unsigned i = 1; i < loops.size(); ++i) {
            if (loops[i].area > 0)
                loops[i].wire.Reverse();
            make.Add(loops[i].wire);
        }
        if (!make.IsDone())
            return {};
        auto face = make.Face();
        if (reversed)
            face.Reverse();
        faces.push_back(face);
    }
    return faces;
}
Boundaries shared_boundaries(const Mesh &mesh, const std::vector<Patch> &patches, double tolerance,
                             bool intersect_analytic = true) {
    trace("surface intersections");
    std::vector<unsigned> owner(mesh.triangles.size());
    for (unsigned i = 0; i < patches.size(); ++i)
        for (auto t : patches[i].triangles)
            owner[t] = i;
    std::map<std::pair<unsigned, unsigned>, std::vector<Handle(Geom_Curve)>> intersections;
    Boundaries out;
    auto top = topology(mesh);
    auto analytic = [](const Patch &p) {
        return p.kind == "plane" || p.kind == "cylinder" || p.kind == "cone" ||
               p.kind == "sphere" || p.kind == "torus";
    };
    // Fit a complete junction, rather than projecting each mesh chord in
    // isolation. The interpolation parameters preserve every mesh endpoint,
    // so both incident faces can still use the same edge partition.
    std::map<Edge, Handle(Geom_Curve)> smooth_boundaries;
    std::map<std::pair<unsigned, unsigned>, std::set<Edge>> fitted_runs;
    for (const auto &item : top.edges)
        if (item.second.size() == 2) {
            auto a = owner[item.second[0]], b = owner[item.second[1]];
            if (a != b && (!analytic(patches[a]) || !analytic(patches[b])))
                fitted_runs[std::minmax(a, b)].insert(item.first);
        }
    for (auto &item : fitted_runs) {
        auto &remaining = item.second;
        std::map<unsigned, std::vector<unsigned>> adjacent;
        for (const auto &edge : remaining) {
            adjacent[edge.a].push_back(edge.b);
            adjacent[edge.b].push_back(edge.a);
        }
        while (!remaining.empty()) {
            unsigned first = remaining.begin()->a;
            for (const auto &edge : remaining) {
                if (adjacent[edge.a].size() != 2) {
                    first = edge.a;
                    break;
                }
                if (adjacent[edge.b].size() != 2) {
                    first = edge.b;
                    break;
                }
            }
            std::vector<unsigned> chain{first};
            unsigned current = first;
            for (;;) {
                unsigned next = ~0u;
                for (auto v : adjacent[current])
                    if (remaining.count(Edge(current, v))) {
                        next = v;
                        break;
                    }
                if (next == ~0u)
                    break;
                remaining.erase(Edge(current, next));
                chain.push_back(next);
                current = next;
                if (current == first || adjacent[current].size() != 2)
                    break;
            }
            if (chain.size() < 3)
                continue;
            try {
                bool periodic = chain.front() == chain.back();
                std::vector<double> lengths;
                for (unsigned i = 1; i < chain.size(); ++i)
                    lengths.push_back((mesh.points[chain[i]] - mesh.points[chain[i - 1]]).norm());
                auto sorted = lengths;
                std::sort(sorted.begin(), sorted.end());
                double typical = lengths.size() > 3 ? sorted[sorted.size() / 2] : sorted.front();
                if (typical <= 1e-8)
                    continue;
                std::vector<Vec> samples;
                for (unsigned i = 1; i < chain.size(); ++i) {
                    unsigned parts =
                        lengths[i - 1] > 3 * typical
                            ? std::min(16u, unsigned(std::ceil(lengths[i - 1] / typical)))
                            : 1u;
                    for (unsigned j = 0; j < parts; ++j) {
                        double t = double(j) / parts;
                        samples.push_back((1 - t) * mesh.points[chain[i - 1]] +
                                          t * mesh.points[chain[i]]);
                    }
                }
                samples.push_back(mesh.points[chain.back()]);
                unsigned count = unsigned(samples.size()) - (periodic ? 1u : 0u);
                Handle(PointHandleArray) points = new PointHandleArray(1, count);
                Handle(RealHandleArray) parameters = new RealHandleArray(1, samples.size());
                double distance = 0;
                for (unsigned i = 0; i < samples.size(); ++i) {
                    if (i)
                        distance += (samples[i] - samples[i - 1]).norm();
                    parameters->SetValue(i + 1, distance);
                    if (i < count)
                        points->SetValue(i + 1, point(samples[i]));
                }
                GeomAPI_Interpolate interpolate(points, parameters, periodic, 1e-8);
                interpolate.Perform();
                if (interpolate.IsDone())
                    for (unsigned i = 1; i < chain.size(); ++i)
                        smooth_boundaries.emplace(Edge(chain[i - 1], chain[i]),
                                                  interpolate.Curve());
            } catch (const Standard_Failure &) { /* A tiny or branching run retains local curves. */
            }
        }
    }
    std::vector<std::unique_ptr<PointProjector>> projectors(patches.size());
    auto projector = [&](unsigned i) -> PointProjector & {
        if (!projectors[i])
            projectors[i] = std::make_unique<PointProjector>(patches[i].surface);
        return *projectors[i];
    };
    for (const auto &kv : top.edges) {
        if (kv.second.size() != 2)
            continue;
        auto a = owner[kv.second[0]], b = owner[kv.second[1]];
        if (a == b)
            continue;
        if (a > b)
            std::swap(a, b);
        const auto &first = patches[a];
        const auto &second = patches[b];
        if (first.kind == "plane" && second.kind == "plane" &&
            std::abs(first.axis.dot(second.axis)) > 0.999999)
            continue;
        if (first.kind == "cylinder" && second.kind == "cylinder" &&
            std::abs(first.axis.dot(second.axis)) > 0.999999) {
            Vec shift = first.origin - second.origin;
            double distance = (shift - shift.dot(first.axis) * first.axis).norm();
            if (distance < tolerance * 0.01)
                continue;
        }
        auto key = std::make_pair(a, b);
        bool fitted_boundary = !analytic(first) || !analytic(second);
        if (!fitted_boundary && !intersect_analytic)
            continue;
        std::vector<Handle(Geom_Curve)> edge_candidates;
        if (fitted_boundary) {
            // At a tangent blend junction a global surface intersection is
            // ill-conditioned. Reconstruct its shared trim from the adjacent
            // support instead, checking the whole edge against both surfaces.
            try {
                auto smooth = smooth_boundaries.find(kv.first);
                if (smooth != smooth_boundaries.end())
                    edge_candidates.push_back(smooth->second);
                int preferred = analytic(first) ? int(a) : (analytic(second) ? int(b) : -1);
                if (preferred >= 0) {
                    auto &project = projector(preferred);
                    const auto &surface = patches[preferred].surface;
                    const auto &p = project.perform(mesh.points[kv.first.a]);
                    double u0, v0;
                    p.LowerDistanceParameters(u0, v0);
                    auto start = p.NearestPoint();
                    const auto &q = project.perform(mesh.points[kv.first.b]);
                    double u1, v1;
                    q.LowerDistanceParameters(u1, v1);
                    auto finish = q.NearestPoint();
                    if (patches[preferred].kind == "plane") {
                        gp_Vec direction(start, finish);
                        if (direction.Magnitude() > 1e-10)
                            edge_candidates.push_back(new Geom_Line(start, gp_Dir(direction)));
                    } else {
                        if (surface->IsUPeriodic())
                            u1 += std::round((u0 - u1) / surface->UPeriod()) * surface->UPeriod();
                        if (surface->IsVPeriodic())
                            v1 += std::round((v0 - v1) / surface->VPeriod()) * surface->VPeriod();
                        edge_candidates.push_back(surface->UIso((u0 + u1) / 2));
                        edge_candidates.push_back(surface->VIso((v0 + v1) / 2));
                    }
                }
                PointArray samples(1, 5);
                for (unsigned i = 0; i < 5; ++i) {
                    double t = double(i) / 4;
                    Vec p = (1 - t) * mesh.points[kv.first.a] + t * mesh.points[kv.first.b];
                    gp_Pnt on;
                    if (preferred >= 0)
                        on = projector(preferred).perform(p).NearestPoint();
                    else {
                        auto p1 = projector(a).perform(p).NearestPoint(),
                             p2 = projector(b).perform(p).NearestPoint();
                        on = gp_Pnt((p1.X() + p2.X()) / 2, (p1.Y() + p2.Y()) / 2,
                                    (p1.Z() + p2.Z()) / 2);
                    }
                    if (i == 0)
                        on = point(mesh.points[kv.first.a]);
                    if (i == 4)
                        on = point(mesh.points[kv.first.b]);
                    if (on.Distance(point(p)) >
                        .1 * (mesh.points[kv.first.a] - mesh.points[kv.first.b]).norm())
                        on = point(p);
                    samples.SetValue(i + 1, on);
                }
                Handle(PointHandleArray) points = new PointHandleArray(samples);
                GeomAPI_Interpolate fit(points, false, 1e-8);
                fit.Perform();
                if (fit.IsDone())
                    edge_candidates.push_back(fit.Curve());
            } catch (const Standard_Failure
                         &) { /* Retain the mesh trim if no common curve can be fitted. */
            }
        } else {
            auto found = intersections.find(key);
            if (found == intersections.end()) {
                std::vector<Handle(Geom_Curve)> curves;
                try {
                    GeomAPI_IntSS intersection(patches[a].surface, patches[b].surface,
                                               std::min(1e-6, tolerance * 0.01));
                    if (intersection.IsDone())
                        for (int i = 1; i <= intersection.NbLines(); ++i)
                            curves.push_back(intersection.Line(i));
                } catch (const Standard_Failure
                             &) { /* Tangent supports may not have a stable intersection. */
                }
                found = intersections.emplace(key, std::move(curves)).first;
            }
            edge_candidates = found->second;
        }
        double best = 1e100, pa = 0, pb = 0;
        Handle(Geom_Curve) selected;
        struct ChartSpan {
            Eigen::Vector2d origin, delta;
            double pu = 0, pv = 0;
        };
        std::array<ChartSpan, 2> charts;
        bool charts_valid = true;
        if (fitted_boundary)
            for (unsigned side = 0; side < 2; ++side) {
                unsigned index = side ? b : a;
                auto &chart = charts[side];
                auto &project = projector(index);
                const auto &start = project.perform(mesh.points[kv.first.a]);
                if (!start.NbPoints()) {
                    charts_valid = false;
                    break;
                }
                start.LowerDistanceParameters(chart.origin.x(), chart.origin.y());
                const auto &end = project.perform(mesh.points[kv.first.b]);
                if (!end.NbPoints()) {
                    charts_valid = false;
                    break;
                }
                Eigen::Vector2d finish;
                end.LowerDistanceParameters(finish.x(), finish.y());
                const auto &support = patches[index].surface;
                if (support->IsUPeriodic()) {
                    chart.pu = support->UPeriod();
                    finish.x() += std::round((chart.origin.x() - finish.x()) / chart.pu) * chart.pu;
                }
                if (support->IsVPeriodic()) {
                    chart.pv = support->VPeriod();
                    finish.y() += std::round((chart.origin.y() - finish.y()) / chart.pv) * chart.pv;
                }
                chart.delta = finish - chart.origin;
                if (chart.delta.squaredNorm() < 1e-24)
                    charts_valid = false;
            }
        for (const auto &c : edge_candidates) {
            if (fitted_boundary && !charts_valid)
                continue;
            GeomAPI_ProjectPointOnCurve p(point(mesh.points[kv.first.a]), c),
                q(point(mesh.points[kv.first.b]), c);
            if (!p.NbPoints() || !q.NbPoints())
                continue;
            double distance = p.LowerDistance() + q.LowerDistance();
            double p0 = p.LowerDistanceParameter(), p1 = q.LowerDistanceParameter();
            if (c->IsPeriodic())
                p1 += std::round((p0 - p1) / c->Period()) * c->Period();
            if (fitted_boundary) {
                if (distance >
                    std::min(tolerance * .001,
                             (mesh.points[kv.first.a] - mesh.points[kv.first.b]).norm() * .001))
                    continue;
                bool on_both = true;
                const Vec start = mesh.points[kv.first.a], chord = mesh.points[kv.first.b] - start;
                const double length2 = chord.squaredNorm();
                if (length2 < 1e-24)
                    continue;
                double previous = 0;
                std::array<double, 2> previous_uv{};
                auto regular_projection = [&](const SurfaceProjection &projection, unsigned side,
                                              double t) {
                    if (!projection.NbPoints())
                        return false;
                    const auto &chart = charts[side];
                    Eigen::Vector2d uv;
                    projection.LowerDistanceParameters(uv.x(), uv.y());
                    auto expected = chart.origin + t * chart.delta;
                    if (chart.pu)
                        uv.x() += std::round((expected.x() - uv.x()) / chart.pu) * chart.pu;
                    if (chart.pv)
                        uv.y() += std::round((expected.y() - uv.y()) / chart.pv) * chart.pv;
                    auto offset = uv - chart.origin;
                    double progress = offset.dot(chart.delta) / chart.delta.squaredNorm();
                    bool regular = progress >= previous_uv[side] - .01 && progress <= 1.01 &&
                                   (offset - progress * chart.delta).squaredNorm() <=
                                       .0625 * chart.delta.squaredNorm() + 1e-14;
                    previous_uv[side] = progress;
                    return regular;
                };
                for (double t : {.125, .25, .5, .75, .875}) {
                    auto v = c->Value((1 - t) * p0 + t * p1);
                    Vec point(v.X(), v.Y(), v.Z());
                    double progress = (point - start).dot(chord) / length2;
                    if (progress < previous - 1e-6 || progress > 1 + 1e-6 ||
                        (point - start - progress * chord).squaredNorm() > .0625 * length2 + 1e-16)
                        on_both = false;
                    previous = progress;
                    const auto &pa = projector(a).perform(point);
                    if (!pa.NbPoints() || pa.LowerDistance() > tolerance ||
                        !regular_projection(pa, 0, t))
                        on_both = false;
                    const auto &pb = projector(b).perform(point);
                    if (!pb.NbPoints() || pb.LowerDistance() > tolerance ||
                        !regular_projection(pb, 1, t))
                        on_both = false;
                }
                if (!on_both)
                    continue;
            }
            if (distance < best) {
                best = distance;
                pa = p0;
                pb = p1;
                selected = c;
                if (fitted_boundary)
                    break;
            }
        }
        if (selected.IsNull() || best > tolerance * 0.2)
            continue;
        if (selected->IsPeriodic()) {
            double period = selected->Period();
            pb += std::round((pa - pb) / period) * period;
        }
        if (std::abs(pa - pb) < 1e-12)
            continue;
        try {
            auto edge =
                BRepBuilderAPI_MakeEdge(selected, std::min(pa, pb), std::max(pa, pb)).Edge();
            if (pb < pa)
                edge.Reverse();
            BRep_Builder builder;
            for (TopExp_Explorer v(edge, TopAbs_VERTEX); v.More(); v.Next())
                builder.UpdateVertex(TopoDS::Vertex(v.Current()), std::max(1e-6, tolerance * 0.1));
            out.emplace(kv.first, edge);
        } catch (const Standard_Failure &) { /* Keep the original boundary for this candidate. */
        }
    }
    trace("surface intersections complete");
    return out;
}
TopoDS_Face curved_triangle(const Mesh &mesh, const Tri &tri, const Handle(Geom_Surface) & surface,
                            double tolerance, const Boundaries &boundaries,
                            PointProjector &projector) {
    gp_Pnt2d uv[3];
    for (unsigned k = 0; k < 3; ++k) {
        const auto &projection = projector.perform(mesh.points[tri[k]]);
        if (projection.NbPoints() == 0)
            throw std::runtime_error("Surface projection failed");
        double u, v;
        projection.LowerDistanceParameters(u, v);
        uv[k] = gp_Pnt2d(u, v);
    }
    // Unwrap only local triangle spans; the seam is sewn back together by OCCT.
    for (unsigned k = 1; k < 3; ++k) {
        if (surface->IsUPeriodic()) {
            double period = surface->UPeriod();
            uv[k].SetX(uv[k].X() + std::round((uv[0].X() - uv[k].X()) / period) * period);
        }
        if (surface->IsVPeriodic()) {
            double period = surface->VPeriod();
            uv[k].SetY(uv[k].Y() + std::round((uv[0].Y() - uv[k].Y()) / period) * period);
        }
    }
    double area = (uv[1].X() - uv[0].X()) * (uv[2].Y() - uv[0].Y()) -
                  (uv[1].Y() - uv[0].Y()) * (uv[2].X() - uv[0].X());
    if (std::abs(area) < 1e-18)
        throw std::runtime_error("Singular surface triangle");
    BRepBuilderAPI_MakeWire wire;
    for (unsigned k = 0; k < 3; ++k) {
        auto boundary = boundaries.find(Edge(tri[k], tri[(k + 1) % 3]));
        if (boundary != boundaries.end()) {
            auto edge = boundary->second;
            if (tri[k] > tri[(k + 1) % 3])
                edge.Reverse();
            ShapeFix_Edge fix;
            fix.FixAddPCurve(edge, surface, TopLoc_Location(), false, tolerance * 0.01);
            wire.Add(edge);
            continue;
        }
        gp_Vec2d d(uv[k], uv[(k + 1) % 3]);
        double length = d.Magnitude();
        Handle(Geom2d_Line) line = new Geom2d_Line(uv[k], gp_Dir2d(d));
        auto edge = BRepBuilderAPI_MakeEdge(line, surface, 0, length).Edge();
        BRepLib::BuildCurve3d(edge, tolerance * 0.1);
        wire.Add(edge);
    }
    if (!wire.IsDone())
        throw std::runtime_error("Cannot build surface wire");
    BRepBuilderAPI_MakeFace make(surface, wire.Wire(), true);
    if (!make.IsDone())
        throw std::runtime_error("Cannot trim surface face");
    auto face = make.Face();
    if (area < 0)
        face.Reverse();
    return face;
}
TopoDS_Face flat_triangle(const Mesh &mesh, const Tri &tri, const Boundaries &boundaries = {}) {
    bool curved = false;
    for (unsigned k = 0; k < 3; ++k)
        curved = curved || boundaries.count(Edge(tri[k], tri[(k + 1) % 3]));
    if (curved) {
        BRepBuilderAPI_MakeWire wire;
        for (unsigned k = 0; k < 3; ++k) {
            auto found = boundaries.find(Edge(tri[k], tri[(k + 1) % 3]));
            if (found == boundaries.end())
                wire.Add(BRepBuilderAPI_MakeEdge(point(mesh.points[tri[k]]),
                                                 point(mesh.points[tri[(k + 1) % 3]]))
                             .Edge());
            else {
                auto edge = found->second;
                if (tri[k] > tri[(k + 1) % 3])
                    edge.Reverse();
                wire.Add(edge);
            }
        }
        if (!wire.IsDone())
            throw std::runtime_error("Cannot connect analytic boundary to plane");
        BRepBuilderAPI_MakeFace make(wire.Wire(), true);
        if (!make.IsDone())
            throw std::runtime_error("Cannot trim planar face with analytic boundary");
        auto face = make.Face();
        BRepAdaptor_Surface support(face);
        auto n = support.Plane().Axis().Direction();
        Vec original = (mesh.points[tri[1]] - mesh.points[tri[0]])
                           .cross(mesh.points[tri[2]] - mesh.points[tri[0]]);
        double dot = original.dot(Vec(n.X(), n.Y(), n.Z()));
        if (face.Orientation() == TopAbs_REVERSED)
            dot = -dot;
        if (dot < 0)
            face.Reverse();
        return face;
    }
    BRepBuilderAPI_MakePolygon poly;
    for (auto v : tri)
        poly.Add(point(mesh.points[v]));
    poly.Close();
    BRepBuilderAPI_MakeFace face(poly.Wire());
    if (!face.IsDone())
        throw std::runtime_error("Cannot build mesh facet");
    return face.Face();
}
TopoDS_Shape build(const Mesh &mesh, const std::vector<Patch> &patches, double tolerance,
                   bool unify, unsigned boundary_mode = 2) {
    BRepBuilderAPI_Sewing sewing(tolerance, true, true, true, false);
    Boundaries boundaries;
    if (!patches.empty() && boundary_mode)
        try {
            boundaries = shared_boundaries(mesh, patches, tolerance, boundary_mode == 2);
        } catch (const Standard_Failure &error) {
            throw std::runtime_error(std::string("Boundaries: ") + failure_message(error));
        }
    if (patches.empty())
        for (const auto &t : mesh.triangles)
            sewing.Add(flat_triangle(mesh, t));
    else
        for (const auto &p : patches) {
            if (p.kind == "plane") {
                for (const auto &face : planar_patch(mesh, p, boundaries))
                    sewing.Add(face);
                continue;
            }
            auto natural = natural_patch(mesh, p, tolerance);
            if (!natural.IsNull()) {
                sewing.Add(natural);
                continue;
            }
            try {
                auto trimmed = trimmed_patch(mesh, p, tolerance, boundaries);
                if (!trimmed.empty()) {
                    for (const auto &face : trimmed)
                        sewing.Add(face);
                    continue;
                }
            } catch (
                const Standard_Failure &) { /* Retain triangulated trims for an unsuitable chart. */
            }
            // Closed natural surfaces require no artificial triangulation seams or
            // parameter-space triangles at their singular poles.
            if ((p.kind == "sphere" || p.kind == "torus") &&
                p.triangles.size() == mesh.triangles.size()) {
                sewing.Add(BRepBuilderAPI_MakeFace(p.surface, tolerance).Face());
                continue;
            }
            PointProjector projector(p.surface);
            for (auto id : p.triangles) {
                const auto &t = mesh.triangles[id];
                // Planes keep original corners and exact shared straight edges.
                try {
                    if (p.kind == "plane")
                        sewing.Add(flat_triangle(mesh, t, boundaries));
                    else
                        sewing.Add(
                            curved_triangle(mesh, t, p.surface, tolerance, boundaries, projector));
                } catch (const Standard_Failure &error) {
                    throw std::runtime_error(p.kind + " face " + std::to_string(id) + ": " +
                                             failure_message(error));
                }
            }
        }
    trace("sewing");
    try {
        sewing.Perform();
    } catch (const Standard_Failure &error) {
        throw std::runtime_error(std::string("Sewing: ") + failure_message(error));
    }
    auto sewn = sewing.SewedShape();
    if (sewn.IsNull())
        throw std::runtime_error("Sewing produced an empty shape");
    TopoDS_Shape shape;
    trace("shell closure");
    try {
        shape = close_shells(sewn, tolerance);
        BRepLib::SameParameter(shape, std::min(1e-5, tolerance * .001), true);
    } catch (const Standard_Failure &error) {
        throw std::runtime_error(std::string("Shell closure: ") + failure_message(error));
    }
    if (!patches.empty() && boundary_mode)
        trace_invalid_faces(shape, &patches);
    if (!patches.empty() && boundary_mode == 2 && !BRepCheck_Analyzer(shape).IsValid()) {
        ShapeFix_Shape heal(shape);
        heal.SetPrecision(tolerance * .1);
        heal.SetMaxTolerance(tolerance);
        heal.Perform();
        auto fixed = heal.Shape();
        if (fixed.IsNull() || !BRepCheck_Analyzer(fixed).IsValid())
            throw std::runtime_error("Intersection-boundary topology is invalid");
        shape = fixed;
    }
    trace("same-domain reduction");
    if (unify) {
        try {
            double before = volume(shape);
            BRepBuilderAPI_Copy copy(shape);
            ShapeUpgrade_UnifySameDomain merge(copy.Shape(), true, true, true);
            merge.SetLinearTolerance(tolerance);
            merge.Build();
            auto reduced = merge.Shape();
            if (!reduced.IsNull() && BRepCheck_Analyzer(reduced).IsValid() &&
                std::abs(volume(reduced) - before) < std::max(.001, std::abs(before) * 1e-7))
                shape = reduced;
        } catch (const Standard_Failure
                     &) { /* Same-domain reduction is optional; retain the sewn shape. */
        }
    }
    return shape;
}
} // namespace
TopoDS_Shape faceted(const Mesh &mesh, double tolerance) {
    return build(mesh, {}, tolerance, true);
}
TopoDS_Shape reconstruct(const Mesh &mesh, const std::vector<Patch> &patches, Context &ctx) {
    // Try exact analytic intersections, then smooth fitted junctions with
    // projected analytic trims, then independent projected trims. Every mode
    // must have valid topology; conversion applies the corrected-volume gate
    // before accepting a candidate or invoking faceted recovery.
    for (unsigned mode : {2u, 1u, 0u}) {
        try {
            auto shape = build(mesh, patches, ctx.options.tolerance, true, mode);
            Diagnostics info;
            if (inspect(shape, info) && (mode == 0 || volume_agrees(shape, mesh, patches)))
                return shape;
            ShapeFix_Shape fix(shape);
            fix.SetPrecision(ctx.options.tolerance);
            fix.SetMaxTolerance(ctx.options.tolerance);
            fix.Perform();
            auto fixed = fix.Shape();
            if (!fixed.IsNull() && inspect(fixed, info) &&
                (mode == 0 || volume_agrees(fixed, mesh, patches)))
                return fixed;
        } catch (const Standard_Failure &error) {
            trace(failure_message(error, "Shared-boundary reconstruction failed"));
        } catch (const std::runtime_error &error) {
            trace(error.what());
        }
    }
    throw std::runtime_error("No valid surface reconstruction passed the corrected-volume checks");
}
TopoDS_Shape recover(const Mesh &mesh, const std::vector<Patch> &patches, Context &ctx) {
    auto best = faceted(mesh, ctx.options.tolerance);
    Diagnostics bare;
    if (!inspect(best, bare))
        throw std::runtime_error("Faceted recovery is not a valid closed solid");
    auto top = topology(mesh);
    std::vector<unsigned> curved;
    std::set<unsigned> kept;
    for (unsigned i = 0; i < patches.size(); ++i)
        if (patches[i].kind != "plane")
            curved.push_back(i);
    unsigned attempts = 0;
    auto assemble = [&](const std::set<unsigned> &selected) {
        std::vector<Patch> candidate;
        for (unsigned i = 0; i < patches.size(); ++i) {
            if (patches[i].kind == "plane" || selected.count(i))
                candidate.push_back(patches[i]);
            else
                for (auto t : patches[i].triangles) {
                    Patch p;
                    p.kind = "plane";
                    p.triangles = {t};
                    p.origin = mesh.points[mesh.triangles[t][0]];
                    p.axis = top.normals[t];
                    p.surface = new Geom_Plane(
                        gp_Ax3(point(p.origin), gp_Dir(p.axis.x(), p.axis.y(), p.axis.z())));
                    candidate.push_back(std::move(p));
                }
        }
        return candidate;
    };
    std::function<void(const std::vector<unsigned> &)> restore;
    restore = [&](const std::vector<unsigned> &group) {
        if (group.empty() || attempts >= 32 || ctx.expired())
            return;
        ++attempts;
        auto selected = kept;
        selected.insert(group.begin(), group.end());
        auto candidate = assemble(selected);
        try {
            auto shape = build(mesh, candidate, ctx.options.tolerance, true, 0);
            Diagnostics info;
            if (inspect(shape, info) && volume_agrees(shape, mesh, candidate)) {
                best = shape;
                kept = std::move(selected);
                return;
            }
        } catch (const Standard_Failure &) {
        } catch (const std::runtime_error &) {
        }
        if (group.size() > 1) {
            auto middle = group.begin() + group.size() / 2;
            restore(std::vector<unsigned>(group.begin(), middle));
            restore(std::vector<unsigned>(middle, group.end()));
        }
    };
    restore(curved);
    ctx.diagnostics.fallbacks.push_back("Recovery retained " + std::to_string(kept.size()) +
                                        " of " + std::to_string(curved.size()) + " curved patches");
    return best;
}
TopoDS_Shape unite(const std::vector<TopoDS_Shape> &shapes, Context &ctx) {
    if (shapes.empty())
        throw std::runtime_error("No shapes to combine");
    if (!ctx.options.fuse)
        return compound(shapes);
    std::vector<TopoDS_Shape> bodies;
    for (const auto &s : shapes)
        for (TopExp_Explorer ex(s, TopAbs_SOLID); ex.More(); ex.Next())
            bodies.push_back(ex.Current());
    if (bodies.empty())
        throw std::runtime_error("No solids to fuse");
    auto shape = bodies.front();
    for (std::size_t i = 1; i < bodies.size(); ++i) {
        // OCCT's persistent default pool would compete with the bounded fitting
        // pool for Emscripten pthread slots on the next conversion.
        BRepAlgoAPI_Fuse fuse(shape, bodies[i]);
        fuse.SetRunParallel(false);
        fuse.SetFuzzyValue(ctx.options.tolerance);
        fuse.Build();
        if (!fuse.IsDone() || fuse.Shape().IsNull() || !BRepCheck_Analyzer(fuse.Shape()).IsValid())
            throw std::runtime_error("Boolean fusion failed validation");
        shape = fuse.Shape();
    }
    return shape;
}
bool inspect(const TopoDS_Shape &shape, Diagnostics &info) {
    trace("solid validation and volume");
    if (shape.IsNull())
        return false;
    info.faces = 0;
    info.solids = 0;
    info.surface_counts.clear();
    ShapeMap faces, solids;
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    TopExp::MapShapes(shape, TopAbs_SOLID, solids);
    info.faces = faces.Extent();
    info.solids = solids.Extent();
    for (int i = 1; i <= faces.Extent(); ++i) {
        BRepAdaptor_Surface surface(TopoDS::Face(faces(i)));
        std::string kind;
        switch (surface.GetType()) {
        case GeomAbs_Plane:
            kind = "plane";
            break;
        case GeomAbs_Cylinder:
            kind = "cylinder";
            break;
        case GeomAbs_Cone:
            kind = "cone";
            break;
        case GeomAbs_Sphere:
            kind = "sphere";
            break;
        case GeomAbs_Torus:
            kind = "torus";
            break;
        case GeomAbs_BSplineSurface:
            kind = "bspline";
            break;
        case GeomAbs_SurfaceOfRevolution:
            kind = "revolution";
            break;
        case GeomAbs_SurfaceOfExtrusion:
            kind = "extrusion";
            break;
        default:
            kind = "other";
        }
        ++info.surface_counts[kind];
    }
    info.volume = volume(shape);
    bool checked = BRepCheck_Analyzer(shape).IsValid();
    if (!checked)
        trace_invalid_faces(shape);
    info.valid = info.solids > 0 && info.volume > 0 && std::isfinite(info.volume) && checked;
    // The validity analyzer permits some open shells inside a solid; count actual
    // face incidences, treating a periodic seam as two uses of the same edge.
    unsigned bad_edges = 0;
    ShapeMap solid_faces;
    for (int solid = 1; solid <= solids.Extent(); ++solid) {
        ShapeAncestors incidence;
        TopExp::MapShapesAndAncestors(solids(solid), TopAbs_EDGE, TopAbs_FACE, incidence);
        TopExp::MapShapes(solids(solid), TopAbs_FACE, solid_faces);
        for (int i = 1; i <= incidence.Extent(); ++i)
            if (!BRep_Tool::Degenerated(TopoDS::Edge(incidence.FindKey(i))) &&
                incidence.FindFromIndex(i).Extent() != 2) {
                info.valid = false;
                ++bad_edges;
            }
    }
    if (solid_faces.Extent() != faces.Extent())
        info.valid = false;
    trace(("validation: BRepCheck=" + std::to_string(checked) +
           " edge-incidence-errors=" + std::to_string(bad_edges) +
           " faces=" + std::to_string(info.faces) + " volume=" + std::to_string(info.volume))
              .c_str());
    trace("solid validation complete");
    return info.valid;
}
bool volume_agrees(const TopoDS_Shape &shape, const Mesh &mesh, const std::vector<Patch> &patches) {
    auto top = topology(mesh);
    double meshvol = signed_volume(mesh), change = 0, gain = 0;
    unsigned curved = 0;
    for (const auto &p : patches)
        if (p.kind != "plane") {
            ++curved;
            PointProjector projector(p.surface);
            for (auto i : p.triangles) {
                const auto &tri = mesh.triangles[i];
                double delta = 0;
                for (unsigned k = 0; k < 3; ++k) {
                    Vec midpoint = (mesh.points[tri[k]] + mesh.points[tri[(k + 1) % 3]]) / 2;
                    const auto &projection = projector.perform(midpoint);
                    if (!projection.NbPoints())
                        return false;
                    gp_Pnt q = projection.NearestPoint();
                    Vec offset(q.X() - midpoint.x(), q.Y() - midpoint.y(), q.Z() - midpoint.z());
                    delta += offset.dot(top.normals[i]) / 3;
                }
                change += delta * top.areas[i];
                gain += std::abs(delta) * top.areas[i];
            }
        }
    double expected = std::abs(meshvol) + (meshvol < 0 ? -change : change);
    double allowed = 0.5 * gain + curved * 1e-3 + std::abs(meshvol) * 1e-6 + 1e-3;
    return std::abs(volume(shape) - expected) <= allowed;
}
} // namespace stl2curves::detail

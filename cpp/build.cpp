#include "internal.hpp"
#include "occt_error.hpp"
#include "projection.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
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
#include <BRepClass_FaceClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRepTools.hxx>
#include <BRepTools_History.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Eigen/Geometry>
#include <GCPnts_AbscissaPoint.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GProp_GProps.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <GeomAPI_IntSS.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomProjLib.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <ShapeAnalysis_ShapeTolerance.hxx>
#include <ShapeFix_Edge.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_ShapeTolerance.hxx>
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
#include <nanoflann.hpp>
#include <set>
#include <tuple>

namespace stl2curves::detail {
TopoDS_Shape reference_surface_face(const Patch &patch, bool generous);
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
        if (std::getenv("STL2CURVES_TRACE_GEOMETRY")) {
            unsigned samples = 0;
            for (TopExp_Explorer edges(face, TopAbs_EDGE); edges.More() && samples < 3;
                 edges.Next()) {
                auto edge = TopoDS::Edge(edges.Current());
                double first, last, pfirst, plast;
                auto curve = BRep_Tool::Curve(edge, first, last);
                auto pc = BRep_Tool::CurveOnSurface(edge, face, pfirst, plast);
                if (curve.IsNull() || pc.IsNull())
                    continue;
                double gap = 0;
                for (unsigned i = 0; i <= 16; ++i) {
                    double t = first + (last - first) * i / 16;
                    auto uv = pc->Value(t);
                    gap = std::max(gap, curve->Value(t).Distance(support->Value(uv.X(), uv.Y())));
                }
                if (gap > BRep_Tool::Tolerance(edge)) {
                    message +=
                        " pcurve-gap/tol=" + std::to_string(gap / BRep_Tool::Tolerance(edge)) +
                        " gap=" + std::to_string(gap) + " range=" + std::to_string(first) + ":" +
                        std::to_string(last) + "/" + std::to_string(pfirst) + ":" +
                        std::to_string(plast);
                    ++samples;
                }
            }
        }
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
bool repair_planar_pcurves(const TopoDS_Shape &shape, double tolerance) {
    BRep_Builder builder;
    ShapeFix_Edge fix;
    bool bounded = true;
    for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next()) {
        auto face = TopoDS::Face(it.Current());
        if (BRepAdaptor_Surface(face).GetType() != GeomAbs_Plane ||
            BRepCheck_Analyzer(face).IsValid())
            continue;
        // Sewing can replace a planar chord by its curved neighbor's edge.
        // Its old pcurve then no longer describes the retained 3D curve.
        for (TopExp_Explorer edges(face, TopAbs_EDGE); edges.More(); edges.Next()) {
            auto edge = TopoDS::Edge(edges.Current());
            if (BRep_Tool::Degenerated(edge))
                continue;
            double first, last;
            auto curve = BRep_Tool::Curve(edge, first, last);
            if (curve.IsNull())
                return false;
            double error = std::min(1e-7, tolerance * .01);
            auto pc = GeomProjLib::Curve2d(curve, first, last, BRep_Tool::Surface(face), error);
            if (pc.IsNull())
                return false;
            // A planar pcurve can be computed implicitly by BRep_Tool, in
            // which case FixAddPCurve considers it present and adds nothing.
            // Store it explicitly so SameParameter measures the real gap.
            builder.UpdateEdge(edge, pc, face, std::max(error, BRep_Tool::Tolerance(edge)));
            builder.Range(edge, face, first, last);
            builder.SameRange(edge, false);
            builder.SameParameter(edge, false);
            fix.FixSameParameter(edge, face, tolerance * .01);
            bounded = bounded && BRep_Tool::Tolerance(edge) <= tolerance;
        }
    }
    return bounded;
}
using SurfaceOwners = std::map<const Geom_Surface *, unsigned>;
void failed_patches(const TopoDS_Shape &shape, const std::vector<Patch> &patches, double tolerance,
                    std::set<unsigned> &failed, const SurfaceOwners &owners) {
    auto blame = [&](const TopoDS_Face &face) {
        auto surface = BRep_Tool::Surface(face);
        auto owner = owners.find(surface.get());
        if (owner != owners.end()) {
            failed.insert(owner->second);
            return;
        }
        for (unsigned i = 0; i < patches.size(); ++i)
            if (surface == patches[i].surface)
                failed.insert(i);
    };
    for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next()) {
        auto face = TopoDS::Face(it.Current());
        bool invalid = !BRepCheck_Analyzer(face).IsValid();
        for (TopExp_Explorer edges(face, TopAbs_EDGE); edges.More(); edges.Next())
            invalid = invalid || BRep_Tool::Tolerance(TopoDS::Edge(edges.Current())) > tolerance;
        if (invalid)
            blame(face);
    }
    ShapeAncestors incidence;
    TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, incidence);
    for (int i = 1; i <= incidence.Extent(); ++i)
        if (!BRep_Tool::Degenerated(TopoDS::Edge(incidence.FindKey(i))) &&
            incidence.FindFromIndex(i).Extent() != 2)
            for (const auto &face : incidence.FindFromIndex(i))
                blame(TopoDS::Face(face));
}
TopoDS_Shape heal_candidate(const TopoDS_Shape &shape, const Mesh &mesh,
                            const std::vector<Patch> &patches, double tolerance) {
    // Precision is a working resolution, not permission to move every corner
    // by the sewing tolerance. Coarse precision collapses short STL edges.
    for (double precision : {std::min(1e-7, tolerance), std::min(1e-5, tolerance * .001),
                             std::min(1e-4, tolerance * .01)}) {
        BRepBuilderAPI_Copy copy(shape, true, false);
        ShapeFix_Shape fix(copy.Shape());
        fix.SetPrecision(precision);
        fix.SetMaxTolerance(tolerance);
        fix.Perform();
        auto repaired = fix.Shape();
        if (repaired.IsNull())
            continue;
        if (!BRepCheck_Analyzer(repaired).IsValid()) {
            if (!repair_planar_pcurves(repaired, tolerance))
                continue;
            ShapeFix_Shape again(repaired);
            again.SetPrecision(precision);
            again.SetMaxTolerance(tolerance);
            again.Perform();
            repaired = again.Shape();
        }
        ShapeAnalysis_ShapeTolerance analysis;
        if (!repaired.IsNull() && BRepCheck_Analyzer(repaired).IsValid() &&
            analysis.Tolerance(repaired, 1) <= tolerance && volume_agrees(repaired, mesh, patches))
            return repaired;
    }
    return {};
}
TopoDS_Shape without_collapsed_edges(const TopoDS_Shape &input, double tolerance) {
    // Upstream _without_collapsed: a short open curve whose vertices sewing
    // merged is not a genuine closed outline. Leaving it pinches the wire.
    if (input.IsNull())
        return input;
    ShapeMap edges;
    TopExp::MapShapes(input, TopAbs_EDGE, edges);
    Handle(BRepTools_ReShape) reshape;
    for (int i = 1; i <= edges.Extent(); ++i) {
        auto edge = TopoDS::Edge(edges(i));
        if (BRep_Tool::Degenerated(edge) ||
            !TopExp::FirstVertex(edge).IsSame(TopExp::LastVertex(edge)))
            continue;
        BRepAdaptor_Curve curve(edge);
        if (curve.Value(curve.FirstParameter()).Distance(curve.Value(curve.LastParameter())) >
                1e-6 &&
            GCPnts_AbscissaPoint::Length(curve) <= 2 * tolerance) {
            if (reshape.IsNull())
                reshape = new BRepTools_ReShape;
            reshape->Remove(edge);
        }
    }
    if (reshape.IsNull())
        return input;
    auto shape = reshape->Apply(input);
    Handle(BRepTools_ReShape) bare;
    for (TopExp_Explorer faces(shape, TopAbs_FACE); faces.More(); faces.Next())
        if (!TopExp_Explorer(faces.Current(), TopAbs_EDGE).More()) {
            if (bare.IsNull())
                bare = new BRepTools_ReShape;
            bare->Remove(faces.Current());
        }
    return bare.IsNull() ? shape : bare->Apply(shape);
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
struct BoundaryCloud {
    const std::vector<Vec> &points;
    std::size_t kdtree_get_point_count() const { return points.size(); }
    double kdtree_get_pt(std::size_t i, std::size_t d) const { return points[i][d]; }
    template <class T> bool kdtree_get_bbox(T &) const { return false; }
};
TopoDS_Shape rejoin_exact(TopoDS_Shape shape, const Mesh &mesh, double tolerance) {
    // Upstream _rejoin: join exact copies before sewing, but leave every
    // corner meeting a short edge to sewing so that it can collapse together.
    constexpr double snap = 1e-7;
    BoundaryCloud cloud{mesh.points};
    using Tree =
        nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, BoundaryCloud>,
                                            BoundaryCloud, 3, std::size_t>;
    Tree tree(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(16));
    tree.buildIndex();
    ShapeMap vertices;
    TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
    std::set<int> short_vertices;
    for (TopExp_Explorer ex(shape, TopAbs_EDGE); ex.More(); ex.Next()) {
        auto edge = TopoDS::Edge(ex.Current());
        if (BRep_Tool::Degenerated(edge))
            continue;
        auto a = TopExp::FirstVertex(edge), b = TopExp::LastVertex(edge);
        if (!a.IsNull() && !b.IsNull() &&
            BRep_Tool::Pnt(a).Distance(BRep_Tool::Pnt(b)) <= 2 * tolerance) {
            short_vertices.insert(vertices.FindIndex(a));
            short_vertices.insert(vertices.FindIndex(b));
        }
    }
    std::map<std::size_t, std::vector<int>> groups;
    for (int i = 1; i <= vertices.Extent(); ++i) {
        auto p = BRep_Tool::Pnt(TopoDS::Vertex(vertices(i)));
        double query[]{p.X(), p.Y(), p.Z()}, distance;
        std::size_t nearest;
        if (tree.knnSearch(query, 1, &nearest, &distance) && distance <= snap * snap)
            groups[nearest].push_back(i);
    }
    BRep_Builder builder;
    Handle(BRepTools_ReShape) reshape = new BRepTools_ReShape;
    for (const auto &item : groups) {
        const auto &ids = item.second;
        if (ids.size() < 2 ||
            std::any_of(ids.begin(), ids.end(), [&](int i) { return short_vertices.count(i); }))
            continue;
        auto keep = TopoDS::Vertex(vertices(ids.front()));
        double reach = BRep_Tool::Tolerance(keep);
        for (int i : ids) {
            auto vertex = TopoDS::Vertex(vertices(i));
            reach = std::max(reach, BRep_Tool::Tolerance(vertex) +
                                        BRep_Tool::Pnt(vertex).Distance(BRep_Tool::Pnt(keep)));
        }
        builder.UpdateVertex(keep, reach);
        for (std::size_t i = 1; i < ids.size(); ++i)
            reshape->Replace(vertices(ids[i]), keep.Oriented(vertices(ids[i]).Orientation()));
    }
    shape = reshape->Apply(shape);
    ShapeAncestors incidence;
    TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, incidence);
    vertices.Clear();
    TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
    std::map<int, std::pair<int, int>> ends;
    short_vertices.clear();
    for (int i = 1; i <= incidence.Extent(); ++i) {
        auto edge = TopoDS::Edge(incidence.FindKey(i));
        if (BRep_Tool::Degenerated(edge))
            continue;
        auto a = TopExp::FirstVertex(edge), b = TopExp::LastVertex(edge);
        if (a.IsNull() || b.IsNull())
            continue;
        int ia = vertices.FindIndex(a), ib = vertices.FindIndex(b);
        ends[i] = {ia, ib};
        if (BRep_Tool::Pnt(a).Distance(BRep_Tool::Pnt(b)) <= 2 * tolerance) {
            short_vertices.insert(ia);
            short_vertices.insert(ib);
        }
    }
    std::map<std::pair<int, int>, std::vector<int>> pairs;
    for (const auto &item : ends)
        if (incidence.FindFromIndex(item.first).Extent() == 1 &&
            !short_vertices.count(item.second.first) && !short_vertices.count(item.second.second))
            pairs[std::minmax(item.second.first, item.second.second)].push_back(item.first);
    reshape = new BRepTools_ReShape;
    auto face_of = [&](int i) { return TopoDS::Face(incidence.FindFromIndex(i).First()); };
    auto flat = [&](int i) { return BRepAdaptor_Surface(face_of(i)).GetType() == GeomAbs_Plane; };
    for (const auto &item : pairs) {
        std::set<int> used;
        for (int x : item.second)
            for (int y : item.second) {
                if (used.count(x))
                    break;
                if (x == y || used.count(y) || face_of(x).IsSame(face_of(y)))
                    continue;
                int k = flat(x) && !flat(y) ? y : x, r = k == x ? y : x;
                auto keep = TopoDS::Edge(incidence.FindKey(k)),
                     replace = TopoDS::Edge(incidence.FindKey(r));
                BRepAdaptor_Curve ck(keep), cr(replace);
                if (ck.Value((ck.FirstParameter() + ck.LastParameter()) / 2)
                        .Distance(cr.Value((cr.FirstParameter() + cr.LastParameter()) / 2)) >
                    10 * snap)
                    continue;
                bool same_way = TopExp::FirstVertex(keep).IsSame(TopExp::FirstVertex(replace));
                bool same_params = std::abs(ck.FirstParameter() - cr.FirstParameter()) <= 1e-9 &&
                                   std::abs(ck.LastParameter() - cr.LastParameter()) <= 1e-9;
                if (!(same_way && same_params) &&
                    (!flat(r) || ck.GetType() != cr.GetType() ||
                     (ck.GetType() != GeomAbs_Line && ck.GetType() != GeomAbs_Circle)))
                    continue;
                double tol = std::max(BRep_Tool::Tolerance(keep), BRep_Tool::Tolerance(replace));
                if (!flat(r)) {
                    auto face = face_of(r);
                    double first, last;
                    auto pc = BRep_Tool::CurveOnSurface(replace, face, first, last);
                    if (pc.IsNull())
                        continue;
                    TopLoc_Location location;
                    auto surface = BRep_Tool::Surface(face, location);
                    builder.UpdateEdge(keep, pc, surface, location, tol);
                } else
                    builder.UpdateEdge(keep, tol);
                reshape->Replace(replace.Oriented(TopAbs_FORWARD),
                                 keep.Oriented(same_way ? TopAbs_FORWARD : TopAbs_REVERSED));
                used.insert(x);
                used.insert(y);
            }
    }
    return reshape->Apply(shape);
}
struct Boundaries {
    std::set<Edge> singles;
    mutable std::map<std::vector<unsigned>, TopoDS_Edge> runs;
    std::vector<unsigned> owners;
    std::vector<bool> curved;
    std::map<Edge, std::pair<unsigned, unsigned>> sides;
    std::map<Edge, TopoDS_Edge> edges;
    std::map<unsigned, TopoDS_Vertex> vertices;
    std::vector<bool> crowded;
    bool separate(const Edge &edge) const {
        return !crowded.empty() && (crowded[edge.a] || crowded[edge.b]);
    }
};
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
TopoDS_Edge mesh_chord(const Handle(Geom_Surface) & surface, const Vec &a, const Vec &b,
                       const gp_Pnt2d &ua, const gp_Pnt2d &ub, TopoDS_Vertex va, TopoDS_Vertex vb,
                       double tolerance) try {
    Vec d = b - a;
    double length = d.norm();
    if (length < 1e-10)
        return {};
    Handle(Geom_Curve) curve = new Geom_Line(point(a), gp_Dir(d.x(), d.y(), d.z()));
    double error = tolerance * .01;
    auto pc = GeomProjLib::Curve2d(curve, 0, length, surface, error);
    if (pc.IsNull())
        return {};
    auto mid = pc->Value(length * .5);
    pc->Translate(
        gp_Vec2d(surface->IsUPeriodic()
                     ? std::round(((ua.X() + ub.X()) * .5 - mid.X()) / surface->UPeriod()) *
                           surface->UPeriod()
                     : 0,
                 surface->IsVPeriodic()
                     ? std::round(((ua.Y() + ub.Y()) * .5 - mid.Y()) / surface->VPeriod()) *
                           surface->VPeriod()
                     : 0));
    double gap = 0;
    for (double t : {0., .125, .25, .5, .75, .875, 1.}) {
        auto uv = pc->Value(t * length);
        gap = std::max(gap, curve->Value(t * length).Distance(surface->Value(uv.X(), uv.Y())));
    }
    if (gap > tolerance)
        return {};
    BRep_Builder builder;
    double ga = BRep_Tool::Pnt(va).Distance(point(a)), gb = BRep_Tool::Pnt(vb).Distance(point(b));
    if (std::max(ga, gb) > tolerance)
        return {};
    builder.UpdateVertex(va, std::max(1e-7, ga * 1.01));
    builder.UpdateVertex(vb, std::max(1e-7, gb * 1.01));
    BRepBuilderAPI_MakeEdge make(curve, va, vb, 0, length);
    if (!make.IsDone())
        return {};
    auto edge = make.Edge();
    builder.UpdateEdge(edge, pc, surface, TopLoc_Location(), std::max({error, 1e-7, gap * 1.01}));
    builder.Range(edge, 0, length);
    return edge;
} catch (const Standard_Failure &) {
    // A failed optional chord projection keeps the existing trim path.
    return {};
}

TopoDS_Wire run_wire(const Mesh &mesh, const Patch &patch, const std::vector<unsigned> &loop,
                     const Boundaries &boundaries, double tolerance) {
    if (loop.size() < 3 || boundaries.owners.empty())
        return {};
    const auto none = ~0u;
    unsigned own = boundaries.owners[patch.triangles.front()];
    if (own >= boundaries.curved.size())
        return {};
    std::vector<unsigned> tags;
    for (unsigned i = 0; i < loop.size(); ++i) {
        auto found = boundaries.sides.find(Edge(loop[i], loop[(i + 1) % loop.size()]));
        unsigned other = none;
        if (found != boundaries.sides.end())
            other = found->second.first == own ? found->second.second : found->second.first;
        bool curved = boundaries.curved[own] ||
                      (other < boundaries.curved.size() && boundaries.curved[other]);
        tags.push_back(curved ? other : none - 1);
    }
    std::vector<unsigned> cuts;
    for (unsigned i = 0; i < loop.size(); ++i) {
        unsigned before = (i + loop.size() - 1) % loop.size();
        Vec u = mesh.points[loop[i]] - mesh.points[loop[before]],
            v = mesh.points[loop[(i + 1) % loop.size()]] - mesh.points[loop[i]];
        if (tags[i] == none - 1 || tags[i] != tags[before] ||
            u.dot(v) < .8660254037844386 * u.norm() * v.norm())
            cuts.push_back(i);
    }
    std::vector<std::vector<unsigned>> runs;
    if (cuts.empty()) {
        runs.push_back(loop);
        runs.back().push_back(loop.front());
    } else
        for (unsigned j = 0; j < cuts.size(); ++j) {
            unsigned begin = cuts[j], end = cuts[(j + 1) % cuts.size()];
            std::vector<unsigned> run{loop[begin]};
            unsigned k = (begin + 1) % loop.size();
            while (k != end) {
                run.push_back(loop[k]);
                k = (k + 1) % loop.size();
            }
            run.push_back(loop[end]);
            runs.push_back(std::move(run));
        }
    BRepBuilderAPI_MakeWire wire;
    for (auto ids : runs) {
        bool closed = ids.front() == ids.back(), flip;
        if (closed) {
            ids.pop_back();
            auto start = std::min_element(ids.begin(), ids.end());
            std::rotate(ids.begin(), start, ids.end());
            flip = ids[1] > ids.back();
            if (flip)
                std::reverse(ids.begin() + 1, ids.end());
            ids.push_back(ids.front());
        } else {
            flip = ids.front() > ids.back();
            if (flip)
                std::reverse(ids.begin(), ids.end());
        }
        auto cached = boundaries.runs.find(ids);
        TopoDS_Edge edge;
        if (cached != boundaries.runs.end())
            edge = cached->second;
        else
            try {
                auto existing = boundaries.edges.find(Edge(ids[0], ids[1]));
                if (ids.size() == 2 && existing != boundaries.edges.end() && patch.kind != "sphere")
                    edge = existing->second;
                else {
                    Handle(Geom_Curve) curve;
                    double first = 0, last = 0;
                    if (existing != boundaries.edges.end()) {
                        auto candidate = BRep_Tool::Curve(existing->second, first, last);
                        bool on = !candidate.IsNull() &&
                                  !candidate->IsKind(STANDARD_TYPE(Geom_BSplineCurve)) &&
                                  !(patch.kind == "sphere" && ids.size() == 2);
                        double previous = 0, sign = 0;
                        for (unsigned i = 0; i < ids.size() && on; ++i) {
                            GeomAPI_ProjectPointOnCurve projection(point(mesh.points[ids[i]]),
                                                                   candidate);
                            if (!projection.NbPoints() ||
                                projection.LowerDistance() > std::max(1e-6, tolerance * .001)) {
                                on = false;
                                break;
                            }
                            double t = projection.LowerDistanceParameter();
                            if (i && candidate->IsPeriodic())
                                t += std::round((previous - t) / candidate->Period()) *
                                     candidate->Period();
                            if (!i)
                                first = t;
                            else if (i == 1)
                                sign = t - previous;
                            else if ((t - previous) * sign < -1e-12)
                                on = false;
                            previous = t;
                            last = t;
                        }
                        if (on && std::abs(last - first) > 1e-10)
                            curve = candidate;
                    }
                    if (curve.IsNull()) {
                        if (ids.size() == 2) {
                            auto delta = mesh.points[ids.back()] - mesh.points[ids.front()];
                            if (delta.norm() < 1e-10)
                                return {};
                            curve = new Geom_Line(point(mesh.points[ids.front()]),
                                                  gp_Dir(delta.x(), delta.y(), delta.z()));
                            first = 0;
                            last = delta.norm();
                            if (patch.kind == "sphere") {
                                Vec a = mesh.points[ids.front()], b = mesh.points[ids.back()],
                                    mid = (a + b) * .5, radial = mid - patch.origin;
                                double bow = radial.norm() - patch.radius;
                                if (radial.norm() > 1e-10 &&
                                    std::abs(bow) <=
                                        delta.squaredNorm() / (4 * patch.radius) + 1e-4 &&
                                    std::abs(bow) > 1e-10) {
                                    Vec at = mid - bow * radial.normalized();
                                    GC_MakeArcOfCircle arc(point(a), point(at), point(b));
                                    if (arc.IsDone()) {
                                        curve = Handle(Geom_Curve)(arc.Value().get());
                                        first = curve->FirstParameter();
                                        last = curve->LastParameter();
                                    }
                                }
                            }

                        } else {
                            std::vector<Vec> points;
                            for (unsigned i = 0; i < ids.size(); ++i)
                                if (points.empty() ||
                                    (mesh.points[ids[i]] - points.back()).norm() > 1e-6)
                                    points.push_back(mesh.points[ids[i]]);
                            if (points.size() < 3)
                                return {};
                            std::vector<double> lengths;
                            for (unsigned i = 1; i < points.size(); ++i)
                                lengths.push_back((points[i] - points[i - 1]).norm());
                            auto sorted = lengths;
                            std::sort(sorted.begin(), sorted.end());
                            double typical = sorted.size() > 3 ? (sorted[(sorted.size() - 1) / 2] +
                                                                  sorted[sorted.size() / 2]) *
                                                                     .5
                                                               : sorted.front();
                            if (typical < 1e-10)
                                return {};
                            std::vector<Vec> samples;
                            for (unsigned i = 1; i < points.size(); ++i) {
                                unsigned count =
                                    lengths[i - 1] > 3 * typical
                                        ? std::min(16u,
                                                   unsigned(std::ceil(lengths[i - 1] / typical)))
                                        : 1;
                                for (unsigned j = 0; j < count; ++j)
                                    samples.push_back(points[i - 1] + (points[i] - points[i - 1]) *
                                                                          (double(j) / count));
                            }
                            if (!closed)
                                samples.push_back(points.back());
                            Handle(PointHandleArray) data = new PointHandleArray(1, samples.size());
                            for (unsigned i = 0; i < samples.size(); ++i)
                                data->SetValue(i + 1, point(samples[i]));
                            GeomAPI_Interpolate interpolate(data, closed, 1e-7);
                            interpolate.Perform();
                            if (!interpolate.IsDone())
                                return {};
                            curve = interpolate.Curve();
                            first = curve->FirstParameter();
                            last = curve->LastParameter();
                        }
                    }
                    auto corner = [&](unsigned id, const gp_Pnt &on) {
                        auto found = boundaries.vertices.find(id);
                        auto vertex =
                            found == boundaries.vertices.end()
                                ? BRepBuilderAPI_MakeVertex(point(mesh.points[id])).Vertex()
                                : found->second;
                        double gap = BRep_Tool::Pnt(vertex).Distance(on);
                        if (gap > tolerance)
                            throw std::runtime_error("Run endpoint off its curve");
                        BRep_Builder().UpdateVertex(vertex, std::max(1e-7, gap * 1.01));
                        return vertex;
                    };
                    auto a = corner(ids.front(), curve->Value(first)),
                         b = closed ? a : corner(ids.back(), curve->Value(last));
                    BRepBuilderAPI_MakeEdge make(curve, last < first ? b : a, last < first ? a : b,
                                                 std::min(first, last), std::max(first, last));
                    if (!make.IsDone())
                        return {};
                    edge = make.Edge();
                    if (last < first)
                        edge.Reverse();
                }
                boundaries.runs.emplace(ids, edge);
            } catch (const Standard_Failure &) {
                return {};
            }
        if (edge.IsNull())
            return {};
        if (flip)
            edge.Reverse();
        wire.Add(edge);
    }
    return wire.IsDone() ? wire.Wire() : TopoDS_Wire{};
}

std::vector<std::vector<unsigned>> unpinched_groups(const std::vector<Tri> &triangles) {
    // Python _unpinched: grow across shared edges without adding a triangle
    // that also touches the group at an isolated corner. Such a pinch cannot
    // form one valid trimming wire, but usually needs only two or three faces.
    std::map<Edge, std::vector<unsigned>> by_edge;
    for (unsigned i = 0; i < triangles.size(); ++i)
        for (unsigned k = 0; k < 3; ++k)
            by_edge[Edge(triangles[i][k], triangles[i][(k + 1) % 3])].push_back(i);
    std::vector<int> owner(triangles.size(), -1);
    std::vector<std::vector<unsigned>> groups;
    for (unsigned start = 0; start < triangles.size(); ++start) {
        if (owner[start] >= 0)
            continue;
        int g = int(groups.size());
        owner[start] = g;
        std::vector<unsigned> members{start}, stack{start};
        std::set<unsigned> vertices(triangles[start].begin(), triangles[start].end());
        while (!stack.empty()) {
            auto t = triangles[stack.back()];
            stack.pop_back();
            for (unsigned k = 0; k < 3; ++k) {
                Edge edge(t[k], t[(k + 1) % 3]);
                for (auto u : by_edge[edge]) {
                    if (owner[u] >= 0)
                        continue;
                    unsigned third = triangles[u][0];
                    for (auto v : triangles[u])
                        if (v != edge.a && v != edge.b)
                            third = v;
                    if (vertices.count(third)) {
                        bool linked = false;
                        for (auto v : {edge.a, edge.b}) {
                            auto found = by_edge.find(Edge(v, third));
                            if (found != by_edge.end())
                                for (auto w : found->second)
                                    linked = linked || (w != u && owner[w] == g);
                        }
                        if (!linked)
                            continue;
                    }
                    owner[u] = g;
                    members.push_back(u);
                    vertices.insert(third);
                    stack.push_back(u);
                }
            }
        }
        groups.push_back(std::move(members));
    }
    return groups;
}
std::vector<TopoDS_Face> planar_patch(const Mesh &mesh, const Patch &patch,
                                      const Boundaries &boundaries, double tolerance) {
    auto top = topology(mesh);
    std::vector<Tri> triangles;
    for (auto t : patch.triangles)
        triangles.push_back(mesh.triangles[t]);
    std::vector<TopoDS_Face> faces;
    for (const auto &local : unpinched_groups(triangles)) {
        std::vector<unsigned> group;
        for (auto t : local)
            group.push_back(patch.triangles[t]);
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
                auto boundary = boundaries.edges.find(Edge(current, to));
                if (boundary == boundaries.edges.end())
                    wire.Add(BRepBuilderAPI_MakeEdge(point(a), point(b)).Edge());
                else {
                    auto edge = boundary->second;
                    if (boundaries.separate(Edge(current, to)))
                        edge = TopoDS::Edge(BRepBuilderAPI_Copy(edge, true, false).Shape());
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
        // Match upstream _plane_face: repair rounding errors before sewing,
        // probing on a copy before a repair can modify shared corners/edges.
        if (!BRepCheck_Analyzer(face).IsValid()) {
            BRepBuilderAPI_Copy copy(face, true, false);
            ShapeFix_Face probe(TopoDS::Face(copy.Shape()));
            probe.SetMaxTolerance(tolerance);
            probe.Perform();
            ShapeMap before, after;
            TopExp::MapShapes(face, TopAbs_EDGE, before);
            TopExp::MapShapes(probe.Face(), TopAbs_EDGE, after);
            ShapeAnalysis_ShapeTolerance analysis;
            bool own = analysis.Tolerance(probe.Face(), 1) > std::min(.005, tolerance) ||
                       before.Extent() != after.Extent();
            if (own)
                face = probe.Face();
            else {
                ShapeFix_Face fix(face);
                fix.SetMaxTolerance(tolerance);
                fix.Perform();
                face = fix.Face();
            }
            if (!BRepCheck_Analyzer(face).IsValid()) {
                if (!own)
                    face = TopoDS::Face(BRepBuilderAPI_Copy(face, true, false).Shape());
                repair_planar_pcurves(face, tolerance);
                ShapeFix_Shape fix(face);
                double noise = std::max(1e-7, patch.deviation);
                fix.SetPrecision(std::min(noise, tolerance));
                fix.SetMaxTolerance(std::min(10 * noise, tolerance));
                fix.Perform();
                ShapeMap faces;
                TopExp::MapShapes(fix.Shape(), TopAbs_FACE, faces);
                if (faces.Extent() == 1)
                    face = TopoDS::Face(faces(1));
            }
        }
        faces.push_back(face);
    }
    return faces;
}
TopoDS_Face natural_patch(const Mesh &mesh, const Patch &patch, double tolerance) {
    if (patch.trimmed)
        return {};
    if (patch.reference_feature) {
        auto shape = reference_surface_face(patch, false);
        if (!shape.IsNull())
            for (TopExp_Explorer face(shape, TopAbs_FACE); face.More(); face.Next())
                return TopoDS::Face(face.Current());
    }
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
        auto normalized = [period](double t) {
            t = std::fmod(t, period);
            return t < 0 ? t + period : t;
        };
        std::vector<double> values;
        for (const auto &item : uv) {
            if (u && poles.count(item.first))
                continue;
            values.push_back(normalized(u ? item.second.X() : item.second.Y()));
        }
        if (values.empty())
            return;
        std::sort(values.begin(), values.end());
        double gap = -1, cut = values.front();
        for (unsigned i = 0; i < values.size(); ++i) {
            double next = values[(i + 1) % values.size()] + (i + 1 == values.size() ? period : 0);
            if (next - values[i] > gap) {
                gap = next - values[i];
                // Keep the actual sample: adding and removing a period can
                // round the cut above its endpoint and wrap it a full turn.
                cut = values[(i + 1) % values.size()];
            }
        }
        for (auto &item : uv) {
            double t = normalized(u ? item.second.X() : item.second.Y());
            if (t < cut)
                t += period;
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
    gp_Pnt middle;
    gp_Vec middle_du, middle_dv;
    patch.surface->D1((ulo + uhi) / 2, (vlo + vhi) / 2, middle, middle_du, middle_dv);
    // Surface parameters have different units (angle, length, or both).
    // A fixed UV epsilon rejects otherwise rectangular cylinder/cone rims
    // merely because STL rounding moves their axial coordinates by microns.
    const double rim_tolerance = std::max(1e-6, tolerance * .1);
    const double eps_u = rim_tolerance / std::max(middle_du.Magnitude(), 1e-12),
                 eps_v = rim_tolerance / std::max(middle_dv.Magnitude(), 1e-12);
    bool all_v = true, all_u = true;
    unsigned edges = 0;
    for (const auto &kv : rim)
        if (kv.second.count == 1) {
            ++edges;
            auto a = uv.at(kv.first.a), b = uv.at(kv.first.b);
            bool v = (std::abs(a.Y() - vlo) < eps_v && std::abs(b.Y() - vlo) < eps_v) ||
                     (std::abs(a.Y() - vhi) < eps_v && std::abs(b.Y() - vhi) < eps_v);
            bool u = ((poles.count(kv.first.a) || std::abs(a.X() - ulo) < eps_u) &&
                      (poles.count(kv.first.b) || std::abs(b.X() - ulo) < eps_u)) ||
                     ((poles.count(kv.first.a) || std::abs(a.X() - uhi) < eps_u) &&
                      (poles.count(kv.first.b) || std::abs(b.X() - uhi) < eps_u));
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
        double u, v;
        if (!projector.chart_parameters(mesh.points[id], u, v))
            return {};
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
                auto at = patch.surface->Value(uv[k].X(), uv[k].Y());
                auto shared = boundaries.vertices.find(t[k]);
                auto v = shared == boundaries.vertices.end()
                             ? BRepBuilderAPI_MakeVertex(at).Vertex()
                             : shared->second;
                double error = at.Distance(BRep_Tool::Pnt(v));
                if (error > tolerance)
                    return {};
                BRep_Builder().UpdateVertex(v, std::max(1e-7, error * 1.01));
                chart.push_back({uv[k], t[k], v});
            }
            mapped[k] = found.first->second;
        }
        triangles.push_back(mapped);
    }
    std::map<Edge, Rim> edges;
    std::vector<TopoDS_Face> faces;
    for (const auto &group : unpinched_groups(triangles)) {
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
                if ((patch.kind == "freeform" || patch.kind == "blend") &&
                    boundaries.singles.count(Edge(a.original, b.original)))
                    edge =
                        mesh_chord(patch.surface, mesh.points[a.original], mesh.points[b.original],
                                   a.uv, b.uv, a.shape, b.shape, tolerance);
                auto exact = boundaries.edges.find(Edge(a.original, b.original));
                if (edge.IsNull() && exact != boundaries.edges.end()) {
                    double first, last;
                    auto curve = BRep_Tool::Curve(exact->second, first, last);
                    double error = tolerance * .01;
                    Handle(Geom2d_Curve) pc;
                    try {
                        pc = GeomProjLib::Curve2d(curve, first, last, patch.surface, error);
                    } catch (const Standard_Failure &failure) {
                        // One ill-conditioned boundary projection must not
                        // discard a valid chart spanning hundreds of facets.
                        trace(failure_message(failure, "Shared trim projection failed"));
                    }
                    if (!pc.IsNull()) {
                        auto middle = pc->Value((first + last) / 2);
                        double du = pu ? std::round(((a.uv.X() + b.uv.X()) / 2 - middle.X()) / pu) *
                                             pu
                                       : 0,
                               dv = pv ? std::round(((a.uv.Y() + b.uv.Y()) / 2 - middle.Y()) / pv) *
                                             pv
                                       : 0;
                        pc->Translate(gp_Vec2d(du, dv));
                        bool reverse = (a.original > b.original) !=
                                       (exact->second.Orientation() == TopAbs_REVERSED);
                        if (boundaries.separate(Edge(a.original, b.original))) {
                            auto at_first = reverse ? b.shape : a.shape,
                                 at_last = reverse ? a.shape : b.shape;
                            double efirst = BRep_Tool::Pnt(at_first).Distance(curve->Value(first)),
                                   elast = BRep_Tool::Pnt(at_last).Distance(curve->Value(last));
                            if (std::max(efirst, elast) > tolerance)
                                return {};
                            builder.UpdateVertex(at_first, std::max(1e-7, efirst * 1.01));
                            builder.UpdateVertex(at_last, std::max(1e-7, elast * 1.01));
                            BRepBuilderAPI_MakeEdge make(curve, at_first, at_last, first, last);
                            if (make.IsDone()) {
                                edge = make.Edge();
                                if (reverse)
                                    edge.Reverse();
                            }
                        } else {
                            // Join both the edge and its corners before sewing.
                            // Sharing only the corners leaves partially joined
                            // topology when sewing later collapses a sliver.
                            edge = exact->second;
                            if (a.original > b.original)
                                edge.Reverse();
                        }
                        if (!edge.IsNull()) {
                            builder.UpdateEdge(
                                edge, pc, patch.surface, TopLoc_Location(),
                                std::max(error, BRep_Tool::Tolerance(exact->second)));
                            builder.Range(edge, first, last);
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
std::vector<TopoDS_Face> outline_patch(const Mesh &mesh, const Patch &patch,
                                       const Boundaries &boundaries, double tolerance,
                                       bool polygon = false) {
    // Python _outline_face/_polygon_face: lay a 3D outline on the support
    // and repair a private copy, without widening a neighbor's shared edges.
    auto rim = patch_rim(mesh, patch.triangles);
    std::map<unsigned, unsigned> next;
    for (const auto &item : rim)
        if (item.second.count == 1 && !next.emplace(item.second.from, item.second.to).second)
            return {};
    auto top = topology(mesh);
    Vec normal = Vec::Zero();
    double target = 0;
    for (auto t : patch.triangles) {
        normal += top.areas[t] * top.normals[t];
        target += top.areas[t];
    }
    struct Loop {
        TopoDS_Wire wire;
        double area;
    };
    std::vector<Loop> loops;
    while (!next.empty()) {
        unsigned start = next.begin()->first, current = start;
        BRepBuilderAPI_MakeWire wire;
        double area = 0;
        std::vector<unsigned> ordered;
        do {
            ordered.push_back(current);
            auto found = next.find(current);
            if (found == next.end())
                return {};
            unsigned to = found->second;
            next.erase(found);
            area += mesh.points[current].cross(mesh.points[to]).dot(normal);
            auto edge = boundaries.edges.find(Edge(current, to));
            if (!polygon && edge != boundaries.edges.end()) {
                auto e = edge->second;
                if (current > to)
                    e.Reverse();
                wire.Add(e);
            } else
                wire.Add(
                    BRepBuilderAPI_MakeEdge(point(mesh.points[current]), point(mesh.points[to]))
                        .Edge());
            current = to;
        } while (current != start);
        if (!wire.IsDone())
            return {};
        auto whole = polygon || patch.kind != "sphere"
                         ? TopoDS_Wire{}
                         : run_wire(mesh, patch, ordered, boundaries, tolerance);
        loops.push_back(
            {TopoDS::Wire(
                 BRepBuilderAPI_Copy(whole.IsNull() ? wire.Wire() : whole, true, false).Shape()),
             area});
    }
    if (loops.empty())
        return {};
    std::sort(loops.begin(), loops.end(),
              [](const Loop &a, const Loop &b) { return std::abs(a.area) > std::abs(b.area); });
    for (bool flip : {false, true}) {
        auto outer = loops.front().wire;
        if (flip)
            outer.Reverse();
        BRepBuilderAPI_MakeFace make(patch.surface, outer, false);
        for (std::size_t i = 1; i < loops.size(); ++i) {
            auto hole = loops[i].wire;
            if (flip)
                hole.Reverse();
            make.Add(hole);
        }
        if (!make.IsDone())
            continue;
        ShapeFix_Face fix(make.Face());
        fix.SetMaxTolerance(std::max(1., tolerance));
        fix.Perform();
        auto face = fix.Face();
        ShapeAnalysis_ShapeTolerance analysis;
        if (!BRepCheck_Analyzer(face).IsValid() ||
            analysis.Tolerance(face, 1) > std::max(1., tolerance))
            continue;
        GProp_GProps properties;
        BRepGProp::SurfaceProperties(face, properties);
        if (std::abs(properties.Mass() - target) > .05 * target + 1e-3)
            continue;
        PointProjector projector(patch.surface);
        const auto &p = projector.perform(top.centers[patch.triangles.front()]);
        if (!p.NbPoints())
            continue;
        double u, v;
        p.LowerDistanceParameters(u, v);
        gp_Pnt position;
        gp_Vec du, dv;
        patch.surface->D1(u, v, position, du, dv);
        auto n = du.Crossed(dv);
        double dot = normal.dot(Vec(n.X(), n.Y(), n.Z()));
        if (face.Orientation() == TopAbs_REVERSED)
            dot = -dot;
        if (dot < 0)
            face.Reverse();
        return {face};
    }
    return {};
}
std::vector<TopoDS_Face> split_patch(const Mesh &mesh, const Patch &patch,
                                     const Boundaries &boundaries, double tolerance) {
    // Python _split_face: cut the generous fitted support with shared 3D
    // outline edges, rather than requiring one manually stitched UV wire.
    // This also handles outlines that pinch into several surface pieces.
    if (Handle(Geom_BSplineSurface)::DownCast(patch.surface).IsNull())
        return {};
    auto rim = patch_rim(mesh, patch.triangles);
    std::set<unsigned> rim_vertices;
    ShapeList tools, arguments;
    for (const auto &item : rim)
        if (item.second.count == 1) {
            auto edge = boundaries.edges.find(item.first);
            if (edge != boundaries.edges.end())
                tools.Append(edge->second);
            else
                tools.Append(BRepBuilderAPI_MakeEdge(point(mesh.points[item.first.a]),
                                                     point(mesh.points[item.first.b]))
                                 .Edge());
            rim_vertices.insert(item.first.a);
            rim_vertices.insert(item.first.b);
        }
    if (tools.IsEmpty())
        return {};
    arguments.Append(BRepBuilderAPI_MakeFace(patch.surface, tolerance * .01).Face());
    BRepAlgoAPI_Splitter splitter;
    splitter.SetArguments(arguments);
    splitter.SetTools(tools);
    splitter.SetNonDestructive(true);
    splitter.SetFuzzyValue(std::max(1e-3, 2 * patch.deviation));
    splitter.Build();
    if (!splitter.IsDone() || splitter.Shape().IsNull())
        return {};
    auto top = topology(mesh);
    double target = 0;
    std::vector<unsigned> inner;
    for (auto t : patch.triangles) {
        target += top.areas[t];
        const auto &tri = mesh.triangles[t];
        if (!rim_vertices.count(tri[0]) && !rim_vertices.count(tri[1]) &&
            !rim_vertices.count(tri[2]))
            inner.push_back(t);
    }
    PointProjector projector(patch.surface);
    struct Piece {
        TopoDS_Face face;
        double area;
    };
    auto pick = [&](const std::vector<unsigned> &samples) {
        std::vector<gp_Pnt2d> probes;
        for (unsigned i = 0; i < std::min<std::size_t>(60, samples.size()); ++i) {
            auto t = samples[i * samples.size() / std::min<std::size_t>(60, samples.size())];
            const auto &p = projector.perform(top.centers[t]);
            if (p.NbPoints()) {
                double u, v;
                p.LowerDistanceParameters(u, v);
                probes.emplace_back(u, v);
            }
        }
        std::vector<Piece> pieces;
        for (TopExp_Explorer ex(splitter.Shape(), TopAbs_FACE); ex.More(); ex.Next()) {
            auto face = TopoDS::Face(ex.Current());
            bool hit = false;
            for (const auto &uv : probes) {
                BRepClass_FaceClassifier classifier(face, uv, 1e-7);
                if (classifier.State() == TopAbs_IN || classifier.State() == TopAbs_ON) {
                    hit = true;
                    break;
                }
            }
            if (hit) {
                GProp_GProps properties;
                BRepGProp::SurfaceProperties(face, properties);
                pieces.push_back({face, properties.Mass()});
            }
        }
        return pieces;
    };
    auto pieces = pick(inner.empty() ? patch.triangles : inner);
    auto total_area = [](const auto &items) {
        double area = 0;
        for (const auto &p : items)
            area += p.area;
        return area;
    };
    if (!inner.empty() && std::abs(total_area(pieces) - target) > .1 * target)
        pieces = pick(patch.triangles);
    if (pieces.empty() || pieces.size() > 10)
        return {};
    unsigned best = 0;
    double error = .1 * target + 1e-3;
    for (unsigned mask = 1; mask < (1u << pieces.size()); ++mask) {
        double area = 0;
        for (unsigned i = 0; i < pieces.size(); ++i)
            if (mask & (1u << i))
                area += pieces[i].area;
        if (std::abs(area - target) <= error) {
            error = std::abs(area - target);
            best = mask;
        }
    }
    std::vector<TopoDS_Face> faces;
    ShapeAnalysis_ShapeTolerance analysis;
    for (unsigned i = 0; i < pieces.size(); ++i)
        if (best & (1u << i)) {
            auto face = pieces[i].face;
            if (!BRepCheck_Analyzer(face).IsValid()) {
                face = TopoDS::Face(BRepBuilderAPI_Copy(face, true, false).Shape());
                ShapeFix_Face fix(face);
                fix.SetMaxTolerance(tolerance);
                fix.Perform();
                face = fix.Face();
            }
            if (!BRepCheck_Analyzer(face).IsValid() || analysis.Tolerance(face, 1) > tolerance)
                return {};
            const auto &p = projector.perform(top.centers[patch.triangles.front()]);
            if (!p.NbPoints())
                return {};
            double u, v;
            p.LowerDistanceParameters(u, v);
            gp_Pnt position;
            gp_Vec du, dv;
            patch.surface->D1(u, v, position, du, dv);
            auto normal = du.Crossed(dv);
            double dot =
                top.normals[patch.triangles.front()].dot(Vec(normal.X(), normal.Y(), normal.Z()));
            if (face.Orientation() == TopAbs_REVERSED)
                dot = -dot;
            if (dot < 0)
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
    // Like upstream SHARED_CROWD: let sewing own close corners rather than
    // pre-joining only some of the faces meeting a short edge.
    // This is a fixed mesh-space distance, independent of adaptive sewing:
    // tying it to sewing would expand the exclusion from 0.02 to 0.4 mm.
    constexpr double shared_crowd = .02;
    BoundaryCloud cloud{mesh.points};
    using Tree =
        nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, BoundaryCloud>,
                                            BoundaryCloud, 3, std::size_t>;
    Tree tree(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(16));
    tree.buildIndex();
    out.crowded.resize(mesh.points.size());
    for (unsigned i = 0; i < mesh.points.size(); ++i) {
        std::size_t nearest[2];
        double distances[2];
        out.crowded[i] = tree.knnSearch(mesh.points[i].data(), 2, nearest, distances) == 2 &&
                         distances[1] <= shared_crowd * shared_crowd;
    }
    auto corner = [&](unsigned id, const gp_Pnt &on_curve) {
        if (out.crowded[id]) {
            auto vertex = BRepBuilderAPI_MakeVertex(point(mesh.points[id])).Vertex();
            BRep_Builder().UpdateVertex(
                vertex, std::max(1e-7, on_curve.Distance(point(mesh.points[id])) * 1.01));
            return vertex;
        }
        auto found = out.vertices.emplace(id, TopoDS_Vertex{});
        if (found.second)
            found.first->second = BRepBuilderAPI_MakeVertex(point(mesh.points[id])).Vertex();
        BRep_Builder().UpdateVertex(
            found.first->second, std::max(1e-7, on_curve.Distance(point(mesh.points[id])) * 1.01));
        return found.first->second;
    };
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
            auto a_vertex = corner(kv.first.a, selected->Value(pa)),
                 b_vertex = corner(kv.first.b, selected->Value(pb));
            BRepBuilderAPI_MakeEdge make(selected, pb < pa ? b_vertex : a_vertex,
                                         pb < pa ? a_vertex : b_vertex, std::min(pa, pb),
                                         std::max(pa, pb));
            if (!make.IsDone())
                continue;
            auto edge = make.Edge();
            if (pb < pa)
                edge.Reverse();
            out.edges.emplace(kv.first, edge);
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
        double u, v;
        if (!projector.chart_parameters(mesh.points[tri[k]], u, v))
            throw std::runtime_error("Surface projection failed");
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
        auto boundary = boundaries.edges.find(Edge(tri[k], tri[(k + 1) % 3]));
        if (boundary != boundaries.edges.end()) {
            auto edge = boundary->second;
            if (boundaries.separate(Edge(tri[k], tri[(k + 1) % 3])))
                edge = TopoDS::Edge(BRepBuilderAPI_Copy(edge, true, false).Shape());
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
        curved = curved || boundaries.edges.count(Edge(tri[k], tri[(k + 1) % 3]));
    if (curved) {
        BRepBuilderAPI_MakeWire wire;
        for (unsigned k = 0; k < 3; ++k) {
            auto found = boundaries.edges.find(Edge(tri[k], tri[(k + 1) % 3]));
            if (found == boundaries.edges.end())
                wire.Add(BRepBuilderAPI_MakeEdge(point(mesh.points[tri[k]]),
                                                 point(mesh.points[tri[(k + 1) % 3]]))
                             .Edge());
            else {
                auto edge = found->second;
                if (boundaries.separate(Edge(tri[k], tri[(k + 1) % 3])))
                    edge = TopoDS::Edge(BRepBuilderAPI_Copy(edge, true, false).Shape());
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
std::vector<TopoDS_Shape> bad_merge_edges(const TopoDS_Shape &source,
                                          const ShapeUpgrade_UnifySameDomain &unify) {
    // Python _bad_merges: keep just the seams of invalid merged faces,
    // allowing independent healthy merges to survive the next attempt.
    auto history = unify.History();
    if (history.IsNull())
        return {};
    ShapeMap faces;
    TopExp::MapShapes(unify.Shape(), TopAbs_FACE, faces);
    BRepCheck_Analyzer check(unify.Shape());
    std::map<int, std::vector<TopoDS_Shape>> sources;
    for (TopExp_Explorer ex(source, TopAbs_FACE); ex.More(); ex.Next())
        for (const auto &image : history->Modified(ex.Current())) {
            int i = faces.FindIndex(image);
            if (i)
                sources[i].push_back(ex.Current());
        }
    std::vector<TopoDS_Shape> apart;
    for (const auto &item : sources) {
        if (item.second.size() < 2 || check.IsValid(faces(item.first)))
            continue;
        ShapeMap seen, repeated;
        for (const auto &face : item.second)
            for (TopExp_Explorer ex(face, TopAbs_EDGE); ex.More(); ex.Next()) {
                if (seen.Contains(ex.Current()))
                    repeated.Add(ex.Current());
                else
                    seen.Add(ex.Current());
            }
        for (int i = 1; i <= repeated.Extent(); ++i)
            apart.push_back(repeated(i));
    }
    return apart;
}
} // namespace
TopoDS_Shape tidy_shape(const TopoDS_Shape &shape, const Mesh &mesh,
                        const std::vector<Patch> &patches, double tolerance) {
    // Copy topology while keeping surface identity. Copying geometry makes
    // one shared B-spline support into a different surface for each triangle,
    // so UnifySameDomain can no longer recognize the faces as one domain.
    BRepBuilderAPI_Copy copy(shape, false, false);
    auto checked = copy.Shape();
    const double target = volume(shape), allowed = std::abs(target) * 1e-3 + 1e-3;
    // Match Python _tidy: a failed face merge does not discard all the
    // successful ones. Last, try joining split edges without merging faces.
    for (bool merge_faces : {true, false}) {
        ShapeMap apart;
        for (unsigned attempt = 0; attempt < (merge_faces ? 3u : 1u); ++attempt) {
            try {
                ShapeUpgrade_UnifySameDomain merge(checked, true, merge_faces, false);
                merge.SetLinearTolerance(std::min(.001, tolerance));
                merge.SetAngularTolerance(.0017453292519943296);
                for (int i = 1; i <= apart.Extent(); ++i)
                    merge.KeepShape(apart(i));
                merge.Build();
                auto reduced = merge.Shape();
                if (reduced.IsNull())
                    break;
                ShapeFix_Shape fix(reduced);
                fix.SetPrecision(std::min(1e-7, tolerance));
                fix.SetMaxTolerance(tolerance);
                fix.Perform();
                reduced = fix.Shape();
                Diagnostics info;
                ShapeAnalysis_ShapeTolerance analysis;
                double delta = reduced.IsNull() ? 1e100 : volume(reduced) - target,
                       widest = reduced.IsNull() ? 1e100 : analysis.Tolerance(reduced, 1);
                trace(("same-domain candidate volume_delta=" + std::to_string(delta) +
                       " max_tolerance=" + std::to_string(widest))
                          .c_str());
                if (!reduced.IsNull() && std::abs(delta) <= allowed &&
                    widest <= std::max(1., tolerance) && volume_agrees(reduced, mesh, patches) &&
                    inspect(reduced, info))
                    return reduced;
                auto more = bad_merge_edges(checked, merge);
                if (!merge_faces || more.empty())
                    break;
                for (const auto &edge : more)
                    apart.Add(edge);
            } catch (const Standard_Failure &) {
                break;
            }
        }
    }
    return shape;
}
TopoDS_Shape join_boundary_runs(const TopoDS_Shape &shape, const Mesh &mesh,
                                const std::vector<Patch> &patches, double tolerance) {
    // The reference facet path already builds one edge per reference run.
    // Concatenating those runs changes Python's deliberately retained seams.
    if (std::any_of(patches.begin(), patches.end(),
                    [](const Patch &p) { return bool(p.reference_feature); }))
        return shape;
    if (std::none_of(patches.begin(), patches.end(), [](const Patch &p) {
            return p.kind == "freeform" || p.kind == "blend" || p.kind == "pipe" ||
                   p.kind == "extrusion";
        }))
        return shape;
    // Python _wire makes one edge per shared run. The native UV-trim fallback
    // can represent that same curve as several mesh-sized spans. Join its C1
    // spline spans only after reconstruction has passed: doing this during
    // face merging changes which reconstruction attempt wins the volume check.
    try {
        Diagnostics before, after;
        if (!inspect(shape, before))
            return shape;
        BRepBuilderAPI_Copy copy(shape, false, false);
        auto working = copy.Shape();
        ShapeMap vertices, source_edges;
        TopExp::MapShapes(working, TopAbs_VERTEX, vertices);
        TopExp::MapShapes(working, TopAbs_EDGE, source_edges);
        std::vector<double> vertex_tolerances;
        const double joining_tolerance = std::max(1e-7, std::min(.005, tolerance * .5));
        ShapeFix_ShapeTolerance set_tolerance;
        for (int i = 1; i <= vertices.Extent(); ++i) {
            auto vertex = TopoDS::Vertex(vertices(i));
            double original = BRep_Tool::Tolerance(vertex);
            vertex_tolerances.push_back(original);
            // ConcatC1 uses vertex tolerances for knot removal. Sewing gaps
            // must not become permission to reshape an interpolated rim.
            set_tolerance.SetTolerance(vertex, joining_tolerance, TopAbs_VERTEX);
        }
        ShapeUpgrade_UnifySameDomain merge(working, true, false, true);
        merge.SetSafeInputMode(true);
        merge.SetLinearTolerance(std::min(.001, tolerance));
        merge.SetAngularTolerance(.0017453292519943296);
        merge.Build();
        auto joined = merge.Shape();
        BRep_Builder builder;
        auto history = merge.History();
        for (int i = 1; i <= vertices.Extent(); ++i) {
            double original = vertex_tolerances[i - 1];
            builder.UpdateVertex(TopoDS::Vertex(vertices(i)), original);
            for (const auto &replacement : history->Modified(vertices(i)))
                if (replacement.ShapeType() == TopAbs_VERTEX)
                    builder.UpdateVertex(TopoDS::Vertex(replacement), original);
        }
        double displacement = 0;
        for (int i = 1; i <= source_edges.Extent(); ++i) {
            auto edge = TopoDS::Edge(source_edges(i));
            double first, last;
            auto curve = BRep_Tool::Curve(edge, first, last);
            if (curve.IsNull() || BRep_Tool::Degenerated(edge))
                continue;
            for (const auto &replacement : history->Modified(edge)) {
                if (replacement.ShapeType() != TopAbs_EDGE)
                    continue;
                double a, b;
                auto target = BRep_Tool::Curve(TopoDS::Edge(replacement), a, b);
                if (target.IsNull())
                    continue;
                for (unsigned j = 0; j <= 16; ++j) {
                    GeomAPI_ProjectPointOnCurve project(
                        curve->Value(first + (last - first) * (double(j) / 16)), target, a, b);
                    displacement = std::max(displacement,
                                            project.NbPoints() ? project.LowerDistance() : 1e100);
                }
            }
        }
        ShapeAnalysis_ShapeTolerance analysis;
        const double widest = joined.IsNull() ? 1e100 : analysis.Tolerance(joined, 1);
        const double previous_widest = analysis.Tolerance(shape, 1);
        bool accepted =
            inspect(joined, after) && after.faces == before.faces &&
            after.solids == before.solids && after.surface_counts == before.surface_counts &&
            displacement <= 2 * joining_tolerance &&
            widest <= std::max(previous_widest, tolerance) &&
            std::abs(after.volume - before.volume) <= .001 * std::abs(before.volume) + .001 &&
            volume_agrees(joined, mesh, patches);
        trace(("boundary run joining accepted=" + std::to_string(accepted) +
               " volume_delta=" + std::to_string(after.volume - before.volume) + " displacement=" +
               std::to_string(displacement) + " max_tolerance=" + std::to_string(widest))
                  .c_str());
        if (accepted)
            return joined;
    } catch (const Standard_Failure &e) {
        trace(failure_message(e, "Boundary run joining failed"));
    }
    return shape;
}
namespace {
TopoDS_Shape build(const Mesh &mesh, const std::vector<Patch> &patches, double tolerance,
                   bool unify, unsigned boundary_mode = 2, std::set<unsigned> *failed = nullptr) {
    BRepBuilderAPI_Sewing sewing(tolerance, true, true, true, false);
    std::vector<TopoDS_Shape> faces;
    auto add = [&](const TopoDS_Shape &face) { faces.push_back(face); };
    Boundaries boundaries;
    SurfaceOwners owners;
    if (!patches.empty() && boundary_mode)
        try {
            boundaries = shared_boundaries(mesh, patches, tolerance, boundary_mode == 2);
        } catch (const Standard_Failure &error) {
            throw std::runtime_error(std::string("Boundaries: ") + failure_message(error));
        }
    if (!patches.empty()) {
        boundaries.owners.assign(mesh.triangles.size(), ~0u);
        for (unsigned i = 0; i < patches.size(); ++i) {
            boundaries.curved.push_back(patches[i].kind != "plane");
            for (auto t : patches[i].triangles)
                boundaries.owners[t] = i;
        }
        for (const auto &item : topology(mesh).edges) {
            auto a = boundaries.owners[item.second[0]],
                 b = item.second.size() > 1 ? boundaries.owners[item.second[1]] : ~0u;
            boundaries.sides.emplace(item.first, std::make_pair(a, b));
        }
    }
    if (!patches.empty()) {
        const auto &owner = boundaries.owners;
        std::map<std::pair<unsigned, unsigned>, std::vector<Edge>> groups;
        for (const auto &item : topology(mesh).edges)
            if (item.second.size() == 2) {
                auto a = owner[item.second[0]], b = owner[item.second[1]];
                if (a != b && a < patches.size() && b < patches.size())
                    groups[std::minmax(a, b)].push_back(item.first);
            }
        for (const auto &group : groups) {
            std::map<unsigned, std::vector<unsigned>> adjacent;
            for (auto e : group.second) {
                adjacent[e.a].push_back(e.b);
                adjacent[e.b].push_back(e.a);
            }
            auto corner = [&](unsigned v) {
                auto &neighbors = adjacent[v];
                if (neighbors.size() != 2)
                    return true;
                Vec a = mesh.points[v] - mesh.points[neighbors[0]],
                    b = mesh.points[neighbors[1]] - mesh.points[v];
                return a.dot(b) < .8660254037844386 * a.norm() * b.norm();
            };
            if (patches[group.first.first].kind == "sphere" ||
                patches[group.first.second].kind == "sphere")
                for (auto e : group.second)
                    if (corner(e.a) && corner(e.b))
                        boundaries.singles.insert(e);
        }
    }
    if (patches.empty())
        for (const auto &t : mesh.triangles)
            add(flat_triangle(mesh, t));
    else
        for (const auto &p : patches) {
            try {
                if (p.kind == "plane") {
                    for (const auto &face : planar_patch(mesh, p, boundaries, tolerance)) {
                        owners[BRep_Tool::Surface(face).get()] = unsigned(&p - patches.data());
                        add(face);
                    }
                    continue;
                }
                auto natural = natural_patch(mesh, p, tolerance);
                if (!natural.IsNull()) {
                    add(natural);
                    continue;
                }
                try {
                    // Prefer the bounded outline work in-process. Python also
                    // uses this path for complex supports; its worker timeout
                    // for expensive kernel cuts has no C++ call-level equivalent.
                    std::vector<TopoDS_Face> trimmed;
                    try {
                        if (p.kind == "sphere")
                            trimmed = outline_patch(mesh, p, boundaries, tolerance);
                        if (trimmed.empty())
                            trimmed = trimmed_patch(mesh, p, tolerance, boundaries);
                    } catch (const Standard_Failure &failure) {
                        trace(failure_message(failure, "Trim chart failed"));
                    }
                    for (auto &face : trimmed)
                        if (!BRepCheck_Analyzer(face).IsValid()) {
                            // Python _outline_face repairs the face before testing
                            // it. In particular, nested UV wires on a fitted
                            // support can need orientation repair even though
                            // every edge and parameter curve is already valid.
                            auto own =
                                TopoDS::Face(BRepBuilderAPI_Copy(face, false, false).Shape());
                            BRepLib::SameParameter(own, std::min(1e-5, tolerance * .001), true);
                            ShapeFix_Face fix(own);
                            fix.SetPrecision(std::min(1e-7, tolerance));
                            fix.SetMaxTolerance(tolerance);
                            fix.Perform();
                            ShapeAnalysis_ShapeTolerance analysis;
                            if (BRepCheck_Analyzer(fix.Face()).IsValid() &&
                                analysis.Tolerance(fix.Face(), 1) <= tolerance)
                                face = fix.Face();
                        }
                    if ((p.kind == "freeform" || p.kind == "blend") &&
                        (trimmed.empty() ||
                         !std::all_of(trimmed.begin(), trimmed.end(), [](const TopoDS_Face &face) {
                             return BRepCheck_Analyzer(face).IsValid();
                         }))) {
                        for (bool polygon : {false, true}) {
                            std::vector<TopoDS_Face> outline;
                            try {
                                outline = outline_patch(mesh, p, boundaries, tolerance, polygon);
                            } catch (const Standard_Failure &) { /* Try the next outline. */
                            }
                            if (!outline.empty()) {
                                trimmed = std::move(outline);
                                break;
                            }
                        }
                    }
                    auto valid_trims = [&] {
                        return !trimmed.empty() &&
                               std::all_of(trimmed.begin(), trimmed.end(),
                                           [](const TopoDS_Face &face) {
                                               return BRepCheck_Analyzer(face).IsValid();
                                           });
                    };
                    if (!valid_trims()) {
                        trimmed.clear();
                        auto spline = Handle(Geom_BSplineSurface)::DownCast(p.surface);
                        if (!spline.IsNull() && spline->NbUPoles() * spline->NbVPoles() <= 64 &&
                            p.triangles.size() < 300) {
                            try {
                                trimmed = split_patch(mesh, p, boundaries, tolerance);
                            } catch (const Standard_Failure &) { /* Keep the triangulated chart. */
                            }
                        }
                    }
                    if (!trimmed.empty()) {
                        for (const auto &face : trimmed) {
                            owners[BRep_Tool::Surface(face).get()] = unsigned(&p - patches.data());
                            add(face);
                        }
                        continue;
                    }
                } catch (const Standard_Failure
                             &) { /* Retain triangulated trims for an unsuitable chart. */
                }
                trace(("triangulated trim fallback patch=" + std::to_string(&p - patches.data()) +
                       " kind=" + p.kind + " triangles=" + std::to_string(p.triangles.size()))
                          .c_str());
                // Closed natural surfaces require no artificial triangulation seams or
                // parameter-space triangles at their singular poles.
                if ((p.kind == "sphere" || p.kind == "torus") &&
                    p.triangles.size() == mesh.triangles.size()) {
                    add(BRepBuilderAPI_MakeFace(p.surface, tolerance).Face());
                    continue;
                }
                PointProjector projector(p.surface);
                for (auto id : p.triangles) {
                    const auto &t = mesh.triangles[id];
                    // Planes keep original corners and exact shared straight edges.
                    try {
                        if (p.kind == "plane")
                            add(flat_triangle(mesh, t, boundaries));
                        else
                            add(curved_triangle(mesh, t, p.surface, tolerance, boundaries,
                                                projector));
                    } catch (const Standard_Failure &error) {
                        throw std::runtime_error(p.kind + " face " + std::to_string(id) + ": " +
                                                 failure_message(error));
                    }
                }
            } catch (const Standard_Failure &) {
                if (failed)
                    failed->insert(unsigned(&p - patches.data()));
                throw;
            } catch (const std::runtime_error &) {
                if (failed)
                    failed->insert(unsigned(&p - patches.data()));
                throw;
            }
        }
    sewing.Add(patches.empty() ? compound(faces) : rejoin_exact(compound(faces), mesh, tolerance));
    trace("sewing");
    try {
        sewing.Perform();
    } catch (const Standard_Failure &error) {
        throw std::runtime_error(std::string("Sewing: ") + failure_message(error));
    }
    auto sewn = without_collapsed_edges(sewing.SewedShape(), tolerance);
    if (sewn.IsNull())
        throw std::runtime_error("Sewing produced an empty shape");
    TopoDS_Shape shape;
    trace("shell closure");
    try {
        shape = close_shells(sewn, tolerance);
        BRepLib::SameParameter(shape, std::min(1e-5, tolerance * .001), true);
    } catch (const Standard_Failure &error) {
        if (failed)
            failed_patches(sewn, patches, tolerance, *failed, owners);
        throw std::runtime_error(std::string("Shell closure: ") + failure_message(error));
    } catch (const std::runtime_error &) {
        if (failed)
            failed_patches(sewn, patches, tolerance, *failed, owners);
        throw;
    }
    if (!patches.empty() && boundary_mode)
        trace_invalid_faces(shape, &patches);
    if (!patches.empty() && boundary_mode == 2 && !BRepCheck_Analyzer(shape).IsValid()) {
        auto fixed = heal_candidate(shape, mesh, patches, tolerance);
        if (fixed.IsNull()) {
            if (failed)
                failed_patches(shape, patches, tolerance, *failed, owners);
            throw std::runtime_error("Intersection-boundary topology is invalid");
        }
        shape = fixed;
    }
    trace("same-domain reduction");
    if (unify)
        shape = tidy_shape(shape, mesh, patches, tolerance);
    if (failed)
        failed_patches(shape, patches, tolerance, *failed, owners);
    return shape;
}
} // namespace
TopoDS_Shape faceted(const Mesh &mesh, double tolerance) {
    return build(mesh, {}, tolerance, true);
}
namespace {
double patch_sewing_tolerance(const Mesh &mesh, const std::vector<Patch> &patches, double minimum) {
    // Python attempt()/sewing_cap(): allow the measured sag of a fitted
    // patch, bounded by 0.2 mm and twice the mesh's median edge length.
    double worst = 0;
    for (const auto &patch : patches) {
        if (patch.kind == "plane")
            continue;
        worst = std::max(worst, patch.deviation);
        PointProjector projector(patch.surface);
        for (auto t : patch.triangles) {
            const auto &tri = mesh.triangles[t];
            for (unsigned k = 0; k < 3; ++k) {
                Vec p = (mesh.points[tri[k]] + mesh.points[tri[(k + 1) % 3]]) / 2;
                const auto &projection = projector.perform(p);
                if (projection.NbPoints())
                    worst = std::max(worst, projection.LowerDistance());
            }
        }
    }
    if (1.5 * worst <= minimum)
        return minimum;
    std::vector<double> lengths;
    lengths.reserve(3 * mesh.triangles.size());
    for (const auto &tri : mesh.triangles)
        for (unsigned k = 0; k < 3; ++k)
            lengths.push_back((mesh.points[tri[k]] - mesh.points[tri[(k + 1) % 3]]).norm());
    if (lengths.empty())
        return minimum;
    auto middle = lengths.begin() + lengths.size() / 2;
    std::nth_element(lengths.begin(), middle, lengths.end());
    double median = *middle;
    if (lengths.size() % 2 == 0)
        median = (median + *std::max_element(lengths.begin(), middle)) / 2;
    double cap = std::max(minimum, std::min(.2, 2 * median));
    return std::max(minimum, std::min(1.5 * worst, cap));
}
} // namespace
#include "build_boundary_runs.hpp"
TopoDS_Shape reconstruct(const Mesh &mesh, const std::vector<Patch> &patches, Context &ctx) {
    // Python's default uses shared runs of original mesh points on fitted
    // outlines. Keep the established analytic path for entirely analytic
    // models; use its alternatives if a fitted outline cannot build safely.
    if (std::any_of(patches.begin(), patches.end(), [](const Patch &p) {
            return p.kind == "freeform" || p.kind == "blend" || p.kind == "pipe" ||
                   p.kind == "extrusion";
        })) {
        try {
            auto shape = mesh_boundary_shape(mesh, patches, ctx);
            Diagnostics info;
            if (inspect(shape, info) && volume_agrees(shape, mesh, patches))
                return shape;
        } catch (const Standard_Failure &error) {
            trace(failure_message(error, "Mesh-boundary reconstruction failed"));
        } catch (const std::runtime_error &error) {
            trace(error.what());
        }
    }

    // Try exact analytic intersections, then smooth fitted junctions with
    // projected analytic trims, then independent projected trims. Every mode
    // must have valid topology; conversion applies the corrected-volume gate
    // before accepting a candidate or invoking faceted recovery.
    double adaptive = patch_sewing_tolerance(mesh, patches, ctx.options.tolerance);
    std::vector<double> tolerances{ctx.options.tolerance};
    if (adaptive > ctx.options.tolerance)
        tolerances.push_back(adaptive);
    for (unsigned mode : {2u, 1u, 0u})
        for (double tolerance : tolerances) {
            try {
                trace(("reconstruction sewing tolerance=" + std::to_string(tolerance)).c_str());
                auto shape = build(mesh, patches, tolerance, true, mode);
                Diagnostics info;
                if (inspect(shape, info) && (mode == 0 || volume_agrees(shape, mesh, patches)))
                    return shape;
                auto fixed = heal_candidate(shape, mesh, patches, tolerance);
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
    std::set<unsigned> faceted_planes;
    std::vector<unsigned> owner(mesh.triangles.size());
    for (unsigned i = 0; i < patches.size(); ++i)
        for (auto t : patches[i].triangles)
            owner[t] = i;
    std::vector<std::set<unsigned>> adjacent(patches.size());
    for (const auto &edge : top.edges)
        if (edge.second.size() == 2) {
            auto a = owner[edge.second[0]], b = owner[edge.second[1]];
            if (a != b) {
                adjacent[a].insert(b);
                adjacent[b].insert(a);
            }
        }
    auto assemble = [&](const std::set<unsigned> &selected,
                        std::vector<unsigned> *sources = nullptr) {
        std::vector<Patch> candidate;
        for (unsigned i = 0; i < patches.size(); ++i) {
            if ((patches[i].kind == "plane" && !faceted_planes.count(i)) || selected.count(i)) {
                candidate.push_back(patches[i]);
                if (sources)
                    sources->push_back(i);
            } else
                for (auto t : patches[i].triangles) {
                    Patch p;
                    p.kind = "plane";
                    p.triangles = {t};
                    p.origin = mesh.points[mesh.triangles[t][0]];
                    p.axis = top.normals[t];
                    p.surface = new Geom_Plane(
                        gp_Ax3(point(p.origin), gp_Dir(p.axis.x(), p.axis.y(), p.axis.z())));
                    candidate.push_back(std::move(p));
                    if (sources)
                        sources->push_back(i);
                }
        }
        return candidate;
    };
    // Keep healthy patches together, dropping only the supports implicated by
    // invalid faces or open seams. A failed fitted plane also needs to revert
    // to its actual triangles; retaining that bad plane poisons every subset.
    std::set<unsigned> selected(curved.begin(), curved.end());
    for (unsigned round = 0; round < 8 && !selected.empty() && !ctx.expired(); ++round) {
        std::vector<unsigned> sources;
        auto candidate = assemble(selected, &sources);
        std::set<unsigned> failed;
        ++attempts;
        try {
            double tolerance = patch_sewing_tolerance(mesh, candidate, ctx.options.tolerance);
            auto shape = build(mesh, candidate, tolerance, true, 2, &failed);
            Diagnostics info;
            if (failed.empty() && inspect(shape, info) && volume_agrees(shape, mesh, candidate)) {
                best = shape;
                kept = selected;
                break;
            }
        } catch (const Standard_Failure &) {
        } catch (const std::runtime_error &) {
        }
        bool changed = false;
        std::set<unsigned> failed_curves;
        std::set<unsigned> failed_planes, failed_blends;
        for (auto index : failed) {
            auto original = sources.at(index);
            if (patches[original].kind == "plane")
                failed_planes.insert(original);
            else
                failed_curves.insert(original);
            // Upstream culprits gives nearby blends the first chance to revert:
            // their overshooting junctions can invalidate healthy neighboring
            // planes and analytic surfaces. Dropping both sides loses curves.
            auto blame_blend = [&](unsigned i) {
                if (selected.count(i) && patches[i].kind == "blend")
                    failed_blends.insert(i);
            };
            blame_blend(original);
            for (auto i : adjacent[original])
                blame_blend(i);
        }
        for (auto original : failed_blends)
            changed = selected.erase(original) != 0 || changed;
        if (!changed)
            for (auto original : failed_planes)
                if (patches[original].triangles.size() > 1)
                    changed = faceted_planes.insert(original).second || changed;
        if (!changed)
            for (auto original : failed_curves)
                changed = selected.erase(original) != 0 || changed;
        if (!changed)
            break;
    }
    std::function<void(const std::vector<unsigned> &)> restore;
    restore = [&](const std::vector<unsigned> &group) {
        if (group.empty() || attempts >= 32 || ctx.expired())
            return;
        ++attempts;
        auto selected = kept;
        selected.insert(group.begin(), group.end());
        std::vector<unsigned> sources;
        auto candidate = assemble(selected, &sources);
        std::set<unsigned> failed;
        try {
            // Restored patches must also restore their boundaries on adjacent
            // planar faces. Independent trims leave mesh chords facing arcs
            // and can pass sewing while changing the reconstructed volume.
            double tolerance = patch_sewing_tolerance(mesh, candidate, ctx.options.tolerance);
            auto shape = build(mesh, candidate, tolerance, true, 2, &failed);
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
            // Python again(): retry the unblamed patches together before
            // splitting blindly, then give each named patch its own chance.
            std::set<unsigned> named;
            for (auto index : failed)
                named.insert(sources.at(index));
            std::vector<unsigned> rest;
            for (auto i : group)
                if (!named.count(i))
                    rest.push_back(i);
            if (!rest.empty() && rest.size() < group.size()) {
                restore(rest);
                for (auto i : group)
                    if (named.count(i))
                        restore({i});
                return;
            }
            auto middle = group.begin() + group.size() / 2;
            restore(std::vector<unsigned>(group.begin(), middle));
            restore(std::vector<unsigned>(middle, group.end()));
        }
    };
    std::vector<unsigned> dropped;
    for (auto i : curved)
        if (!kept.count(i))
            dropped.push_back(i);
    restore(dropped);
    if (!faceted_planes.empty())
        ctx.diagnostics.fallbacks.push_back("Rebuilt " + std::to_string(faceted_planes.size()) +
                                            " invalid planar patches from their triangles");
    ctx.diagnostics.fallbacks.push_back("Recovery retained " + std::to_string(kept.size()) +
                                        " of " + std::to_string(curved.size()) + " curved patches");
    return best;
}
TopoDS_Shape unite(const std::vector<TopoDS_Shape> &shapes, Context &ctx, bool overlapping_parts) {
    if (shapes.empty())
        throw std::runtime_error("No shapes to combine");
    if (!ctx.options.fuse)
        return compound(shapes);
    std::vector<TopoDS_Shape> checked = shapes;
    if (shapes.size() > 1)
        for (auto &s : checked)
            // finish() first calls fuse_overlapping inside each reconstructed
            // group. Only the subsequent join of distinct STL groups must
            // preserve their summed volume; 3MF parts may overlap at that join.
            s = unite({s}, ctx, true);
    std::vector<TopoDS_Shape> bodies;
    for (const auto &s : checked)
        for (TopExp_Explorer ex(s, TopAbs_SOLID); ex.More(); ex.Next())
            bodies.push_back(ex.Current());
    if (bodies.empty())
        throw std::runtime_error("No solids to fuse");
    if (bodies.size() == 1)
        return bodies.front();
    // Python fuse_overlapping/fuse_checked: a failed or material-losing union
    // keeps the already checked bodies side by side. A boolean is optional;
    // throwing here used to turn a successful reconstruction into a failure.
    auto apart = compound(checked);
    double total = 0, biggest = 0;
    for (const auto &body : bodies) {
        double v = volume(body);
        total += v;
        biggest = std::max(biggest, v);
    }
    double low = shapes.size() == 1 ? biggest : total;
    if (shapes.size() > 1 && overlapping_parts) {
        low = 0;
        for (const auto &s : checked)
            low = std::max(low, volume(s));
    }
    const double slack = .001 * std::abs(total) + .001;
    try {
        ShapeList args, tools;
        args.Append(bodies.front());
        for (std::size_t i = 1; i < bodies.size(); ++i)
            tools.Append(bodies[i]);
        BRepAlgoAPI_Fuse fuse;
        fuse.SetArguments(args);
        fuse.SetTools(tools);
        fuse.SetRunParallel(false);
        fuse.SetUseOBB(true);
        fuse.SetFuzzyValue(1e-5);
        fuse.Build();
        if (fuse.IsDone() && !fuse.Shape().IsNull()) {
            Diagnostics info;
            if (inspect(fuse.Shape(), info) && info.volume >= low - slack &&
                info.volume <= total + slack)
                return fuse.Shape();
        }
    } catch (const Standard_Failure &) {
        // Preserve the checked input compound, as upstream does.
    }
    ctx.diagnostics.warnings.push_back("Boolean fusion failed checks; retained separate bodies");
    return apart;
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
    trace(("validation: BRepCheck=" + std::to_string(checked) + " edge-incidence-errors=" +
           std::to_string(bad_edges) + " solids=" + std::to_string(info.solids) +
           " solid-faces=" + std::to_string(solid_faces.Extent()) +
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

// Derived from stl2curves build.py, Copyright (c) 2026 Duvdivan. MIT license.
#pragma once
// Included by build.cpp inside stl2curves::detail. Boundary helpers remain
// private to that translation unit, alongside the sewing/cleanup helpers.
// Python uses original mesh-point runs by default; exact intersections are an
// opt-in experiment. Preserve the runs and their shared corners on fitted faces.

namespace {
struct RunBounds {
    bool analytic = false, full_u = false, full_v = false;
    double u0 = 0, u1 = 0, v0 = 0, v1 = 0, pu = 0, pv = 0, eu = 0, ev = 0;
    std::map<unsigned, gp_Pnt2d> uv;
    static double wrap(double a, double p) {
        a = std::fmod(a, p);
        return a < 0 ? a + p : a;
    }
    static std::tuple<double, double, bool> arc(std::vector<double> a, double period) {
        for (auto &v : a)
            v = wrap(v, period);
        std::sort(a.begin(), a.end());
        std::vector<double> distinct;
        for (unsigned i = 0; i < a.size(); ++i)
            if (i == 0 || a[i] - a[i - 1] > 1e-6)
                distinct.push_back(a[i]);
        a = std::move(distinct);
        if (a.size() < 2)
            return {a.front(), a.front(), false};
        double gap = -1, cut = 0;
        for (unsigned i = 0; i < a.size(); ++i) {
            double d = a[(i + 1) % a.size()] + (i + 1 == a.size() ? period : 0) - a[i];
            if (d > gap) {
                gap = d;
                cut = a[(i + 1) % a.size()];
            }
        }
        if (a.size() >= 6 && gap <= period / 8)
            return {0, period, true};
        return {cut, cut + period - gap, false};
    }
    RunBounds(const Mesh &m, const Patch &p, double tol) {
        analytic = p.kind == "cylinder" || p.kind == "cone" || p.kind == "torus";
        if (!analytic)
            return;
        PointProjector project(p.surface);
        std::vector<double> us, vs;
        for (auto i : vertices(m, p.triangles)) {
            auto &q = project.perform(m.points[i]);
            if (!q.NbPoints()) {
                analytic = false;
                return;
            }
            double u, v;
            q.LowerDistanceParameters(u, v);
            uv.emplace(i, gp_Pnt2d(u, v));
            us.push_back(u);
            vs.push_back(v);
        }
        pu = p.surface->IsUPeriodic() ? p.surface->UPeriod() : 0;
        pv = p.surface->IsVPeriodic() ? p.surface->VPeriod() : 0;
        if (pu)
            std::tie(u0, u1, full_u) = arc(us, pu);
        else {
            u0 = *std::min_element(us.begin(), us.end());
            u1 = *std::max_element(us.begin(), us.end());
        }
        if (pv)
            std::tie(v0, v1, full_v) = arc(vs, pv);
        else {
            v0 = *std::min_element(vs.begin(), vs.end());
            v1 = *std::max_element(vs.begin(), vs.end());
        }
        for (auto &entry : uv) {
            auto &q = entry.second;
            if (pu) {
                double u = wrap(q.X(), pu);
                if (u < u0)
                    u += pu;
                q.SetX(u);
            }
            if (pv) {
                double v = wrap(q.Y(), pv);
                if (v < v0)
                    v += pv;
                q.SetY(v);
            }
        }
        gp_Pnt q;
        gp_Vec du, dv;
        p.surface->D1((u0 + u1) * .5, (v0 + v1) * .5, q, du, dv);
        eu = tol / std::max(du.Magnitude(), 1e-12);
        ev = tol / std::max(dv.Magnitude(), 1e-12);
    }
    unsigned labels(unsigned i) const {
        if (!analytic)
            return 0;
        auto q = uv.at(i);
        unsigned mask = 0;
        auto gap = [](double a, double b, double period) {
            double d = std::abs(a - b);
            return period ? std::abs(std::remainder(a - b, period)) : d;
        };
        if (!full_v) {
            if (gap(q.Y(), v0, pv) <= ev)
                mask |= 1;
            if (gap(q.Y(), v1, pv) <= ev)
                mask |= 2;
        }
        if (!full_u) {
            if (gap(q.X(), u0, pu) <= eu)
                mask |= 4;
            if (gap(q.X(), u1, pu) <= eu)
                mask |= 8;
        }
        return mask;
    }
};
struct RunTag {
    unsigned patch = ~0u, neighbor = ~0u;
    int label = -1;
    bool operator==(const RunTag &o) const {
        return patch == o.patch && neighbor == o.neighbor && label == o.label;
    }
};
struct BoundaryRun {
    RunTag tag;
    std::vector<unsigned> ids;
};
struct RunOutline {
    std::vector<std::vector<unsigned>> loops;
    std::vector<TopoDS_Wire> wires;
    ShapeList tools;
};
class BoundaryRuns {
    const Mesh &m;
    const std::vector<Patch> &patches;
    const Topology &top;
    std::vector<unsigned> owner;
    std::vector<RunBounds> bounds;
    std::vector<bool> crowded;
    std::map<unsigned, TopoDS_Vertex> corners;
    std::map<Edge, TopoDS_Edge> lines;
    std::map<std::vector<unsigned>, TopoDS_Edge> cache;
    TopoDS_Vertex corner(unsigned i) {
        if (crowded[i])
            return BRepBuilderAPI_MakeVertex(point(m.points[i])).Vertex();
        auto it = corners.emplace(i, TopoDS_Vertex{});
        if (it.second)
            it.first->second = BRepBuilderAPI_MakeVertex(point(m.points[i])).Vertex();
        return it.first->second;
    }
    TopoDS_Edge line(unsigned a, unsigned b) {
        if (crowded[a] || crowded[b])
            return BRepBuilderAPI_MakeEdge(point(m.points[a]), point(m.points[b])).Edge();
        Edge key(a, b);
        auto it = lines.find(key);
        if (it == lines.end()) {
            BRepBuilderAPI_MakeEdge make(corner(key.a), corner(key.b));
            auto edge = make.IsDone() ? make.Edge()
                                      : BRepBuilderAPI_MakeEdge(point(m.points[key.a]),
                                                                point(m.points[key.b]))
                                            .Edge();
            it = lines.emplace(key, edge).first;
        }
        auto edge = it->second;
        if (a > b)
            edge.Reverse();
        return edge;
    }
    TopoDS_Edge arc(const Vec &a, const Vec &mid, const Vec &b) {
        GC_MakeArcOfCircle make(point(a), point(mid), point(b));
        if (!make.IsDone())
            return BRepBuilderAPI_MakeEdge(point(a), point(b)).Edge();
        return BRepBuilderAPI_MakeEdge(make.Value()).Edge();
    }
    TopoDS_Edge spline(const std::vector<unsigned> &ids) {
        bool closed = ids.size() > 3 && ids.front() == ids.back();
        std::vector<Vec> pts;
        for (unsigned i = 0; i < ids.size() - (closed ? 1 : 0); ++i)
            if (i == 0 || (m.points[ids[i - 1]] - m.points[ids[i]]).norm() > 1e-6)
                pts.push_back(m.points[ids[i]]);
        if (pts.size() == 2 && !closed)
            return line(ids.front(), ids.back());
        if (pts.size() < 3)
            return {};
        auto ring = pts;
        if (closed)
            ring.push_back(pts.front());
        std::vector<double> lengths;
        for (unsigned i = 1; i < ring.size(); ++i)
            lengths.push_back((ring[i] - ring[i - 1]).norm());
        auto sorted = lengths;
        std::sort(sorted.begin(), sorted.end());
        double typical = sorted.size() > 3
                             ? (sorted[(sorted.size() - 1) / 2] + sorted[sorted.size() / 2]) * .5
                             : sorted.front();
        if (typical <= 1e-12)
            return {};
        std::vector<Vec> samples;
        for (unsigned i = 1; i < ring.size(); ++i) {
            unsigned n = lengths[i - 1] > 3 * typical
                             ? std::min(16u, unsigned(std::ceil(lengths[i - 1] / typical)))
                             : 1;
            for (unsigned j = 0; j < n; ++j)
                samples.push_back(ring[i - 1] + (ring[i] - ring[i - 1]) * (double(j) / n));
        }
        if (!closed)
            samples.push_back(ring.back());
        Handle(PointHandleArray) array = new PointHandleArray(1, samples.size());
        for (unsigned i = 0; i < samples.size(); ++i)
            array->SetValue(i + 1, point(samples[i]));
        GeomAPI_Interpolate fit(array, closed, 1e-7);
        fit.Perform();
        if (!fit.IsDone())
            return {};
        return BRepBuilderAPI_MakeEdge(fit.Curve()).Edge();
    }
    TopoDS_Edge circle(const RunTag &tag, const std::vector<unsigned> &ids) {
        const auto &p = patches[tag.patch];
        const auto &b = bounds[tag.patch];
        Handle(Geom_Curve) c;
        if (tag.label < 2)
            c = p.surface->VIso(tag.label == 0 ? b.v0 : b.v1);
        else
            c = p.surface->UIso(tag.label == 2 ? b.u0 : b.u1);
        BRepBuilderAPI_MakeEdge make(c);
        if (!make.IsDone())
            return {};
        auto edge = make.Edge();
        GeomAdaptor_Curve adaptor(c);
        if (adaptor.GetType() == GeomAbs_Circle) {
            auto cir = adaptor.Circle();
            auto center = cir.Location();
            Vec o(center.X(), center.Y(), center.Z());
            Vec n(cir.Axis().Direction().X(), cir.Axis().Direction().Y(),
                  cir.Axis().Direction().Z());
            double turn = 0;
            for (unsigned i = 1; i < ids.size(); ++i)
                turn += (m.points[ids[i - 1]] - o).cross(m.points[ids[i]] - o).dot(n);
            if (turn < 0)
                edge.Reverse();
        }
        return edge;
    }
    TopoDS_Edge natural(const RunTag &tag, const std::vector<unsigned> &ids) {
        const auto &p = patches[tag.patch];
        Vec a = m.points[ids.front()], b = m.points[ids.back()];
        if (tag.label >= 2 && (p.kind == "cylinder" || p.kind == "cone"))
            return line(ids.front(), ids.back());
        if (ids.size() > 2)
            return arc(a, m.points[ids[ids.size() / 2]], b);
        Vec qa = a - p.origin, qb = b - p.origin;
        double za = qa.dot(p.axis), zb = qb.dot(p.axis);
        Vec ra = qa - za * p.axis, rb = qb - zb * p.axis,
            radial = ra.normalized() + rb.normalized();
        if (radial.norm() < 1e-10)
            return {};
        radial.normalize();
        Vec mid;
        if (tag.label < 2)
            mid = p.origin + (za + zb) * .5 * p.axis + (ra.norm() + rb.norm()) * .5 * radial;
        else {
            double va = std::atan2(za, ra.norm() - p.radius),
                   vb = std::atan2(zb, rb.norm() - p.radius),
                   v = va + std::remainder(vb - va, 2 * 3.14159265358979323846) * .5;
            mid = p.origin + (p.radius + p.secondary * std::cos(v)) * radial +
                  p.secondary * std::sin(v) * p.axis;
        }
        return arc(a, mid, b);
    }
    TopoDS_Edge edge(const RunTag &tag, std::vector<unsigned> ids, bool patch_side) {
        if (tag.patch == ~0u)
            return line(ids.front(), ids.back());
        const auto &p = patches[tag.patch];
        if (tag.label < 0 && ids.size() == 2) {
            if (!patch_side || p.kind == "freeform" || p.kind == "blend" || p.kind == "pipe" ||
                p.kind == "extrusion")
                return line(ids[0], ids[1]);
            Vec a = m.points[ids[0]], b = m.points[ids[1]], mid = (a + b) * .5;
            PointProjector proj(p.surface);
            auto &q = proj.perform(mid);
            if (!q.NbPoints())
                return line(ids[0], ids[1]);
            double radius = p.kind == "sphere" || p.kind == "cylinder" ? p.radius
                            : p.kind == "torus"                        ? p.secondary
                                                                       : 1e-3;
            if (q.LowerDistance() > (a - b).squaredNorm() / (4 * std::max(1e-3, radius)) + 1e-4)
                return line(ids[0], ids[1]);
            auto v = q.NearestPoint();
            return arc(a, Vec(v.X(), v.Y(), v.Z()), b);
        }
        bool closed = ids.size() > 3 && ids.front() == ids.back(), flip;
        if (closed && tag.label >= 0 &&
            ((tag.label < 2 && bounds[tag.patch].full_u) ||
             (tag.label >= 2 && bounds[tag.patch].full_v)))
            return circle(tag, ids);
        if (closed) {
            ids.pop_back();
            auto at = std::min_element(ids.begin(), ids.end());
            std::rotate(ids.begin(), at, ids.end());
            flip = ids[1] > ids.back();
            if (flip)
                std::reverse(ids.begin() + 1, ids.end());
            ids.push_back(ids.front());
        } else {
            flip = ids.front() > ids.back();
            if (flip)
                std::reverse(ids.begin(), ids.end());
        }
        auto found = cache.find(ids);
        if (found == cache.end()) {
            TopoDS_Edge value;
            try {
                value = tag.label < 0 ? spline(ids) : natural(tag, ids);
                if (value.IsNull())
                    value = spline(ids);
            } catch (const Standard_Failure &) {
                value = spline(ids);
            }
            if (!value.IsNull() && !crowded[ids.front()] && !crowded[ids.back()])
                try {
                    double a, b;
                    auto curve = BRep_Tool::Curve(value, a, b);
                    BRepBuilderAPI_MakeEdge make(curve, corner(ids.front()), corner(ids.back()), a,
                                                 b);
                    if (make.IsDone())
                        value = make.Edge();
                } catch (const Standard_Failure &) {
                }
            found = cache.emplace(ids, value).first;
        }
        auto value = found->second;
        if (!value.IsNull() && flip)
            value.Reverse();
        return value;
    }

  public:
    BoundaryRuns(const Mesh &mesh, const std::vector<Patch> &ps, const Topology &t, double mesh_tol)
        : m(mesh), patches(ps), top(t), owner(mesh.triangles.size(), ~0u),
          crowded(mesh.points.size(), false) {
        for (unsigned i = 0; i < ps.size(); ++i) {
            for (auto t : ps[i].triangles)
                owner[t] = i;
            bounds.emplace_back(m, ps[i], 10 * mesh_tol);
        }
        BoundaryCloud cloud{m.points};
        using Tree =
            nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, BoundaryCloud>,
                                                BoundaryCloud, 3, std::size_t>;
        Tree tree(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(16));
        tree.buildIndex();
        for (unsigned i = 0; i < m.points.size(); ++i) {
            std::size_t at[2];
            double ds[2];
            crowded[i] = tree.knnSearch(m.points[i].data(), 2, at, ds) == 2 && ds[1] <= .02 * .02;
        }
    }
    RunOutline outline(unsigned own) {
        const auto &p = patches[own];
        auto rim = patch_rim(m, p.triangles);
        std::map<unsigned, unsigned> next;
        for (auto &entry : rim)
            if (entry.second.count == 1)
                if (!next.emplace(entry.second.from, entry.second.to).second)
                    throw std::runtime_error("mesh boundary branching rim");
        RunOutline result;
        while (!next.empty()) {
            unsigned start = next.begin()->first, at = start;
            std::vector<unsigned> loop;
            do {
                loop.push_back(at);
                auto found = next.find(at);
                if (found == next.end())
                    throw std::runtime_error("mesh boundary open rim");
                at = found->second;
                next.erase(found);
            } while (at != start);
            std::vector<RunTag> tags;
            for (unsigned i = 0; i < loop.size(); ++i) {
                Edge e(loop[i], loop[(i + 1) % loop.size()]);
                unsigned other = ~0u;
                for (auto t : top.edges.at(e))
                    if (owner[t] != own)
                        other = owner[t];
                unsigned ref = p.kind == "plane" ? other : own;
                RunTag tag;
                if (ref != ~0u && patches[ref].kind != "plane") {
                    tag.patch = ref;
                    tag.neighbor = other;
                    unsigned labels = bounds[ref].labels(e.a) & bounds[ref].labels(e.b);
                    for (int k : {1, 0, 2, 3})
                        if (labels & (1u << k)) {
                            tag.label = k;
                            break;
                        }
                }
                tags.push_back(tag);
            }
            std::vector<unsigned> cuts;
            for (unsigned i = 0; i < loop.size(); ++i) {
                unsigned prev = (i + loop.size() - 1) % loop.size();
                bool cut = !(tags[i] == tags[prev]) || tags[i].patch == ~0u;
                Vec a = m.points[loop[i]] - m.points[loop[prev]],
                    b = m.points[loop[(i + 1) % loop.size()]] - m.points[loop[i]];
                if (tags[i].label < 0 && a.dot(b) < .8660254037844386 * a.norm() * b.norm())
                    cut = true;
                if (cut)
                    cuts.push_back(i);
            }
            std::vector<BoundaryRun> runs;
            if (cuts.empty()) {
                runs.push_back({tags[0], loop});
                runs.back().ids.push_back(loop.front());
            } else
                for (unsigned i = 0; i < cuts.size(); ++i) {
                    unsigned begin = cuts[i], end = cuts[(i + 1) % cuts.size()];
                    BoundaryRun run;
                    run.tag = tags[begin];
                    run.ids.push_back(loop[begin]);
                    unsigned k = (begin + 1) % loop.size();
                    while (k != end) {
                        run.ids.push_back(loop[k]);
                        k = (k + 1) % loop.size();
                    }
                    run.ids.push_back(loop[end]);
                    runs.push_back(std::move(run));
                }
            BRepBuilderAPI_MakeWire wire;
            for (auto &run : runs) {
                auto e = edge(run.tag, run.ids, p.kind != "plane");
                if (e.IsNull())
                    throw std::runtime_error("mesh boundary failed run");
                wire.Add(e);
                result.tools.Append(e);
            }
            if (!wire.IsDone())
                throw std::runtime_error("mesh boundary wire disconnected");
            result.loops.push_back(loop);
            result.wires.push_back(wire.Wire());
        }
        return result;
    }
};
void orient_run_face(TopoDS_Face &face, const Patch &p, const Topology &top) {
    PointProjector proj(p.surface);
    auto &q = proj.perform(top.centers[p.triangles.front()]);
    if (!q.NbPoints())
        return;
    double u, v;
    q.LowerDistanceParameters(u, v);
    gp_Pnt at;
    gp_Vec du, dv;
    p.surface->D1(u, v, at, du, dv);
    auto normal = du.Crossed(dv);
    double dot = top.normals[p.triangles.front()].dot(Vec(normal.X(), normal.Y(), normal.Z()));
    if (face.Orientation() == TopAbs_REVERSED)
        dot = -dot;
    if (dot < 0)
        face.Reverse();
}
std::vector<TopoDS_Face> run_outline_face(const Mesh &m, const Patch &p, const Topology &top,
                                          const RunOutline &o) {
    if (o.wires.empty())
        return {};
    Vec normal = Vec::Zero();
    double target = 0;
    for (auto t : p.triangles) {
        normal += top.areas[t] * top.normals[t];
        target += top.areas[t];
    }
    std::vector<std::pair<double, unsigned>> order;
    for (unsigned i = 0; i < o.loops.size(); ++i) {
        double area = 0;
        auto &loop = o.loops[i];
        for (unsigned j = 0; j < loop.size(); ++j)
            area += m.points[loop[j]].cross(m.points[loop[(j + 1) % loop.size()]]).dot(normal);
        order.push_back({std::abs(area), i});
    }
    std::sort(order.rbegin(), order.rend());
    for (bool flip : {false, true}) {
        std::vector<TopoDS_Wire> wires;
        for (auto entry : order) {
            auto w = TopoDS::Wire(BRepBuilderAPI_Copy(o.wires[entry.second], true, false).Shape());
            if (flip)
                w.Reverse();
            wires.push_back(w);
        }
        BRepBuilderAPI_MakeFace make(p.surface, wires[0], p.kind == "plane");
        for (unsigned i = 1; i < wires.size(); ++i)
            make.Add(wires[i]);
        if (!make.IsDone())
            continue;
        ShapeFix_Face fix(make.Face());
        fix.Perform();
        auto face = fix.Face();
        ShapeAnalysis_ShapeTolerance analysis;
        GProp_GProps props;
        BRepGProp::SurfaceProperties(face, props);
        if (BRepCheck_Analyzer(face).IsValid() && analysis.Tolerance(face, 1) <= 1. &&
            std::abs(props.Mass() - target) <= .05 * target + 1e-3) {
            orient_run_face(face, p, top);
            return {face};
        }
    }
    return {};
}
std::vector<TopoDS_Face> split_run_face(const Mesh &m, const Patch &p, const Topology &top,
                                        const RunOutline &o, double mesh_tol, double fuzz = 0) {
    ShapeList args;
    auto base = reference_surface_face(p, true);
    if (base.IsNull())
        base = BRepBuilderAPI_MakeFace(p.surface, 1e-7).Face();
    args.Append(base);
    BRepAlgoAPI_Splitter splitter;
    splitter.SetArguments(args);
    splitter.SetTools(o.tools);
    splitter.SetNonDestructive(true);
    // build._split_face uses the height field's rim fit error, independently
    // of the analytic feature's midpoint sag used by sewing.
    splitter.SetFuzzyValue(std::max({10 * mesh_tol, 1e-4, 2 * p.rim_deviation, fuzz}));
    splitter.Build();
    if (!splitter.IsDone())
        return {};
    std::set<unsigned> rim;
    for (auto &loop : o.loops)
        rim.insert(loop.begin(), loop.end());
    std::vector<unsigned> inner;
    double target = 0;
    for (auto t : p.triangles) {
        target += top.areas[t];
        auto tri = m.triangles[t];
        if (!rim.count(tri[0]) && !rim.count(tri[1]) && !rim.count(tri[2]))
            inner.push_back(t);
    }
    PointProjector proj(p.surface);
    struct Piece {
        TopoDS_Face face;
        double area;
    };
    auto pick = [&](const std::vector<unsigned> &samples) {
        std::vector<TopoDS_Vertex> probes;
        unsigned n = std::min<std::size_t>(60, samples.size());
        for (unsigned i = 0; i < n; ++i) {
            auto t = samples[n == 1 ? 0 : std::size_t(i) * (samples.size() - 1) / (n - 1)];
            auto &q = proj.perform(top.centers[t]);
            if (q.NbPoints())
                probes.push_back(BRepBuilderAPI_MakeVertex(q.NearestPoint()).Vertex());
        }
        std::vector<Piece> found;
        for (TopExp_Explorer ex(splitter.Shape(), TopAbs_FACE); ex.More(); ex.Next()) {
            auto face = TopoDS::Face(ex.Current());
            bool hit = false;
            for (auto &probe : probes) {
                BRepExtrema_DistShapeShape dist(probe, face);
                if (dist.IsDone() && dist.Value() < 1e-3) {
                    hit = true;
                    break;
                }
            }
            if (hit) {
                GProp_GProps props;
                BRepGProp::SurfaceProperties(face, props);
                found.push_back({face, props.Mass()});
            }
        }
        return found;
    };
    auto total = [](const std::vector<Piece> &items) {
        double sum = 0;
        for (auto &i : items)
            sum += i.area;
        return sum;
    };
    auto found = pick(inner.empty() ? p.triangles : inner);
    if (!inner.empty() && std::abs(total(found) - target) > .1 * target)
        found = pick(p.triangles);
    if (found.empty() || found.size() > 10)
        return {};
    unsigned mask = (1u << found.size()) - 1;
    if (std::abs(total(found) - target) > .1 * target) {
        double best = .1 * target;
        mask = 0;
        for (unsigned bits = 1; bits < (1u << found.size()); ++bits) {
            double area = 0;
            for (unsigned i = 0; i < found.size(); ++i)
                if (bits & (1u << i))
                    area += found[i].area;
            if (std::abs(area - target) <= best) {
                best = std::abs(area - target);
                mask = bits;
            }
        }
    }
    if (!mask)
        return {};
    std::vector<TopoDS_Face> result;
    for (unsigned i = 0; i < found.size(); ++i)
        if (mask & (1u << i)) {
            auto face = found[i].face;
            orient_run_face(face, p, top);
            result.push_back(face);
        }
    return result;
}
} // namespace
TopoDS_Shape mesh_boundary_shape(const Mesh &m, const std::vector<Patch> &input, Context &ctx) {
    std::vector<Patch> patches;
    for (auto &p : input) {
        std::vector<Tri> tris;
        for (auto t : p.triangles)
            tris.push_back(m.triangles[t]);
        for (auto &group : unpinched_groups(tris)) {
            auto piece = p;
            piece.triangles.clear();
            for (auto i : group)
                piece.triangles.push_back(p.triangles[i]);
            patches.push_back(piece);
        }
    }
    std::stable_partition(patches.begin(), patches.end(),
                          [](const Patch &p) { return p.kind != "plane"; });
    auto top = topology(m);
    double magnitude = 0;
    for (auto &q : m.points)
        magnitude = std::max(magnitude, q.cwiseAbs().maxCoeff());
    double mesh_tol = .001 + 2e-7 * magnitude;
    BoundaryRuns pool(m, patches, top, mesh_tol);
    std::vector<TopoDS_Face> natural(patches.size());
    std::vector<RunOutline> outlines(patches.size());
    for (unsigned i = 0; i < patches.size(); ++i)
        if (patches[i].kind != "plane")
            natural[i] = natural_patch(m, patches[i], ctx.options.tolerance);
    for (unsigned i = 0; i < patches.size(); ++i)
        if (patches[i].kind != "plane" && natural[i].IsNull()) {
            trace(("mesh boundary outline=" + std::to_string(i)).c_str());
            outlines[i] = pool.outline(i);
        }
    std::vector<TopoDS_Shape> faces;
    for (unsigned i = 0; i < patches.size(); ++i) {
        auto &p = patches[i];
        trace(("mesh boundary patch=" + std::to_string(i) + " " + p.kind).c_str());
        if (!natural[i].IsNull()) {
            faces.push_back(natural[i]);
            continue;
        }
        if (p.kind == "plane")
            outlines[i] = pool.outline(i);
        std::vector<TopoDS_Face> made;
        bool outline_first =
            p.kind == "plane" || p.kind == "blend" || p.kind == "pipe" || p.triangles.size() >= 600;
        if (outline_first)
            made = run_outline_face(m, p, top, outlines[i]);
        if (made.empty()) {
            auto spline = Handle(Geom_BSplineSurface)::DownCast(p.surface);
            if (spline.IsNull() || spline->NbUPoles() * spline->NbVPoles() <= 256)
                made = split_run_face(m, p, top, outlines[i], mesh_tol);
        }
        if (made.empty() && !outline_first)
            made = run_outline_face(m, p, top, outlines[i]);
        if (made.empty()) {
            RunOutline polygon;
            polygon.loops = outlines[i].loops;
            for (auto &loop : polygon.loops) {
                BRepBuilderAPI_MakeWire wire;
                for (unsigned j = 0; j < loop.size(); ++j)
                    wire.Add(BRepBuilderAPI_MakeEdge(point(m.points[loop[j]]),
                                                     point(m.points[loop[(j + 1) % loop.size()]]))
                                 .Edge());
                if (wire.IsDone())
                    polygon.wires.push_back(wire.Wire());
            }
            if (polygon.wires.size() == polygon.loops.size())
                made = run_outline_face(m, p, top, polygon);
        }
        if (made.empty()) {
            PointProjector project(p.surface);
            double gap = 0;
            for (const auto &tool : outlines[i].tools) {
                BRepAdaptor_Curve c(TopoDS::Edge(tool));
                for (unsigned j = 0; j < 9; ++j) {
                    auto q = c.Value(c.FirstParameter() +
                                     (c.LastParameter() - c.FirstParameter()) * j / 8.);
                    auto &at = project.perform(Vec(q.X(), q.Y(), q.Z()));
                    if (at.NbPoints())
                        gap = std::max(gap, at.LowerDistance());
                }
            }
            if (gap <= .1 && 1.25 * gap > 10 * mesh_tol)
                made = split_run_face(m, p, top, outlines[i], mesh_tol, 1.25 * gap);
        }
        if (made.empty())
            throw std::runtime_error("mesh boundary face failed " + std::to_string(i) + " " +
                                     p.kind);
        for (auto &f : made)
            faces.push_back(f);
    }
    double tol = patch_sewing_tolerance(m, patches, ctx.options.tolerance);
    auto comp = rejoin_exact(compound(faces), m, tol);
    BRepBuilderAPI_Sewing sewing(tol);
    sewing.Add(comp);
    sewing.Perform();
    trace(("mesh boundary sewing free=" + std::to_string(sewing.NbFreeEdges())).c_str());
    auto shape = close_shells(without_collapsed_edges(sewing.SewedShape(), tol), tol);
    return tidy_shape(shape, m, patches, tol);
}

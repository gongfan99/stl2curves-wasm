// Integration of the facet-based fillet/transition path from the reference.
namespace {
namespace RF = reference_features;
Patch reference_patch(const RF::Pass &pass, const RF::Found &feature) {
    Patch patch;
    patch.kind = feature.model.kind;
    patch.origin = feature.model.center();
    patch.axis = feature.model.d;
    patch.radius = feature.model.r;
    patch.secondary = patch.kind == "torus" ? feature.model.tube : feature.model.k;
    patch.deviation = feature.worst;
    patch.trimmed = feature.trimmed;
    for (auto facet : feature.facets)
        patch.triangles.insert(patch.triangles.end(), pass.fs.triangles[facet].begin(),
                               pass.fs.triangles[facet].end());
    patch.reference_feature = std::make_shared<RF::Found>(feature);
    return patch;
}
std::vector<RF::Found> reference_features_of(const std::vector<Patch> &patches) {
    std::vector<RF::Found> features;
    for (const auto &patch : patches)
        features.push_back(*patch.reference_feature);
    return features;
}
bool reference_analysis(const std::vector<Patch> &patches) {
    return !patches.empty() && std::all_of(patches.begin(), patches.end(), [](const Patch &p) {
        return bool(p.reference_feature);
    });
}
std::optional<std::vector<Patch>> reference_analysis(const Mesh &mesh, const Topology &top,
                                                     Context &context) {
    // This port handles fillet transitions with planar height-field bases.
    // Other feature families continue through the established recognizer.
    if (!context.options.curves || !context.options.blends || context.options.true_size ||
        context.expired())
        return {};
    RF::Pass pass(mesh, top, {});
    RF::FilletPort fillets(pass);
    fillets.fillet_tol = .02 * context.options.accuracy / .02;
    fillets.grow_tol = .05 * context.options.accuracy / .02;
    auto rolled = fillets.find();
    if (rolled.empty())
        return {};
    std::vector<Patch> known;
    for (const auto &feature : rolled)
        known.push_back(reference_patch(pass, feature));
    RF::Pass loose(mesh, top, known);
    auto found = loose.run();
    rolled.insert(rolled.end(), found.begin(), found.end());
    auto features = RF::fill_pinholes(loose, RF::merge_same_surface(loose, std::move(rolled)));

    // Only select this currently supported path when every loose curved facet
    // is explained by these features or by a whole transition. This is a
    // capability check; it does not change the reference's ownership decisions.
    auto snapped = RF::snap(loose, features);
    auto areas = RF::transition_areas(loose, snapped, std::min(.01, .1 * context.options.accuracy));
    if (areas.empty())
        return {};
    auto owner = RF::owners(loose, snapped);
    for (const auto &area : areas)
        for (auto facet : area.facets)
            owner[facet] = 0;
    double small = .02 * *std::max_element(loose.fs.areas.begin(), loose.fs.areas.end());
    for (unsigned facet = 0; facet < owner.size(); ++facet)
        if (owner[facet] < 0 && loose.fs.areas[facet] <= small &&
            !loose.fs.neighbors[facet].empty())
            return {};
    std::vector<Patch> patches;
    for (const auto &feature : features)
        patches.push_back(reference_patch(loose, feature));
    trace("reference facet fillet analysis");
    return patches;
}
void reference_snap(const Mesh &mesh, const Topology &top, std::vector<Patch> &patches) {
    RF::Pass pass(mesh, top, {});
    auto features = RF::snap(pass, reference_features_of(patches));
    patches.clear();
    for (const auto &feature : features)
        patches.push_back(reference_patch(pass, feature));
}
std::optional<std::vector<Patch>> reference_finish(const Mesh &mesh, const Topology &top,
                                                   const std::vector<Patch> &patches,
                                                   Context &context) {
    RF::Pass pass(mesh, top, {});
    auto features = reference_features_of(patches);
    std::set<unsigned> used;
    std::vector<Patch> transitions;
    if (context.options.blends)
        for (const auto &area :
             RF::transition_areas(pass, features, std::min(.01, .1 * context.options.accuracy))) {
            Patch patch;
            auto exact = RF::corner_surface(pass, area.facets, context.options.accuracy);
            if (exact) {
                patch = reference_patch(pass, exact);
                surface(patch);
            } else
                patch = reference_freeform_patch(mesh, top, area.facets, context);
            if (patch.surface.IsNull())
                return {};
            transitions.push_back(std::move(patch));
            used.insert(area.parts.begin(), area.parts.end());
        }
    std::vector<Patch> out;
    for (unsigned k = 0; k < patches.size(); ++k)
        if (!used.count(k)) {
            auto patch = patches[k];
            surface(patch);
            out.push_back(std::move(patch));
        }
    out.insert(out.end(), transitions.begin(), transitions.end());
    if (context.options.blends)
        absorb_feature_islands(mesh, top, out, context.options.accuracy);
    std::vector<bool> claimed(mesh.triangles.size(), false);
    for (const auto &patch : out)
        for (auto triangle : patch.triangles)
            claimed[triangle] = true;
    for (const auto &facet : pass.fs.triangles) {
        std::vector<unsigned> left;
        for (auto triangle : facet)
            if (!claimed[triangle])
                left.push_back(triangle);
        if (!left.empty()) {
            auto plane = fit_plane(mesh, left);
            surface(plane);
            out.push_back(std::move(plane));
        }
    }
    std::vector<unsigned> bodies(mesh.triangles.size());
    unsigned body = 0;
    for (const auto &group : components(top.neighbors)) {
        for (auto triangle : group)
            bodies[triangle] = body;
        ++body;
    }
    for (auto &patch : out)
        patch.body = bodies[patch.triangles.front()];
    trace("reference whole-unit transitions");
    return out;
}
} // namespace
TopoDS_Shape reference_surface_face(const Patch &patch, bool generous) {
    if (!patch.reference_feature || (!generous && patch.trimmed))
        return {};
    const auto &feature = *patch.reference_feature;
    const auto &model = feature.model;
    double lo = feature.lo, hi = feature.hi, u0 = feature.u0, span = feature.span;
    if (generous) {
        if (model.kind == "cylinder" || model.kind == "cone") {
            double margin = .2 * (hi - lo) + .5, new_lo = lo - margin, new_hi = hi + margin;
            if (model.k != 0) {
                double apex = -model.r / model.k;
                if (new_lo < apex && apex <= lo)
                    new_lo = apex;
                if (hi <= apex && apex < new_hi)
                    new_hi = apex;
            }
            lo = new_lo;
            hi = new_hi;
        } else if (model.kind == "torus" || model.sphere_cap) {
            double radius = model.sphere_cap ? model.r : model.tube;
            double center = model.sphere_cap ? 0 : model.r;
            double margin = hi - lo >= 2 * RF::pi - 1e-9 ? 0 : 10 * RF::pi / 180;
            auto crosses = [&](double a, double b) {
                for (unsigned i = 0; i < 60; ++i)
                    if (center + radius * std::cos(a + (b - a) * i / 59) < 0)
                        return true;
                return false;
            };
            double grow_lo = 0, grow_hi = 0;
            for (double fraction : {1., .5, .25, .1})
                if (!crosses(lo - margin * fraction, lo)) {
                    grow_lo = margin * fraction;
                    break;
                }
            for (double fraction : {1., .5, .25, .1})
                if (!crosses(hi, hi + margin * fraction)) {
                    grow_hi = margin * fraction;
                    break;
                }
            lo -= grow_lo;
            hi += grow_hi;
            if ((margin && hi - lo >= 2 * RF::pi - 1e-6) || crosses(lo, hi))
                return {};
        }
        if (span < 2 * RF::pi) {
            double grow = 10 * RF::pi / 180;
            if (span + 2 * grow >= 2 * RF::pi) {
                u0 = 0;
                span = 2 * RF::pi;
            } else {
                u0 -= grow;
                span += 2 * grow;
            }
        }
    }
    TopoDS_Shape shape;
    if (model.kind == "torus" && model.r > 0 &&
        (hi - lo >= 2 * RF::pi - 1e-9 || model.r < model.tube)) {
        auto torus = new Geom_ToroidalSurface(
            gp_Ax3(point(model.center()), direction(model.d), direction(model.e1)), model.r,
            model.tube);
        BRepBuilderAPI_MakeFace make(torus, u0, u0 + span, lo, hi, 1e-7);
        if (make.IsDone())
            shape = make.Face();
    } else if (model.kind == "cylinder" || model.kind == "cone" || model.kind == "torus" ||
               model.sphere_cap) {
        bool circular = model.kind == "torus" || model.sphere_cap;
        auto profile_point = [&](double height) -> gp_Pnt {
            double radius = model.r + model.k * height;
            double z = height;
            if (circular) {
                double r = model.sphere_cap ? model.r : model.tube;
                radius = (model.sphere_cap ? 0 : model.r) + r * std::cos(height);
                z = model.zc + r * std::sin(height);
            }
            if (std::abs(radius) < 1e-4)
                radius = 0;
            return point(model.a + z * model.d +
                         radius * (std::cos(u0) * model.e1 + std::sin(u0) * model.e2));
        };
        TopoDS_Edge edge;
        if (!circular)
            edge = BRepBuilderAPI_MakeEdge(profile_point(lo), profile_point(hi)).Edge();
        else if (hi - lo >= 2 * RF::pi - 1e-9) {
            double radius = model.sphere_cap ? model.r : model.tube;
            Vec out = std::cos(u0) * model.e1 + std::sin(u0) * model.e2;
            Vec center = model.a + model.zc * model.d + (model.sphere_cap ? 0 : model.r) * out;
            gp_Circ circle(gp_Ax2(point(center), direction(out.cross(model.d)), direction(out)),
                           radius);
            edge = BRepBuilderAPI_MakeEdge(circle).Edge();
        } else {
            GC_MakeArcOfCircle arc(profile_point(lo), profile_point((lo + hi) / 2),
                                   profile_point(hi));
            if (arc.IsDone())
                edge = BRepBuilderAPI_MakeEdge(arc.Value()).Edge();
        }
        if (!edge.IsNull()) {
            BRepPrimAPI_MakeRevol revol(edge, gp_Ax1(point(model.a), direction(model.d)), span);
            if (revol.IsDone())
                shape = revol.Shape();
        }
    }
    if (!generous && !shape.IsNull() && !feature.convex)
        shape.Reverse();
    return shape;
}

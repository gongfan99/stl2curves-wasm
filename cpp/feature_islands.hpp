// blends._islands and Freeform.signed from stl2curves (MIT).
#pragma once
#include "feature_facets.hpp"
#include "projection.hpp"
#include <gp_Vec.hxx>
#include <limits>
#include <memory>

namespace stl2curves::detail {
class FeatureSignedDistance {
    const Patch &patch;
    std::unique_ptr<PointProjector> projector;

  public:
    explicit FeatureSignedDistance(const Patch &p) : patch(p) {
        if (p.kind == "freeform" || p.kind == "blend" || p.kind == "pipe" || p.kind == "extrusion")
            projector = std::make_unique<PointProjector>(p.surface);
    }
    double operator()(const Vec &p) {
        const Vec q = p - patch.origin;
        if (patch.kind == "sphere")
            return std::abs(q.norm() - patch.radius);
        if (!projector) {
            double z = q.dot(patch.axis), rho = (q - z * patch.axis).norm();
            if (patch.kind == "cylinder")
                return std::abs(rho - patch.radius);
            if (patch.kind == "cone")
                return std::abs(rho - patch.radius - z * patch.secondary) /
                       std::sqrt(1 + patch.secondary * patch.secondary);
            if (patch.kind == "torus")
                return std::abs(std::hypot(rho - patch.radius, z) - patch.secondary);
            return std::abs(q.dot(patch.axis));
        }
        // A plane-based height field's uvh() keeps its base coordinates.
        // Its signed residual is normalized by the implicit gradient, rather
        // than being a nearest-point projection onto the fitted surface.
        double u, v;
        if (!projector->chart_parameters(p, u, v))
            return std::numeric_limits<double>::infinity();
        gp_Pnt at;
        gp_Vec du, dv;
        patch.surface->D1(u, v, at, du, dv);
        auto normal = du.Crossed(dv);
        if (normal.SquareMagnitude() <= 1e-24)
            return std::numeric_limits<double>::infinity();
        return std::abs(Vec(p.x() - at.X(), p.y() - at.Y(), p.z() - at.Z())
                            .dot(Vec(normal.X(), normal.Y(), normal.Z()).normalized()));
    }
};

inline unsigned absorb_feature_islands(const Mesh &mesh, const Topology &top,
                                       std::vector<Patch> &patches,
                                       double max_blend_deviation = .1) {
    FeatureFacets facets(mesh, top);
    std::vector<int> triangle_owner(mesh.triangles.size(), -1);
    for (unsigned k = 0; k < patches.size(); ++k)
        if (patches[k].kind != "plane")
            for (auto t : patches[k].triangles)
                triangle_owner[t] = int(k);
    std::vector<int> owner(facets.triangles.size(), -1);
    for (unsigned f = 0; f < owner.size(); ++f) {
        owner[f] = triangle_owner[facets.triangles[f].front()];
        for (auto t : facets.triangles[f])
            if (triangle_owner[t] != owner[f]) {
                owner[f] = -2;
                break;
            }
    }
    std::vector<bool> seen(owner.size(), false);
    std::map<unsigned, std::vector<unsigned>> grown;
    std::vector<std::unique_ptr<FeatureSignedDistance>> distance(patches.size());
    std::vector<double> near(patches.size(), -1);
    for (unsigned seed = 0; seed < owner.size(); ++seed) {
        if (owner[seed] != -1 || seen[seed])
            continue;
        std::vector<unsigned> group, stack{seed};
        std::set<int> around;
        seen[seed] = true;
        while (!stack.empty() && group.size() <= 4) {
            auto g = stack.back();
            stack.pop_back();
            group.push_back(g);
            for (auto h : facets.beside[g]) {
                if (owner[h] != -1)
                    around.insert(owner[h]);
                else if (!seen[h]) {
                    seen[h] = true;
                    stack.push_back(h);
                }
            }
        }
        if (!stack.empty() || group.size() > 4 || around.size() != 1 || *around.begin() < 0)
            continue;
        unsigned k = unsigned(*around.begin());
        const auto &patch = patches[k];
        if (patch.kind == "blend") {
            bool smooth = true, thin = true;
            for (auto f : group) {
                for (auto h : facets.beside[f])
                    if (owner[h] == int(k) && facets.normals[f].dot(facets.normals[h]) <=
                                                  std::cos(35 * 3.14159265358979323846 / 180))
                        smooth = false;
                for (auto t : facets.triangles[f]) {
                    auto tri = mesh.triangles[t];
                    double longest = 0;
                    for (unsigned j = 0; j < 3; ++j)
                        longest = std::max(
                            longest, (mesh.points[tri[j]] - mesh.points[tri[(j + 1) % 3]]).norm());
                    if (4 * top.areas[t] / std::max(longest, 1e-12) > max_blend_deviation)
                        thin = false;
                }
            }
            if (smooth || thin)
                grown[k].insert(grown[k].end(), group.begin(), group.end());
            continue;
        }
        if (!distance[k])
            distance[k] = std::make_unique<FeatureSignedDistance>(patch);
        if (near[k] < 0) {
            double worst = patch.kind == "freeform" ? patch.deviation : 0;
            for (auto t : patch.triangles) {
                auto tri = mesh.triangles[t];
                for (unsigned j = 0; j < 3; ++j)
                    worst = std::max(
                        worst,
                        (*distance[k])((mesh.points[tri[j]] + mesh.points[tri[(j + 1) % 3]]) * .5));
            }
            near[k] = std::max(worst, facets.tolerance);
        }
        bool on = true;
        for (auto f : group) {
            for (auto v : facets.verts[f])
                if ((*distance[k])(mesh.points[v]) > near[k])
                    on = false;
            if ((*distance[k])(facets.centers[f]) > near[k])
                on = false;
        }
        if (on)
            grown[k].insert(grown[k].end(), group.begin(), group.end());
    }
    unsigned absorbed = 0;
    for (const auto &entry : grown)
        for (auto f : entry.second)
            for (auto t : facets.triangles[f]) {
                patches[entry.first].triangles.push_back(t);
                triangle_owner[t] = int(entry.first);
                ++absorbed;
            }
    for (auto &p : patches)
        if (p.kind == "plane")
            p.triangles.erase(std::remove_if(p.triangles.begin(), p.triangles.end(),
                                             [&](unsigned t) { return triangle_owner[t] >= 0; }),
                              p.triangles.end());
    patches.erase(std::remove_if(patches.begin(), patches.end(),
                                 [](const Patch &p) { return p.triangles.empty(); }),
                  patches.end());
    return absorbed;
}
} // namespace stl2curves::detail

// Facet-space operations from stl2curves/features.py (MIT).
#pragma once
#include "feature_order.hpp"
#include "internal.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <map>
#include <numeric>
#include <set>

namespace stl2curves::detail {
struct FeatureFacets {
    std::vector<unsigned> facet_of;
    std::vector<std::vector<unsigned>> triangles, verts, neighbors, beside;
    std::vector<Vec> normals, centers;
    std::vector<double> areas;
    std::map<std::pair<unsigned, unsigned>, std::vector<Edge>> shared;
    double noise = 0, tolerance = .001;

    FeatureFacets(const Mesh &mesh, const Topology &connections) {
        constexpr double pi = 3.14159265358979323846;
        struct TriangleGeometry {
            const EdgeFaces &edges;
            std::vector<Vec> normals, centers;
            std::vector<double> areas;
        } top{connections.edges, {}, {}, {}};
        top.normals.resize(mesh.triangles.size());
        top.centers.resize(mesh.triangles.size());
        top.areas.resize(mesh.triangles.size());
        for (unsigned t = 0; t < mesh.triangles.size(); ++t) {
            auto tri = mesh.triangles[t];
            const auto &a = mesh.points[tri[0]], &b = mesh.points[tri[1]], &c = mesh.points[tri[2]];
            Vec cross = (b - a).cross(c - a);
            double length =
                std::sqrt(cross.x() * cross.x() + cross.y() * cross.y() + cross.z() * cross.z());
            top.areas[t] = length / 2;
            for (unsigned k = 0; k < 3; ++k) {
                top.normals[t][k] = length > 0 ? cross[k] / length : 0;
                top.centers[t][k] = (a[k] + b[k] + c[k]) / 3;
            }
        }
        double magnitude = 0;
        for (const auto &p : mesh.points)
            magnitude = std::max(magnitude, p.cwiseAbs().maxCoeff());
        tolerance += 2e-7 * magnitude;
        for (double step : {.01, .001, .0001}) {
            std::size_t near = 0;
            for (const auto &p : mesh.points)
                for (unsigned k = 0; k < 3; ++k)
                    near += std::abs(p[k] - std::nearbyint(p[k] / step) * step) <=
                            2.5e-7 * std::abs(p[k]) + 1e-7;
            if (near >= .75 * mesh.points.size() * 3) {
                noise = step * std::sqrt(3.) / 2;
                break;
            }
        }
        std::vector<unsigned> parent(mesh.triangles.size());
        std::iota(parent.begin(), parent.end(), 0);
        auto root = [&](unsigned i) {
            while (parent[i] != i) {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            return i;
        };
        struct Pair {
            unsigned a, b;
            Edge edge;
            double dot;
        };
        std::vector<Pair> pairs;
        std::map<Edge, std::vector<unsigned>> edges(top.edges.begin(), top.edges.end());
        for (const auto &entry : edges)
            for (unsigned i = 1; i < entry.second.size(); ++i) {
                unsigned a = entry.second[i - 1], b = entry.second[i];
                double dot = top.normals[a].dot(top.normals[b]);
                if (dot > std::cos(.05 * pi / 180))
                    parent[root(a)] = root(b);
                else
                    pairs.push_back({a, b, entry.first, std::clamp(dot, -1., 1.)});
            }
        if (noise > 0) {
            std::vector<double> width(mesh.triangles.size());
            std::map<unsigned, std::set<unsigned>> vertices;
            for (unsigned t = 0; t < mesh.triangles.size(); ++t) {
                auto tri = mesh.triangles[t];
                double longest = 0;
                for (unsigned k = 0; k < 3; ++k) {
                    longest = std::max(
                        longest, (mesh.points[tri[k]] - mesh.points[tri[(k + 1) % 3]]).norm());
                    vertices[root(t)].insert(tri[k]);
                }
                width[t] = 2 * top.areas[t] / std::max(longest, 1e-12);
            }
            auto candidates = pairs;
            std::stable_sort(candidates.begin(), candidates.end(),
                             [](const Pair &a, const Pair &b) { return a.dot > b.dot; });
            for (const auto &pair : candidates) {
                unsigned a = root(pair.a), b = root(pair.b);
                double allowed = std::min(noise / std::max(width[pair.a], 1e-12) +
                                              noise / std::max(width[pair.b], 1e-12),
                                          pi / 180);
                if (a == b || std::acos(pair.dot) > allowed)
                    continue;
                auto both = vertices[a];
                both.insert(vertices[b].begin(), vertices[b].end());
                Vec center = Vec::Zero();
                for (auto v : both)
                    center += mesh.points[v];
                center /= double(both.size());
                Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
                for (auto v : both) {
                    Vec q = mesh.points[v] - center;
                    covariance += q * q.transpose();
                }
                Vec normal = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(covariance)
                                 .eigenvectors()
                                 .col(0);
                double worst = 0;
                for (auto v : both)
                    worst = std::max(worst, std::abs((mesh.points[v] - center).dot(normal)));
                if (worst <= tolerance + noise) {
                    parent[a] = b;
                    vertices[b] = std::move(both);
                    vertices.erase(a);
                }
            }
        }
        std::map<unsigned, unsigned> labels;
        for (unsigned t = 0; t < parent.size(); ++t)
            labels.emplace(root(t), 0);
        unsigned next = 0;
        for (auto &entry : labels)
            entry.second = next++;
        triangles.resize(next);
        verts.resize(next);
        normals.resize(next, Vec::Zero());
        centers.resize(next, Vec::Zero());
        areas.resize(next, 0);
        neighbors.resize(next);
        beside.resize(next);
        for (unsigned t = 0; t < parent.size(); ++t) {
            unsigned f = labels.at(root(t));
            facet_of.push_back(f);
            triangles[f].push_back(t);
            normals[f] += top.areas[t] * top.normals[t];
            centers[f] += top.areas[t] * top.centers[t];
            areas[f] += top.areas[t];
        }
        for (unsigned f = 0; f < next; ++f) {
            verts[f] = vertices(mesh, triangles[f]);
            double length =
                std::sqrt(normals[f].x() * normals[f].x() + normals[f].y() * normals[f].y() +
                          normals[f].z() * normals[f].z());
            // NumPy divides each component. Eigen's vector/scalar operator
            // multiplies by a reciprocal, which can flip an exact <= test.
            for (unsigned k = 0; k < 3; ++k) {
                normals[f][k] = normals[f][k] / length;
                centers[f][k] = centers[f][k] / areas[f];
            }
        }
        for (const auto &entry : edges)
            for (auto a : entry.second)
                for (auto b : entry.second) {
                    unsigned x = facet_of[a], y = facet_of[b];
                    if (x == y)
                        continue;
                    beside[x].push_back(y);
                }
        for (const auto &pair : pairs) {
            unsigned a = facet_of[pair.a], b = facet_of[pair.b];
            if (a == b || normals[a].dot(normals[b]) <= std::cos(40 * pi / 180))
                continue;
            neighbors[a].push_back(b);
            neighbors[b].push_back(a);
            shared[std::minmax(a, b)].push_back(pair.edge);
        }
        for (auto &list : neighbors) {
            FeatureNeighborOrder order;
            for (auto neighbor : list)
                order.add(neighbor);
            list = order.values();
        }
        // blends._islands walks edge_facets in edge order, including repeats.
        // Its bounded traversal makes that order observable on large islands.
    }
};
} // namespace stl2curves::detail

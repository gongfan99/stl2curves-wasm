#include "internal.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <nanoflann.hpp>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <tuple>

namespace stl2curves::detail {
namespace {
std::uint32_t u32(const std::uint8_t *p) {
    return std::uint32_t(p[0]) | (std::uint32_t(p[1]) << 8) | (std::uint32_t(p[2]) << 16) |
           (std::uint32_t(p[3]) << 24);
}
double f32(const std::uint8_t *p) {
    auto bits = u32(p);
    float f;
    std::memcpy(&f, &bits, 4);
    return f;
}
bool forward(const Tri &t, const Edge &e) {
    for (int k = 0; k < 3; ++k)
        if (t[k] == e.a && t[(k + 1) % 3] == e.b)
            return true;
    return false;
}
struct Cloud {
    const std::vector<Vec> &p;
    std::size_t kdtree_get_point_count() const { return p.size(); }
    double kdtree_get_pt(std::size_t i, std::size_t d) const { return p[i][d]; }
    template <class T> bool kdtree_get_bbox(T &) const { return false; }
};
} // namespace
Mesh read_stl(const Bytes &b) {
    if (b.empty())
        throw std::runtime_error("Empty STL input");
    Mesh m;
    std::size_t n = b.size() >= 84 ? u32(b.data() + 80) : 0;
    if (b.size() >= 84 && n <= (b.size() - 84) / 50 && 84 + n * 50 == b.size()) {
        m.points.reserve(n * 3);
        m.triangles.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            const auto *p = b.data() + 84 + i * 50 + 12;
            for (int k = 0; k < 3; ++k)
                m.points.emplace_back(f32(p + k * 12), f32(p + k * 12 + 4), f32(p + k * 12 + 8));
            auto v = std::uint32_t(i * 3);
            m.triangles.push_back({v, v + 1, v + 2});
        }
    } else {
        std::istringstream in(std::string(b.begin(), b.end()));
        std::string line;
        while (std::getline(in, line)) {
            std::istringstream row(line);
            std::string token;
            row >> token;
            if (token != "vertex")
                continue;
            double x, y, z;
            if (!(row >> x >> y >> z))
                throw std::runtime_error("Malformed ASCII STL vertex");
            m.points.emplace_back(x, y, z);
        }
        if (m.points.size() % 3)
            throw std::runtime_error("Incomplete ASCII STL triangle");
        for (std::uint32_t i = 0; i < m.points.size(); i += 3)
            m.triangles.push_back({i, i + 1, i + 2});
    }
    if (m.triangles.empty())
        throw std::runtime_error("STL contains no triangles");
    for (const auto &p : m.points)
        if (!p.allFinite() || p.cwiseAbs().maxCoeff() > 1e10)
            throw std::runtime_error("Nonfinite or out-of-range STL coordinate");
    weld(m);
    return m;
}
void weld(Mesh &m) {
    // Lexicographic, 0.0001 mm quantization, matching upstream load_stl.
    using Key = std::array<std::int64_t, 3>;
    std::map<Key, unsigned> keys;
    std::vector<Key> original;
    original.reserve(m.points.size());
    for (const auto &p : m.points) {
        if (!p.allFinite() || p.cwiseAbs().maxCoeff() > 1e10)
            throw std::runtime_error("Invalid mesh coordinate");
        Key k{std::int64_t(std::nearbyint(p.x() * 1e4)), std::int64_t(std::nearbyint(p.y() * 1e4)),
              std::int64_t(std::nearbyint(p.z() * 1e4))};
        original.push_back(k);
        keys.emplace(k, 0);
    }
    unsigned i = 0;
    for (auto &kv : keys)
        kv.second = i++;
    std::vector<Vec> points(keys.size());
    std::vector<unsigned> remap;
    for (std::size_t j = 0; j < original.size(); ++j) {
        auto v = keys.at(original[j]);
        points[v] = m.points[j];
        remap.push_back(v);
    }
    std::vector<Tri> tris;
    std::vector<unsigned> parts;
    for (std::size_t j = 0; j < m.triangles.size(); ++j) {
        auto t = m.triangles[j];
        for (auto &v : t) {
            if (v >= remap.size())
                throw std::runtime_error("Triangle index out of range");
            v = remap[v];
        }
        if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2])
            continue;
        tris.push_back(t);
        if (!m.parts.empty())
            parts.push_back(m.parts.at(j));
    }
    m.points = std::move(points);
    m.triangles = std::move(tris);
    m.parts = std::move(parts);
}
Topology topology(const Mesh &m) {
    Topology top;
    auto n = m.triangles.size();
    top.normals.resize(n);
    top.centers.resize(n);
    top.areas.resize(n);
    top.neighbors.resize(n);
    for (unsigned i = 0; i < n; ++i) {
        auto t = m.triangles[i];
        const Vec &a = m.points[t[0]];
        const Vec &b = m.points[t[1]];
        const Vec &c = m.points[t[2]];
        Vec cross = (b - a).cross(c - a);
        double norm = cross.norm();
        top.normals[i] = norm > 0 ? Vec(cross / norm) : Vec::Zero();
        top.areas[i] = norm / 2;
        top.centers[i] = (a + b + c) / 3;
        for (int k = 0; k < 3; ++k)
            top.edges[Edge(t[k], t[(k + 1) % 3])].push_back(i);
    }
    for (const auto &kv : top.edges)
        if (kv.second.size() == 2) {
            auto a = kv.second[0], b = kv.second[1];
            if (!m.parts.empty() && m.parts[a] != m.parts[b])
                continue;
            top.neighbors[a].push_back(b);
            top.neighbors[b].push_back(a);
        }
    for (auto &v : top.neighbors)
        std::sort(v.begin(), v.end());
    return top;
}
std::vector<std::pair<Vec, double>> symmetry_planes(const Mesh &mesh) {
    if (mesh.points.size() < 8)
        return {};
    Cloud cloud{mesh.points};
    using Tree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, Cloud>,
                                                     Cloud, 3, std::size_t>;
    Tree tree(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(16));
    tree.buildIndex();
    Vec center = Vec::Zero(), lo = mesh.points.front(), hi = lo;
    for (const auto &p : mesh.points) {
        center += p;
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
    }
    center /= double(mesh.points.size());
    Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
    for (const auto &p : mesh.points) {
        Vec q = p - center;
        covariance += q * q.transpose();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(covariance);
    std::vector<Vec> candidates{Vec::UnitX(),
                                Vec::UnitY(),
                                Vec::UnitZ(),
                                eigen.eigenvectors().col(0),
                                eigen.eigenvectors().col(1),
                                eigen.eigenvectors().col(2)};
    std::vector<std::pair<Vec, double>> out;
    const double tolerance = std::max(1e-4, (hi - lo).norm() * 0.0001);
    for (const auto &n : candidates) {
        double low = 1e100, high = -1e100;
        for (const auto &p : mesh.points) {
            low = std::min(low, p.dot(n));
            high = std::max(high, p.dot(n));
        }
        double offset = (low + high) / 2;
        unsigned good = 0, total = 0;
        std::size_t stride = std::max<std::size_t>(1, mesh.points.size() / 1024);
        for (std::size_t i = 0; i < mesh.points.size(); i += stride) {
            Vec reflected = mesh.points[i] - 2 * (n.dot(mesh.points[i]) - offset) * n;
            std::size_t nearest;
            double distance;
            if (tree.knnSearch(reflected.data(), 1, &nearest, &distance) == 1 &&
                distance <= tolerance * tolerance)
                ++good;
            ++total;
        }
        if (good * 100 >= total * 98)
            out.emplace_back(n, offset);
    }
    return out;
}
double signed_volume(const Mesh &m) {
    if (m.points.empty())
        return 0;
    const Vec origin = m.points.front();
    double sum = 0, correction = 0;
    for (const auto &t : m.triangles) {
        double value = (m.points[t[0]] - origin)
                           .dot((m.points[t[1]] - origin).cross(m.points[t[2]] - origin)) /
                       6;
        double y = value - correction, next = sum + y;
        correction = (next - sum) - y;
        sum = next;
    }
    return sum;
}
void repair_mesh(Mesh &m, Context &ctx) {
    auto before = m.triangles.size();
    std::set<std::pair<unsigned, Tri>> seen;
    std::vector<Tri> triangles;
    std::vector<unsigned> parts;
    for (unsigned i = 0; i < m.triangles.size(); ++i) {
        auto t = m.triangles[i], key = t;
        std::sort(key.begin(), key.end());
        unsigned part = m.parts.empty() ? 0 : m.parts[i];
        if (!seen.emplace(part, key).second)
            continue;
        if ((m.points[t[1]] - m.points[t[0]]).cross(m.points[t[2]] - m.points[t[0]]).squaredNorm() <
            1e-24)
            continue;
        triangles.push_back(t);
        if (!m.parts.empty())
            parts.push_back(part);
    }
    m.triangles = std::move(triangles);
    m.parts = std::move(parts);
    if (before != m.triangles.size())
        ctx.diagnostics.repairs.push_back("Removed " + std::to_string(before - m.triangles.size()) +
                                          " duplicate or zero-area triangles");
    auto top = topology(m);
    // Propagate orientation without crossing edges shared by more than two facets.
    std::vector<int> flip(m.triangles.size(), -1);
    std::size_t flipped = 0;
    for (unsigned root = 0; root < m.triangles.size(); ++root) {
        if (flip[root] >= 0)
            continue;
        std::vector<unsigned> queue{root};
        flip[root] = 0;
        for (std::size_t q = 0; q < queue.size(); ++q) {
            auto i = queue[q];
            const auto &t = m.triangles[i];
            for (int k = 0; k < 3; ++k) {
                Edge e(t[k], t[(k + 1) % 3]);
                const auto &fs = top.edges.at(e);
                if (fs.size() != 2)
                    continue;
                auto j = fs[0] == i ? fs[1] : fs[0];
                if (!m.parts.empty() && m.parts[i] != m.parts[j])
                    continue;
                int wanted = flip[i] ^ (forward(t, e) == forward(m.triangles[j], e));
                if (flip[j] < 0) {
                    flip[j] = wanted;
                    queue.push_back(j);
                }
            }
        }
    }
    for (unsigned i = 0; i < flip.size(); ++i)
        if (flip[i]) {
            std::swap(m.triangles[i][1], m.triangles[i][2]);
            ++flipped;
        }
    if (flipped)
        ctx.diagnostics.repairs.push_back("Reoriented " + std::to_string(flipped) + " triangles");
    top = topology(m);
    // Fill simple planar holes, traversing each rim opposite to the adjacent
    // facets. Ear clipping supports concave rims without a fan outside the hole.
    std::map<unsigned, std::vector<unsigned>> outgoing;
    for (const auto &kv : top.edges)
        if (kv.second.size() == 1) {
            const auto &t = m.triangles[kv.second[0]];
            unsigned a = kv.first.a, b = kv.first.b;
            if (forward(t, kv.first))
                outgoing[b].push_back(a);
            else
                outgoing[a].push_back(b);
        }
    std::set<unsigned> visited;
    unsigned filled = 0;
    for (const auto &item : outgoing) {
        unsigned start = item.first;
        if (visited.count(start))
            continue;
        std::vector<unsigned> loop;
        unsigned current = start;
        bool closed = false;
        for (unsigned step = 0; step < 129; ++step) {
            if (visited.count(current)) {
                closed = current == start;
                break;
            }
            visited.insert(current);
            loop.push_back(current);
            auto next = outgoing.find(current);
            if (next == outgoing.end() || next->second.size() != 1)
                break;
            current = next->second.front();
        }
        if (!closed || loop.size() < 3 || loop.size() > 128)
            continue;
        Vec normal = Vec::Zero();
        for (unsigned i = 0; i < loop.size(); ++i)
            normal += m.points[loop[i]].cross(m.points[loop[(i + 1) % loop.size()]]);
        if (normal.norm() < 1e-12)
            continue;
        normal.normalize();
        Vec origin = m.points[loop[0]];
        bool planar = true, nonflat = false;
        for (auto v : loop)
            if (std::abs((m.points[v] - origin).dot(normal)) >
                std::min(0.002, ctx.options.accuracy))
                planar = false;
        for (unsigned i = 0; i < loop.size(); ++i) {
            const auto &rim = top.edges.at(Edge(loop[i], loop[(i + 1) % loop.size()]));
            if (std::abs(top.normals[rim[0]].dot(normal)) < 0.95)
                nonflat = true;
        }
        if (!planar || !nonflat)
            continue;
        Vec e = normal.unitOrthogonal(), f = normal.cross(e);
        auto cross2 = [&](unsigned a, unsigned b, unsigned c) {
            Vec ab = m.points[b] - m.points[a], ac = m.points[c] - m.points[a];
            return ab.dot(e) * ac.dot(f) - ab.dot(f) * ac.dot(e);
        };
        std::vector<Tri> cap;
        auto rim = loop;
        while (rim.size() > 3) {
            bool found = false;
            for (unsigned k = 0; k < rim.size(); ++k) {
                auto a = rim[(k + rim.size() - 1) % rim.size()], b = rim[k],
                     c = rim[(k + 1) % rim.size()];
                if (cross2(a, b, c) < 1e-14)
                    continue;
                bool inside = false;
                for (auto v : rim)
                    if (v != a && v != b && v != c && cross2(a, b, v) >= -1e-14 &&
                        cross2(b, c, v) >= -1e-14 && cross2(c, a, v) >= -1e-14)
                        inside = true;
                if (inside)
                    continue;
                cap.push_back({a, b, c});
                rim.erase(rim.begin() + k);
                found = true;
                break;
            }
            if (!found)
                break;
        }
        if (rim.size() != 3)
            continue;
        cap.push_back({rim[0], rim[1], rim[2]});
        unsigned part = 0;
        if (!m.parts.empty())
            part = m.parts[top.edges.at(Edge(loop[0], loop[1]))[0]];
        m.triangles.insert(m.triangles.end(), cap.begin(), cap.end());
        if (!m.parts.empty())
            m.parts.insert(m.parts.end(), cap.size(), part);
        ++filled;
    }
    if (filled) {
        ctx.diagnostics.repairs.push_back("Filled " + std::to_string(filled) +
                                          " planar mesh holes");
        top = topology(m);
    }
    for (const auto &kv : top.edges) {
        if (kv.second.size() == 1)
            ++ctx.diagnostics.open_edges;
        if (kv.second.size() > 2)
            ++ctx.diagnostics.nonmanifold_edges;
    }
    if (ctx.diagnostics.open_edges)
        ctx.diagnostics.warnings.push_back(
            "Mesh has " + std::to_string(ctx.diagnostics.open_edges) + " boundary edges");
    if (ctx.diagnostics.nonmanifold_edges)
        ctx.diagnostics.warnings.push_back(
            "Mesh has " + std::to_string(ctx.diagnostics.nonmanifold_edges) + " nonmanifold edges");
}
std::vector<Mesh> split_bodies(const Mesh &m) {
    // Do not invert each component: inward cavity shells must retain their sign.
    auto top = topology(m);
    auto groups = components(top.neighbors);
    std::vector<Mesh> out;
    for (const auto &ids : groups) {
        Mesh part;
        std::map<unsigned, unsigned> local;
        for (auto i : ids) {
            Tri t = m.triangles[i];
            for (auto &v : t) {
                auto inserted = local.emplace(v, unsigned(part.points.size()));
                if (inserted.second)
                    part.points.push_back(m.points[v]);
                v = inserted.first->second;
            }
            part.triangles.push_back(t);
            if (!m.parts.empty())
                part.parts.push_back(m.parts[i]);
        }
        out.push_back(std::move(part));
    }
    return out;
}
void simplify_mesh(Mesh &m, double tolerance, Context &ctx) {
    if (tolerance <= 0)
        return;
    // Conservative vertex removal: retain original vertices, bound every removed
    // corner against surviving local triangles, preserve topology and sharp edges.
    auto top = topology(m);
    std::vector<std::vector<unsigned>> star(m.points.size());
    for (unsigned i = 0; i < m.triangles.size(); ++i)
        for (auto v : m.triangles[i])
            star[v].push_back(i);
    std::vector<bool> locked(m.points.size(), false), dead(m.triangles.size(), false);
    for (const auto &kv : top.edges) {
        bool crease = kv.second.size() != 2;
        if (!crease)
            crease = top.normals[kv.second[0]].dot(top.normals[kv.second[1]]) < 0.98;
        if (crease)
            locked[kv.first.a] = locked[kv.first.b] = true;
    }
    std::size_t removed = 0;
    for (const auto &kv : top.edges) {
        auto a = kv.first.a, b = kv.first.b;
        if (locked[a] || locked[b])
            continue;
        if ((m.points[a] - m.points[b]).norm() > tolerance)
            continue;
        // Link condition for a triangular two-manifold edge.
        std::set<unsigned> na, nb;
        for (auto i : star[a])
            for (auto v : m.triangles[i])
                if (v != a)
                    na.insert(v);
        for (auto i : star[b])
            for (auto v : m.triangles[i])
                if (v != b)
                    nb.insert(v);
        std::vector<unsigned> common;
        std::set_intersection(na.begin(), na.end(), nb.begin(), nb.end(),
                              std::back_inserter(common));
        if (common.size() != 2)
            continue;
        bool good = true;
        for (auto i : star[b]) {
            auto t = m.triangles[i];
            for (auto &v : t)
                if (v == b)
                    v = a;
            if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2])
                continue;
            Vec normal = (m.points[t[1]] - m.points[t[0]]).cross(m.points[t[2]] - m.points[t[0]]);
            if (normal.norm() < 1e-12 || normal.normalized().dot(top.normals[i]) < 0.99)
                good = false;
        }
        if (!good)
            continue;
        for (auto i : star[b]) {
            auto &t = m.triangles[i];
            for (auto &v : t)
                if (v == b)
                    v = a;
            if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2])
                dead[i] = true;
            for (auto v : t)
                locked[v] = true;
        }
        locked[a] = locked[b] = true;
        ++removed;
    }
    std::vector<Tri> kept;
    std::vector<unsigned> parts;
    for (unsigned i = 0; i < m.triangles.size(); ++i)
        if (!dead[i]) {
            kept.push_back(m.triangles[i]);
            if (!m.parts.empty())
                parts.push_back(m.parts[i]);
        }
    m.triangles = std::move(kept);
    m.parts = std::move(parts);
    ctx.diagnostics.repairs.push_back("Simplification removed " + std::to_string(removed) +
                                      " short-edge vertices");
}
} // namespace stl2curves::detail

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
#include "mesh_repair.hpp"
void repair_mesh(Mesh &m, Context &ctx) {
    const auto duplicates = remove_duplicates(m);
    if (duplicates.first)
        ctx.diagnostics.repairs.push_back(std::to_string(duplicates.first) +
                                          " duplicate triangles removed");
    if (duplicates.second)
        ctx.diagnostics.repairs.push_back(std::to_string(duplicates.second) +
                                          " zero-thickness wall triangles removed");
    if (m.triangles.empty())
        return;
    auto slivers = split_slivers(m);
    if (slivers)
        ctx.diagnostics.repairs.push_back(std::to_string(slivers) + " sliver triangles split away");
    Vec lo = m.points.front(), hi = lo;
    for (const auto &p : m.points) {
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
    }
    auto joined = close_cracks(m, std::max(.001, 1e-5 * (hi - lo).maxCoeff()));
    if (joined)
        ctx.diagnostics.repairs.push_back(std::to_string(joined) + " corners along cracks joined");
    const auto before_orientation = topology(m);
    bool has_open = std::any_of(before_orientation.edges.begin(), before_orientation.edges.end(),
                                [](const auto &e) { return e.second.size() == 1; });
    auto flaps = has_open ? remove_flaps(m) : 0;
    if (flaps)
        ctx.diagnostics.repairs.push_back(std::to_string(flaps) + " stray flap triangles removed");
    auto knots = untie_knots(m);
    if (knots)
        ctx.diagnostics.repairs.push_back(std::to_string(knots) +
                                          " knots of tiny triangles pulled to a point");
    auto flipped = orient_mesh(m);
    if (flipped)
        ctx.diagnostics.repairs.push_back(std::to_string(flipped) +
                                          " triangles turned to face the right way");
    auto filled = fill_holes(m);
    if (filled)
        ctx.diagnostics.repairs.push_back(std::to_string(filled) + " holes filled");
    if (flipped || filled)
        orient_mesh(m);
    const auto top = topology(m);
    unsigned open = 0, nonmanifold = 0;
    for (const auto &kv : top.edges) {
        open += kv.second.size() == 1;
        nonmanifold += kv.second.size() > 2;
    }
    ctx.diagnostics.open_edges += open;
    ctx.diagnostics.nonmanifold_edges += nonmanifold;
    if (open)
        ctx.diagnostics.warnings.push_back("Mesh has " + std::to_string(open) + " boundary edges");
    if (nonmanifold)
        ctx.diagnostics.warnings.push_back("Mesh has " + std::to_string(nonmanifold) +
                                           " nonmanifold edges");
}
#include "mesh_bodies.hpp"
std::vector<Mesh> split_bodies(const Mesh &m) {
    // Keep inward cavity shells with their host, and colour bodies that meet
    // along overused edges apart before sewing, as bodies.split_bodies does.
    auto groups = mesh_body_groups(m);
    std::vector<Mesh> out;
    for (const auto &ids : groups) {
        Mesh part;
        std::map<unsigned, unsigned> local;
        // Preserve the reference's vertex order. First-encounter compaction
        // changes the order of sorted corner samples used by fitting and can
        // change which numerically tied proposal wins on another body.
        std::set<unsigned> used;
        for (auto i : ids)
            for (auto v : m.triangles[i])
                used.insert(v);
        for (auto v : used) {
            local.emplace(v, unsigned(part.points.size()));
            part.points.push_back(m.points[v]);
        }
        for (auto i : ids) {
            Tri t = m.triangles[i];
            for (auto &v : t)
                v = local.at(v);
            part.triangles.push_back(t);
            if (!m.parts.empty())
                part.parts.push_back(m.parts[i]);
        }
        out.push_back(std::move(part));
    }
    return out;
}
#include "mesh_simplify.hpp"
void simplify_mesh(Mesh &m, double tolerance, Context &ctx) {
    if (tolerance <= 0 || m.triangles.empty())
        return;
    simplify_quadrics(m, tolerance, ctx);
}
} // namespace stl2curves::detail

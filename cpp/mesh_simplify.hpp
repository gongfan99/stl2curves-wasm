// Derived from stl2curves simplify.py, Copyright (c) 2026 Duvdivan. MIT license.
#pragma once
namespace {
double simplify_triangle_distance(const Vec &p, const Vec &a, const Vec &b, const Vec &c) {
    Vec ab = b - a, ac = c - a, ap = p - a, bp = p - b, cp = p - c;
    double d1 = ab.dot(ap), d2 = ac.dot(ap), d3 = ab.dot(bp), d4 = ac.dot(bp), d5 = ab.dot(cp),
           d6 = ac.dot(cp);
    double va = d3 * d6 - d5 * d4, vb = d5 * d2 - d1 * d6, vc = d1 * d4 - d3 * d2;
    double denominator = va + vb + vc;
    if (std::abs(denominator) <= 1e-300)
        denominator = 1e-300;
    double v = vb / denominator, w = vc / denominator;
    double distance = v >= 0 && w >= 0 && v + w <= 1 ? (p - (a + ab * v + ac * w)).norm()
                                                     : std::numeric_limits<double>::infinity();
    const Vec corners[]{a, b, c};
    for (unsigned k = 0; k < 3; ++k) {
        Vec e = corners[(k + 1) % 3] - corners[k];
        double t = std::clamp((p - corners[k]).dot(e) / std::max(e.squaredNorm(), 1e-300), 0., 1.);
        distance = std::min(distance, (p - (corners[k] + e * t)).norm());
    }
    return distance;
}
using MeshQuadric = std::array<double, 10>;
double evaluate_quadric(const MeshQuadric &q, const Vec &p) {
    double x = p.x(), y = p.y(), z = p.z();
    return q[0] * x * x + 2 * q[1] * x * y + 2 * q[2] * x * z + 2 * q[3] * x + q[4] * y * y +
           2 * q[5] * y * z + 2 * q[6] * y + q[7] * z * z + 2 * q[8] * z + q[9];
}
void simplify_quadrics(Mesh &m, double error, Context &ctx) {
    std::vector<MeshQuadric> Q(m.points.size(), MeshQuadric{});
    std::vector<double> W(m.points.size(), 0);
    // Accumulate in the reference's corner-column order.
    std::vector<MeshQuadric> coefficients(m.triangles.size());
    std::vector<double> area(m.triangles.size());
    for (unsigned i = 0; i < m.triangles.size(); ++i) {
        auto t = m.triangles[i];
        const auto &a = m.points[t[0]], &b = m.points[t[1]], &c = m.points[t[2]];
        Vec cross = (b - a).cross(c - a);
        double norm = cross.norm();
        Vec n = cross / std::max(norm, 1e-300);
        double x = n.x(), y = n.y(), z = n.z(), d = -n.dot(a);
        area[i] = norm / 2;
        coefficients[i] = {x * x, x * y, x * z, x * d, y * y, y * z, y * d, z * z, z * d, d * d};
        for (auto &value : coefficients[i])
            value *= area[i];
    }
    for (unsigned k = 0; k < 3; ++k)
        for (unsigned i = 0; i < m.triangles.size(); ++i) {
            auto v = m.triangles[i][k];
            for (unsigned j = 0; j < 10; ++j)
                Q[v][j] += coefficients[i][j];
            W[v] += area[i];
        }
    coefficients = std::vector<MeshQuadric>{};
    area = std::vector<double>{};
    struct Attached {
        unsigned vertex, triangle;
    };
    std::vector<Attached> gone;
    unsigned total = 0, passes = 0;
    for (; passes < 200; ++passes) {
        const auto top = topology(m);
        std::map<Edge, std::vector<unsigned>> edges(top.edges.begin(), top.edges.end());
        std::vector<bool> locked(m.points.size(), false);
        std::vector<std::vector<unsigned>> star(m.points.size()), neighbors(m.points.size()),
            attached(m.triangles.size());
        for (unsigned i = 0; i < m.triangles.size(); ++i)
            for (auto v : m.triangles[i])
                star[v].push_back(i);
        for (unsigned i = 0; i < gone.size(); ++i)
            attached[gone[i].triangle].push_back(i);
        for (const auto &entry : edges) {
            auto a = entry.first.a, b = entry.first.b;
            neighbors[a].push_back(b);
            neighbors[b].push_back(a);
            if (entry.second.size() != 2 ||
                (!m.parts.empty() && m.parts[entry.second[0]] != m.parts[entry.second[1]]))
                locked[a] = locked[b] = true;
        }
        for (auto &nb : neighbors)
            std::sort(nb.begin(), nb.end());
        struct Collapse {
            unsigned src, dst;
            double cost;
        };
        std::vector<Collapse> candidates;
        for (bool reverse : {false, true})
            for (const auto &entry : edges) {
                auto a = entry.first.a, b = entry.first.b;
                if (entry.second.size() != 2 || locked[a] || locked[b])
                    continue;
                if (reverse)
                    std::swap(a, b);
                MeshQuadric q{};
                for (unsigned k = 0; k < 10; ++k)
                    q[k] = Q[a][k] + Q[b][k];
                double cost =
                    std::max(evaluate_quadric(q, m.points[b]) / std::max(W[a] + W[b], 1e-300), 0.);
                if (cost <= error * error)
                    candidates.push_back({a, b, cost});
            }
        std::stable_sort(candidates.begin(), candidates.end(),
                         [](const Collapse &a, const Collapse &b) { return a.cost < b.cost; });
        if (candidates.empty())
            break;
        struct Chosen {
            Collapse edge;
            std::vector<Attached> points;
        };
        std::vector<Chosen> chosen;
        std::vector<bool> touched(m.points.size(), false), reattached(gone.size(), false);
        for (auto candidate : candidates) {
            auto src = candidate.src, dst = candidate.dst;
            if (touched[src] || touched[dst])
                continue;
            const auto &ns = neighbors[src], &nd = neighbors[dst];
            if (std::any_of(ns.begin(), ns.end(), [&](unsigned v) { return touched[v]; }))
                continue;
            unsigned common = 0;
            auto a = ns.begin(), b = nd.begin();
            while (a != ns.end() && b != nd.end()) {
                if (*a == *b) {
                    ++common;
                    ++a;
                    ++b;
                } else if (*a < *b)
                    ++a;
                else
                    ++b;
            }
            if (common != 2)
                continue;
            std::vector<std::pair<unsigned, Tri>> surviving;
            std::vector<unsigned> points{src};
            bool good = true;
            for (auto i : star[src]) {
                for (auto index : attached[i])
                    points.push_back(gone[index].vertex);
                auto t = m.triangles[i];
                if (std::find(t.begin(), t.end(), dst) != t.end())
                    continue;
                for (auto &v : t)
                    if (v == src)
                        v = dst;
                const auto &p = m.points[t[0]], &q = m.points[t[1]], &r = m.points[t[2]];
                Vec n = (q - p).cross(r - p);
                double norm = n.norm(),
                       longest = std::max({(q - p).norm(), (r - q).norm(), (p - r).norm()});
                double turn = top.normals[i].dot(n) / std::max(norm, 1e-300);
                double shape = norm / 2 / std::max(longest * longest, 1e-300);
                if (turn < .9396926207859084 || shape < .02)
                    good = false;
                surviving.emplace_back(i, t);
            }
            if (!good || surviving.empty())
                continue;
            Chosen move{candidate, {}};
            for (auto v : points) {
                double nearest = std::numeric_limits<double>::infinity();
                unsigned best = 0;
                for (const auto &entry : surviving) {
                    const auto &t = entry.second;
                    double d = simplify_triangle_distance(m.points[v], m.points[t[0]],
                                                          m.points[t[1]], m.points[t[2]]);
                    if (d < nearest) {
                        nearest = d;
                        best = entry.first;
                    }
                }
                if (nearest > error) {
                    good = false;
                    break;
                }
                move.points.push_back({v, best});
            }
            if (!good)
                continue;
            touched[src] = true;
            for (auto v : ns)
                touched[v] = true;
            chosen.push_back(std::move(move));
        }
        if (chosen.empty())
            break;
        for (const auto &move : chosen) {
            unsigned src = move.edge.src, dst = move.edge.dst;
            for (auto i : star[src]) {
                for (auto &v : m.triangles[i])
                    if (v == src)
                        v = dst;
                for (auto index : attached[i])
                    reattached[index] = true;
            }
            for (unsigned k = 0; k < 10; ++k)
                Q[dst][k] += Q[src][k];
            W[dst] += W[src];
        }
        std::vector<unsigned> renumber(m.triangles.size(), ~0u), labels;
        std::vector<Tri> kept;
        for (unsigned i = 0; i < m.triangles.size(); ++i) {
            auto t = m.triangles[i];
            if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2])
                continue;
            renumber[i] = unsigned(kept.size());
            kept.push_back(t);
            if (!m.parts.empty())
                labels.push_back(m.parts[i]);
        }
        std::vector<Attached> next;
        for (unsigned i = 0; i < gone.size(); ++i)
            if (!reattached[i])
                next.push_back({gone[i].vertex, renumber[gone[i].triangle]});
        for (const auto &move : chosen)
            for (auto p : move.points)
                next.push_back({p.vertex, renumber[p.triangle]});
        gone = std::move(next);
        m.triangles = std::move(kept);
        m.parts = std::move(labels);
        total += unsigned(chosen.size());
        if (chosen.size() < .002 * m.triangles.size()) {
            ++passes;
            break;
        }
    }
    ctx.diagnostics.repairs.push_back("Simplification removed " + std::to_string(total) +
                                      " corners in " + std::to_string(passes) + " passes");
}
} // namespace

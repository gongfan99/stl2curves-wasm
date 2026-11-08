// Derived from stl2curves repair.py, Copyright (c) 2026 Duvdivan. MIT license.
// Private mesh repair helpers; included by mesh.cpp inside detail's namespace.
#pragma once
namespace {
bool repair_edge_forward(const Tri &t, const Edge &e) {
    for (unsigned k = 0; k < 3; ++k)
        if (t[k] == e.a && t[(k + 1) % 3] == e.b)
            return true;
    return false;
}
void keep_triangles(Mesh &m, const std::vector<bool> &keep) {
    std::vector<Tri> tris;
    std::vector<unsigned> labels;
    for (unsigned i = 0; i < m.triangles.size(); ++i)
        if (keep[i]) {
            tris.push_back(m.triangles[i]);
            if (!m.parts.empty())
                labels.push_back(m.parts[i]);
        }
    m.triangles = std::move(tris);
    m.parts = std::move(labels);
}
std::pair<unsigned, unsigned> remove_duplicates(Mesh &m) {
    using Key = std::pair<unsigned, Tri>;
    std::map<Key, std::array<std::vector<unsigned>, 2>> groups;
    std::vector<bool> keep(m.triangles.size(), false);
    for (unsigned i = 0; i < m.triangles.size(); ++i) {
        auto t = m.triangles[i], key = t;
        if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2])
            continue;
        auto first = unsigned(std::min_element(t.begin(), t.end()) - t.begin());
        bool parity = t[(first + 1) % 3] < t[(first + 2) % 3];
        std::sort(key.begin(), key.end());
        groups[{m.parts.empty() ? 0 : m.parts[i], key}][parity].push_back(i);
    }
    unsigned same = 0, opposite = 0;
    for (const auto &entry : groups) {
        const auto &up = entry.second[1], &down = entry.second[0];
        auto pairs = std::min(up.size(), down.size());
        opposite += unsigned(2 * pairs);
        const auto &rest = up.size() > down.size() ? up : down;
        if (rest.size() > pairs) {
            keep[rest[pairs]] = true;
            same += unsigned(rest.size() - pairs - 1);
        }
    }
    keep_triangles(m, keep);
    return {same, opposite};
}
unsigned close_cracks(Mesh &m, double tolerance) {
    auto top = topology(m);
    std::set<unsigned> ends;
    for (const auto &entry : top.edges)
        if (entry.second.size() == 1) {
            ends.insert(entry.first.a);
            ends.insert(entry.first.b);
        }
    std::vector<unsigned> ids(ends.begin(), ends.end()), parent(ids.size());
    std::iota(parent.begin(), parent.end(), 0);
    auto root = [&](unsigned i) {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    };
    using Cell = std::array<std::int64_t, 3>;
    std::map<Cell, std::vector<unsigned>> cells;
    for (unsigned i = 0; i < ids.size(); ++i) {
        Vec q = m.points[ids[i]] / tolerance;
        cells[{std::int64_t(std::floor(q.x())), std::int64_t(std::floor(q.y())),
               std::int64_t(std::floor(q.z()))}]
            .push_back(i);
    }
    // Match Python's lexicographic cell traversal and representative selection.
    for (const auto &entry : cells)
        for (int x = -1; x <= 1; ++x)
            for (int y = -1; y <= 1; ++y)
                for (int z = -1; z <= 1; ++z) {
                    auto c = entry.first;
                    c[0] += x;
                    c[1] += y;
                    c[2] += z;
                    auto found = cells.find(c);
                    if (found == cells.end())
                        continue;
                    for (auto j : found->second)
                        for (auto i : entry.second)
                            if (i < j &&
                                (m.points[ids[i]] - m.points[ids[j]]).norm() <= tolerance) {
                                auto a = root(i), b = root(j);
                                if (a != b)
                                    parent[a] = b;
                            }
                }
    unsigned joined = 0;
    std::vector<unsigned> remap(m.points.size());
    std::iota(remap.begin(), remap.end(), 0);
    for (unsigned i = 0; i < ids.size(); ++i) {
        auto r = root(i);
        joined += r != i;
        remap[ids[i]] = ids[r];
    }
    if (joined) {
        for (auto &t : m.triangles)
            for (auto &v : t)
                v = remap[v];
        remove_duplicates(m);
    }
    return joined;
}
bool ray_inside(const Mesh &m, const std::vector<unsigned> &triangles, const Vec &q) {
    const Vec ray = Vec(.5773, .5774, .5775).normalized();
    unsigned hits = 0;
    for (auto i : triangles) {
        auto t = m.triangles[i];
        Vec a = m.points[t[0]], b = m.points[t[1]], c = m.points[t[2]];
        if ((a.cwiseMax(b).cwiseMax(c).array() < q.array() - 1e-6).any())
            continue;
        Vec e1 = b - a, e2 = c - a, p = ray.cross(e2);
        double det = e1.dot(p);
        if (std::abs(det) <= 1e-15)
            continue;
        Vec s = q - a, cross = s.cross(e1);
        double u = s.dot(p) / det, v = ray.dot(cross) / det;
        if (u >= 0 && v >= 0 && u + v <= 1 && e2.dot(cross) / det > 1e-9)
            ++hits;
    }
    return hits % 2 != 0;
}
unsigned orient_mesh(Mesh &m) {
    auto top = topology(m);
    std::vector<int> flip(m.triangles.size(), -1);
    std::vector<std::vector<unsigned>> pieces;
    for (unsigned start = 0; start < m.triangles.size(); ++start) {
        if (flip[start] >= 0)
            continue;
        std::vector<unsigned> piece, stack{start};
        flip[start] = 0;
        while (!stack.empty()) {
            auto i = stack.back();
            stack.pop_back();
            piece.push_back(i);
            for (auto j : top.neighbors[i])
                if (flip[j] < 0) {
                    Edge shared(0, 0);
                    const auto &t = m.triangles[i];
                    for (unsigned k = 0; k < 3; ++k) {
                        Edge e(t[k], t[(k + 1) % 3]);
                        const auto &fs = top.edges.at(e);
                        if (fs.size() == 2 && (fs[0] == j || fs[1] == j)) {
                            shared = e;
                            break;
                        }
                    }
                    flip[j] = flip[i] ^ (repair_edge_forward(t, shared) ==
                                         repair_edge_forward(m.triangles[j], shared));
                    stack.push_back(j);
                }
        }
        std::sort(piece.begin(), piece.end());
        unsigned changed = 0;
        for (auto i : piece)
            changed += flip[i] != 0;
        if (2 * changed > piece.size())
            for (auto i : piece)
                flip[i] ^= 1;
        pieces.push_back(std::move(piece));
    }
    unsigned count = 0;
    for (unsigned i = 0; i < flip.size(); ++i)
        if (flip[i]) {
            std::swap(m.triangles[i][0], m.triangles[i][2]);
            ++count;
        }
    struct Closed {
        std::vector<unsigned> triangles, corners;
        Vec lo, hi;
        double volume = 0;
    };
    std::vector<Closed> closed;
    for (const auto &piece : pieces) {
        if (piece.size() < 4)
            continue;
        bool open = false;
        std::set<unsigned> corners;
        double v = 0;
        for (auto i : piece) {
            auto t = m.triangles[i];
            for (unsigned k = 0; k < 3; ++k) {
                open = open || top.edges.at(Edge(t[k], t[(k + 1) % 3])).size() == 1;
                corners.insert(t[k]);
            }
            v += m.points[t[0]].dot(m.points[t[1]].cross(m.points[t[2]]));
        }
        if (open)
            continue;
        Closed c;
        c.triangles = piece;
        c.corners.assign(corners.begin(), corners.end());
        c.lo = c.hi = m.points[c.corners.front()];
        for (auto i : c.corners) {
            c.lo = c.lo.cwiseMin(m.points[i]);
            c.hi = c.hi.cwiseMax(m.points[i]);
        }
        c.volume = v;
        closed.push_back(std::move(c));
    }
    // Only an inward piece not contained by another closed piece is inverted.
    // Keep the file's inward cavity facing; depth alone cannot distinguish a
    // cavity from a separate body buried in another body.
    for (unsigned k = 0; k < closed.size(); ++k) {
        auto &c = closed[k];
        if (c.volume > 0)
            continue;
        bool contained = false;
        for (unsigned j = 0; j < closed.size() && !contained; ++j) {
            if (j == k)
                continue;
            const auto &host = closed[j];
            bool inside = true;
            auto stride = std::max<std::size_t>(1, c.corners.size() / 12);
            for (std::size_t i = 0; i < c.corners.size() && inside; i += stride) {
                const auto &q = m.points[c.corners[i]];
                inside = (q.array() >= host.lo.array()).all() &&
                         (q.array() <= host.hi.array()).all() && ray_inside(m, host.triangles, q);
            }
            contained = inside;
        }
        if (!contained)
            for (auto i : c.triangles) {
                std::swap(m.triangles[i][0], m.triangles[i][2]);
                ++count;
            }
    }
    return count;
}
// Liepa-style minimum-area triangulation. For small holes this also works
// in 3D, where projection and a planar-only test reject valid twisted rims.
std::vector<Tri> triangulate_hole(const Mesh &m, const std::vector<unsigned> &loop,
                                  const std::set<Edge> &taken) {
    auto n = loop.size();
    if (n < 3)
        return {};
    if (n == 3)
        return {{loop[0], loop[1], loop[2]}};
    auto free = [&](unsigned i, unsigned j) { return !taken.count(Edge(loop[i], loop[j])); };
    if (n <= 60) {
        std::vector<std::vector<double>> best(n, std::vector<double>(n, 0));
        std::vector<std::vector<unsigned>> split(n, std::vector<unsigned>(n));
        for (unsigned gap = 2; gap < n; ++gap)
            for (unsigned i = 0; i + gap < n; ++i) {
                auto j = i + gap;
                best[i][j] = std::numeric_limits<double>::infinity();
                for (unsigned k = i + 1; k < j; ++k) {
                    if ((k != i + 1 && !free(i, k)) || (k != j - 1 && !free(k, j)))
                        continue;
                    double cost = best[i][k] + best[k][j] +
                                  .5 * (m.points[loop[k]] - m.points[loop[i]])
                                           .cross(m.points[loop[j]] - m.points[loop[i]])
                                           .norm();
                    if (cost < best[i][j]) {
                        best[i][j] = cost;
                        split[i][j] = k;
                    }
                }
            }
        if (!std::isfinite(best[0][n - 1]))
            return {};
        std::vector<Tri> out;
        std::vector<std::pair<unsigned, unsigned>> stack{{0, unsigned(n - 1)}};
        while (!stack.empty()) {
            auto ij = stack.back();
            stack.pop_back();
            auto i = ij.first, j = ij.second;
            if (j - i < 2)
                continue;
            auto k = split[i][j];
            out.push_back({loop[i], loop[k], loop[j]});
            stack.push_back({i, k});
            stack.push_back({k, j});
        }
        return out;
    }
    Vec mid = Vec::Zero();
    for (auto i : loop)
        mid += m.points[i];
    mid /= double(n);
    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    for (auto i : loop) {
        Vec d = m.points[i] - mid;
        cov += d * d.transpose();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(cov);
    Vec e = eigen.eigenvectors().col(2), f = eigen.eigenvectors().col(1);
    auto cross = [&](unsigned a, unsigned b, unsigned c) {
        Vec ab = m.points[loop[b]] - m.points[loop[a]], ac = m.points[loop[c]] - m.points[loop[a]];
        return ab.dot(e) * ac.dot(f) - ab.dot(f) * ac.dot(e);
    };
    double signed_area = 0;
    for (unsigned i = 0; i < n; ++i) {
        Vec a = m.points[loop[i]] - mid, b = m.points[loop[(i + 1) % n]] - mid;
        signed_area += a.dot(e) * b.dot(f) - b.dot(e) * a.dot(f);
    }
    double sign = signed_area > 0 ? 1 : -1;
    std::vector<unsigned> rim(n);
    std::iota(rim.begin(), rim.end(), 0);
    if (sign < 0)
        std::reverse(rim.begin(), rim.end());
    std::vector<Tri> out;
    while (rim.size() > 3) {
        bool clipped = false;
        for (unsigned k = 0; k < rim.size(); ++k) {
            auto a = rim[(k + rim.size() - 1) % rim.size()], b = rim[k],
                 c = rim[(k + 1) % rim.size()];
            if (cross(a, b, c) <= 0 || !free(a, c))
                continue;
            bool inside = false;
            for (auto v : rim)
                if (v != a && v != b && v != c && cross(a, b, v) >= 0 && cross(b, c, v) >= 0 &&
                    cross(c, a, v) >= 0)
                    inside = true;
            if (inside)
                continue;
            out.push_back(sign > 0 ? Tri{loop[a], loop[b], loop[c]}
                                   : Tri{loop[c], loop[b], loop[a]});
            rim.erase(rim.begin() + k);
            clipped = true;
            break;
        }
        if (!clipped)
            return {};
    }
    out.push_back(sign > 0 ? Tri{loop[rim[0]], loop[rim[1]], loop[rim[2]]}
                           : Tri{loop[rim[2]], loop[rim[1]], loop[rim[0]]});
    return out;
}
unsigned fill_holes(Mesh &m) {
    auto top = topology(m);
    std::map<unsigned, std::vector<unsigned>> next;
    std::vector<std::pair<unsigned, unsigned>> open;
    std::set<Edge> taken;
    for (const auto &entry : top.edges)
        taken.insert(entry.first);
    // _open_edges preserves _edges order: the first edge of every triangle,
    // then the second and third. Loop start changes minimum-area tie breaks.
    for (unsigned k = 0; k < 3; ++k)
        for (const auto &t : m.triangles) {
            auto a = t[k], b = t[(k + 1) % 3];
            if (top.edges.at(Edge(a, b)).size() == 1) {
                open.emplace_back(a, b);
                next[a].push_back(b);
            }
        }
    std::set<std::pair<unsigned, unsigned>> used;
    unsigned filled = 0;
    for (auto edge : open) {
        if (used.count(edge))
            continue;
        std::vector<unsigned> loop{edge.first};
        unsigned at = edge.second;
        used.insert(edge);
        while (at != loop.front() && loop.size() < 100000) {
            loop.push_back(at);
            auto &choices = next[at];
            auto found = std::find_if(choices.begin(), choices.end(),
                                      [&](unsigned b) { return !used.count({at, b}); });
            if (found == choices.end())
                break;
            used.emplace(at, *found);
            at = *found;
        }
        if (at != loop.front())
            continue;
        unsigned part = m.parts.empty() ? 0 : m.parts[top.edges.at(Edge(loop[0], loop[1]))[0]];
        std::reverse(loop.begin(), loop.end());
        auto cap = triangulate_hole(m, loop, taken);
        if (cap.empty())
            cap = triangulate_hole(m, loop, {});
        if (cap.empty())
            continue;
        m.triangles.insert(m.triangles.end(), cap.begin(), cap.end());
        if (!m.parts.empty())
            m.parts.insert(m.parts.end(), cap.size(), part);
        for (auto t : cap)
            for (unsigned k = 0; k < 3; ++k)
                taken.emplace(t[k], t[(k + 1) % 3]);
        ++filled;
    }
    return filled;
}
unsigned split_slivers(Mesh &m) {
    unsigned total = 0;
    for (unsigned round = 0; round < 20; ++round) {
        std::map<std::pair<unsigned, unsigned>, unsigned> directed;
        for (unsigned i = 0; i < m.triangles.size(); ++i)
            for (unsigned k = 0; k < 3; ++k) {
                const auto &t = m.triangles[i];
                directed[{t[k], t[(k + 1) % 3]}] = i;
            }
        std::vector<bool> keep(m.triangles.size(), true), touched(m.triangles.size(), false);
        std::set<Edge> made;
        std::vector<Tri> add;
        std::vector<unsigned> labels;
        unsigned count = 0;
        // Python visits flags by corner, then by triangle.
        for (unsigned k = 0; k < 3; ++k)
            for (unsigned i = 0; i < m.triangles.size(); ++i) {
                auto t = m.triangles[i];
                unsigned v = t[k], u = t[(k + 1) % 3], w = t[(k + 2) % 3];
                Vec e = m.points[w] - m.points[u], q = m.points[v] - m.points[u];
                double length = std::max(e.norm(), 1e-12), s = q.dot(e) / (length * length);
                double h = q.cross(e).norm() / length;
                if (h > 1e-5 * length || h >= 1e-4 || s <= 0 || s >= 1)
                    continue;
                auto across = directed.find({w, u});
                if (across == directed.end())
                    continue;
                unsigned j = across->second;
                if (touched[i] || touched[j] || (!m.parts.empty() && m.parts[i] != m.parts[j]))
                    continue;
                const auto &other = m.triangles[j];
                unsigned x = other[0] + other[1] + other[2] - u - w;
                if (x != v &&
                    (directed.count({v, x}) || directed.count({x, v}) || made.count(Edge(v, x))))
                    continue;
                made.emplace(v, x);
                touched[i] = touched[j] = true;
                keep[i] = keep[j] = false;
                ++count;
                if (x != v) {
                    add.push_back({w, v, x});
                    add.push_back({v, u, x});
                    if (!m.parts.empty()) {
                        labels.push_back(m.parts[i]);
                        labels.push_back(m.parts[i]);
                    }
                }
            }
        if (!count)
            break;
        keep_triangles(m, keep);
        m.triangles.insert(m.triangles.end(), add.begin(), add.end());
        m.parts.insert(m.parts.end(), labels.begin(), labels.end());
        total += count;
    }
    return total;
}
unsigned remove_flaps(Mesh &m) {
    auto top = topology(m);
    std::vector<bool> keep(m.triangles.size(), true);
    unsigned dropped = 0;
    std::map<Edge, std::vector<unsigned>> edges(top.edges.begin(), top.edges.end());
    for (const auto &entry : edges)
        if (entry.second.size() >= 3) {
            // _edges lists edge 0 for every triangle, then edge 1, then edge 2.
            auto members = entry.second;
            auto position = [&](unsigned i) {
                const auto &t = m.triangles[i];
                for (unsigned k = 0; k < 3; ++k)
                    if (Edge(t[k], t[(k + 1) % 3]) == entry.first)
                        return std::make_pair(k, i);
                return std::make_pair(3u, i);
            };
            std::sort(members.begin(), members.end(),
                      [&](unsigned a, unsigned b) { return position(a) < position(b); });
            for (unsigned i = 0; i < members.size(); ++i)
                for (unsigned j = i + 1; j < members.size(); ++j) {
                    auto a = members[i], b = members[j];
                    if (keep[a] && keep[b] &&
                        repair_edge_forward(m.triangles[a], entry.first) ==
                            repair_edge_forward(m.triangles[b], entry.first) &&
                        std::abs(top.normals[a].dot(top.normals[b])) > .9999) {
                        keep[top.areas[a] < top.areas[b] ? a : b] = false;
                        ++dropped;
                    }
                }
        }
    if (dropped)
        keep_triangles(m, keep);
    return dropped;
}
unsigned untie_knots(Mesh &m) {
    const auto top = topology(m);
    std::set<unsigned> over;
    unsigned over_count = 0;
    std::map<unsigned, std::set<unsigned>> adjacent;
    for (const auto &e : top.edges) {
        if (e.second.size() > 2) {
            over.insert(e.first.a);
            over.insert(e.first.b);
            ++over_count;
        }
        if ((m.points[e.first.a] - m.points[e.first.b]).norm() <= .005) {
            adjacent[e.first.a].insert(e.first.b);
            adjacent[e.first.b].insert(e.first.a);
        }
    }
    if (!over_count)
        return 0;
    std::set<unsigned> seen;
    std::vector<std::vector<unsigned>> knots;
    for (auto start : over) {
        if (seen.count(start) || !adjacent.count(start))
            continue;
        std::vector<unsigned> knot, stack{start};
        seen.insert(start);
        while (!stack.empty()) {
            auto v = stack.back();
            stack.pop_back();
            knot.push_back(v);
            for (auto w : adjacent[v])
                if (!seen.count(w) && (m.points[w] - m.points[start]).norm() <= .02) {
                    seen.insert(w);
                    stack.push_back(w);
                }
        }
        if (knot.size() > 1)
            knots.push_back(std::move(knot));
    }
    if (knots.empty())
        return 0;
    Mesh original = m;
    std::vector<double> weight(m.points.size(), 0), width(m.triangles.size());
    for (unsigned i = 0; i < m.triangles.size(); ++i) {
        auto t = m.triangles[i];
        double longest = 0;
        for (unsigned k = 0; k < 3; ++k) {
            weight[t[k]] += top.areas[i];
            longest = std::max(longest, (m.points[t[k]] - m.points[t[(k + 1) % 3]]).norm());
        }
        width[i] = 2 * top.areas[i] / std::max(longest, 1e-300);
    }
    auto counts = [&](const std::vector<Tri> &tris, const std::set<unsigned> &at) {
        std::map<Edge, unsigned> edges;
        for (auto t : tris)
            for (unsigned k = 0; k < 3; ++k) {
                auto a = t[k], b = t[(k + 1) % 3];
                if (at.count(a) || at.count(b))
                    ++edges[Edge(a, b)];
            }
        unsigned over = 0, open = 0;
        for (const auto &e : edges) {
            over += e.second > 2;
            open += e.second == 1;
        }
        return std::make_pair(over, open);
    };
    std::vector<bool> alive(m.triangles.size(), true);
    unsigned pulled = 0;
    for (const auto &knot : knots) {
        auto keep = *std::max_element(knot.begin(), knot.end(), [&](unsigned a, unsigned b) {
            return weight[a] < weight[b];
        });
        std::set<unsigned> mine(knot.begin(), knot.end());
        std::vector<unsigned> near;
        Mesh after;
        after.points = m.points;
        std::vector<Tri> before;
        bool turned = false;
        for (unsigned i = 0; i < m.triangles.size(); ++i) {
            auto t = m.triangles[i];
            if (!alive[i] || (!mine.count(t[0]) && !mine.count(t[1]) && !mine.count(t[2])))
                continue;
            near.push_back(i);
            before.push_back(t);
            for (auto &v : t)
                if (mine.count(v))
                    v = keep;
            after.triangles.push_back(t);
            if (t[0] != t[1] && t[1] != t[2] && t[0] != t[2] && t != m.triangles[i]) {
                Vec n = (m.points[t[1]] - m.points[t[0]]).cross(m.points[t[2]] - m.points[t[0]]);
                turned = turned || (n.dot(top.normals[i]) <= 0 && width[i] > .005);
            }
        }
        auto changed = after.triangles;
        remove_duplicates(after);
        auto was = counts(before, mine), now = counts(after.triangles, {keep});
        if (turned || now.first >= was.first || now.second > was.second)
            continue;
        std::set<Tri> survivors(after.triangles.begin(), after.triangles.end());
        for (unsigned i = 0; i < near.size(); ++i) {
            if (survivors.erase(changed[i]))
                m.triangles[near[i]] = changed[i];
            else
                alive[near[i]] = false;
        }
        ++pulled;
    }
    keep_triangles(m, alive);
    unsigned remaining = 0;
    for (const auto &e : topology(m).edges)
        remaining += e.second.size() > 2;
    if (remaining > .25 * over_count) {
        m = std::move(original);
        return 0;
    }
    return pulled;
}
} // namespace

// Derived from stl2curves bodies.py, Copyright (c) 2026 Duvdivan. MIT license.
#pragma once
namespace {
struct BodyUnion {
    std::vector<unsigned> parent;
    explicit BodyUnion(unsigned n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
    unsigned find(unsigned i) {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    }
    void join(unsigned a, unsigned b) {
        a = find(a);
        b = find(b);
        if (a != b)
            parent[a] = b;
    }
};
std::vector<std::pair<unsigned, unsigned>>
material_partners(const Mesh &m, const std::vector<unsigned> &triangles, const Edge &edge) {
    Vec d = m.points[edge.b] - m.points[edge.a];
    if (d.norm() == 0)
        return {};
    d.normalize();
    struct Side {
        Vec w;
        double turn, angle;
        bool valid, forward;
    };
    std::vector<Side> sides;
    for (auto i : triangles) {
        auto t = m.triangles[i];
        Vec a = m.points[t[0]], b = m.points[t[1]], c = m.points[t[2]];
        Vec n = (b - a).cross(c - a), third = Vec::Zero();
        for (auto v : t)
            if (v != edge.a && v != edge.b)
                third = m.points[v] - m.points[edge.a];
        Vec w = third - third.dot(d) * d;
        bool valid = n.norm() >= 1e-12 && w.norm() >= 1e-9;
        if (!valid) {
            n = Vec::UnitX();
            w = d.cross(Vec::UnitX()) + Vec::Constant(.001);
        }
        n.normalize();
        w.normalize();
        bool forward = false;
        for (unsigned k = 0; k < 3; ++k)
            if (t[k] == edge.a)
                forward = t[(k + 1) % 3] == edge.b;
        sides.push_back({w, -n.dot(d.cross(w)) > 0 ? 1. : -1., 0, valid, forward});
    }
    Vec ref = sides.front().w, ref2 = d.cross(ref);
    for (auto &s : sides)
        s.angle = std::atan2(s.w.dot(ref2), s.w.dot(ref));
    std::vector<std::tuple<double, unsigned, unsigned>> options;
    constexpr double turn = 6.28318530717958647692;
    for (unsigned i = 0; i < sides.size(); ++i)
        for (unsigned j = 0; j < sides.size(); ++j)
            if (i != j && sides[i].valid && sides[j].valid &&
                sides[i].forward != sides[j].forward) {
                double g = std::fmod(sides[i].turn * (sides[j].angle - sides[i].angle), turn);
                if (g < 0)
                    g += turn;
                if (g > 1e-4 && g < turn - 1e-4)
                    options.emplace_back(g, i, j);
            }
    std::sort(options.begin(), options.end());
    std::vector<bool> used(sides.size(), false);
    std::vector<std::pair<unsigned, unsigned>> pairs;
    for (const auto &entry : options) {
        auto i = std::get<1>(entry), j = std::get<2>(entry);
        if (!used[i] && !used[j]) {
            used[i] = used[j] = true;
            pairs.emplace_back(triangles[i], triangles[j]);
        }
    }
    return pairs;
}
std::vector<std::vector<unsigned>> mesh_body_groups(const Mesh &m) {
    if (m.triangles.empty())
        return {};
    // First keep disconnected shells together only where their boxes overlap
    // in all three dimensions: an inner cavity belongs with its outer shell.
    BodyUnion vertices(unsigned(m.points.size()));
    for (auto t : m.triangles) {
        vertices.join(t[0], t[1]);
        vertices.join(t[1], t[2]);
    }
    std::map<unsigned, std::vector<unsigned>> pieces;
    for (unsigned i = 0; i < m.triangles.size(); ++i)
        pieces[vertices.find(m.triangles[i][0])].push_back(i);
    std::vector<std::vector<unsigned>> separate;
    std::vector<std::pair<Vec, Vec>> boxes;
    double size = 0;
    for (const auto &entry : pieces) {
        Vec lo = m.points[m.triangles[entry.second[0]][0]], hi = lo;
        for (auto i : entry.second)
            for (auto v : m.triangles[i]) {
                lo = lo.cwiseMin(m.points[v]);
                hi = hi.cwiseMax(m.points[v]);
            }
        size = std::max(size, (hi - lo).maxCoeff());
        separate.push_back(entry.second);
        boxes.emplace_back(lo, hi);
    }
    BodyUnion clusters(unsigned(separate.size()));
    for (unsigned i = 0; i < boxes.size(); ++i)
        for (unsigned j = i + 1; j < boxes.size(); ++j)
            if ((boxes[i].second.cwiseMin(boxes[j].second) -
                 boxes[i].first.cwiseMax(boxes[j].first))
                    .minCoeff() > 1e-4 * size)
                clusters.join(i, j);
    std::map<unsigned, std::vector<unsigned>> apart;
    for (unsigned i = 0; i < separate.size(); ++i) {
        auto &group = apart[clusters.find(i)];
        group.insert(group.end(), separate[i].begin(), separate[i].end());
    }
    std::vector<std::vector<unsigned>> result;
    for (auto &entry : apart) {
        auto &group = entry.second;
        std::sort(group.begin(), group.end());
        std::map<Edge, std::vector<unsigned>> edges;
        // Stable Python edge ordering matters to partner tie breaks.
        for (unsigned k = 0; k < 3; ++k)
            for (auto i : group) {
                auto t = m.triangles[i];
                edges[Edge(t[k], t[(k + 1) % 3])].push_back(i);
            }
        bool touching = std::any_of(edges.begin(), edges.end(),
                                    [](const auto &e) { return e.second.size() > 2; });
        if (!touching) {
            result.push_back(group);
            continue;
        }
        BodyUnion triangles(unsigned(m.triangles.size()));
        std::vector<std::vector<unsigned>> seams;
        for (const auto &e : edges) {
            if (e.second.size() == 2)
                triangles.join(e.second[0], e.second[1]);
            else if (e.second.size() > 2) {
                for (auto pair : material_partners(m, e.second, e.first))
                    triangles.join(pair.first, pair.second);
                seams.push_back(e.second);
            }
        }
        std::map<unsigned, std::vector<unsigned>> components;
        std::map<unsigned, double> volumes;
        for (auto i : group) {
            unsigned r = triangles.find(i);
            components[r].push_back(i);
            auto t = m.triangles[i];
            volumes[r] += m.points[t[0]].dot(m.points[t[1]].cross(m.points[t[2]])) / 6;
        }
        double biggest = 0;
        for (auto &v : volumes) {
            v.second = std::abs(v.second);
            biggest = std::max(biggest, v.second);
        }
        for (auto it = components.begin(); it != components.end();)
            if (it->second.size() < 4 || volumes[it->first] <= 1e-6 * biggest)
                it = components.erase(it);
            else
                ++it;
        std::map<unsigned, std::set<unsigned>> clash;
        for (const auto &seam : seams) {
            std::set<unsigned> roots;
            for (auto i : seam) {
                auto r = triangles.find(i);
                if (components.count(r))
                    roots.insert(r);
            }
            for (auto a : roots)
                for (auto b : roots)
                    if (a != b)
                        clash[a].insert(b);
        }
        std::vector<unsigned> order;
        for (const auto &c : components)
            order.push_back(c.first);
        std::stable_sort(order.begin(), order.end(), [&](unsigned a, unsigned b) {
            return components[a].size() > components[b].size();
        });
        std::map<unsigned, unsigned> colours;
        for (auto r : order) {
            std::set<unsigned> used;
            for (auto q : clash[r])
                if (colours.count(q))
                    used.insert(colours[q]);
            unsigned colour = 0;
            while (used.count(colour))
                ++colour;
            colours[r] = colour;
        }
        std::map<unsigned, std::vector<unsigned>> groups;
        for (const auto &c : components) {
            auto &g = groups[colours[c.first]];
            g.insert(g.end(), c.second.begin(), c.second.end());
        }
        for (auto &g : groups)
            result.push_back(std::move(g.second));
    }
    return result;
}
} // namespace

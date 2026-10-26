#pragma once
#include "occt_collections.hpp"
#include <Eigen/Core>
#include <Geom_Surface.hxx>
#include <TopoDS_Shape.hxx>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <stl2curves/converter.hpp>
#include <unordered_map>

namespace stl2curves::detail {
inline constexpr double volume_epsilon = 1e-7;
using Vec = Eigen::Vector3d;
using Tri = std::array<std::uint32_t, 3>;
using Clock = std::chrono::steady_clock;
struct Mesh {
    std::vector<Vec> points;
    std::vector<Tri> triangles;
    std::vector<unsigned> parts;
};
struct Object {
    std::string name;
    Mesh mesh;
};
struct Edge {
    std::uint32_t a, b;
    Edge(std::uint32_t x, std::uint32_t y) : a(std::min(x, y)), b(std::max(x, y)) {}
    bool operator==(const Edge &o) const { return a == o.a && b == o.b; }
    bool operator<(const Edge &o) const { return a < o.a || (a == o.a && b < o.b); }
};
struct EdgeHash {
    std::size_t operator()(const Edge &e) const { return (std::size_t(e.a) * 0x9e3779b1u) ^ e.b; }
};
using EdgeFaces = std::unordered_map<Edge, std::vector<unsigned>, EdgeHash>;
struct Topology {
    EdgeFaces edges;
    std::vector<std::vector<unsigned>> neighbors;
    std::vector<Vec> normals, centers;
    std::vector<double> areas;
};
struct Context {
    Options options;
    Diagnostics diagnostics;
    Clock::time_point deadline = Clock::time_point::max();
    bool expired() {
        if (options.time_limit == 0)
            return false;
        if (Clock::now() < deadline)
            return false;
        diagnostics.time_limited = true;
        return true;
    }
};
struct Patch {
    std::string kind;
    unsigned part = 0, body = 0;
    std::vector<unsigned> triangles;
    Handle(Geom_Surface) surface;
    Vec origin = Vec::Zero(), axis = Vec::UnitZ();
    double radius = 0, secondary = 0, deviation = 0;
};
Mesh read_stl(const Bytes &bytes);
std::vector<Object> read_3mf(const Bytes &bytes);
void weld(Mesh &mesh);
Topology topology(const Mesh &mesh);
std::vector<std::pair<Vec, double>> symmetry_planes(const Mesh &mesh);
void repair_mesh(Mesh &mesh, Context &context);
void simplify_mesh(Mesh &mesh, double tolerance, Context &context);
std::vector<Mesh> split_bodies(const Mesh &mesh);
double signed_volume(const Mesh &mesh);
std::vector<Patch> detect(const Mesh &mesh, Context &context);
std::vector<Patch> advanced_patches(const Mesh &mesh, const Topology &top,
                                    const std::vector<int> &owner,
                                    const std::vector<Patch> &recognized, Context &context);
Patch thread_patch(const Mesh &mesh, const Topology &top, const std::vector<unsigned> &triangles,
                   const Vec &axis, double tolerance);
void design_size(Mesh &mesh, std::vector<Patch> &patches, Context &context);
TopoDS_Shape reconstruct(const Mesh &mesh, const std::vector<Patch> &patches, Context &context);
TopoDS_Shape faceted(const Mesh &mesh, double tolerance);
TopoDS_Shape recover(const Mesh &mesh, const std::vector<Patch> &patches, Context &context);
TopoDS_Shape unite(const std::vector<TopoDS_Shape> &shapes, Context &context);
bool inspect(const TopoDS_Shape &shape, Diagnostics &info);
bool volume_agrees(const TopoDS_Shape &shape, const Mesh &mesh, const std::vector<Patch> &patches);
Bytes write_step(const TopoDS_Shape &shape, double tolerance, bool surface_curves = true);
bool verify_step(const Bytes &bytes, const TopoDS_Shape &shape, double tolerance,
                 Diagnostics *exported_info = nullptr, TopoDS_Shape *exported_shape = nullptr);
// Stable indexed results: workers never mutate the mesh or merge topology.
unsigned parallel_for(std::size_t count, const std::function<void(std::size_t)> &task);
unsigned worker_count();
void trace(const char *stage);
std::vector<unsigned> vertices(const Mesh &mesh, const std::vector<unsigned> &triangles);
std::vector<std::vector<unsigned>> components(const std::vector<std::vector<unsigned>> &adjacency);
struct Fit {
    Eigen::VectorXd x;
    double cost = 0;
    bool converged = false;
};
Fit least_squares(const std::function<Eigen::VectorXd(const Eigen::VectorXd &)> &residual,
                  Eigen::VectorXd initial, unsigned iterations = 100, double loss_scale = 0.01);
std::string safe_name(std::string name);
} // namespace stl2curves::detail

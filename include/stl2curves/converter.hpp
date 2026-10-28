// Derived from stl2curves, Copyright (c) 2026 Duvdivan. MIT license.
#pragma once
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace stl2curves {
using Bytes = std::vector<std::uint8_t>;
struct InputFile {
    std::string name;
    Bytes data;
};
struct Options {
    double tolerance = 0.01;
    double accuracy = 0.1;
    double time_limit = 600;
    // -1: automatic (0.005 mm above 1,000,000 triangles); 0: disabled.
    double simplify = -1;
    bool fuse = true;
    bool curves = true;
    bool blends = true;
    bool repair = true;
    bool true_size = false;
    bool details = false;
    std::string merge;
};
struct FeatureInfo {
    std::string kind;
    std::size_t triangles = 0;
    double deviation = 0;
    double radius = 0;
    double pitch = 0;
};
struct Diagnostics {
    std::size_t input_triangles = 0, triangles = 0, faces = 0, solids = 0;
    std::size_t open_edges = 0, nonmanifold_edges = 0, parallel_tasks = 0;
    double mesh_volume = 0, volume = 0, scale = 1;
    std::string design_unit, design_label;
    bool valid = false, fallback = false, time_limited = false;
    std::vector<std::string> repairs, warnings, fallbacks;
    std::vector<FeatureInfo> features;
    std::map<std::string, std::size_t> surface_counts;
};
struct OutputFile {
    std::string name;
    Bytes data;
    Diagnostics diagnostics;
};
struct GraphId {
    int kind = 0;
    std::uint32_t index = 0;
    bool operator==(const GraphId &other) const {
        return kind == other.kind && index == other.index;
    }
};
struct UsagePathStep {
    GraphId node;
    std::optional<GraphId> ref;
    int step_index = -1;
    bool operator==(const UsagePathStep &other) const {
        return node == other.node && ref == other.ref && step_index == other.step_index;
    }
};
struct ShapeRef {
    std::uint32_t model_id = 0;
    std::optional<GraphId> brepgraph_id;
    std::optional<std::vector<UsagePathStep>> usage_path;
};
struct GltfExportOptions {
    double linear_deflection = 0.1, angular_deflection = 0.5;
    bool faces = true, edges = true, vertices = true;
};
struct GltfFile {
    std::string name;
    Bytes data;
    std::vector<std::string> warnings;
};
struct ClosestPoints {
    std::vector<double> p1, p2;
    double distance = 0;
};
struct ShapeProperty {
    std::string name;
    std::vector<double> value;
};
namespace detail {
struct ModelAccess;
}
class Model {
  public:
    ~Model();
    Model(const Model &) = delete;
    Model &operator=(const Model &) = delete;
    std::uint32_t id() const;
    const std::string &name() const;
    const Diagnostics &diagnostics() const;
    ShapeRef root() const;
    OutputFile export_step() const;
    GltfFile export_gltf(const GltfExportOptions &options = {}) const;
    std::vector<ShapeProperty> get_shape_properties(const ShapeRef &selection) const;
    ClosestPoints get_closest_points(const ShapeRef &selection, const Model &other,
                                     const ShapeRef &other_selection) const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    explicit Model(std::unique_ptr<Impl> impl);
    friend struct detail::ModelAccess;
};
struct Failure {
    std::string input, object, code, message;
    std::optional<Diagnostics> diagnostics = std::nullopt;
};
struct Result {
    std::vector<std::shared_ptr<Model>> models;
    std::vector<Failure> errors;
};
struct BuildInfo {
    std::string version, occt_version;
    bool threaded;
    unsigned workers;
};
BuildInfo build_info();
// Partial success is represented by retained models and per-object errors together.
// Configuration errors throw std::invalid_argument before processing inputs.
Result convert(const std::vector<InputFile> &inputs, const Options &options = {});
} // namespace stl2curves

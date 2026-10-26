// Derived from stl2curves, Copyright (c) 2026 Duvdivan. MIT license.
#pragma once
#include <cstddef>
#include <cstdint>
#include <map>
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
struct Failure {
    std::string input, object, code, message;
    std::optional<Diagnostics> diagnostics = std::nullopt;
};
struct Result {
    std::vector<OutputFile> files;
    std::vector<Failure> errors;
};
struct BuildInfo {
    std::string version, occt_version;
    bool threaded;
    unsigned workers;
};
BuildInfo build_info();
// Partial success is represented by files and per-object errors together.
// Configuration errors throw std::invalid_argument before processing inputs.
Result convert(const std::vector<InputFile> &inputs, const Options &options = {});
} // namespace stl2curves

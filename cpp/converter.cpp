#include "model_internal.hpp"
#include "occt_error.hpp"
#include <ShapeFix_Shape.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <limits>
#include <mutex>
#include <set>

namespace stl2curves {
namespace {
std::map<unsigned, detail::Mesh> split_parts(detail::Mesh &mesh) {
    std::map<unsigned, detail::Mesh> parts;
    if (mesh.parts.empty()) {
        parts.emplace(0, std::move(mesh));
        return parts;
    }
    for (std::size_t i = 0; i < mesh.triangles.size(); ++i)
        parts[mesh.parts.at(i)].triangles.push_back(mesh.triangles[i]);

    // Reuse one assembly-sized index map. Each part owns only its referenced
    // vertices, in their original order, rather than a copy of every part's
    // points. Reset only touched entries between parts.
    constexpr auto unset = std::numeric_limits<std::uint32_t>::max();
    std::vector<std::uint32_t> remap(mesh.points.size(), unset), used;
    for (auto &entry : parts) {
        auto &part = entry.second;
        used.clear();
        for (const auto &triangle : part.triangles)
            for (auto vertex : triangle)
                if (remap.at(vertex) == unset) {
                    remap[vertex] = 0;
                    used.push_back(vertex);
                }
        std::sort(used.begin(), used.end());
        part.points.reserve(used.size());
        for (auto vertex : used) {
            remap[vertex] = std::uint32_t(part.points.size());
            part.points.push_back(mesh.points[vertex]);
        }
        for (auto &triangle : part.triangles)
            for (auto &vertex : triangle)
                vertex = remap[vertex];
        for (auto vertex : used)
            remap[vertex] = unset;
    }
    mesh = detail::Mesh{};
    return parts;
}

} // namespace
BuildInfo build_info() {
    return {"0.2.0", OCC_VERSION_COMPLETE, STL2CURVES_THREADS != 0, detail::worker_count()};
}
Result convert(const std::vector<InputFile> &inputs, const Options &options) {
    using namespace detail;
    if (!std::isfinite(options.tolerance) || options.tolerance <= 0 ||
        !std::isfinite(options.accuracy) || options.accuracy <= 0 ||
        !std::isfinite(options.time_limit) || options.time_limit < 0 ||
        !std::isfinite(options.simplify) || (options.simplify < 0 && options.simplify != -1))
        throw std::invalid_argument(
            "Invalid tolerance, accuracy, time limit, or simplification option");
    // Serialize conversion calls while internal fitting tasks remain bounded.
    static std::mutex conversion_mutex;
    std::lock_guard<std::mutex> guard(conversion_mutex);
    Result result;
    std::vector<TopoDS_Shape> successful;
    std::set<std::string> used_names;
    auto unique = [&](std::string name) {
        std::string base = name;
        unsigned n = 1;
        while (!used_names.insert(name).second)
            name = base + " (" + std::to_string(++n) + ")";
        return name;
    };
    for (const auto &input : inputs) {
        std::vector<Object> objects;
        std::string stem = safe_name(std::filesystem::path(input.name).stem().string());
        try {
            auto extension = std::filesystem::path(input.name).extension().string();
            std::transform(extension.begin(), extension.end(), extension.begin(),
                           [](unsigned char c) { return char(std::tolower(c)); });
            if (extension == ".stl")
                objects.push_back({stem, read_stl(input.data)});
            else if (extension == ".3mf") {
                objects = read_3mf(input.data);
                for (auto &o : objects)
                    o.name = stem + " - " + o.name;
            } else
                throw std::runtime_error("Expected an STL or 3MF input filename");
        } catch (const std::exception &e) {
            result.errors.push_back({input.name, "", "input", e.what()});
            continue;
        }
        for (auto &object : objects) {
            Context ctx;
            ctx.options = options;
            ctx.diagnostics.input_triangles = object.mesh.triangles.size();
            ctx.deadline = Clock::time_point::max();
            try {
                if (object.mesh.triangles.empty())
                    throw std::runtime_error("Mesh contains no nondegenerate triangles");
                // Per-part reconstruction avoids treating overlapping 3MF parts as
                // a self-intersecting shell; each part may include cavity shells.
                auto part_meshes = split_parts(object.mesh);
                std::map<unsigned, std::vector<Patch>> analyses;
                std::vector<Patch> combined;
                Mesh measured;
                std::map<unsigned, unsigned> triangle_offsets;
                std::map<unsigned, std::pair<std::size_t, std::size_t>> vertex_ranges;
                for (auto &kv : part_meshes) {
                    auto &mesh = kv.second;
                    if (options.repair)
                        repair_mesh(mesh, ctx);
                    else {
                        auto top = topology(mesh);
                        for (const auto &edge : top.edges) {
                            if (edge.second.size() == 1)
                                ++ctx.diagnostics.open_edges;
                            if (edge.second.size() > 2)
                                ++ctx.diagnostics.nonmanifold_edges;
                        }
                    }
                    double simplify = options.simplify < 0
                                          ? (mesh.triangles.size() > 1000000 ? 0.005 : 0)
                                          : options.simplify;
                    if (simplify)
                        simplify_mesh(mesh, simplify, ctx);
                    ctx.diagnostics.triangles += mesh.triangles.size();
                    auto patches = detect(mesh, ctx);
                    unsigned offset = unsigned(measured.triangles.size());
                    triangle_offsets[kv.first] = offset;
                    const auto vertex_offset = measured.points.size();
                    vertex_ranges[kv.first] = {vertex_offset, mesh.points.size()};
                    if (vertex_offset + mesh.points.size() >
                        std::numeric_limits<std::uint32_t>::max())
                        throw std::runtime_error("Assembly has too many vertices");
                    if (part_meshes.size() == 1)
                        measured.points = std::move(mesh.points);
                    else
                        measured.points.insert(measured.points.end(), mesh.points.begin(),
                                               mesh.points.end());
                    for (auto triangle : mesh.triangles) {
                        for (auto &vertex : triangle)
                            vertex += std::uint32_t(vertex_offset);
                        measured.triangles.push_back(triangle);
                    }
                    for (auto p : patches) {
                        p.part = kv.first;
                        for (auto &t : p.triangles)
                            t += offset;
                        combined.push_back(std::move(p));
                    }
                    analyses.emplace(kv.first, std::move(patches));
                }
                design_size(measured, combined, ctx);
                std::size_t patch_offset = 0;
                for (auto &kv : part_meshes) {
                    const auto range = vertex_ranges.at(kv.first);
                    if (part_meshes.size() == 1)
                        kv.second.points = std::move(measured.points);
                    else
                        kv.second.points.assign(measured.points.begin() + range.first,
                                                measured.points.begin() + range.first +
                                                    range.second);
                    auto &patches = analyses.at(kv.first);
                    for (auto &patch : patches) {
                        patch = std::move(combined[patch_offset++]);
                        for (auto &t : patch.triangles)
                            t -= triangle_offsets[kv.first];
                    }
                }
                // Release assembly measurement storage before allocating BRep
                // topology. Each part and its patches now own the final data.
                measured = Mesh{};
                combined = std::vector<Patch>{};
                std::vector<TopoDS_Shape> shapes;
                for (auto &kv : part_meshes) {
                    auto &mesh = kv.second;
                    auto &patches = analyses.at(kv.first);
                    ctx.diagnostics.mesh_volume += std::abs(signed_volume(mesh));
                    if (options.details)
                        for (const auto &p : patches)
                            ctx.diagnostics.features.push_back(
                                {p.kind, p.triangles.size(), p.deviation, p.radius,
                                 p.kind.rfind("thread_", 0) == 0 ? std::abs(p.secondary) : 0});
                    TopoDS_Shape shape;
                    try {
                        shape = reconstruct(mesh, patches, ctx);
                        Diagnostics info;
                        if (!inspect(shape, info))
                            throw std::runtime_error("Reconstructed solid is invalid");
                        if (!volume_agrees(shape, mesh, patches))
                            throw std::runtime_error("Reconstructed volume is inconsistent with "
                                                     "the fitted surface corrections");
                    } catch (const Standard_Failure &e) {
                        ctx.diagnostics.fallback = true;
                        ctx.diagnostics.fallbacks.push_back(
                            failure_message(e, "OpenCascade reconstruction failure"));
                        shape = recover(mesh, patches, ctx);
                    } catch (const std::exception &e) {
                        ctx.diagnostics.fallback = true;
                        ctx.diagnostics.fallbacks.push_back(e.what());
                        shape = recover(mesh, patches, ctx);
                    }
                    shapes.push_back(shape);
                    ctx.expired();
                }
                auto shape = unite(shapes, ctx);
                if (!inspect(shape, ctx.diagnostics))
                    throw std::runtime_error("Final solid failed validation");
                auto model = ModelAccess::create(unique(safe_name(object.name)), shape,
                                                 ctx.diagnostics, options.tolerance);
                result.models.push_back(std::move(model));
                if (!options.merge.empty())
                    successful.push_back(shape);
            } catch (const Standard_Failure &e) {
                result.errors.push_back({input.name, object.name, "geometry", failure_message(e),
                                         std::move(ctx.diagnostics)});
            } catch (const std::exception &e) {
                result.errors.push_back(
                    {input.name, object.name, "conversion", e.what(), std::move(ctx.diagnostics)});
            }
        }
    }
    if (!options.merge.empty() && !successful.empty()) {
        Context ctx;
        ctx.options = options;
        ctx.options.fuse = false;
        try {
            auto shape = unite(successful, ctx);
            inspect(shape, ctx.diagnostics);
            result.models.push_back(ModelAccess::create(
                unique(safe_name(std::filesystem::path(options.merge).stem().string())), shape,
                ctx.diagnostics, options.tolerance));
        } catch (const Standard_Failure &e) {
            result.errors.push_back({"", "merge", "geometry", failure_message(e, "Merge failure"),
                                     std::move(ctx.diagnostics)});
        } catch (const std::exception &e) {
            result.errors.push_back(
                {"", "merge", "geometry", e.what(), std::move(ctx.diagnostics)});
        }
    }
    return result;
}
} // namespace stl2curves

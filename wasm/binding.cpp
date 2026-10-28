#include <cmath>
#include <emscripten/bind.h>
#include <limits>
#include <map>
#include <stdexcept>
#include <stl2curves/converter.hpp>

namespace {
using emscripten::val;
using namespace stl2curves;
std::map<std::uint32_t, std::shared_ptr<Model>> models;
val strings(const std::vector<std::string> &values) {
    auto a = val::array();
    for (const auto &v : values)
        a.call<void>("push", v);
    return a;
}
val diagnostics(const Diagnostics &d) {
    auto o = val::object();
    o.set("inputTriangles", double(d.input_triangles));
    o.set("triangles", double(d.triangles));
    o.set("faces", double(d.faces));
    o.set("solids", double(d.solids));
    o.set("openEdges", double(d.open_edges));
    o.set("nonmanifoldEdges", double(d.nonmanifold_edges));
    o.set("parallelTasks", double(d.parallel_tasks));
    o.set("meshVolume", d.mesh_volume);
    o.set("volume", d.volume);
    o.set("scale", d.scale);
    o.set("valid", d.valid);
    o.set("fallback", d.fallback);
    o.set("timeLimited", d.time_limited);
    o.set("designUnit", d.design_unit);
    o.set("designLabel", d.design_label);
    o.set("repairs", strings(d.repairs));
    o.set("warnings", strings(d.warnings));
    o.set("fallbacks", strings(d.fallbacks));
    auto counts = val::object();
    for (const auto &kv : d.surface_counts)
        counts.set(kv.first, double(kv.second));
    o.set("surfaceCounts", counts);
    auto features = val::array();
    for (const auto &f : d.features) {
        auto item = val::object();
        item.set("kind", f.kind);
        item.set("triangles", double(f.triangles));
        item.set("deviation", f.deviation);
        item.set("radius", f.radius);
        item.set("pitch", f.pitch);
        features.call<void>("push", item);
    }
    o.set("features", features);
    return o;
}
[[noreturn]] void type_error(const std::string &message) {
    throw std::invalid_argument(message);
}
bool has(const val &object, const char *key) {
    return !object[key].isUndefined();
}
double number(const val &o, const char *key, double fallback) {
    if (!has(o, key))
        return fallback;
    if (o[key].typeOf().as<std::string>() != "number")
        type_error(std::string(key) + " must be a number");
    double v = o[key].as<double>();
    if (!std::isfinite(v))
        type_error(std::string(key) + " must be finite");
    return v;
}
bool boolean(const val &o, const char *key, bool fallback) {
    if (!has(o, key))
        return fallback;
    if (o[key].typeOf().as<std::string>() != "boolean")
        type_error(std::string(key) + " must be a boolean");
    return o[key].as<bool>();
}
std::uint32_t integer(val value, const char *label, bool positive = false) {
    if (value.typeOf().as<std::string>() != "number")
        type_error(std::string(label) + " must be an integer");
    double n = value.as<double>();
    if (!std::isfinite(n) || n != std::floor(n) || n < (positive ? 1 : 0) || n > 4294967295.)
        type_error(std::string(label) + " is out of range");
    return std::uint32_t(n);
}
val graph_id(const GraphId &id) {
    auto result = val::object();
    result.set("kind", id.kind);
    result.set("index", double(id.index));
    return result;
}
val selection(const ShapeRef &ref) {
    auto result = val::object();
    result.set("modelId", double(ref.model_id));
    if (ref.brepgraph_id)
        result.set("brepgraph_id", graph_id(*ref.brepgraph_id));
    if (ref.usage_path) {
        auto path = val::array();
        for (const auto &step : *ref.usage_path) {
            auto item = val::object();
            item.set("node", graph_id(step.node));
            item.set("stepIndex", step.step_index);
            if (step.ref)
                item.set("ref", graph_id(*step.ref));
            path.call<void>("push", item);
        }
        result.set("usagePath", path);
    }
    return result;
}
GraphId parse_id(val value) {
    if (value.isNull() || value.isUndefined() || value.typeOf().as<std::string>() != "object")
        type_error("Graph ID must be an object");
    const auto kind = integer(value["kind"], "kind");
    if (kind > 11 || kind == 9)
        type_error("Unknown graph kind");
    return {int(kind), integer(value["index"], "index")};
}
ShapeRef parse_selection(val value) {
    if (value.isNull() || value.isUndefined() || value.typeOf().as<std::string>() != "object")
        type_error("Selection must be an object");
    ShapeRef ref;
    ref.model_id = integer(value["modelId"], "modelId", true);
    if (has(value, "brepgraph_id"))
        ref.brepgraph_id = parse_id(value["brepgraph_id"]);
    if (has(value, "usagePath")) {
        auto path = value["usagePath"];
        if (!val::global("Array").call<bool>("isArray", path))
            type_error("usagePath must be an array");
        ref.usage_path.emplace();
        for (unsigned i = 0; i < path["length"].as<unsigned>(); ++i) {
            auto step = path[i];
            if (step.isNull() || step.isUndefined())
                type_error("Invalid usagePath step");
            double index = number(step, "stepIndex", -2);
            if (index != std::floor(index) || index < -1 || index > 2147483647.)
                type_error("Invalid stepIndex");
            UsagePathStep item{parse_id(step["node"]), std::nullopt, int(index)};
            if (has(step, "ref"))
                item.ref = parse_id(step["ref"]);
            ref.usage_path->push_back(item);
        }
    }
    return ref;
}
std::shared_ptr<Model> model(std::uint32_t id) {
    auto found = models.find(id);
    if (found == models.end())
        throw std::invalid_argument("Unknown or released modelId");
    return found->second;
}
val bytes(const Bytes &data) {
    return val::global("Uint8Array")
        .new_(val(emscripten::typed_memory_view(data.size(), data.data())));
}
template <class F> val invoke(F &&operation) {
    // Return errors to post.js instead of throwing JS through active C++ frames.
    // This lets vectors, shared_ptrs and caught exceptions unwind normally.
    const auto error_result = [](const char *name, const char *message) {
        auto result = val::object(), info = val::object();
        info.set("name", std::string(name));
        info.set("message", std::string(message));
        result.set("__stl2curvesError", info);
        return result;
    };
    try {
        return operation();
    } catch (const std::invalid_argument &error) {
        return error_result("TypeError", error.what());
    } catch (const std::exception &error) {
        return error_result("Error", error.what());
    }
}
val convert_files(val input, val config) {
    if (!val::global("Array").call<bool>("isArray", input))
        type_error("inputs must be an array of {name, data}");
    Options options;
    if (!config.isUndefined() && !config.isNull()) {
        if (config.typeOf().as<std::string>() != "object")
            type_error("options must be an object");
        options.tolerance = number(config, "tolerance", options.tolerance);
        options.accuracy = number(config, "accuracy", options.accuracy);
        options.time_limit = number(config, "timeLimit", options.time_limit);
        options.simplify = number(config, "simplify", options.simplify);
        options.fuse = boolean(config, "fuse", options.fuse);
        options.curves = boolean(config, "curves", options.curves);
        options.blends = boolean(config, "blends", options.blends);
        options.repair = boolean(config, "repair", options.repair);
        options.true_size = boolean(config, "trueSize", options.true_size);
        options.details = boolean(config, "details", options.details);
        if (has(config, "merge")) {
            if (config["merge"].typeOf().as<std::string>() != "string")
                type_error("merge must be a string");
            options.merge = config["merge"].as<std::string>();
        }
    }
    if (options.tolerance <= 0 || options.accuracy <= 0 || options.time_limit < 0 ||
        (options.simplify < 0 && options.simplify != -1))
        type_error("Invalid conversion options");
    std::vector<InputFile> files;
    auto count = input["length"].as<unsigned>();
    files.reserve(count);
    for (unsigned i = 0; i < count; ++i) {
        auto file = input[i];
        if (file.isNull() || file.isUndefined() ||
            file["name"].typeOf().as<std::string>() != "string")
            type_error("Each input needs a name");
        auto data = file["data"];
        if (data.isNull() || data.isUndefined() ||
            !val::global("ArrayBuffer").call<bool>("isView", data) ||
            data["BYTES_PER_ELEMENT"].isUndefined() ||
            data["BYTES_PER_ELEMENT"].as<unsigned>() != 1)
            type_error("Input data must be a Uint8Array");
        auto length = data["byteLength"].as<std::size_t>();
        InputFile f{file["name"].as<std::string>(), Bytes(length)};
        if (length)
            val(emscripten::typed_memory_view(length, f.data.data())).call<void>("set", data);
        files.push_back(std::move(f));
    }
    auto result = convert(files, options);
    auto o = val::object();
    auto outputs = val::array(), errors = val::array();
    for (const auto &entry : result.models) {
        auto f = val::object();
        f.set("modelId", double(entry->id()));
        f.set("name", entry->name());
        f.set("root", selection(entry->root()));
        f.set("diagnostics", diagnostics(entry->diagnostics()));
        outputs.call<void>("push", f);
    }
    for (const auto &error : result.errors) {
        auto e = val::object();
        e.set("input", error.input);
        e.set("object", error.object);
        e.set("code", error.code);
        e.set("message", error.message);
        if (error.diagnostics)
            e.set("diagnostics", diagnostics(*error.diagnostics));
        errors.call<void>("push", e);
    }
    o.set("models", outputs);
    o.set("errors", errors);
    try {
        for (const auto &entry : result.models)
            models.emplace(entry->id(), entry);
    } catch (...) {
        for (const auto &entry : result.models)
            models.erase(entry->id());
        throw;
    }
    return o;
}
val export_step(val id) {
    return invoke([&]() {
        auto file = model(integer(id, "modelId", true))->export_step();
        auto result = val::object();
        result.set("name", file.name);
        result.set("data", bytes(file.data));
        result.set("diagnostics", diagnostics(file.diagnostics));
        return result;
    });
}
val export_gltf(val id, val config) {
    return invoke([&]() {
        GltfExportOptions options;
        if (!config.isUndefined() && !config.isNull()) {
            if (config.typeOf().as<std::string>() != "object")
                type_error("GLB options must be an object");
            options.linear_deflection =
                number(config, "linearDeflection", options.linear_deflection);
            options.angular_deflection =
                number(config, "angularDeflection", options.angular_deflection);
            options.faces = boolean(config, "faces", options.faces);
            options.edges = boolean(config, "edges", options.edges);
            options.vertices = boolean(config, "vertices", options.vertices);
        }
        auto file = model(integer(id, "modelId", true))->export_gltf(options);
        auto result = val::object();
        result.set("name", file.name);
        result.set("data", bytes(file.data));
        result.set("warnings", strings(file.warnings));
        return result;
    });
}
val closest(val a, val b) {
    return invoke([&]() {
        auto first = parse_selection(a), second = parse_selection(b);
        auto value =
            model(first.model_id)->get_closest_points(first, *model(second.model_id), second);
        auto result = val::object(), p1 = val::array(), p2 = val::array();
        for (double n : value.p1)
            p1.call<void>("push", n);
        for (double n : value.p2)
            p2.call<void>("push", n);
        result.set("p1", p1);
        result.set("p2", p2);
        result.set("distance", value.distance);
        return result;
    });
}
val properties(val input) {
    return invoke([&]() {
        auto ref = parse_selection(input);
        auto result = val::object(), items = val::array();
        for (const auto &property : model(ref.model_id)->get_shape_properties(ref)) {
            auto item = val::object(), value = val::array();
            item.set("name", property.name);
            for (double n : property.value)
                value.call<void>("push", n);
            item.set("value", value);
            items.call<void>("push", item);
        }
        result.set("properties", items);
        return result;
    });
}
val info() {
    auto b = build_info();
    auto v = val::object();
    v.set("version", b.version);
    v.set("occtVersion", b.occt_version);
    v.set("threaded", b.threaded);
    v.set("workers", b.workers);
    return v;
}
} // namespace
EMSCRIPTEN_BINDINGS(stl2curves_bindings) {
    emscripten::function("convert", emscripten::optional_override([](val input, val config) {
                             return invoke([&]() { return convert_files(input, config); });
                         }));
    emscripten::function("exportStep", &export_step);
    emscripten::function("exportGltf", &export_gltf);
    emscripten::function("getClosestPoints", &closest);
    emscripten::function("getShapeProperties", &properties);
    emscripten::function(
        "releaseModel", emscripten::optional_override([](val id) {
            return invoke([&]() { return val(models.erase(integer(id, "modelId", true)) != 0); });
        }));
    emscripten::function("_clearModels", emscripten::optional_override([]() { models.clear(); }));
    emscripten::function("buildInfo", &info);
}

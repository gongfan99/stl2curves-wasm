#include <cmath>
#include <emscripten/bind.h>
#include <limits>
#include <stl2curves/converter.hpp>

namespace {
using emscripten::val;
using namespace stl2curves;
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
    val::global("TypeError").new_(message).throw_();
    __builtin_unreachable();
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
    for (const auto &file : result.files) {
        auto f = val::object();
        f.set("name", file.name);
        // Constructor copies the view before the native vector is destroyed.
        f.set("data",
              val::global("Uint8Array")
                  .new_(val(emscripten::typed_memory_view(file.data.size(), file.data.data()))));
        f.set("diagnostics", diagnostics(file.diagnostics));
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
    o.set("files", outputs);
    o.set("errors", errors);
    return o;
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
    emscripten::function("convert", &convert_files);
    emscripten::function("buildInfo", &info);
}

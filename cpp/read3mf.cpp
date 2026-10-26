#include "internal.hpp"
#include <Eigen/Geometry>
#include <cctype>
#include <filesystem>
#include <memory>
#include <pugixml.hpp>
#include <set>
#include <sstream>
#include <zip.h>

namespace stl2curves::detail {
namespace {
std::string local(const char *name) {
    std::string s(name);
    auto p = s.find(':');
    return p == s.npos ? s : s.substr(p + 1);
}
pugi::xml_node child(pugi::xml_node node, const char *name) {
    for (auto c : node.children())
        if (local(c.name()) == name)
            return c;
    return {};
}
pugi::xml_attribute attr(pugi::xml_node node, const char *name) {
    for (auto a : node.attributes())
        if (local(a.name()) == name)
            return a;
    return {};
}
std::string path(std::string value) {
    std::replace(value.begin(), value.end(), '\\', '/');
    while (!value.empty() && value.front() == '/')
        value.erase(value.begin());
    auto normalized = std::filesystem::path(value).lexically_normal().generic_string();
    if (normalized == ".." || normalized.rfind("../", 0) == 0)
        throw std::runtime_error("3MF reference escapes archive");
    return normalized;
}
Eigen::Matrix4d matrix(const char *text, double unit) {
    Eigen::Matrix4d m = Eigen::Matrix4d::Identity();
    if (!text || !*text)
        return m;
    std::istringstream in(text);
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 3; ++col) {
            if (!(in >> m(col, row)) || !std::isfinite(m(col, row)))
                throw std::runtime_error("Invalid 3MF transform");
        }
    std::string extra;
    if (in >> extra)
        throw std::runtime_error("Invalid 3MF transform length");
    // Coordinates have already been converted to millimeters. Translation
    // entries are locations in the owning model's unit; linear entries are
    // dimensionless and must not receive a second unit conversion.
    m.block<3, 1>(0, 3) *= unit;
    return m;
}
double coordinate(pugi::xml_node node, const char *key) {
    auto a = node.attribute(key);
    if (!a)
        throw std::runtime_error("Missing 3MF coordinate");
    std::string text = a.value();
    std::size_t end = 0;
    double value = std::stod(text, &end);
    if (end != text.size() || !std::isfinite(value) || std::abs(value) > 1e10)
        throw std::runtime_error("Invalid 3MF coordinate");
    return value;
}
} // namespace
std::vector<Object> read_3mf(const Bytes &bytes) {
    zip_error_t error;
    zip_error_init(&error);
    auto *source = zip_source_buffer_create(bytes.data(), bytes.size(), 0, &error);
    if (!source) {
        std::string message = zip_error_strerror(&error);
        zip_error_fini(&error);
        throw std::runtime_error(message);
    }
    zip_t *raw = zip_open_from_source(source, ZIP_RDONLY, &error);
    if (!raw) {
        zip_source_free(source);
        std::string message = zip_error_strerror(&error);
        zip_error_fini(&error);
        throw std::runtime_error(message);
    }
    zip_error_fini(&error);
    std::unique_ptr<zip_t, decltype(&zip_discard)> archive(raw, zip_discard);
    auto exists = [&](const std::string &name) {
        return zip_name_locate(raw, name.c_str(), 0) >= 0;
    };
    auto read = [&](const std::string &name) {
        zip_stat_t st;
        zip_stat_init(&st);
        if (zip_stat(raw, name.c_str(), 0, &st) < 0)
            throw std::runtime_error("Missing 3MF archive entry: " + name);
        if (st.size > 512ull * 1024 * 1024)
            throw std::runtime_error("3MF XML entry exceeds 512 MiB");
        std::unique_ptr<zip_file_t, decltype(&zip_fclose)> file(zip_fopen(raw, name.c_str(), 0),
                                                                zip_fclose);
        if (!file)
            throw std::runtime_error("Cannot open 3MF XML: " + name);
        std::string text(std::size_t(st.size), '\0');
        std::size_t offset = 0;
        while (offset < text.size()) {
            auto n = zip_fread(file.get(), text.data() + offset, text.size() - offset);
            if (n <= 0)
                throw std::runtime_error("Truncated 3MF XML");
            offset += std::size_t(n);
        }
        return text;
    };
    auto xml = [&](const std::string &name) {
        auto document = std::make_unique<pugi::xml_document>();
        auto text = read(name);
        auto result = document->load_buffer(text.data(), text.size());
        if (!result)
            throw std::runtime_error("Invalid XML in " + name + ": " + result.description());
        return document;
    };
    std::string root = "3D/3dmodel.model";
    if (exists("_rels/.rels")) {
        auto rels = xml("_rels/.rels");
        for (auto rel : rels->document_element().children()) {
            std::string type = rel.attribute("Type").value();
            if (type.size() >= 8 && type.substr(type.size() - 8) == "/3dmodel")
                root = path(rel.attribute("Target").value());
        }
    }
    struct Model {
        std::unique_ptr<pugi::xml_document> document;
        double scale;
        std::map<std::string, pugi::xml_node> objects;
    };
    std::map<std::string, Model> models;
    const std::map<std::string, double> units{{"micron", .001},   {"millimeter", 1},
                                              {"centimeter", 10}, {"inch", 25.4},
                                              {"foot", 304.8},    {"meter", 1000}};
    auto model = [&](const std::string &file) -> Model & {
        auto it = models.find(file);
        if (it != models.end())
            return it->second;
        Model m;
        m.document = xml(file);
        auto node = m.document->document_element();
        std::string unit = node.attribute("unit").as_string("millimeter");
        auto u = units.find(unit);
        if (u == units.end())
            throw std::runtime_error("Unsupported 3MF unit: " + unit);
        m.scale = u->second;
        for (auto object : child(node, "resources").children())
            if (local(object.name()) == "object")
                m.objects.emplace(object.attribute("id").value(), object);
        return models.emplace(file, std::move(m)).first->second;
    };
    std::map<std::string, std::string> titles;
    std::map<std::pair<std::string, std::string>, std::string> kinds;
    if (exists("Metadata/model_settings.config")) {
        auto settings = xml("Metadata/model_settings.config");
        for (auto object : settings->document_element().children("object")) {
            std::string oid = object.attribute("id").value();
            for (auto md : object.children("metadata"))
                if (std::string(md.attribute("key").value()) == "name")
                    titles[oid] = md.attribute("value").value();
            for (auto part : object.children("part"))
                kinds[{oid, part.attribute("id").value()}] =
                    part.attribute("subtype").as_string("normal_part");
        }
    }
    std::set<std::pair<std::string, std::string>> active;
    unsigned part_label = 0;
    std::function<void(const std::string &, const std::string &, const Eigen::Matrix4d &,
                       const std::string &, Mesh &, unsigned)>
        collect;
    collect = [&](const std::string &file, const std::string &id, const Eigen::Matrix4d &transform,
                  const std::string &top, Mesh &mesh, unsigned depth) {
        if (depth > 128 || !active.emplace(file, id).second)
            throw std::runtime_error("Cyclic or excessively nested 3MF components");
        auto &m = model(file);
        auto it = m.objects.find(id);
        if (it == m.objects.end())
            throw std::runtime_error("Missing 3MF object " + id);
        auto obj = it->second;
        if (std::string(obj.attribute("type").as_string("model")) != "model") {
            active.erase({file, id});
            return;
        }
        if (auto data = child(obj, "mesh")) {
            auto offset = mesh.points.size();
            for (auto v : child(data, "vertices").children())
                if (local(v.name()) == "vertex") {
                    Eigen::Vector4d p(coordinate(v, "x") * m.scale, coordinate(v, "y") * m.scale,
                                      coordinate(v, "z") * m.scale, 1);
                    mesh.points.push_back((transform * p).head<3>());
                }
            unsigned label = part_label++;
            for (auto t : child(data, "triangles").children())
                if (local(t.name()) == "triangle") {
                    Tri tri;
                    for (int k = 0; k < 3; ++k) {
                        auto key = "v" + std::to_string(k + 1);
                        auto a = t.attribute(key.c_str());
                        std::string value = a.value();
                        if (value.empty() ||
                            !std::all_of(value.begin(), value.end(),
                                         [](unsigned char c) { return std::isdigit(c); }))
                            throw std::runtime_error("Invalid 3MF triangle index");
                        auto index = std::stoull(value);
                        if (index >= mesh.points.size() - offset || index + offset > UINT32_MAX)
                            throw std::runtime_error("3MF triangle index out of range");
                        tri[k] = std::uint32_t(offset + index);
                    }
                    if (transform.block<3, 3>(0, 0).determinant() < 0)
                        std::swap(tri[1], tri[2]);
                    mesh.triangles.push_back(tri);
                    mesh.parts.push_back(label);
                }
        }
        for (auto c : child(obj, "components").children())
            if (local(c.name()) == "component") {
                std::string cid = c.attribute("objectid").value();
                auto kind = kinds.find({top, cid});
                if (kind != kinds.end() && kind->second != "normal_part")
                    continue;
                auto ref = attr(c, "path");
                std::string sub = ref ? path(ref.value()) : file;
                collect(sub, cid, transform * matrix(c.attribute("transform").value(), m.scale),
                        top, mesh, depth + 1);
            }
        active.erase({file, id});
    };
    auto &main = model(root);
    std::vector<Object> out;
    std::map<std::string, unsigned> names;
    for (auto item : child(main.document->document_element(), "build").children())
        if (local(item.name()) == "item") {
            std::string id = item.attribute("objectid").value();
            Object object;
            auto found = main.objects.find(id);
            if (found == main.objects.end())
                throw std::runtime_error("Missing 3MF build object");
            auto title = titles.find(id);
            object.name =
                safe_name(title == titles.end()
                              ? found->second.attribute("name").as_string(("object " + id).c_str())
                              : title->second);
            auto instance = ++names[object.name];
            if (instance > 1)
                object.name += " (" + std::to_string(instance) + ")";
            collect(root, id, matrix(item.attribute("transform").value(), main.scale), id,
                    object.mesh, 0);
            if (!object.mesh.triangles.empty()) {
                weld(object.mesh);
                out.push_back(std::move(object));
            }
        }
    if (out.empty())
        throw std::runtime_error("3MF contains no printable objects");
    return out;
}
} // namespace stl2curves::detail

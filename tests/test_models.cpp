#include "model_internal.hpp"
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRep_Builder.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <TopoDS_Compound.hxx>
#include <cstring>
#include <fstream>
#include <gp_Circ.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <iostream>
#include <nlohmann/json.hpp>
#include <set>

using namespace stl2curves;
using namespace stl2curves::detail;
using Json = nlohmann::json;
namespace {
void check(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
template <class F> void rejects(F &&f, const char *message) {
    bool failed = false;
    try {
        f();
    } catch (const std::invalid_argument &) {
        failed = true;
    }
    check(failed, message);
}
std::uint32_t u32(const Bytes &bytes, std::size_t pos) {
    check(pos + 4 <= bytes.size(), "truncated GLB");
    return bytes[pos] | (std::uint32_t(bytes[pos + 1]) << 8) |
           (std::uint32_t(bytes[pos + 2]) << 16) | (std::uint32_t(bytes[pos + 3]) << 24);
}
Json json(const GltfFile &file) {
    check(u32(file.data, 0) == 0x46546c67 && u32(file.data, 4) == 2 &&
              u32(file.data, 8) == file.data.size(),
          "GLB header");
    check(u32(file.data, 16) == 0x4e4f534a, "GLB JSON chunk");
    return Json::parse(file.data.begin() + 20, file.data.begin() + 20 + u32(file.data, 12));
}
ShapeRef selection_ref(const Json &value) {
    ShapeRef out;
    out.model_id = value.at("modelId");
    out.brepgraph_id = GraphId{value["brepgraph_id"]["kind"], value["brepgraph_id"]["index"]};
    out.usage_path.emplace();
    for (const auto &step : value["usagePath"]) {
        UsagePathStep entry{
            {step["node"]["kind"], step["node"]["index"]}, std::nullopt, step["stepIndex"]};
        if (step.contains("ref"))
            entry.ref = GraphId{step["ref"]["kind"], step["ref"]["index"]};
        out.usage_path->push_back(entry);
    }
    return out;
}
double property(const std::vector<ShapeProperty> &properties, const std::string &name) {
    for (const auto &p : properties)
        if (p.name == name)
            return p.value.at(0);
    throw std::runtime_error("missing property: " + name);
}
std::shared_ptr<Model> retain(const TopoDS_Shape &shape, std::string name = "fixture") {
    Diagnostics info;
    inspect(shape, info);
    return ModelAccess::create(std::move(name), shape, info, .01);
}
void cube(const char *output, const char *picked) {
    auto model = retain(BRepPrimAPI_MakeBox(1, 1, 1).Shape(), "cube");
    auto first = model->export_gltf();
    if (output) {
        std::ofstream stream(output, std::ios::binary);
        stream.write(reinterpret_cast<const char *>(first.data.data()), first.data.size());
        check(bool(stream), "write native GLB fixture");
    }
    if (picked) {
        std::ifstream stream(picked);
        Json selection;
        stream >> selection;
        auto selected = selection_ref(selection);
        check(property(model->get_shape_properties(selected), "shape_type") == TopAbs_FACE,
              "Three.js picked face resolves to native OCCT");
        check(model->get_closest_points(selected, *model, selected).distance < 1e-8,
              "Three.js selection closest points");
    }
    auto parsed = json(first);
    check(first.name == "cube.glb" && first.warnings.empty(), "cube export warnings/name");
    check(parsed["nodes"][0]["matrix"] ==
              Json::array({.001, 0, 0, 0, 0, 0, -.001, 0, 0, .001, 0, 0, 0, 0, 0, 1}),
          "GLB axis/unit transform");
    check(parsed["nodes"][0]["extras"]["coordinates"]["modelUnits"] == "mm",
          "GLB coordinate metadata");
    std::map<int, int> counts;
    std::vector<ShapeRef> edges;
    for (const auto &node : parsed["nodes"]) {
        auto selected = selection_ref(node["extras"]);
        auto properties = model->get_shape_properties(selected);
        check(!properties.empty(), "exported selection must resolve");
        if (!node.contains("mesh"))
            continue;
        int kind = selected.brepgraph_id->kind;
        ++counts[kind];
        const auto &primitive = parsed["meshes"][node["mesh"].get<int>()]["primitives"][0];
        check(primitive["mode"] == (kind == 2 ? 4 : kind == 4 ? 1 : 0), "GLB primitive mode");
        if (kind == 4) {
            check(std::abs(property(properties, "length") - 1) < 1e-9, "cube edge length");
            edges.push_back(selected);
        }
        if (kind == 2)
            check(primitive["attributes"].contains("NORMAL"), "face normals");
    }
    check(counts[2] == 6 && counts[4] == 12 && counts[5] == 8,
          "cube selectable face/edge/vertex counts");
    auto step = model->export_step();
    check(step.diagnostics.valid && step.diagnostics.faces == 6 &&
              std::abs(step.diagnostics.volume - 1) < 1e-8,
          "explicit STEP validation");
    auto second = json(model->export_gltf());
    check(parsed["nodes"] == second["nodes"], "STEP export changes graph identifiers");
    check(model->diagnostics().faces == 6 && std::abs(model->diagnostics().volume - 1) < 1e-9,
          "exports mutate model diagnostics");
    auto distance = model->get_closest_points(edges.front(), *model, edges.front());
    check(distance.distance < 1e-8 && distance.p1.size() == 3, "same edge distance");
    auto invalid = edges.front();
    invalid.usage_path->back().node.index += 100;
    rejects([&] { model->get_shape_properties(invalid); }, "inconsistent usage path accepted");
    invalid = edges.front();
    invalid.model_id++;
    rejects([&] { model->get_shape_properties(invalid); }, "wrong model selection accepted");
    GltfExportOptions options;
    options.faces = false;
    options.vertices = false;
    auto lines = json(model->export_gltf(options));
    check(lines["meshes"].size() == 12, "edge-only export");
    options.linear_deflection = 0;
    rejects([&] { model->export_gltf(options); }, "invalid tessellation accepted");
}
void placed() {
    auto box = BRepPrimAPI_MakeBox(1, 1, 1).Shape();
    gp_Trsf transform;
    transform.SetTranslation(gp_Vec(5, 0, 0));
    auto moved = box.Moved(TopLoc_Location(transform));
    auto first = retain(box), second = retain(moved);
    auto distance = first->get_closest_points(first->root(), *second, second->root());
    check(std::abs(distance.distance - 4) < 1e-8, "cross-model placed solid distance");
    auto moved_glb = json(second->export_gltf());
    for (const auto &node : moved_glb["nodes"]) {
        auto selected = selection_ref(node["extras"]);
        if (selected.brepgraph_id->kind == 2) {
            auto separation = second->get_closest_points(selected, *first, first->root());
            check(separation.distance >= 4 - 1e-8 && separation.distance <= 5 + 1e-8,
                  "root placement preserved in face selection");
        }
    }
    TopoDS_Compound compound;
    BRep_Builder builder;
    builder.MakeCompound(compound);
    builder.Add(compound, box);
    builder.Add(compound, moved);
    auto assembly = retain(compound);
    auto glb = json(assembly->export_gltf());
    std::map<std::pair<int, unsigned>, std::vector<ShapeRef>> by_id;
    unsigned solids = 0, edge_count = 0;
    for (const auto &node : glb["nodes"]) {
        auto selected = selection_ref(node["extras"]);
        int kind = selected.brepgraph_id->kind;
        if (kind == 0)
            solids++;
        if (kind == 4) {
            edge_count++;
            auto properties = assembly->get_shape_properties(selected);
            check(std::abs(property(properties, "length") - 1) < 1e-9, "placed edge query");
            const auto &primitive = glb["meshes"][node["mesh"].get<int>()]["primitives"][0];
            const auto &accessor = glb["accessors"][primitive["attributes"]["POSITION"].get<int>()];
            check(accessor["min"][0].get<double>() >= -1e-8 &&
                      accessor["max"][0].get<double>() <= 6 + 1e-8,
                  "repeated instance placement applied more than once");
        }
        by_id[{kind, selected.brepgraph_id->index}].push_back(selected);
    }
    check(solids == 2 && edge_count == 24, "compound grouping and repeated placements");
    bool ambiguity_tested = false;
    for (const auto &entry : by_id)
        if (entry.second.size() > 1) {
            auto selected = entry.second.front();
            selected.usage_path.reset();
            rejects([&] { assembly->get_shape_properties(selected); },
                    "ambiguous definition accepted");
            ambiguity_tested = true;
            break;
        }
    check(ambiguity_tested, "fixture must exercise shared graph definitions");
}
void curves() {
    auto circle =
        retain(BRepBuilderAPI_MakeEdge(gp_Circ(gp_Ax2(gp_Pnt(2, 3, 4), gp_Dir(0, 0, 1)), 2), 0,
                                       std::acos(-1.) / 2)
                   .Shape());
    auto properties = circle->get_shape_properties(circle->root());
    check(std::abs(property(properties, "radius") - 2) < 1e-10 &&
              std::abs(property(properties, "length") - std::acos(-1.)) < 1e-9 &&
              std::abs(property(properties, "angle_degrees") - 90) < 1e-9,
          "arc properties");
    bool step_failed = false;
    try {
        circle->export_step();
    } catch (const std::runtime_error &) {
        step_failed = true;
    }
    check(step_failed && property(circle->get_shape_properties(circle->root()), "radius") == 2,
          "failed solid STEP validation must leave the model queryable");
    auto unbounded =
        retain(BRepBuilderAPI_MakeEdge(gp_Lin(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0))).Shape());
    check(unbounded->get_shape_properties(unbounded->root()).size() == 2,
          "unbounded edge metadata");
    auto sphere = retain(BRepPrimAPI_MakeSphere(2).Shape());
    auto before = sphere->root();
    auto sphere_file = sphere->export_gltf();
    check(!sphere_file.warnings.empty(), "degenerate edges retain metadata and produce warnings");
    auto glb = json(sphere_file);
    check(!glb["meshes"].empty(), "sphere tessellation");
    check(property(sphere->get_shape_properties(before), "shape_type") == TopAbs_SOLID,
          "sphere root identity");
    auto plane = retain(
        BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), -1, 1, -1, 1).Shape());
    auto distance = plane->get_closest_points(plane->root(), *sphere, sphere->root());
    check(distance.distance < 1e-7, "intersecting shape distance");
    for (const auto &surface : std::vector<Handle(Geom_Surface)>{
             new Geom_SphericalSurface(gp_Ax3(), 2), new Geom_ToroidalSurface(gp_Ax3(), 3, 1)}) {
        auto natural = retain(BRepBuilderAPI_MakeFace(surface, .01).Face());
        auto natural_glb = json(natural->export_gltf());
        check(!natural_glb["meshes"].empty(), "natural-bound periodic face tessellation");
        check(property(natural->get_shape_properties(natural->root()), "shape_type") == TopAbs_FACE,
              "natural-bound face identity");
    }
}
} // namespace
int main(int argc, char **argv) {
    try {
        cube(argc > 1 ? argv[1] : nullptr, argc > 2 ? argv[2] : nullptr);
        placed();
        curves();
        std::cout << "Model/export/query tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

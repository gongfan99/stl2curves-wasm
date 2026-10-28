#include "model_internal.hpp"
#include "occt_error.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepLib_ToolTriangulatedShape.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <GCPnts_TangentialDeflection.hxx>
#include <IMeshTools_Parameters.hxx>
#include <Poly_Triangulation.hxx>
#include <Precision.hxx>
#include <TopoDS.hxx>
#include <cstring>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#define TINYGLTF_IMPLEMENTATION
#define TINYGLTF_NO_STB_IMAGE
#define TINYGLTF_NO_STB_IMAGE_WRITE
#define TINYGLTF_NO_EXTERNAL_IMAGE
#define TINYGLTF_NO_FS
#include <tiny_gltf.h>

namespace stl2curves {
namespace {
using Value = tinygltf::Value;
Value id_value(const GraphId &id) {
    return Value(Value::Object{{"kind", Value(id.kind)}, {"index", Value(double(id.index))}});
}
Value metadata(const ShapeRef &ref) {
    Value::Array steps;
    for (const auto &step : *ref.usage_path) {
        Value::Object value{{"node", id_value(step.node)}, {"stepIndex", Value(step.step_index)}};
        if (step.ref)
            value["ref"] = id_value(*step.ref);
        steps.emplace_back(value);
    }
    return Value(Value::Object{{"modelId", Value(double(ref.model_id))},
                               {"brepgraph_id", id_value(*ref.brepgraph_id)},
                               {"usagePath", Value(steps)}});
}
std::string path_key(const std::vector<UsagePathStep> &steps, std::size_t count) {
    std::ostringstream s;
    for (std::size_t i = 0; i < count; ++i) {
        const auto &step = steps[i];
        s << step.node.kind << ':' << step.node.index << ':' << step.step_index << ':';
        if (step.ref)
            s << step.ref->kind << ':' << step.ref->index;
        s << '/';
    }
    return s.str();
}
bool container(int kind) {
    return kind == 0 || kind == 1 || kind == 6 || kind == 7;
}
struct Writer {
    tinygltf::Model model;
    std::vector<unsigned char> binary;
    template <class T>
    int accessor(const std::vector<T> &values, int component, int type, int target) {
        while (binary.size() % 4)
            binary.push_back(0);
        tinygltf::BufferView view;
        view.buffer = 0;
        view.byteOffset = binary.size();
        view.byteLength = values.size() * sizeof(T);
        view.target = target;
        const auto *bytes = reinterpret_cast<const unsigned char *>(values.data());
        binary.insert(binary.end(), bytes, bytes + view.byteLength);
        tinygltf::Accessor accessor;
        accessor.bufferView = int(model.bufferViews.size());
        accessor.componentType = component;
        accessor.type = type;
        accessor.count = values.size() / (type == TINYGLTF_TYPE_VEC3 ? 3 : 1);
        if (type == TINYGLTF_TYPE_VEC3) {
            accessor.minValues.assign(3, std::numeric_limits<double>::infinity());
            accessor.maxValues.assign(3, -std::numeric_limits<double>::infinity());
            for (std::size_t i = 0; i < values.size(); ++i) {
                if (!std::isfinite(double(values[i])))
                    throw std::runtime_error("Non-finite GLB geometry");
                accessor.minValues[i % 3] = std::min(accessor.minValues[i % 3], double(values[i]));
                accessor.maxValues[i % 3] = std::max(accessor.maxValues[i % 3], double(values[i]));
            }
        }
        model.bufferViews.push_back(view);
        model.accessors.push_back(accessor);
        return int(model.accessors.size()) - 1;
    }
    void mesh(int node, const ShapeRef &ref, int mode, const std::vector<float> &positions,
              const std::vector<float> &normals, const std::vector<std::uint32_t> &indices) {
        tinygltf::Primitive primitive;
        primitive.mode = mode;
        primitive.material = mode == TINYGLTF_MODE_TRIANGLES ? 0 : 1;
        primitive.attributes["POSITION"] =
            accessor(positions, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC3,
                     TINYGLTF_TARGET_ARRAY_BUFFER);
        if (!normals.empty())
            primitive.attributes["NORMAL"] =
                accessor(normals, TINYGLTF_COMPONENT_TYPE_FLOAT, TINYGLTF_TYPE_VEC3,
                         TINYGLTF_TARGET_ARRAY_BUFFER);
        if (!indices.empty())
            primitive.indices =
                accessor(indices, TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT, TINYGLTF_TYPE_SCALAR,
                         TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER);
        tinygltf::Mesh mesh;
        mesh.extras = metadata(ref);
        mesh.primitives.push_back(primitive);
        model.nodes[node].mesh = int(model.meshes.size());
        model.meshes.push_back(std::move(mesh));
    }
};
void point(std::vector<float> &values, const gp_Pnt &p) {
    values.insert(values.end(), {float(p.X()), float(p.Y()), float(p.Z())});
}
} // namespace

GltfFile Model::export_gltf(const GltfExportOptions &options) const {
    if (!std::isfinite(options.linear_deflection) || options.linear_deflection <= 0 ||
        !std::isfinite(options.angular_deflection) || options.angular_deflection <= 0 ||
        options.angular_deflection > std::acos(-1.) ||
        !(options.faces || options.edges || options.vertices))
        throw std::invalid_argument("Invalid GLB meshing options");
    std::lock_guard<std::recursive_mutex> guard(impl_->mutex);
    try {
        Writer writer;
        auto &model = writer.model;
        GltfFile result{impl_->name + ".glb", {}, {}};
        model.asset.version = "2.0";
        model.asset.generator = "stl2curves-wasm 0.2.0";
        model.defaultScene = 0;
        model.scenes.emplace_back();
        model.scenes[0].nodes = {0};
        tinygltf::Material surface, line;
        surface.pbrMetallicRoughness.baseColorFactor = {0.7, 0.7, 0.7, 1};
        surface.pbrMetallicRoughness.metallicFactor = 0;
        surface.pbrMetallicRoughness.roughnessFactor = 0.8;
        line.pbrMetallicRoughness.baseColorFactor = {0.1, 0.1, 0.1, 1};
        model.materials = {surface, line};
        tinygltf::Node root;
        root.name = impl_->name;
        auto root_extras = metadata(impl_->occurrences.front().ref).Get<Value::Object>();
        root_extras["coordinates"] = Value(Value::Object{{"modelUnits", Value(std::string("mm"))},
                                                         {"modelUp", Value(std::string("Z"))},
                                                         {"gltfUnits", Value(std::string("m"))},
                                                         {"gltfUp", Value(std::string("Y"))}});
        root.extras = Value(root_extras);
        root.matrix = {.001, 0, 0, 0, 0, 0, -.001, 0, 0, .001, 0, 0, 0, 0, 0, 1};
        model.nodes.push_back(root);
        std::map<std::string, int> groups;
        const auto &root_path = *impl_->occurrences.front().ref.usage_path;
        groups[path_key(root_path, root_path.size())] = 0;
        std::set<std::string> rendered;
        for (std::size_t i = 0; i < impl_->occurrences.size(); ++i) {
            const auto &entry = impl_->occurrences[i];
            const auto &ref = entry.ref;
            const int kind = ref.brepgraph_id->kind;
            const bool group = container(kind);
            if (!group && kind != 2 && kind != 4 && kind != 5)
                continue;
            if (!group && ((kind == 2 && !options.faces) || (kind == 4 && !options.edges) ||
                           (kind == 5 && !options.vertices)))
                continue;
            const auto &path = *ref.usage_path;
            int parent = 0;
            for (std::size_t n = path.size(); n > 1; --n) {
                // Shared edges/vertices are attached to the owning solid, not a face's wire.
                if ((kind == 4 || kind == 5) && path[n - 2].node.kind != 0 &&
                    path[n - 2].node.kind != 6 && path[n - 2].node.kind != 7)
                    continue;
                auto found = groups.find(path_key(path, n - 1));
                if (found != groups.end()) {
                    parent = found->second;
                    break;
                }
            }
            if (kind == 4 || kind == 5) {
                std::ostringstream key;
                key << parent << ':' << kind << ':' << ref.brepgraph_id->index
                    << std::setprecision(17);
                const auto transform = entry.shape.Location().Transformation();
                for (int row = 1; row <= 3; ++row)
                    for (int col = 1; col <= 4; ++col)
                        key << ':' << transform.Value(row, col);
                if (!rendered.insert(key.str()).second)
                    continue;
            }
            int node = 0;
            if (i != 0) {
                tinygltf::Node value;
                value.name = std::string(group       ? "shape_"
                                         : kind == 2 ? "face_"
                                         : kind == 4 ? "edge_"
                                                     : "vertex_") +
                             std::to_string(ref.brepgraph_id->index);
                value.extras = metadata(ref);
                node = int(model.nodes.size());
                model.nodes.push_back(value);
                model.nodes[parent].children.push_back(node);
            }
            if (group) {
                groups[path_key(path, path.size())] = node;
                continue;
            }
            std::vector<float> positions, normals;
            std::vector<std::uint32_t> indices;
            if (kind == 2) {
                BRepBuilderAPI_Copy copy(entry.shape, true, false);
                auto face = TopoDS::Face(copy.Shape());
                IMeshTools_Parameters parameters;
                parameters.Deflection = options.linear_deflection;
                parameters.Angle = options.angular_deflection;
                parameters.InParallel = false;
                BRepMesh_IncrementalMesh mesh(face, parameters);
                if (!mesh.IsDone())
                    throw std::runtime_error("Face tessellation failed");
                TopLoc_Location location;
                auto triangulation = BRep_Tool::Triangulation(face, location);
                if (triangulation.IsNull() || !triangulation->NbTriangles())
                    throw std::runtime_error("Face tessellation produced no triangles");
                BRepLib_ToolTriangulatedShape::ComputeNormals(face, triangulation);
                if (!triangulation->HasNormals())
                    triangulation->ComputeNormals();
                const auto transform = location.Transformation();
                const bool reversed = face.Orientation() == TopAbs_REVERSED;
                for (int n = 1; n <= triangulation->NbNodes(); ++n) {
                    point(positions, triangulation->Node(n).Transformed(transform));
                    auto normal = triangulation->Normal(n).Transformed(transform);
                    if (reversed)
                        normal.Reverse();
                    normals.insert(normals.end(),
                                   {float(normal.X()), float(normal.Y()), float(normal.Z())});
                }
                for (int t = 1; t <= triangulation->NbTriangles(); ++t) {
                    int a, b, c;
                    triangulation->Triangle(t).Get(a, b, c);
                    if (reversed != transform.IsNegative())
                        std::swap(b, c);
                    indices.insert(indices.end(), {std::uint32_t(a - 1), std::uint32_t(b - 1),
                                                   std::uint32_t(c - 1)});
                }
                writer.mesh(node, ref, TINYGLTF_MODE_TRIANGLES, positions, normals, indices);
            } else if (kind == 4) {
                auto edge = TopoDS::Edge(entry.shape);
                if (!BRep_Tool::Degenerated(edge)) {
                    BRepAdaptor_Curve curve(edge);
                    if (std::isfinite(curve.FirstParameter()) &&
                        std::isfinite(curve.LastParameter()) &&
                        !Precision::IsInfinite(curve.FirstParameter()) &&
                        !Precision::IsInfinite(curve.LastParameter())) {
                        GCPnts_TangentialDeflection points(curve, options.angular_deflection,
                                                           options.linear_deflection);
                        for (int p = 1; p <= points.NbPoints(); ++p)
                            point(positions, points.Value(p));
                        for (std::size_t p = 1; p < positions.size() / 3; ++p)
                            indices.insert(indices.end(), {std::uint32_t(p - 1), std::uint32_t(p)});
                    }
                }
                if (indices.empty())
                    result.warnings.push_back("Edge " + std::to_string(ref.brepgraph_id->index) +
                                              " has no drawable geometry");
                else
                    writer.mesh(node, ref, TINYGLTF_MODE_LINE, positions, normals, indices);
            } else {
                point(positions, BRep_Tool::Pnt(TopoDS::Vertex(entry.shape)));
                writer.mesh(node, ref, TINYGLTF_MODE_POINTS, positions, {}, {});
            }
        }
        if (!writer.binary.empty()) {
            tinygltf::Buffer buffer;
            buffer.data = std::move(writer.binary);
            model.buffers.push_back(std::move(buffer));
        }
        tinygltf::TinyGLTF exporter;
        std::ostringstream stream(std::ios::out | std::ios::binary);
        if (!exporter.WriteGltfSceneToStream(&model, stream, false, true))
            throw std::runtime_error("GLB serialization failed");
        const auto data = stream.str();
        result.data.assign(data.begin(), data.end());
        return result;
    } catch (const Standard_Failure &error) {
        throw std::runtime_error(detail::failure_message(error));
    }
}
} // namespace stl2curves

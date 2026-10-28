#include "model_internal.hpp"
#include "occt_error.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGraph_ChildExplorer.hxx>
#include <BRepGraph_ShapesView.hxx>
#include <BRep_Tool.hxx>
#include <Precision.hxx>
#include <TopoDS.hxx>
#include <limits>

namespace stl2curves {
namespace {
GraphId graph_id(BRepGraph_NodeId id) {
    return {int(id.NodeKind), id.Index};
}
std::vector<UsagePathStep> path_value(const BRepGraph_UsagePath &path) {
    std::vector<UsagePathStep> result;
    for (std::size_t i = 0; i < path.Size(); ++i) {
        const auto &step = path.Value(i);
        UsagePathStep value{graph_id(step.Node), std::nullopt, step.StepIndex};
        if (step.Ref.IsValid())
            value.ref = GraphId{int(step.Ref.RefKind), step.Ref.Index};
        result.push_back(value);
    }
    return result;
}
std::vector<double> coordinates(const gp_Pnt &p) {
    return {p.X(), p.Y(), p.Z()};
}
} // namespace

std::shared_ptr<Model> detail::ModelAccess::create(std::string name, const TopoDS_Shape &shape,
                                                   Diagnostics diagnostics, double tolerance) {
    static std::mutex id_mutex;
    static std::uint64_t next_id = 1;
    std::uint32_t id;
    {
        std::lock_guard<std::mutex> guard(id_mutex);
        if (next_id > std::numeric_limits<std::uint32_t>::max())
            throw std::runtime_error("Model handle space exhausted");
        id = std::uint32_t(next_id++);
    }
    if (shape.IsNull())
        throw std::invalid_argument("Cannot retain a null shape");
    auto impl = std::make_unique<Model::Impl>();
    impl->id = id;
    impl->name = std::move(name);
    impl->shape = shape;
    impl->diagnostics = std::move(diagnostics);
    impl->tolerance = tolerance;
    BRepGraph::ShapesView::Options options;
    options.CreateAutoProduct = false;
    options.Parallel = false;
    const auto added = impl->graph.Shapes().Add(shape, options);
    if (!added.IsOk() || !added.TopologyRoot.IsValid())
        throw std::runtime_error("BRepGraph construction failed");
    impl->root = added.TopologyRoot;
    ShapeRef root{impl->id, graph_id(impl->root),
                  std::vector<UsagePathStep>{{graph_id(impl->root), std::nullopt, -1}}};
    impl->occurrences.push_back({root, shape});
    for (BRepGraph_ChildExplorer it(impl->graph, impl->root); it.More(); it.Next()) {
        const auto instance = it.Current();
        const auto kind = instance.DefId.NodeKind;
        // Coedges have no independent TopoDS counterpart; retain them in paths.
        if (kind == BRepGraph_NodeId::Kind::CoEdge || kind == BRepGraph_NodeId::Kind::Product ||
            kind == BRepGraph_NodeId::Kind::Occurrence)
            continue;
        auto path = path_value(it.CurrentUsagePath());
        if (path.size() == 1 && instance.DefId == impl->root)
            continue;
        auto placed = impl->graph.Shapes().Shape(instance.DefId);
        if (placed.IsNull())
            throw std::runtime_error("BRepGraph subshape reconstruction failed");
        placed.Move(instance.Location);
        placed.Orientation(instance.Orientation);
        impl->occurrences.push_back(
            {{impl->id, graph_id(instance.DefId), std::move(path)}, placed});
    }
    return std::shared_ptr<Model>(new Model(std::move(impl)));
}

const TopoDS_Shape &Model::Impl::resolve(const ShapeRef &ref) const {
    if (ref.model_id != id)
        throw std::invalid_argument("Selection belongs to a different model");
    if (!ref.brepgraph_id) {
        if (ref.usage_path)
            throw std::invalid_argument("A usage path requires brepgraph_id");
        return shape;
    }
    const TopoDS_Shape *found = nullptr;
    for (const auto &entry : occurrences) {
        if (!(*entry.ref.brepgraph_id == *ref.brepgraph_id))
            continue;
        if (ref.usage_path && entry.ref.usage_path != ref.usage_path)
            continue;
        if (found && !found->IsSame(entry.shape))
            throw std::invalid_argument("Ambiguous graph ID; supply its usagePath");
        found = &entry.shape;
    }
    if (!found)
        throw std::invalid_argument("Unknown graph ID or inconsistent usagePath");
    return *found;
}

Model::Model(std::unique_ptr<Impl> impl) : impl_(std::move(impl)) {
}
Model::~Model() = default;
std::uint32_t Model::id() const {
    return impl_->id;
}
const std::string &Model::name() const {
    return impl_->name;
}
const Diagnostics &Model::diagnostics() const {
    return impl_->diagnostics;
}
ShapeRef Model::root() const {
    return impl_->occurrences.front().ref;
}
OutputFile Model::export_step() const {
    // OCCT STEP configuration is process-wide; model state stays read-only.
    static std::mutex step_mutex;
    std::lock_guard<std::mutex> step_guard(step_mutex);
    std::lock_guard<std::recursive_mutex> guard(impl_->mutex);
    BRepBuilderAPI_Copy copy(impl_->shape, true, false);
    auto diagnostics = impl_->diagnostics;
    try {
        auto bytes = detail::validated_step(copy.Shape(), impl_->tolerance, diagnostics);
        return {impl_->name + ".step", std::move(bytes), std::move(diagnostics)};
    } catch (const Standard_Failure &error) {
        throw std::runtime_error(detail::failure_message(error));
    }
}
ClosestPoints Model::get_closest_points(const ShapeRef &selection, const Model &other,
                                        const ShapeRef &other_selection) const {
    const auto compute = [&]() {
        try {
            BRepExtrema_DistShapeShape distance;
            distance.SetMultiThread(false);
            distance.LoadS1(impl_->resolve(selection));
            distance.LoadS2(other.impl_->resolve(other_selection));
            if (!distance.Perform() || !distance.IsDone() || distance.NbSolution() < 1)
                throw std::runtime_error("Closest-point computation failed");
            auto result = ClosestPoints{coordinates(distance.PointOnShape1(1)),
                                        coordinates(distance.PointOnShape2(1)), distance.Value()};
            if (!std::isfinite(result.distance))
                throw std::runtime_error("Non-finite distance");
            return result;
        } catch (const Standard_Failure &error) {
            throw std::runtime_error(detail::failure_message(error));
        }
    };
    if (this == &other) {
        std::lock_guard<std::recursive_mutex> guard(impl_->mutex);
        return compute();
    }
    std::scoped_lock guard(impl_->mutex, other.impl_->mutex);
    return compute();
}
std::vector<ShapeProperty> Model::get_shape_properties(const ShapeRef &selection) const {
    std::lock_guard<std::recursive_mutex> guard(impl_->mutex);
    const auto &shape = impl_->resolve(selection);
    std::vector<ShapeProperty> result;
    auto type = GeomAbs_OtherCurve;
    const auto finish = [&]() {
        result.push_back({"shape_type", {double(shape.ShapeType())}});
        if (shape.ShapeType() == TopAbs_EDGE)
            result.push_back({"curve_type", {double(type)}});
        return result;
    };
    if (shape.ShapeType() != TopAbs_EDGE)
        return finish();
    try {
        const auto edge = TopoDS::Edge(shape);
        if (BRep_Tool::Degenerated(edge)) {
            double first, last;
            if (!BRep_Tool::Curve(edge, first, last).IsNull()) {
                try {
                    type = BRepAdaptor_Curve(edge).GetType();
                } catch (const Standard_Failure &) {
                }
            }
            return finish();
        }
        BRepAdaptor_Curve curve(edge);
        type = curve.GetType();
        if (type != GeomAbs_Line && type != GeomAbs_Circle)
            return finish();
        double first = curve.FirstParameter(), last = curve.LastParameter();
        if (!std::isfinite(first) || !std::isfinite(last))
            throw std::runtime_error("Non-finite edge parameters");
        if (Precision::IsInfinite(first) || Precision::IsInfinite(last))
            return finish();
        auto start = curve.Value(first), end = curve.Value(last);
        double length = type == GeomAbs_Circle ? curve.Circle().Radius() * std::abs(last - first)
                                               : start.Distance(end);
        if (!std::isfinite(length))
            throw std::runtime_error("Non-finite edge length");
        if (length == 0)
            return finish();
        result = {{"start_point", coordinates(start)},
                  {"end_point", coordinates(end)},
                  {"length", {length}}};
        if (type == GeomAbs_Circle) {
            auto circle = curve.Circle();
            result.push_back({"center", coordinates(circle.Location())});
            result.push_back({"radius", {circle.Radius()}});
            result.push_back({"angle_degrees", {std::abs(last - first) * 180 / std::acos(-1.)}});
        }
        return finish();
    } catch (const Standard_Failure &error) {
        throw std::runtime_error(detail::failure_message(error));
    }
}
} // namespace stl2curves

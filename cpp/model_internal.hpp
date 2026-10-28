#pragma once
#include "internal.hpp"
#include <BRepGraph.hxx>
#include <BRepGraph_UsagePath.hxx>
#include <mutex>

namespace stl2curves {
struct Model::Impl {
    struct Occurrence {
        ShapeRef ref;
        TopoDS_Shape shape;
    };
    std::uint32_t id;
    std::string name;
    Diagnostics diagnostics;
    double tolerance;
    TopoDS_Shape shape;
    BRepGraph graph;
    BRepGraph_NodeId root;
    std::vector<Occurrence> occurrences;
    mutable std::recursive_mutex mutex;
    const TopoDS_Shape &resolve(const ShapeRef &ref) const;
};
namespace detail {
struct ModelAccess {
    static std::shared_ptr<Model> create(std::string name, const TopoDS_Shape &shape,
                                         Diagnostics diagnostics, double tolerance);
    static const Model::Impl &get(const Model &model) { return *model.impl_; }
};
Bytes validated_step(const TopoDS_Shape &shape, double tolerance, Diagnostics &diagnostics);
} // namespace detail
} // namespace stl2curves

#include "model_internal.hpp"
#include "occt_error.hpp"
#include <BRepAdaptor_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <STEPControl_Reader.hxx>
#include <ShapeAnalysis_ShapeTolerance.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <sstream>

namespace stl2curves {
std::shared_ptr<Model> import_step(const InputFile &input) {
    if (input.name.empty() || input.data.empty())
        throw std::invalid_argument("STEP input needs a name and nonempty data");
    try {
        STEPControl_Reader reader;
        std::istringstream stream(std::string(input.data.begin(), input.data.end()));
        if (reader.ReadStream(input.name.c_str(), stream) != IFSelect_RetDone)
            throw std::runtime_error("Could not read this STEP file");
        // All public properties, distances, and GLB coordinates use millimeters,
        // including STEP files whose source units are inches or meters.
        reader.SetSystemLengthUnit(1.0);
        if (!reader.TransferRoots())
            throw std::runtime_error("STEP file contains no transferable geometry");
        auto shape = reader.OneShape();
        if (shape.IsNull())
            throw std::runtime_error("STEP file contains no geometry");
        detail::ShapeMap faces, solids, vertices;
        TopExp::MapShapes(shape, TopAbs_FACE, faces);
        TopExp::MapShapes(shape, TopAbs_SOLID, solids);
        TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
        if (faces.IsEmpty() && vertices.IsEmpty())
            throw std::runtime_error("STEP file contains no viewable geometry");
        Diagnostics diagnostics;
        diagnostics.faces = faces.Extent();
        diagnostics.solids = solids.Extent();
        diagnostics.valid = BRepCheck_Analyzer(shape).IsValid();
        diagnostics.design_unit = "mm";
        diagnostics.design_label = "STEP";
        static const char *kinds[] = {"plane",     "cylinder", "cone",    "sphere",
                                      "torus",     "bezier",   "bspline", "revolution",
                                      "extrusion", "offset",   "other"};
        for (int i = 1; i <= faces.Extent(); ++i) {
            auto type = unsigned(BRepAdaptor_Surface(TopoDS::Face(faces(i))).GetType());
            ++diagnostics.surface_counts[kinds[std::min(type, 10u)]];
        }
        if (!solids.IsEmpty()) {
            GProp_GProps properties;
            BRepGProp::VolumeProperties(shape, properties);
            diagnostics.volume = properties.Mass();
        }
        if (!diagnostics.valid)
            diagnostics.warnings.push_back("The imported STEP topology has validity issues");
        ShapeAnalysis_ShapeTolerance analysis;
        double tolerance = std::max(1e-7, analysis.Tolerance(shape, 1));
        auto name = input.name.substr(input.name.find_last_of("/\\") + 1);
        auto dot = name.find_last_of('.');
        if (dot != std::string::npos && dot > 0)
            name.resize(dot);
        return detail::ModelAccess::create(std::move(name), shape, std::move(diagnostics),
                                           tolerance);
    } catch (const Standard_Failure &error) {
        throw std::runtime_error(detail::failure_message(error, "STEP import failed"));
    }
}
} // namespace stl2curves

#include "internal.hpp"
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <ShapeAnalysis_ShapeTolerance.hxx>
#include <Standard_Version.hxx>
#if OCC_VERSION_HEX < 0x070800
#include <StepData.hxx>
#include <StepData_Protocol.hxx>
#include <StepData_StepModel.hxx>
#include <StepData_StepWriter.hxx>
#include <XSControl_WorkSession.hxx>
#endif
#include <sstream>

namespace stl2curves::detail {
Bytes write_step(const TopoDS_Shape &shape, double tolerance, bool surface_curves) {
    STEPControl_Writer writer;
    if (!surface_curves) {
        // An average file uncertainty can understate a seam's existing native
        // tolerance. Preserve its greatest tolerance on a recovery export,
        // capped by the caller's sewing tolerance, without changing surfaces.
        ShapeAnalysis_ShapeTolerance analysis;
        writer.SetTolerance(std::min(tolerance, std::max(1e-7, analysis.Tolerance(shape, 1))));
    }
#if OCC_VERSION_HEX >= 0x080000
    DESTEP_Parameters parameters;
    parameters.WriteSurfaceCurMode = surface_curves;
    if (writer.Transfer(shape, STEPControl_AsIs, parameters) != IFSelect_RetDone)
#else
    (void)surface_curves;
    if (writer.Transfer(shape, STEPControl_AsIs) != IFSelect_RetDone)
#endif
        throw std::runtime_error("STEP transfer failed");
    std::ostringstream stream;
#if OCC_VERSION_HEX >= 0x070800
    if (writer.WriteStream(stream) != IFSelect_RetDone)
        throw std::runtime_error("STEP serialization failed");
#else
    // Allows native development checks with Ubuntu's older system OCCT package.
    StepData_StepWriter serializer(writer.Model());
    serializer.SendModel(Handle(StepData_Protocol)::DownCast(writer.WS()->Protocol()));
    if (!serializer.Print(stream))
        throw std::runtime_error("STEP serialization failed");
#endif
    const auto text = stream.str();
    return Bytes(text.begin(), text.end());
}
bool verify_step(const Bytes &bytes, const TopoDS_Shape &shape, double tolerance,
                 Diagnostics *exported_info, TopoDS_Shape *exported_shape) {
    STEPControl_Reader reader;
    std::istringstream stream(std::string(bytes.begin(), bytes.end()));
    if (reader.ReadStream("memory.step", stream) != IFSelect_RetDone || reader.TransferRoots() == 0)
        return false;
    auto restored = reader.OneShape();
    // Retain the imported shape for bounded topology repair even when its
    // validity check fails. It is never returned to a caller as an accepted STEP.
    if (exported_shape)
        *exported_shape = restored;
    Diagnostics info;
    if (!inspect(restored, info))
        return false;
    Diagnostics original;
    if (!inspect(shape, original))
        return false;
    trace(("STEP volume before=" + std::to_string(original.volume) +
           " after=" + std::to_string(info.volume))
              .c_str());
    if (std::abs(original.volume - info.volume) >
        std::max(tolerance * tolerance * tolerance, std::abs(original.volume) * 1e-3 + 1e-3))
        return false;
    if (exported_info) {
        exported_info->faces = info.faces;
        exported_info->solids = info.solids;
        exported_info->volume = info.volume;
        exported_info->valid = info.valid;
        exported_info->surface_counts = std::move(info.surface_counts);
    }
    return true;
}
} // namespace stl2curves::detail

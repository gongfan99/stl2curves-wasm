#include "model_internal.hpp"
#include <BRepAdaptor_Curve.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace stl2curves;
using namespace stl2curves::detail;
namespace {
bool check_rims(const TopoDS_Shape &shape, const char *file, const char *label) {
    std::ifstream in(file);
    if (!in)
        throw std::runtime_error("Cannot open reference rim samples");
    ShapeMap edges;
    TopExp::MapShapes(shape, TopAbs_EDGE, edges);
    unsigned count;
    in >> count;
    bool ok = true;
    for (unsigned i = 0; i < count; ++i) {
        std::string name;
        unsigned n;
        double tolerance;
        in >> name >> n >> tolerance;
        std::vector<gp_Pnt> samples;
        for (unsigned j = 0; j < n; ++j) {
            double x, y, z;
            in >> x >> y >> z;
            samples.emplace_back(x, y, z);
        }
        if (!in || samples.empty())
            throw std::runtime_error("Malformed reference rim samples");
        double best = std::numeric_limits<double>::infinity();
        int best_edge = 0;
        for (int j = 1; j <= edges.Extent(); ++j) {
            auto edge = TopoDS::Edge(edges(j));
            BRepAdaptor_Curve adaptor(edge);
            if (adaptor.GetType() == GeomAbs_Line ||
                (adaptor.GetType() == GeomAbs_BSplineCurve && adaptor.Degree() < 2))
                continue;
            double first, last;
            auto curve = BRep_Tool::Curve(edge, first, last);
            if (curve.IsNull())
                continue;
            double worst = 0;
            for (const auto &p : samples) {
                // Restrict to this edge's trim: the old segmented edges can
                // share an underlying whole curve without being whole edges.
                GeomAPI_ProjectPointOnCurve project(p, curve, first, last);
                worst =
                    std::max(worst, project.NbPoints() ? project.LowerDistance()
                                                       : std::numeric_limits<double>::infinity());
                if (worst > best)
                    break;
            }
            if (worst < best) {
                best = worst;
                best_edge = j;
            }
        }
        std::cout << label << " RIM " << name << " wholeEdge=" << best_edge << " gap=" << best
                  << " allowed=" << tolerance << '\n';
        ok = ok && best <= tolerance;
    }
    return ok;
}
} // namespace
int main(int argc, char **argv) {
    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: inspect_conversion INPUT.stl OUTPUT.step [REFERENCE_RIMS.txt]\n";
        return 2;
    }
    try {
        std::ifstream in(argv[1], std::ios::binary);
        if (!in)
            throw std::runtime_error("Cannot open input mesh");
        Bytes bytes((std::istreambuf_iterator<char>(in)), {});
        Options options;
        options.details = true;
        auto result = convert(
            {{std::filesystem::path(argv[1]).filename().string(), std::move(bytes)}}, options);
        for (const auto &e : result.errors)
            std::cerr << e.message << '\n';
        if (!result.errors.empty() || result.models.size() != 1)
            return 1;
        auto model = result.models.front();
        const auto &shape = ModelAccess::get(*model).shape;
        auto report = [&](const char *label, const Diagnostics &d) {
            std::cout << label << std::setprecision(17) << " valid=" << d.valid
                      << " faces=" << d.faces << " solids=" << d.solids << " volume=" << d.volume
                      << " fallback=" << d.fallback << " timeLimited=" << d.time_limited << '\n';
            for (const auto &v : d.surface_counts)
                std::cout << label << " surface " << v.first << ' ' << v.second << '\n';
        };
        report("BREP", model->diagnostics());
        ShapeMap edges;
        TopExp::MapShapes(shape, TopAbs_EDGE, edges);
        std::map<int, unsigned> curves;
        for (int i = 1; i <= edges.Extent(); ++i)
            ++curves[int(BRepAdaptor_Curve(TopoDS::Edge(edges(i))).GetType())];
        for (const auto &v : curves)
            std::cout << "BREP curve " << v.first << ' ' << v.second << '\n';
        auto brep = std::filesystem::path(argv[2]);
        brep.replace_extension(".brep");
        BRepTools::Write(shape, brep.string().c_str());
        auto step = model->export_step();
        std::ofstream out(argv[2], std::ios::binary);
        out.write(reinterpret_cast<const char *>(step.data.data()),
                  std::streamsize(step.data.size()));
        report("STEP", step.diagnostics);
        for (const auto &f : step.diagnostics.fallbacks)
            std::cout << "FALLBACK " << f << '\n';
        bool rims_ok = true;
        if (argc == 4) {
            rims_ok = check_rims(shape, argv[3], "BREP");
            auto imported = import_step({"export.step", step.data});
            rims_ok = check_rims(ModelAccess::get(*imported).shape, argv[3], "STEP") && rims_ok;
        }
        return out && step.diagnostics.valid && rims_ok ? 0 : 1;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

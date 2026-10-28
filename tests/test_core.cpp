#include "internal.hpp"
#include "projection.hpp"
#include <BRepGProp.hxx>
#include <Eigen/Geometry>
#include <GProp_GProps.hxx>
#include <Geom_SphericalSurface.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <STEPControl_Reader.hxx>
#include <atomic>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <gp_Ax3.hxx>
#include <iostream>
#include <mutex>
#include <set>
#include <sstream>
#include <thread>
#include <zip.h>

using namespace stl2curves;
using namespace stl2curves::detail;
namespace {
void check(bool value, const std::string &message) {
    if (!value)
        throw std::runtime_error(message);
}
Mesh cube() {
    Mesh m;
    m.points = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    m.triangles = {{0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7}, {0, 1, 5}, {0, 5, 4},
                   {1, 2, 6}, {1, 6, 5}, {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
    return m;
}
Bytes ascii(const Mesh &m) {
    std::ostringstream s;
    s.precision(17);
    s << "solid test\n";
    for (auto t : m.triangles) {
        s << "facet normal 0 0 0\nouter loop\n";
        for (auto v : t)
            s << "vertex " << m.points[v].x() << ' ' << m.points[v].y() << ' ' << m.points[v].z()
              << '\n';
        s << "endloop\nendfacet\n";
    }
    s << "endsolid test\n";
    auto t = s.str();
    return {t.begin(), t.end()};
}
Mesh short_cylinder() {
    Mesh m;
    constexpr unsigned n = 32;
    constexpr double pi = 3.14159265358979323846;
    for (unsigned layer = 0; layer < 2; ++layer)
        for (unsigned i = 0; i < n; ++i)
            m.points.emplace_back(std::cos(2 * pi * i / n), std::sin(2 * pi * i / n), .2 * layer);
    m.points.emplace_back(0, 0, 0);
    m.points.emplace_back(0, 0, .2);
    for (unsigned i = 0; i < n; ++i) {
        unsigned j = (i + 1) % n;
        m.triangles.push_back({i, j, j + n});
        m.triangles.push_back({i, j + n, i + n});
        m.triangles.push_back({2 * n, j, i});
        m.triangles.push_back({2 * n + 1, i + n, j + n});
    }
    return m;
}
Bytes archive(const std::map<std::string, std::string> &entries) {
    static unsigned serial = 0;
    auto file = std::filesystem::temp_directory_path() /
                ("stl2curves-test-" +
                 std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + "-" +
                 std::to_string(serial++) + ".3mf");
    int error;
    auto *z = zip_open(file.string().c_str(), ZIP_CREATE | ZIP_TRUNCATE, &error);
    check(z, "zip_open");
    for (const auto &kv : entries) {
        auto *source = zip_source_buffer(z, kv.second.data(), kv.second.size(), 0);
        if (!source || zip_file_add(z, kv.first.c_str(), source, ZIP_FL_OVERWRITE) < 0) {
            if (source)
                zip_source_free(source);
            zip_discard(z);
            throw std::runtime_error("zip_file_add");
        }
    }
    check(zip_close(z) == 0, "zip_close");
    std::ifstream in(file, std::ios::binary);
    Bytes bytes((std::istreambuf_iterator<char>(in)), {});
    in.close();
    std::filesystem::remove(file);
    return bytes;
}
std::string cube_xml() {
    auto m = cube();
    std::ostringstream s;
    s << "<mesh><vertices>";
    for (auto v : m.points)
        s << "<vertex x='" << v.x() << "' y='" << v.y() << "' z='" << v.z() << "'/>";
    s << "</vertices><triangles>";
    for (auto t : m.triangles)
        s << "<triangle v1='" << t[0] << "' v2='" << t[1] << "' v3='" << t[2] << "'/>";
    s << "</triangles></mesh>";
    return s.str();
}
void test_mesh() {
    auto m = read_stl(ascii(cube()));
    check(m.points.size() == 8 && m.triangles.size() == 12, "STL welding");
    check(std::abs(signed_volume(m) - 1) < 1e-12, "mesh volume");
    Context c;
    auto damaged = m;
    damaged.triangles.push_back(damaged.triangles[0]);
    std::swap(damaged.triangles[2][0], damaged.triangles[2][1]);
    repair_mesh(damaged, c);
    check(damaged.triangles.size() == 12 && c.diagnostics.open_edges == 0,
          "repair duplicates/orientation");
    check(!c.diagnostics.repairs.empty(), "repair diagnostics");
    auto parts = split_bodies(m);
    check(parts.size() == 1, "connected body");
    check(symmetry_planes(m).size() >= 3, "mesh reflection symmetry");
    auto open = m;
    open.triangles.erase(open.triangles.begin() + 2, open.triangles.begin() + 4);
    Context holes;
    repair_mesh(open, holes);
    check(open.triangles.size() == 12 && holes.diagnostics.open_edges == 0, "planar hole repair");
    Bytes binary(84 + 12 * 50);
    binary[80] = 12;
    auto store = [&](unsigned offset, double value) {
        float f = float(value);
        std::uint32_t bits;
        std::memcpy(&bits, &f, 4);
        for (int k = 0; k < 4; ++k)
            binary[offset + k] = std::uint8_t(bits >> (8 * k));
    };
    for (unsigned i = 0; i < 12; ++i)
        for (unsigned k = 0; k < 3; ++k)
            for (unsigned d = 0; d < 3; ++d)
                store(84 + i * 50 + 12 + k * 12 + d * 4, cube().points[cube().triangles[i][k]][d]);
    check(read_stl(binary).triangles.size() == 12, "binary STL");
}
void test_numeric() {
    Eigen::VectorXd x(2);
    x << 0, 0;
    auto fit = least_squares(
        [](const Eigen::VectorXd &x) -> Eigen::VectorXd {
            Eigen::VectorXd r(3);
            r << x[0] - 2, x[1] + 3, x[0] + x[1] + 1;
            return r;
        },
        x);
    check((fit.x - Eigen::Vector2d(2, -3)).norm() < 1e-5, "robust nonlinear fit");
    std::atomic<unsigned> sum{0};
    std::mutex mutex;
    std::set<std::thread::id> ids;
    auto parallel = parallel_for(100, [&](std::size_t i) {
        sum += unsigned(i);
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        std::lock_guard<std::mutex> lock(mutex);
        ids.insert(std::this_thread::get_id());
    });
    check(sum == 4950, "parallel task accounting");
    check(worker_count() > 1 ? (parallel > 0 && ids.size() > 1) : parallel == 0,
          "actual parallel worker execution");
    bool caught = false;
    try {
        parallel_for(20, [](std::size_t i) {
            if (i == 3)
                throw std::runtime_error("expected");
        });
    } catch (const std::runtime_error &) {
        caught = true;
    }
    check(caught, "worker exception propagation");
}
void test_projection() {
    Vec normal = Vec(1, 2, 3).normalized(), e = normal.unitOrthogonal(), f = normal.cross(e),
        offset(20, -17, 8);
    PointGrid poles(1, 4, 1, 4);
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j) {
            double x = 2. * i / 3, y = 3. * j / 3;
            Vec p = offset + x * e + y * f + (.1 * x * x + .05 * x * y + .03 * y * y) * normal;
            poles.SetValue(i + 1, j + 1, gp_Pnt(p.x(), p.y(), p.z()));
        }
    RealArray knots(1, 2);
    knots.SetValue(1, 0);
    knots.SetValue(2, 1);
    IntArray multiplicities(1, 2);
    multiplicities.Init(4);
    Handle(Geom_Surface) surface =
        new Geom_BSplineSurface(poles, knots, knots, multiplicities, multiplicities, 3, 3);
    PointProjector projector(surface);
    for (double distance : {.01, -.02}) {
        gp_Pnt at;
        gp_Vec du, dv;
        surface->D1(.37, .63, at, du, dv);
        auto n = du.Crossed(dv).Normalized();
        Vec q(at.X() + distance * n.X(), at.Y() + distance * n.Y(), at.Z() + distance * n.Z());
        const auto &projection = projector.perform(q);
        double u, v;
        projection.LowerDistanceParameters(u, v);
        check(projection.NbPoints() > 0 &&
                  std::abs(projection.LowerDistance() - std::abs(distance)) < 1e-8 &&
                  std::abs(u - .37) < 1e-7 && std::abs(v - .63) < 1e-7,
              "rotated height-field normal projection");
    }
}
void test_surfaces() {
    auto m = short_cylinder();
    Context ctx;
    auto patches = detect(m, ctx);
    unsigned cylinders = 0;
    for (const auto &p : patches)
        if (p.kind == "cylinder") {
            ++cylinders;
            check(std::abs(p.radius - 1) < 1e-6, "cylinder radius");
        }
    check(cylinders == 1, "two-ring cylinder must not be mistaken for a sphere");
    Context fitted;
    Mesh freeform;
    constexpr unsigned n = 7;
    for (unsigned i = 0; i < n; ++i)
        for (unsigned j = 0; j < n; ++j) {
            double x = double(i) / (n - 1), y = double(j) / (n - 1);
            freeform.points.emplace_back(x, y, .1 * x * y + .02 * x * x * y);
        }
    for (unsigned i = 0; i + 1 < n; ++i)
        for (unsigned j = 0; j + 1 < n; ++j) {
            unsigned a = i * n + j, b = (i + 1) * n + j;
            freeform.triangles.push_back({a, b, b + 1});
            freeform.triangles.push_back({a, b + 1, a + 1});
        }
    auto top = topology(freeform);
    std::vector<int> owner(freeform.triangles.size(), -1);
    fitted.deadline = Clock::time_point::max();
    auto smooth = advanced_patches(freeform, top, owner, {}, fitted);
    check(!smooth.empty(), "freeform fitting");
    for (const auto &p : smooth)
        check(!p.surface.IsNull() && p.deviation <= .002, "freeform fitting tolerance");
}
void test_closed_surfaces() {
    auto mesh = cube();
    for (auto &v : mesh.points)
        v *= 10;
    Context ctx;
    auto patches = detect(mesh, ctx);
    auto append = [&](Vec center, double radius, bool inward) {
        unsigned offset = mesh.points.size();
        for (const Vec v : std::vector<Vec>{Vec::UnitX(), -Vec::UnitX(), Vec::UnitY(),
                                            -Vec::UnitY(), Vec::UnitZ(), -Vec::UnitZ()})
            mesh.points.push_back(center + radius * v);
        Patch p;
        p.kind = "sphere";
        p.origin = center;
        p.radius = radius;
        p.surface = new Geom_SphericalSurface(
            gp_Ax3(gp_Pnt(center.x(), center.y(), center.z()), gp_Dir(0, 0, 1)), radius);
        for (auto t : std::vector<Tri>{{0, 2, 4},
                                       {2, 1, 4},
                                       {1, 3, 4},
                                       {3, 0, 4},
                                       {2, 0, 5},
                                       {1, 2, 5},
                                       {3, 1, 5},
                                       {0, 3, 5}}) {
            for (auto &v : t)
                v += offset;
            if (inward)
                std::swap(t[1], t[2]);
            p.triangles.push_back(mesh.triangles.size());
            mesh.triangles.push_back(t);
        }
        patches.push_back(p);
    };
    append(Vec(5, 5, 5), 1, true);
    append(Vec(20, 5, 5), 2, false);
    auto shape = reconstruct(mesh, patches, ctx);
    Diagnostics info;
    check(inspect(shape, info) && info.solids == 2 && info.faces == 8 &&
              std::abs(info.volume - (1000 + 4 * 3.14159265358979323846 * 7 / 3)) < 1e-5,
          "closed single-face sphere bodies and cavities");
}
void test_3mf() {
    const std::string prefix =
        "<model xmlns='http://schemas.microsoft.com/3dmanufacturing/core/2015/02' "
        "unit='inch'><resources>";
    auto xml =
        prefix + "<object id='1' name='cube'>" + cube_xml() +
        "</object><object id='2'><components><component objectid='1' transform='1 0 0 0 1 0 "
        "0 0 1 10 20 30'/></components></object></resources><build><item "
        "objectid='2' transform='1 0 0 0 1 0 0 0 1 2 3 4'/><item objectid='2'/></build></model>";
    auto objects = read_3mf(archive({{"3D/3dmodel.model", xml}}));
    check(objects.size() == 2, "3MF instances");
    check(objects[0].name != objects[1].name, "3MF duplicate names");
    check(std::abs(signed_volume(objects[0].mesh) - std::pow(25.4, 3)) < 1e-7, "3MF units");
    Vec minimum = objects[0].mesh.points.front();
    for (auto v : objects[0].mesh.points)
        minimum = minimum.cwiseMin(v);
    check((minimum - 25.4 * Vec(12, 23, 34)).norm() < 1e-9,
          "3MF component and build translation units");
    minimum = objects[1].mesh.points.front();
    for (auto v : objects[1].mesh.points)
        minimum = minimum.cwiseMin(v);
    check((minimum - 25.4 * Vec(10, 20, 30)).norm() < 1e-9, "3MF component translation units");
    auto reflected_xml = std::string("<model unit='centimeter'><resources><object id='1'>") +
                         cube_xml() +
                         "</object></resources><build><item objectid='1' transform='-1 0 0 0 1 0 0 "
                         "0 1 3 0 0'/></build></model>";
    auto reflected = read_3mf(archive({{"3D/3dmodel.model", reflected_xml}}));
    check(reflected.size() == 1 && std::abs(signed_volume(reflected[0].mesh) - 1000) < 1e-9,
          "3MF reflection preserves positive volume");
    minimum = reflected[0].mesh.points.front();
    for (auto v : reflected[0].mesh.points)
        minimum = minimum.cwiseMin(v);
    check((minimum - Vec(20, 0, 0)).norm() < 1e-9, "3MF reflection and translation coordinates");
    auto external = "<model xmlns='http://schemas.microsoft.com/3dmanufacturing/core/2015/02' "
                    "xmlns:p='http://schemas.microsoft.com/3dmanufacturing/production/2015/"
                    "06'><resources><object id='2'><components><component objectid='1' "
                    "p:path='/3D/sub.model'/><component objectid='3' "
                    "p:path='/3D/sub.model'/></components></object></resources><build><item "
                    "objectid='2'/></build></model>";
    auto sub = "<model><resources><object id='1'>" + cube_xml() + "</object><object id='3'>" +
               cube_xml() + "</object></resources></model>";
    auto settings = "<config><object id='2'><metadata key='name' value='part'/><part id='3' "
                    "subtype='negative_part'/></object></config>";
    auto selected = read_3mf(archive({{"3D/3dmodel.model", external},
                                      {"3D/sub.model", sub},
                                      {"Metadata/model_settings.config", settings}}));
    check(selected[0].name == "part" && selected[0].mesh.triangles.size() == 12,
          "external components/helper filtering");
    bool caught = false;
    try {
        read_3mf(archive(
            {{"3D/3dmodel.model", "<model><resources><object id='1'><components><component "
                                  "objectid='1'/></components></object></resources><build><item "
                                  "objectid='1'/></build></model>"}}));
    } catch (const std::runtime_error &) {
        caught = true;
    }
    check(caught, "component cycle");
}
void test_conversion() {
    auto open = cube();
    open.triangles.erase(open.triangles.begin() + 2, open.triangles.begin() + 4);
    Options no_repair;
    no_repair.repair = false;
    no_repair.curves = false;
    auto failed = convert({{"open.stl", ascii(open)}}, no_repair);
    check(failed.models.empty() && failed.errors.size() == 1 && failed.errors[0].diagnostics &&
              failed.errors[0].diagnostics->open_edges == 4,
          "failed-object diagnostics retain unresolved mesh defects");
    Options o;
    o.details = true;
    auto bytes = ascii(cube());
    auto r = convert({{"cube.stl", bytes}}, o);
    if (!r.errors.empty())
        throw std::runtime_error(r.errors[0].message);
    check(r.models.size() == 1, "cube output");
    auto d = r.models[0]->diagnostics();
    check(d.valid && d.faces == 6 && d.solids == 1 && std::abs(d.volume - 1) < 1e-8,
          "cube validity and planar consolidation");
    auto exported = r.models[0]->export_step();
    auto retained = exported.data;
    o.curves = false;
    o.blends = false;
    o.repair = false;
    o.simplify = 0;
    o.merge = "assembly";
    auto repeated = convert({{"a.stl", bytes}, {"a.stl", bytes}, {"invalid.stl", {1, 2, 3}}}, o);
    check(repeated.models.size() == 3 && repeated.errors.size() == 1,
          "batch partial success and merge");
    check(repeated.models[0]->name() != repeated.models[1]->name(), "unique names");
    check(exported.data == retained, "retained buffer ownership");
    bool caught = false;
    o.accuracy = 0;
    try {
        convert({}, o);
    } catch (const std::invalid_argument &) {
        caught = true;
    }
    check(caught, "invalid options");
    auto container = cube();
    for (auto &p : container.points)
        p *= 3;
    auto inner = cube();
    unsigned offset = container.points.size();
    for (auto &p : inner.points)
        container.points.push_back(p + Vec(1, 1, 1));
    for (auto t : inner.triangles) {
        for (auto &v : t)
            v += offset;
        std::swap(t[1], t[2]);
        container.triangles.push_back(t);
    }
    auto cavity = convert({{"cavity.stl", ascii(container)}}, Options{});
    check(cavity.errors.empty() && cavity.models.size() == 1 &&
              std::abs(cavity.models[0]->diagnostics().volume - 26) < 1e-6,
          "nested cavity shells");
    auto project =
        "<model><resources><object id='1'>" + cube_xml() +
        "</object><object id='2'><components><component objectid='1'/><component objectid='1' "
        "transform='1 0 0 0 1 0 0 0 1 .5 0 0'/></components></object></resources><build><item "
        "objectid='2'/></build></model>";
    auto overlapping =
        convert({{"overlap.3mf", archive({{"3D/3dmodel.model", project}})}}, Options{});
    check(overlapping.errors.empty() && overlapping.models.size() == 1 &&
              std::abs(overlapping.models[0]->diagnostics().volume - 1.5) < 1e-6,
          "overlapping 3MF parts must be repaired separately and fused");
    auto scaled = cube();
    for (auto &point : scaled.points)
        point = point.cwiseProduct(Vec(.99, 1.98, 2.97));
    Options sizing;
    sizing.true_size = true;
    auto designed = convert({{"scaled.stl", ascii(scaled)}}, sizing);
    check(designed.errors.empty() && designed.models.size() == 1 &&
              std::abs(designed.models[0]->diagnostics().volume - 6) < 1e-6,
          "common design-size scale");
    check(std::abs(designed.models[0]->diagnostics().scale - 100. / 99) < 1e-8,
          "99 percent scale hypothesis");
    const auto scaled_parts =
        "<model><resources><object id='1'>" + cube_xml() +
        "</object><object id='2'><components>"
        "<component objectid='1' transform='.99 0 0 0 1.98 0 0 0 2.97 0 0 0'/>"
        "<component objectid='1' transform='.99 0 0 0 1.98 0 0 0 2.97 5.94 0 0'/>"
        "</components></object></resources><build><item objectid='2'/></build></model>";
    sizing.fuse = false;
    auto assembly =
        convert({{"scaled-parts.3mf", archive({{"3D/3dmodel.model", scaled_parts}})}}, sizing);
    check(assembly.errors.empty() && assembly.models.size() == 1 &&
              assembly.models[0]->diagnostics().solids == 2 &&
              std::abs(assembly.models[0]->diagnostics().volume - 12) < 1e-6 &&
              std::abs(assembly.models[0]->diagnostics().scale - 100. / 99) < 1e-8,
          "multipart design sizing preserves part vertex ranges");
    STEPControl_Reader assembly_reader;
    auto assembly_step = assembly.models[0]->export_step();
    const auto &assembly_bytes = assembly_step.data;
    std::istringstream assembly_stream(std::string(assembly_bytes.begin(), assembly_bytes.end()));
    check(assembly_reader.ReadStream("scaled-parts.step", assembly_stream) == IFSelect_RetDone &&
              assembly_reader.TransferRoots() > 0,
          "multipart sizing STEP re-import");
    GProp_GProps assembly_mass;
    BRepGProp::VolumeProperties(assembly_reader.OneShape(), assembly_mass);
    check(assembly_mass.CentreOfMass().Distance(gp_Pnt(3.5, 1, 1.5)) < 1e-6,
          "multipart sizing preserves component positions");
    auto metres = cube();
    for (auto &point : metres.points)
        point = point.cwiseProduct(Vec(.123, .234, .345));
    auto units = convert({{"metres.stl", ascii(metres)}}, sizing);
    check(units.errors.empty() && units.models.size() == 1 &&
              units.models[0]->diagnostics().scale == 1000,
          "metre unit hypothesis");
}
} // namespace
int main() {
    try {
        test_mesh();
        test_numeric();
        test_projection();
        test_surfaces();
        test_closed_surfaces();
        test_3mf();
        test_conversion();
        std::cout << "Core tests passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

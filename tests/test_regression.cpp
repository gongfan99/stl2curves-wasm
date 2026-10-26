#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stl2curves/converter.hpp>
using namespace stl2curves;
int main(int argc, char **argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: test_regression FIXTURE_DIRECTORY [FILENAME]\n";
        return 2;
    }
    struct Case {
        const char *name;
        double volume, error;
        unsigned faces;
    };
    const Case cases[] = {
        {"test_fine.stl", 34615.8107, .01, 29},   {"test_coarse.stl", 34615.8107, .01, 29},
        {"test_hard.stl", 35134.2571, .05, 31},   {"test_lid.stl", 57600.4078, .05, 45},
        {"test_sphere.stl", 34303.804, .01, 34},  {"test_corners.stl", 32796.5606, .1, 30},
        {"test_corners2.stl", 14072.0767, .1, 34}};
    unsigned failures = 0, examined = 0;
    Options options;
    options.fuse = false;
    options.details = true;
    for (auto c : cases) {
        if (argc == 3 && std::string(argv[2]) != c.name)
            continue;
        ++examined;
        std::cout << "Converting " << c.name << std::endl;
        std::ifstream in(std::filesystem::path(argv[1]) / c.name, std::ios::binary);
        if (!in) {
            std::cerr << "Missing fixture " << c.name << '\n';
            ++failures;
            continue;
        }
        Bytes bytes((std::istreambuf_iterator<char>(in)), {});
        auto result = convert({{c.name, std::move(bytes)}}, options);
        if (result.files.size() != 1 || !result.errors.empty()) {
            std::cerr << "FAIL " << c.name;
            for (const auto &e : result.errors)
                std::cerr << " " << e.message;
            std::cerr << '\n';
            ++failures;
            continue;
        }
        const auto &d = result.files[0].diagnostics;
        bool ok = d.valid && std::abs(d.volume - c.volume) <= c.error && d.faces <= c.faces;
        std::cout << (ok ? "PASS " : "FAIL ") << c.name << " volume error " << d.volume - c.volume
                  << " faces " << d.faces << " max " << c.faces << " fallback " << d.fallback
                  << '\n';
        failures += !ok;
        if (!ok) {
            for (const auto &reason : d.fallbacks)
                std::cout << "  fallback: " << reason << '\n';
            for (const auto &count : d.surface_counts)
                std::cout << "  " << count.first << ": " << count.second << '\n';
            for (const auto &feature : d.features)
                if (feature.kind != "plane")
                    std::cout << "  recognized " << feature.kind << " triangles "
                              << feature.triangles << " radius " << feature.radius << '\n';
        }
    }
    return examined == 0 ? 2 : (failures ? 1 : 0);
}

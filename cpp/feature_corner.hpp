// Derived from stl2curves/features.py, Copyright (c) 2026 Duvdivan. MIT license.
#pragma once
#include "feature_geometry.hpp"
#include "reference_least_squares.hpp"
namespace stl2curves::detail::reference_features {
inline Found corner_surface(const Pass &pass, const std::vector<unsigned> &facets,
                            double accuracy) {
    if (facets.size() < 8)
        return {};
    double least = 1;
    Eigen::Matrix3d normal_matrix = Eigen::Matrix3d::Zero();
    Vec outward = Vec::Zero();
    for (auto a : facets) {
        const auto &normal = pass.fs.normals[a];
        normal_matrix += normal * normal.transpose();
        outward += normal;
        for (auto b : facets)
            least = std::min(least, normal.dot(pass.fs.normals[b]));
    }
    if (std::acos(std::clamp(least, -1., 1.)) < 15 * pi / 180)
        return {};
    outward /= double(facets.size());
    outward = divided(outward, std::max(norm3(outward), 1e-12));
    auto P = points(pass.m, pass.fs, facets);
    Vec center = Vec::Zero();
    for (const auto &p : P)
        center += p;
    center = divided(center, double(P.size()));
    Eigen::MatrixXd centered(P.size(), 3);
    for (unsigned i = 0; i < P.size(); ++i)
        centered.row(i) = P[i] - center;
    std::array<Vec, 2> guesses{
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(normal_matrix).eigenvectors().col(0),
        centered.jacobiSvd(Eigen::ComputeThinV).matrixV().col(0)};
    double best_error = std::numeric_limits<double>::infinity();
    Model best;
    for (const auto &guess : guesses) {
        Model frame;
        frame.frame(guess);
        auto residual = [&](const Eigen::VectorXd &x) -> Eigen::VectorXd {
            Vec axis = normalized(frame.d + x[0] * frame.e1 + x[1] * frame.e2);
            Vec origin = center + x[2] * frame.e1 + x[3] * frame.e2;
            Eigen::VectorXd result(P.size());
            for (unsigned i = 0; i < P.size(); ++i) {
                Vec q = P[i] - origin;
                q -= q.dot(axis) * axis;
                result[i] = norm3(q) - x[4];
            }
            return result;
        };
        for (double radius : {5., 20., 80.})
            for (double sign : {1., -1.}) {
                Vec off = sign * radius * outward;
                Eigen::VectorXd initial(5);
                initial << 0, 0, off.dot(frame.e1), off.dot(frame.e2), radius;
                auto fit = trust_region_fit(residual, initial);
                double error = residual(fit.x).cwiseAbs().maxCoeff();
                if (error < best_error) {
                    best_error = error;
                    best.kind = "cylinder";
                    best.frame(frame.d + fit.x[0] * frame.e1 + fit.x[1] * frame.e2);
                    Vec origin = center + fit.x[2] * frame.e1 + fit.x[3] * frame.e2;
                    best.a = origin - origin.dot(best.d) * best.d;
                    best.r = fit.x[4];
                }
            }
    }
    if (!best || best_error > .005 * accuracy / .02 + pass.fs.noise ||
        !(best.r > 0 && best.r < 1000))
        return {};
    bool convex =
        pass.fs.normals[facets.front()].dot(best.normal(pass.fs.centers[facets.front()])) > 0;
    for (auto id : facets)
        if (std::abs(best.signed_at(pass.fs.centers[id])) > .25 ||
            (pass.fs.normals[id].dot(best.normal(pass.fs.centers[id])) > 0) != convex)
            return {};
    auto feature = pass.describe(best, facets, true);
    if (feature)
        feature.convex = convex;
    return feature;
}
} // namespace stl2curves::detail::reference_features

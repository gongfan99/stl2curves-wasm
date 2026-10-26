#pragma once
#include "internal.hpp"
#include <Eigen/QR>
#include <Extrema_ExtAlgo.hxx>
#include <Extrema_ExtFlag.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_BSplineSurface.hxx>

namespace stl2curves::detail {
class SurfaceProjection {
    friend class PointProjector;
    const GeomAPI_ProjectPointOnSurf *general_ = nullptr;
    gp_Pnt point_;
    double u_ = 0, v_ = 0, distance_ = 0;

  public:
    int NbPoints() const { return general_ ? general_->NbPoints() : 1; }
    double LowerDistance() const { return general_ ? general_->LowerDistance() : distance_; }
    gp_Pnt NearestPoint() const { return general_ ? general_->NearestPoint() : point_; }
    void LowerDistanceParameters(double &u, double &v) const {
        if (general_)
            general_->LowerDistanceParameters(u, v);
        else {
            u = u_;
            v = v_;
        }
    }
};

// Reuse search structures across a patch. Height fields have an exact affine
// inverse for their two plane coordinates, providing a reliable Newton seed.
// Other bounded surfaces retain the global tree search for multiple extrema.
class PointProjector {
    Handle(Geom_Surface) surface_;
    GeomAPI_ProjectPointOnSurf projection_;
    SurfaceProjection result_;
    bool bounded_ = false, affine_ = false;
    double a_, b_, c_, d_;
    Vec center_ = Vec::Zero();
    Eigen::Matrix<double, 4, 2> inverse_;

    bool affine_inverse() {
        auto spline = Handle(Geom_BSplineSurface)::DownCast(surface_);
        if (spline.IsNull() || spline->IsUPeriodic() || spline->IsVPeriodic() ||
            spline->IsURational() || spline->IsVRational())
            return false;
        unsigned nu = spline->NbUPoles(), nv = spline->NbVPoles();
        Eigen::MatrixXd positions(nu * nv, 4), parameters(nu * nv, 2);
        for (unsigned i = 1; i <= nu; ++i)
            for (unsigned j = 1; j <= nv; ++j) {
                auto p = spline->Pole(i, j);
                center_ += Vec(p.X(), p.Y(), p.Z());
            }
        center_ /= double(nu * nv);
        const auto &uk = spline->UKnotSequence();
        const auto &vk = spline->VKnotSequence();
        for (unsigned i = 1; i <= nu; ++i)
            for (unsigned j = 1; j <= nv; ++j) {
                auto p = spline->Pole(i, j);
                unsigned row = (i - 1) * nv + j - 1;
                Vec q = Vec(p.X(), p.Y(), p.Z()) - center_;
                positions.row(row) << q.x(), q.y(), q.z(), 1;
                double u = 0, v = 0;
                for (int k = 1; k <= spline->UDegree(); ++k)
                    u += uk.Value(i + k);
                for (int k = 1; k <= spline->VDegree(); ++k)
                    v += vk.Value(j + k);
                parameters.row(row) << u / spline->UDegree(), v / spline->VDegree();
            }
        auto solve = positions.colPivHouseholderQr();
        if (solve.rank() != 4)
            return false;
        inverse_ = solve.solve(parameters);
        return inverse_.allFinite() && (positions * inverse_ - parameters).cwiseAbs().maxCoeff() <
                                           1e-9 * std::max({1., b_ - a_, d_ - c_});
    }

  public:
    explicit PointProjector(const Handle(Geom_Surface) & surface) : surface_(surface) {
        surface->Bounds(a_, b_, c_, d_);
        bounded_ = std::isfinite(a_) && std::isfinite(b_) && std::isfinite(c_) &&
                   std::isfinite(d_) &&
                   std::max({std::abs(a_), std::abs(b_), std::abs(c_), std::abs(d_)}) < 1e12;
        affine_ = bounded_ && affine_inverse();
        if (bounded_ && !affine_) {
            projection_.Init(surface, a_, b_, c_, d_, 1e-8, Extrema_ExtAlgo_Tree);
            projection_.SetExtremaFlag(Extrema_ExtFlag_MIN);
        }
    }
    const SurfaceProjection &perform(const Vec &p) {
        gp_Pnt point(p.x(), p.y(), p.z());
        if (affine_) {
            Vec q = p - center_;
            Eigen::Vector4d position(q.x(), q.y(), q.z(), 1);
            Eigen::Vector2d uv = inverse_.transpose() * position;
            double u = std::clamp(uv[0], a_, b_), v = std::clamp(uv[1], c_, d_);
            for (unsigned iteration = 0; iteration < 8; ++iteration) {
                gp_Pnt at;
                gp_Vec du, dv;
                surface_->D1(u, v, at, du, dv);
                gp_Vec r(point, at);
                double aa = du.SquareMagnitude(), bb = du.Dot(dv), cc = dv.SquareMagnitude(),
                       det = aa * cc - bb * bb;
                if (det <= 1e-24)
                    break;
                double x = (cc * du.Dot(r) - bb * dv.Dot(r)) / det,
                       y = (aa * dv.Dot(r) - bb * du.Dot(r)) / det;
                double nextu = std::clamp(u - x, a_, b_), nextv = std::clamp(v - y, c_, d_);
                bool done = std::abs(nextu - u) + std::abs(nextv - v) < 1e-11;
                u = nextu;
                v = nextv;
                if (done)
                    break;
            }
            result_.general_ = nullptr;
            result_.u_ = u;
            result_.v_ = v;
            result_.point_ = surface_->Value(u, v);
            result_.distance_ = result_.point_.Distance(point);
        } else {
            if (bounded_)
                projection_.Perform(point);
            else
                projection_.Init(point, surface_);
            result_.general_ = &projection_;
        }
        return result_;
    }
};
} // namespace stl2curves::detail

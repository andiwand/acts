// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Helix transport, Jacobians, and RZ intersections in an axial field.

#include "Acts/TrackFinding/Rz/RzTypes.hpp"

#include <cmath>
#include <optional>
#include <span>

namespace Acts::Experimental {

namespace detail {

/// Trigonometric terms for a helix step and its Jacobian.
struct StepTrig {
  double sn{};
  double cs{};
  /// `sin(u) / u`
  double sinc{};
  /// `(1 - cos(u)) / u`
  double versc{};
  /// `d/du (sin(u) / u)`
  double dsinc{};
  /// `d/du ((1 - cos(u)) / u)`
  double dversc{};
};

inline StepTrig stepTrig(double u) {
  StepTrig t;
  const double u2 = u * u;
  // The direct derivative formulas cancel near zero. At this threshold the
  // omitted series terms are at double-precision rounding scale.
  if (std::abs(u) < 1e-2) {
    t.sinc = 1. + u2 * (-1. / 6. + u2 / 120.);
    t.versc = u * (0.5 + u2 * (-1. / 24. + u2 / 720.));
    t.dsinc = u * (-1. / 3. + u2 * (1. / 30. - u2 / 840.));
    t.dversc = 0.5 + u2 * (-1. / 8. + u2 * (1. / 144. - u2 / 5760.));
    t.sn = u * t.sinc;
    t.cs = 1. - u * t.versc;
    return t;
  }
  t.sn = std::sin(u);
  t.cs = std::cos(u);
  // Avoid cancellation except near cos(u) = -1.
  const double oneMinusCos =
      1. + t.cs > 1e-3 ? t.sn * t.sn / (1. + t.cs) : 1. - t.cs;
  t.sinc = t.sn / u;
  t.versc = oneMinusCos / u;
  t.dsinc = (u * t.cs - t.sn) / u2;
  t.dversc = (u * t.sn - oneMinusCos) / u2;
  return t;
}

}  // namespace detail

/// Helix in a constant axial field, with curvature `q/p * Bz`.
struct RzHelix {
  /// Field along z in native units
  double bz{};

  /// Turning rate of the state
  /// @param v the state
  /// @return `q/p * Bz`, in inverse length
  double kappa(const RzVector& v) const { return v[eRzQOverP] * bz; }

  /// Reverse direction and charge to traverse the same helix backwards.
  /// @param v the state
  /// @return the reversed state
  static RzVector reversed(const RzVector& v) {
    RzVector r = v;
    r.segment<3>(eRzDir0) *= -1.;
    r[eRzQOverP] *= -1.;
    return r;
  }

  /// `d(state)/ds`
  /// @param v the state
  /// @return the derivative
  RzVector derivative(const RzVector& v) const {
    const double k = kappa(v);
    RzVector d;
    d[eRzPos0] = v[eRzDir0];
    d[eRzPos1] = v[eRzDir1];
    d[eRzPos2] = v[eRzDir2];
    d[eRzDir0] = k * v[eRzDir1];
    d[eRzDir1] = -k * v[eRzDir0];
    d[eRzDir2] = 0.;
    d[eRzQOverP] = 0.;
    return d;
  }

  /// Transport the state by a path length, in place.
  /// @param v the state
  /// @param s the path length, may be negative
  void step(RzVector& v, double s) const {
    step(v, s, detail::stepTrig(kappa(v) * s));
  }

  /// Step with cached trigonometry, also reusable for the Jacobian.
  /// @param v the state
  /// @param s the path length
  /// @param t `stepTrig(kappa(v) * s)`
  void step(RzVector& v, double s, const detail::StepTrig& t) const {
    const double f1 = s * t.sinc;
    const double f2 = s * t.versc;
    const double sn = t.sn;
    const double cs = t.cs;
    const double dx = v[eRzDir0];
    const double dy = v[eRzDir1];
    v[eRzPos0] += dx * f1 + dy * f2;
    v[eRzPos1] += -dx * f2 + dy * f1;
    v[eRzPos2] += v[eRzDir2] * s;
    v[eRzDir0] = dx * cs + dy * sn;
    v[eRzDir1] = -dx * sn + dy * cs;
  }

  /// Jacobian of `step` at a fixed path length, evaluated on the state before
  /// the step.
  /// @param v0 the state before the step
  /// @param s the path length
  /// @return the 7x7 Jacobian
  RzMatrix stepJacobian(const RzVector& v0, double s) const {
    const double k = kappa(v0);
    const detail::StepTrig t = detail::stepTrig(k * s);
    const double f1 = s * t.sinc;
    const double f2 = s * t.versc;
    const double g1 = s * s * t.dsinc;
    const double g2 = s * s * t.dversc;
    const double sn = t.sn;
    const double cs = t.cs;
    const double dx = v0[eRzDir0];
    const double dy = v0[eRzDir1];
    const double dxs = dx * cs + dy * sn;
    const double dys = -dx * sn + dy * cs;

    RzMatrix j = RzMatrix::Identity();
    j(eRzPos0, eRzDir0) = f1;
    j(eRzPos0, eRzDir1) = f2;
    j(eRzPos0, eRzQOverP) = bz * (dx * g1 + dy * g2);
    j(eRzPos1, eRzDir0) = -f2;
    j(eRzPos1, eRzDir1) = f1;
    j(eRzPos1, eRzQOverP) = bz * (-dx * g2 + dy * g1);
    j(eRzPos2, eRzDir2) = s;
    j(eRzDir0, eRzDir0) = cs;
    j(eRzDir0, eRzDir1) = sn;
    j(eRzDir0, eRzQOverP) = bz * s * dys;
    j(eRzDir1, eRzDir0) = -sn;
    j(eRzDir1, eRzDir1) = cs;
    j(eRzDir1, eRzQOverP) = -bz * s * dxs;
    return j;
  }

  /// Sparse step Jacobian plus the rank-one path-length correction.
  struct StepJacobian {
    /// The helix entries: position from direction and q/p, direction from
    /// direction and q/p, z from dz
    double f1{};
    double f2{};
    double a1{};
    double a2{};
    double s{};
    double cs{};
    double sn{};
    double b1{};
    double b2{};
    /// z and dz from q/p: zero for the helix, what a field off the z axis
    /// adds
    double a3{};
    double b3{};
    /// The derivative of the end state and `ds/d(start)`; the rank-1 term
    /// is their outer product
    RzVector d{RzVector::Zero()};
    RzVector dsdv{RzVector::Zero()};

    /// `J x` for one column
    /// @param x the column
    /// @param y the result
    void applyLeft(std::span<const double, eRzSize> x,
                   std::span<double, eRzSize> y) const {
      const double w = dsdv[eRzPos0] * x[eRzPos0] + dsdv[eRzPos1] * x[eRzPos1] +
                       dsdv[eRzPos2] * x[eRzPos2] + dsdv[eRzDir0] * x[eRzDir0] +
                       dsdv[eRzDir1] * x[eRzDir1] + dsdv[eRzDir2] * x[eRzDir2] +
                       dsdv[eRzQOverP] * x[eRzQOverP];
      y[eRzPos0] = x[eRzPos0] + f1 * x[eRzDir0] + f2 * x[eRzDir1] +
                   a1 * x[eRzQOverP] + d[eRzPos0] * w;
      y[eRzPos1] = x[eRzPos1] - f2 * x[eRzDir0] + f1 * x[eRzDir1] +
                   a2 * x[eRzQOverP] + d[eRzPos1] * w;
      y[eRzPos2] =
          x[eRzPos2] + s * x[eRzDir2] + a3 * x[eRzQOverP] + d[eRzPos2] * w;
      y[eRzDir0] = cs * x[eRzDir0] + sn * x[eRzDir1] + b1 * x[eRzQOverP] +
                   d[eRzDir0] * w;
      y[eRzDir1] = -sn * x[eRzDir0] + cs * x[eRzDir1] + b2 * x[eRzQOverP] +
                   d[eRzDir1] * w;
      y[eRzDir2] = x[eRzDir2] + b3 * x[eRzQOverP];
      y[eRzQOverP] = x[eRzQOverP];
    }

    /// `J X`
    /// @param x the matrix
    /// @return the product
    RzMatrix applyLeft(const RzMatrix& x) const {
      RzMatrix y;
      for (std::uint32_t c = 0; c < eRzSize; ++c) {
        applyLeft(std::span<const double, eRzSize>{x.col(c).data(), eRzSize},
                  std::span<double, eRzSize>{y.col(c).data(), eRzSize});
      }
      return y;
    }

    /// `J C J^T` for a symmetric `C`
    /// @param c the covariance
    /// @return the transported covariance
    RzMatrix transport(const RzMatrix& c) const {
      // J (J C)^T = J C J^T, symmetric, so the second pass fills the lower
      // triangle only and mirrors it
      const RzMatrix jc = applyLeft(c);
      const RzMatrix x = jc.transpose();
      double w[eRzSize];
      for (std::uint32_t col = 0; col < eRzSize; ++col) {
        w[col] = dsdv.dot(x.col(col));
      }
      RzMatrix y;
      auto row = [&](std::uint32_t r, auto&& entry) {
        for (std::uint32_t col = 0; col <= r; ++col) {
          y(r, col) = entry(col) + d[r] * w[col];
        }
      };
      row(eRzPos0, [&](std::uint32_t col) {
        return x(eRzPos0, col) + f1 * x(eRzDir0, col) + f2 * x(eRzDir1, col) +
               a1 * x(eRzQOverP, col);
      });
      row(eRzPos1, [&](std::uint32_t col) {
        return x(eRzPos1, col) - f2 * x(eRzDir0, col) + f1 * x(eRzDir1, col) +
               a2 * x(eRzQOverP, col);
      });
      row(eRzPos2, [&](std::uint32_t col) {
        return x(eRzPos2, col) + s * x(eRzDir2, col) + a3 * x(eRzQOverP, col);
      });
      row(eRzDir0, [&](std::uint32_t col) {
        return cs * x(eRzDir0, col) + sn * x(eRzDir1, col) +
               b1 * x(eRzQOverP, col);
      });
      row(eRzDir1, [&](std::uint32_t col) {
        return -sn * x(eRzDir0, col) + cs * x(eRzDir1, col) +
               b2 * x(eRzQOverP, col);
      });
      for (std::uint32_t col = 0; col <= eRzDir2; ++col) {
        y(eRzDir2, col) = x(eRzDir2, col) + b3 * x(eRzQOverP, col);
      }
      for (std::uint32_t col = 0; col < eRzSize; ++col) {
        y(eRzQOverP, col) = x(eRzQOverP, col);
      }
      for (std::uint32_t r = 0; r < eRzSize; ++r) {
        for (std::uint32_t col = r + 1; col < eRzSize; ++col) {
          y(r, col) = y(col, r);
        }
      }
      return y;
    }

    /// The position rows of the Jacobian
    /// @return the 3x7 block
    Eigen::Matrix<double, 3, eRzSize> positionRows() const {
      Eigen::Matrix<double, 3, eRzSize> r =
          Eigen::Matrix<double, 3, eRzSize>::Zero();
      r(0, eRzPos0) = 1.;
      r(0, eRzDir0) = f1;
      r(0, eRzDir1) = f2;
      r(0, eRzQOverP) = a1;
      r(1, eRzPos1) = 1.;
      r(1, eRzDir0) = -f2;
      r(1, eRzDir1) = f1;
      r(1, eRzQOverP) = a2;
      r(2, eRzPos2) = 1.;
      r(2, eRzDir2) = s;
      r(2, eRzQOverP) = a3;
      for (std::uint32_t i = 0; i < 3; ++i) {
        r.row(i) += d[eRzPos0 + i] * dsdv.transpose();
      }
      return r;
    }

    /// The whole matrix, for checks
    RzMatrix dense() const {
      RzMatrix j = RzMatrix::Identity();
      j.block<3, eRzSize>(eRzPos0, 0) = positionRows();
      j(eRzDir0, eRzDir0) = cs;
      j(eRzDir0, eRzDir1) = sn;
      j(eRzDir0, eRzQOverP) = b1;
      j(eRzDir1, eRzDir0) = -sn;
      j(eRzDir1, eRzDir1) = cs;
      j(eRzDir1, eRzQOverP) = b2;
      j(eRzDir2, eRzQOverP) = b3;
      j.row(eRzDir0) += d[eRzDir0] * dsdv.transpose();
      j.row(eRzDir1) += d[eRzDir1] * dsdv.transpose();
      return j;
    }
  };

  /// The Jacobian of a step of length `s` from `v0` that ends at `end` on a
  /// surface with the given normal, in the form above
  /// @param v0 the state before the step
  /// @param s the path length
  /// @param end the state after the step
  /// @param normal the surface normal at the end
  /// @return the Jacobian
  StepJacobian stepJacobianOnto(const RzVector& v0, double s,
                                const RzVector& end,
                                const Vector3& normal) const {
    return stepJacobianOnto(v0, s, end, normal,
                            detail::stepTrig(kappa(v0) * s));
  }

  /// The same with the trigonometry of the step already taken
  /// @param v0 the state before the step
  /// @param s the path length
  /// @param end the state after the step
  /// @param normal the surface normal at the end
  /// @param t `stepTrig(kappa(v0) * s)`
  /// @return the Jacobian
  /// @param qopPosition added to the position rows' q/p column, for what
  ///        the helix does not model
  /// @param qopDirection the same for the direction rows
  StepJacobian stepJacobianOnto(
      const RzVector& v0, double s, const RzVector& end, const Vector3& normal,
      const detail::StepTrig& t, const Vector3& qopPosition = Vector3::Zero(),
      const Vector3& qopDirection = Vector3::Zero()) const {
    const double dx = v0[eRzDir0];
    const double dy = v0[eRzDir1];
    const double g1 = s * s * t.dsinc;
    const double g2 = s * s * t.dversc;
    const double dxs = dx * t.cs + dy * t.sn;
    const double dys = -dx * t.sn + dy * t.cs;
    StepJacobian j;
    j.f1 = s * t.sinc;
    j.f2 = s * t.versc;
    j.a1 = bz * (dx * g1 + dy * g2);
    j.a2 = bz * (-dx * g2 + dy * g1);
    j.s = s;
    j.cs = t.cs;
    j.sn = t.sn;
    j.b1 = bz * s * dys;
    j.b2 = -bz * s * dxs;
    j.a1 += qopPosition.x();
    j.a2 += qopPosition.y();
    j.a3 = qopPosition.z();
    j.b1 += qopDirection.x();
    j.b2 += qopDirection.y();
    j.b3 = qopDirection.z();
    j.d = derivative(end);
    // ds/d(start) = -(n^T J_pos) / (n . d_end), with the position rows of
    // the free Jacobian written out
    const double along = normal.dot(j.d.segment<3>(eRzPos0));
    const double n0 = normal.x();
    const double n1 = normal.y();
    const double n2 = normal.z();
    j.dsdv[eRzPos0] = n0;
    j.dsdv[eRzPos1] = n1;
    j.dsdv[eRzPos2] = n2;
    j.dsdv[eRzDir0] = n0 * j.f1 - n1 * j.f2;
    j.dsdv[eRzDir1] = n0 * j.f2 + n1 * j.f1;
    j.dsdv[eRzDir2] = n2 * s;
    j.dsdv[eRzQOverP] = n0 * j.a1 + n1 * j.a2 + n2 * j.a3;
    j.dsdv *= -1. / along;
    return j;
  }

  /// Fold into a fixed-path Jacobian that the step ends on a surface, i.e.
  /// that the path length itself depends on the start state.
  /// @param j the Jacobian from `stepJacobian`, updated in place
  /// @param derivativeEnd `derivative` of the end state
  /// @param normal the surface normal at the end point
  static void constrainToSurface(RzMatrix& j, const RzVector& derivativeEnd,
                                 const Vector3& normal) {
    const double along = normal.dot(derivativeEnd.segment<3>(eRzPos0));
    const Eigen::Matrix<double, 1, eRzSize> dsdv =
        -(normal.transpose() * j.block<3, eRzSize>(eRzPos0, 0)) / along;
    // the outer product by hand: the derivative has no dz and no q/p
    // component, and Eigen's general path is a call with a temporary
    for (const std::uint32_t r :
         {eRzPos0, eRzPos1, eRzPos2, eRzDir0, eRzDir1}) {
      const double dr = derivativeEnd[r];
      for (std::uint32_t c = 0; c < eRzSize; ++c) {
        j(r, c) += dr * dsdv[c];
      }
    }
  }

  /// Smallest positive path to a cylinder, or nothing if unreachable.
  /// Uses transverse circle intersections with an angular-solve fallback.
  /// @param v the state
  /// @param radius the cylinder radius
  /// @return the path length
  std::optional<double> pathToCylinder(const RzVector& v, double radius) const;

  /// Angular cylinder-intersection solution, polished by Newton iteration.
  /// @param v the state
  /// @param radius the cylinder radius
  /// @return the path length
  std::optional<double> pathToCylinderClosedForm(const RzVector& v,
                                                 double radius) const;

  /// Path length to a plane perpendicular to the beam axis, or nothing if it
  /// is behind the state or the state runs parallel to it.
  /// @param v the state
  /// @param z the plane position
  /// @return the path length
  std::optional<double> pathToDisc(const RzVector& v, double z) const {
    const double dz = v[eRzDir2];
    if (dz == 0.) {
      return std::nullopt;
    }
    const double s = (z - v[eRzPos2]) / dz;
    if (s <= 0.) {
      return std::nullopt;
    }
    return s;
  }

  /// Plane crossing by Newton iteration from a parabolic estimate.
  /// May be negative for modules behind the current RZ stop.
  /// @param v the state
  /// @param point a point on the plane
  /// @param normal the plane normal
  /// @return the path length, or nothing if parallel or the solve does not converge
  std::optional<double> pathToPlane(const RzVector& v, const Vector3& point,
                                    const Vector3& normal) const;

  /// A plane crossing together with the step that reached it
  struct PlaneStep {
    double s{};
    RzVector state;
    detail::StepTrig trig;
  };

  /// Plane crossing with cached state and trigonometry at the last iterate.
  /// Omits the final path correction, whose magnitude is below 1e-5.
  /// @param v the state
  /// @param point a point on the plane
  /// @param normal the plane normal
  /// @return the step, or nothing if parallel or the solve does not converge
  std::optional<PlaneStep> stepToPlane(const RzVector& v, const Vector3& point,
                                       const Vector3& normal) const;

  /// Signed path to transverse closest approach to the x = y = 0 axis.
  /// For an offset beam line this approximates the perigee only when the
  /// offset is small compared with the curvature radius.
  /// @param v the state
  /// @return the path length
  double pathToPerigee(const RzVector& v) const;
};

/// First-order radial-field correction after a Bz helix step.
/// Integrates q/p Br (dir x r_hat), with Br linear along the path and the
/// starting direction and radial unit vector held fixed. Higher-order field
/// terms are omitted.
/// @param v the state after the helix step, corrected in place
/// @param from the state before the step
/// @param s the path length of the step, may be negative
/// @param br0 the radial field at the start
/// @param br1 the radial field at the end
/// @param qopPosition accumulated `d(position)/d(q/p)` of the kicks since the
///        covariance was last moved, updated: an earlier change of direction
///        carries on over this step as a straight line
/// @param qopDirection the same for the direction
void rzRadialKick(RzVector& v, const RzVector& from, double s, double br0,
                  double br1, Vector3& qopPosition, Vector3& qopDirection);

}  // namespace Acts::Experimental

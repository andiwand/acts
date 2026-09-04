// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Closed-form helix transport of a free track state in a field along z, with
/// the Jacobian the Kalman filter needs and the intersections the RZ navigation
/// needs. Everything here is in double and cancellation-safe for stiff tracks.

#include "Acts/Definitions/Algebra.hpp"

#include <cmath>
#include <optional>

namespace Acts::Experimental {

/// Components of the free state the RZ finder transports: position, unit
/// direction and q/p. Time is not carried.
enum RzIndices : unsigned int {
  eRzPos0 = 0,
  eRzPos1 = 1,
  eRzPos2 = 2,
  eRzDir0 = 3,
  eRzDir1 = 4,
  eRzDir2 = 5,
  eRzQOverP = 6,
  eRzSize = 7,
};

using RzVector = Eigen::Matrix<double, eRzSize, 1>;
using RzMatrix = Eigen::Matrix<double, eRzSize, eRzSize>;

namespace detail {

/// The trigonometry of one helix step, from a single sincos: the sine and
/// cosine of the turning angle `u`, and the four finite ratios the step and
/// its Jacobian are built from, each free of cancellation.
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
  t.sn = std::sin(u);
  t.cs = std::cos(u);
  const double u2 = u * u;
  if (std::abs(u) < 1e-2) {
    t.sinc = 1. + u2 * (-1. / 6. + u2 / 120.);
    t.versc = u * (0.5 + u2 * (-1. / 24. + u2 / 720.));
    t.dsinc = u * (-1. / 3. + u2 * (1. / 30. - u2 / 840.));
    t.dversc = 0.5 + u2 * (-1. / 8. + u2 * (1. / 144. - u2 / 5760.));
    return t;
  }
  // 1 - cos(u) as sin^2 / (1 + cos), which does not cancel for small u; the
  // direct form is fine where the cosine is near -1
  const double oneMinusCos =
      1. + t.cs > 1e-3 ? t.sn * t.sn / (1. + t.cs) : 1. - t.cs;
  t.sinc = t.sn / u;
  t.versc = oneMinusCos / u;
  t.dsinc = (u * t.cs - t.sn) / u2;
  t.dversc = (u * t.sn - oneMinusCos) / u2;
  return t;
}

}  // namespace detail

/// A helix in a constant magnetic field along z. The equation of motion is
/// `d(dir)/ds = q/p * dir x B`, so with `kappa = q/p * Bz` the transverse
/// direction turns by `-kappa * s` over a path length `s` and the transverse
/// position circles the centre `(x + dy / kappa, y - dx / kappa)`.
struct RzHelix {
  /// Field along z in native units
  double bz{};

  /// Turning rate of the state
  /// @param v the state
  /// @return `q/p * Bz`, in inverse length
  double kappa(const RzVector& v) const { return v[eRzQOverP] * bz; }

  /// The state that runs the same helix the other way: direction and charge
  /// flipped. A positive path length on it is a negative one on the original,
  /// which is how the intersections are asked for backwards.
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
    const double k = kappa(v);
    const detail::StepTrig t = detail::stepTrig(k * s);
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

  /// The Jacobian of one step onto a surface kept as what it is: the
  /// identity, twelve helix entries and a rank-1 term for the path length
  /// depending on the start state. Applying it to a matrix column by column
  /// costs half of what the dense product does, and the three position
  /// rows an innovation needs come without the rest.
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
    /// The derivative of the end state and `ds/d(start)`; the rank-1 term
    /// is their outer product
    RzVector d{RzVector::Zero()};
    RzVector dsdv{RzVector::Zero()};

    /// `J x` for one column
    /// @param x the column
    /// @param y the result
    void applyLeft(const double* x, double* y) const {
      const double w = dsdv[eRzPos0] * x[eRzPos0] + dsdv[eRzPos1] * x[eRzPos1] +
                       dsdv[eRzPos2] * x[eRzPos2] + dsdv[eRzDir0] * x[eRzDir0] +
                       dsdv[eRzDir1] * x[eRzDir1] + dsdv[eRzDir2] * x[eRzDir2] +
                       dsdv[eRzQOverP] * x[eRzQOverP];
      y[eRzPos0] = x[eRzPos0] + f1 * x[eRzDir0] + f2 * x[eRzDir1] +
                   a1 * x[eRzQOverP] + d[eRzPos0] * w;
      y[eRzPos1] = x[eRzPos1] - f2 * x[eRzDir0] + f1 * x[eRzDir1] +
                   a2 * x[eRzQOverP] + d[eRzPos1] * w;
      y[eRzPos2] = x[eRzPos2] + s * x[eRzDir2] + d[eRzPos2] * w;
      y[eRzDir0] = cs * x[eRzDir0] + sn * x[eRzDir1] + b1 * x[eRzQOverP] +
                   d[eRzDir0] * w;
      y[eRzDir1] = -sn * x[eRzDir0] + cs * x[eRzDir1] + b2 * x[eRzQOverP] +
                   d[eRzDir1] * w;
      y[eRzDir2] = x[eRzDir2];
      y[eRzQOverP] = x[eRzQOverP];
    }

    /// `J X`
    /// @param x the matrix
    /// @return the product
    RzMatrix applyLeft(const RzMatrix& x) const {
      RzMatrix y;
      for (unsigned int c = 0; c < eRzSize; ++c) {
        applyLeft(x.col(c).data(), y.col(c).data());
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
      for (unsigned int col = 0; col < eRzSize; ++col) {
        w[col] = dsdv.dot(x.col(col));
      }
      RzMatrix y;
      auto row = [&](unsigned int r, auto&& entry) {
        for (unsigned int col = 0; col <= r; ++col) {
          y(r, col) = entry(col) + d[r] * w[col];
        }
      };
      row(eRzPos0, [&](unsigned int col) {
        return x(eRzPos0, col) + f1 * x(eRzDir0, col) + f2 * x(eRzDir1, col) +
               a1 * x(eRzQOverP, col);
      });
      row(eRzPos1, [&](unsigned int col) {
        return x(eRzPos1, col) - f2 * x(eRzDir0, col) + f1 * x(eRzDir1, col) +
               a2 * x(eRzQOverP, col);
      });
      row(eRzPos2, [&](unsigned int col) {
        return x(eRzPos2, col) + s * x(eRzDir2, col);
      });
      row(eRzDir0, [&](unsigned int col) {
        return cs * x(eRzDir0, col) + sn * x(eRzDir1, col) +
               b1 * x(eRzQOverP, col);
      });
      row(eRzDir1, [&](unsigned int col) {
        return -sn * x(eRzDir0, col) + cs * x(eRzDir1, col) +
               b2 * x(eRzQOverP, col);
      });
      for (unsigned int col = 0; col <= eRzDir2; ++col) {
        y(eRzDir2, col) = x(eRzDir2, col);
      }
      for (unsigned int col = 0; col < eRzSize; ++col) {
        y(eRzQOverP, col) = x(eRzQOverP, col);
      }
      for (unsigned int r = 0; r < eRzSize; ++r) {
        for (unsigned int col = r + 1; col < eRzSize; ++col) {
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
      for (unsigned int i = 0; i < 3; ++i) {
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
    const double k = kappa(v0);
    const detail::StepTrig t = detail::stepTrig(k * s);
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
    j.dsdv[eRzQOverP] = n0 * j.a1 + n1 * j.a2;
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
    for (const unsigned int r : {eRzPos0, eRzPos1, eRzPos2, eRzDir0, eRzDir1}) {
      const double dr = derivativeEnd[r];
      for (unsigned int c = 0; c < eRzSize; ++c) {
        j(r, c) += dr * dsdv[c];
      }
    }
  }

  /// Smallest positive path length to a cylinder around the beam axis, or
  /// nothing if the helix never reaches it on its way forward. Newton from
  /// the straight-line crossing for the mildly curved case, which is every
  /// layer gap a tracker has, and the closed form otherwise.
  /// @param v the state
  /// @param radius the cylinder radius
  /// @return the path length
  std::optional<double> pathToCylinder(const RzVector& v, double radius) const;

  /// The same by the closed form alone: the law of cosines for the turning
  /// angle, polished once. Exact for any curvature, dearer in transcendentals.
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

  /// Path length to a plane, by Newton iteration from the straight-line
  /// estimate. May be negative: a module sits on either side of the RZ
  /// surface the state was brought to.
  /// @param v the state
  /// @param point a point on the plane
  /// @param normal the plane normal
  /// @return the path length, or nothing if the state runs parallel
  std::optional<double> pathToPlane(const RzVector& v, const Vector3& point,
                                    const Vector3& normal) const;

  /// Path length to the point of closest approach to the beam axis, in the
  /// transverse plane. Negative when the perigee is behind the state.
  /// @param v the state
  /// @return the path length
  double pathToPerigee(const RzVector& v) const;

  /// Residual of the cylinder equation at a path length, for polishing
  /// @param v the state
  /// @param s the path length
  /// @param radius the cylinder radius
  /// @return `|p_T(s)|^2 - radius^2`
  double cylinderResidual(const RzVector& v, double s, double radius) const {
    RzVector w = v;
    step(w, s);
    return w[eRzPos0] * w[eRzPos0] + w[eRzPos1] * w[eRzPos1] - radius * radius;
  }
};

}  // namespace Acts::Experimental

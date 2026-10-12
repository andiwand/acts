// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"

#include <array>
#include <cstdint>

namespace Acts::detail {

/// Transport the angular and scaled q/p columns of the ATLAS Jacobian.
///
/// This applies the two-lane angular arrangement of Athena's Trk::propJacobian
/// to the ACTS scalar equations and 8-by-6 storage (including time). Eigen's
/// fixed-size row-major rows provide portable SIMD pairs without CxxUtils.
/// Source:
/// https://gitlab.cern.ch/atlas/athena/-/blob/c5982a1efd76fe5528cc550f4dc0a87371ec1c67/Tracking/TrkExtrapolation/TrkExUtils/TrkExUtils/JacobianHelper.h
/// Field gradients are deliberately excluded from this fast path.
/// @param P Parameters and derivatives in AtlasStepper's flat storage
/// @param H0 Half-step field at the start
/// @param H1 Half-step field at the midpoint
/// @param H2 Half-step field at the end
/// @param A Direction before the step
/// @param A0 First half-step bend
/// @param A3 Second intermediate direction
/// @param A4 Third intermediate direction
/// @param A6 Last half-step bend
/// @param S3 One third of the step length
inline void atlasTransportJacobian(double* P, const double* H0,
                                   const double* H1, const double* H2,
                                   const double* A, const double* A0,
                                   const double* A3, const double* A4,
                                   const double* A6, double S3) {
  using Angular = Eigen::Matrix<double, 3, 2, Eigen::RowMajor>;
  const auto bend = [](const Angular& d, const double* H) -> Angular {
    Angular result;
    result.row(0) = d.row(1) * H[2] - d.row(2) * H[1];
    result.row(1) = d.row(2) * H[0] - d.row(0) * H[2];
    result.row(2) = d.row(0) * H[1] - d.row(1) * H[0];
    return result;
  };
  Angular d;
  d << P[28], P[36], P[29], P[37], P[30], P[38];
  const Angular kick0 = bend(d, H0);
  const Angular dir2 = d + kick0;
  const Angular dir3 = d + bend(dir2, H1);
  const Angular dir4 = d + bend(dir3, H1);
  const Angular dir5 = 2. * dir4 - d;
  const Angular kick6 = bend(dir5, H2);
  const Angular position = (dir2 + dir3 + dir4) * S3;
  const Angular direction = (kick0 + 2. * dir3 + dir5 + kick6) * (1. / 3.);
  for (int i = 0; i < 3; ++i) {
    P[24 + i] += position(i, 0);
    P[32 + i] += position(i, 1);
    P[28 + i] = direction(i, 0);
    P[36 + i] = direction(i, 1);
  }

  double* d4A = P + 44;
  double d4A0 = (A0[0] + H0[2] * d4A[1]) - H0[1] * d4A[2];
  double d4B0 = (A0[1] + H0[0] * d4A[2]) - H0[2] * d4A[0];
  double d4C0 = (A0[2] + H0[1] * d4A[0]) - H0[0] * d4A[1];
  double d4A2 = d4A0 + d4A[0];
  double d4B2 = d4B0 + d4A[1];
  double d4C2 = d4C0 + d4A[2];
  double d0 = d4A[0] - A[0];
  double d1 = d4A[1] - A[1];
  double d2 = d4A[2] - A[2];
  double d4A3 = ((A3[0] + d0) + d4B2 * H1[2]) - d4C2 * H1[1];
  double d4B3 = ((A3[1] + d1) + d4C2 * H1[0]) - d4A2 * H1[2];
  double d4C3 = ((A3[2] + d2) + d4A2 * H1[1]) - d4B2 * H1[0];
  double d4A4 = ((A4[0] + d0) + d4B3 * H1[2]) - d4C3 * H1[1];
  double d4B4 = ((A4[1] + d1) + d4C3 * H1[0]) - d4A3 * H1[2];
  double d4C4 = ((A4[2] + d2) + d4A3 * H1[1]) - d4B3 * H1[0];
  double d4A5 = 2. * d4A4 - d4A[0];
  double d4B5 = 2. * d4B4 - d4A[1];
  double d4C5 = 2. * d4C4 - d4A[2];
  double d4A6 = d4B5 * H2[2] - d4C5 * H2[1];
  double d4B6 = d4C5 * H2[0] - d4A5 * H2[2];
  double d4C6 = d4A5 * H2[1] - d4B5 * H2[0];
  P[40] += (d4A2 + d4A3 + d4A4) * S3;
  P[41] += (d4B2 + d4B3 + d4B4) * S3;
  P[42] += (d4C2 + d4C3 + d4C4) * S3;
  d4A[0] = ((d4A0 + 2. * d4A3) + (d4A5 + d4A6 + A6[0])) * (1. / 3.);
  d4A[1] = ((d4B0 + 2. * d4B3) + (d4B5 + d4B6 + A6[1])) * (1. / 3.);
  d4A[2] = ((d4C0 + 2. * d4C3) + (d4C5 + d4C6 + A6[2])) * (1. / 3.);
}

/// Transport all spatial Jacobian columns including field derivatives.
///
/// Adaptation of Athena's rungeKuttaStepWithGradient variational equations to
/// the ACTS 8-by-6 layout and scaled q/p column. Unlike the fast angular-only
/// helper, this also evolves the local-position columns in a varying field.
/// @param P Parameters and derivatives in AtlasStepper's flat storage
/// @param fields Magnetic fields at the three RK sample positions
/// @param gradients Field derivatives, G(i,j) = d B_i / d x_j
/// @param h Signed step length
/// @param qOverP Signed curvature coefficient
inline void atlasTransportJacobianWithGradient(
    double* P, const std::array<Vector3, 3>& fields,
    const std::array<SquareMatrix3, 3>& gradients, double h, double qOverP) {
  const Vector3 A(P[4], P[5], P[6]);
  const double halfStep = 0.5 * h * qOverP;
  const Vector3 H0 = halfStep * fields[0];
  const Vector3 H1 = halfStep * fields[1];
  const Vector3 H2 = halfStep * fields[2];
  const Vector3 A2 = A + A.cross(H0);
  const Vector3 A3 = A + A2.cross(H1);
  const Vector3 A4 = A + A3.cross(H1);
  const Vector3 A5 = 2. * A4 - A;

  // The time column has no spatial components in a static magnetic field.
  for (std::uint32_t column = 0; column < eBoundTime; ++column) {
    double* J = P + eFreeSize * (column + 1);
    const Vector3 dR(J[0], J[1], J[2]);
    const Vector3 dA(J[4], J[5], J[6]);
    const auto fieldDerivative =
        [&](int stage, const Vector3& positionDerivative) -> Vector3 {
      Vector3 result = halfStep * gradients[stage] * positionDerivative;
      if (column == eBoundQOverP) {
        // This column stores q/p times the spatial derivative w.r.t. q/p.
        result += halfStep * fields[stage];
      }
      return result;
    };
    const Vector3 dH0 = fieldDerivative(0, dR);
    const Vector3 kick0 = dA.cross(H0) + A.cross(dH0);
    const Vector3 dir2 = dA + kick0;
    const Vector3 dH1 = fieldDerivative(1, dR + (0.25 * h) * (dA + dir2));
    const Vector3 dir3 = dA + dir2.cross(H1) + A2.cross(dH1);
    const Vector3 dir4 = dA + dir3.cross(H1) + A3.cross(dH1);
    const Vector3 dir5 = 2. * dir4 - dA;
    const Vector3 dH2 = fieldDerivative(2, dR + h * dir4);
    const Vector3 kick6 = dir5.cross(H2) + A5.cross(dH2);
    Eigen::Map<Vector3>{J} = dR + (h / 3.) * (dir2 + dir3 + dir4);
    Eigen::Map<Vector3>(J + eFreeDir0) =
        (kick0 + 2. * dir3 + dir5 + kick6) / 3.;
  }
}

}  // namespace Acts::detail

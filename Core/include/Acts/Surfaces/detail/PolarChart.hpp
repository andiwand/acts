// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cmath>

namespace Acts::detail {

/// Polar coordinates (r, phi) in a surface's local Cartesian frame.
/// These maps describe the coordinate chart, independently of its accepted
/// region, sector orientation, or placement in global coordinates.
struct PolarChart {
  /// Convert local polar coordinates to local Cartesian coordinates.
  /// @param polar Local (r, phi) position
  /// @return Local (x, y) position
  static Vector2 toCartesian(const Vector2& polar) {
    return {polar[0] * std::cos(polar[1]), polar[0] * std::sin(polar[1])};
  }

  /// Convert local Cartesian coordinates to local polar coordinates.
  /// @param cartesian Local (x, y) position
  /// @return Local (r, phi) position
  static Vector2 fromCartesian(const Vector2& cartesian) {
    return {cartesian.norm(), std::atan2(cartesian[1], cartesian[0])};
  }

  /// Derivative of local Cartesian position with respect to (r, phi).
  /// @param polar Local (r, phi) position
  /// @return Polar-to-Cartesian Jacobian
  static SquareMatrix2 toCartesianJacobian(const Vector2& polar) {
    const double cosPhi = std::cos(polar[1]);
    const double sinPhi = std::sin(polar[1]);
    SquareMatrix2 jacobian;
    jacobian << cosPhi, -polar[0] * sinPhi, sinPhi, polar[0] * cosPhi;
    return jacobian;
  }

  /// Derivative of (r, phi) with respect to local Cartesian position.
  /// @param cartesian Nonzero local (x, y) position; undefined at the origin
  /// @return Cartesian-to-polar Jacobian
  static SquareMatrix2 fromCartesianJacobian(const Vector2& cartesian) {
    const double r = cartesian.norm();
    const double cosPhi = cartesian[0] / r;
    const double sinPhi = cartesian[1] / r;
    SquareMatrix2 jacobian;
    jacobian << cosPhi, sinPhi, -sinPhi / r, cosPhi / r;
    return jacobian;
  }

  /// Euclidean local metric expressed in (r, phi) coordinates.
  /// @param polar Local (r, phi) position
  /// @return The metric J.transpose() * J
  static SquareMatrix2 metric(const Vector2& polar) {
    SquareMatrix2 result = SquareMatrix2::Zero();
    result(0, 0) = 1.;
    result(1, 1) = polar[0] * polar[0];
    return result;
  }
};

}  // namespace Acts::detail

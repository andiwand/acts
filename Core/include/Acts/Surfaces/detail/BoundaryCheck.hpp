// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"

#include <cmath>
#include <stdexcept>

namespace Acts::detail {

/// Evaluate a boundary tolerance using an explicitly supplied coordinate chart.
/// The region only needs inside(position) and closestPoint(position, metric);
/// it need not inherit from SurfaceBounds or describe a coordinate system.
///
/// The metric and boundary point are expressed in the input coordinates.
/// Cartesian tolerances use the Jacobian at the query position, i.e. a local
/// linear approximation, not the finite Cartesian displacement to the boundary.
/// The chart is evaluated lazily so exact/infinite checks do not need a regular
/// Jacobian. No inverse chart is required, including at a polar origin.
///
/// @tparam Region Region providing strict containment and boundary projection
/// @tparam Jacobian Callable returning the input-to-Cartesian Jacobian
/// @param region Accepted region
/// @param position Query in surface local coordinates
/// @param tolerance Boundary tolerance
/// @param jacobian Coordinate chart evaluated at the query position
/// @return Whether the position satisfies the boundary tolerance
template <typename Region, typename Jacobian>
bool insideWithTolerance(const Region& region, const Vector2& position,
                         const BoundaryTolerance& tolerance,
                         const Jacobian& jacobian) {
  using enum BoundaryTolerance::ToleranceMode;
  if (tolerance.isInfinite()) {
    return true;
  }
  const auto mode = tolerance.toleranceMode();
  const bool strictlyInside = region.inside(position);
  if (mode == None || (mode == Extend && strictlyInside)) {
    return strictlyInside;
  }

  const SquareMatrix2 boundToCartesian = jacobian(position);
  SquareMatrix2 metric;
  if (tolerance.hasAbsoluteEuclidean()) {
    metric = boundToCartesian.transpose() * boundToCartesian;
  } else if (tolerance.hasChi2Bound()) {
    // Preserve the legacy projection convention during extraction. This
    // extra chart transformation for a bound-coordinate weight needs a
    // separate numerical correction together with the polar region kernels.
    metric = boundToCartesian.transpose() *
             tolerance.asChi2Bound().weightMatrix() * boundToCartesian;
  } else if (tolerance.hasChi2Cartesian()) {
    metric = boundToCartesian.transpose() *
             tolerance.asChi2Cartesian().weightMatrix() * boundToCartesian;
  } else {
    throw std::runtime_error("Unsupported boundary tolerance type.");
  }

  const Vector2 delta = region.closestPoint(position, metric) - position;
  const bool tolerated = tolerance.isTolerated(delta, boundToCartesian);
  return tolerated && (mode != Shrink || strictlyInside);
}

/// Distance to a region boundary in an explicitly supplied local metric.
/// @tparam Region Region providing closestPoint(position, metric)
/// @param region Accepted region
/// @param position Query in surface local coordinates
/// @param metric Metric expressed in those same coordinates, at the query
/// @return Unsigned distance using the local metric (not a geodesic distance)
template <typename Region>
double boundaryDistance(const Region& region, const Vector2& position,
                        const SquareMatrix2& metric) {
  const Vector2 delta = region.closestPoint(position, metric) - position;
  return std::sqrt(delta.dot(metric * delta));
}

}  // namespace Acts::detail

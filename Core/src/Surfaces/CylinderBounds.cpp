// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/CylinderBounds.hpp"

#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/detail/OstreamStateGuard.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <utility>

namespace Acts {

using VectorHelpers::perp;
using VectorHelpers::phi;

std::vector<double> CylinderBounds::values() const {
  return {m_values.begin(), m_values.end()};
}

Vector2 CylinderBounds::shifted(const Vector2& lposition) const {
  return {detail::radian_sym((lposition[0] / get(eR)) - get(eAveragePhi)),
          lposition[1]};
}

bool CylinderBounds::inside(const Vector2& lposition) const {
  double halfLengthZ = get(eHalfLengthZ);
  if (coversFullAzimuth()) {
    return -halfLengthZ <= lposition.y() && lposition.y() < halfLengthZ;
  }
  double halfPhi = get(eHalfPhiSector);
  return detail::VerticesHelper::isInsideRectangle(
      shifted(lposition), Vector2(-halfPhi, -halfLengthZ),
      Vector2(halfPhi, halfLengthZ));
}

Vector2 CylinderBounds::closestPoint(const Vector2& lposition,
                                     const SquareMatrix2& metric) const {
  const double halfZ = get(eHalfLengthZ);
  const double radius = get(eR);
  const double halfRPhi = radius * get(eHalfPhiSector);
  const double period = 2. * std::numbers::pi * radius;
  // Work in an unrolled chart centered on the sector. Return a representative
  // near the query, so callers can subtract coordinates across the phi seam.
  const Vector2 point(radius * shifted(lposition).x(), lposition.y());
  Vector2 closest = point;
  double bestDistance = std::numeric_limits<double>::infinity();
  auto consider = [&](const Vector2& candidate) {
    const Vector2 delta = candidate - point;
    const double distance = delta.dot(metric * delta);
    if (distance < bestDistance) {
      bestDistance = distance;
      closest = candidate;
    }
  };

  // On each z boundary, first project onto the infinite line in the supplied
  // metric. For a sector, restrict that projection to the nearest periodic
  // copy of the horizontal segment.
  for (double z : {-halfZ, halfZ}) {
    double x = point.x() - metric(0, 1) / metric(0, 0) * (z - point.y());
    if (!coversFullAzimuth()) {
      const double center = period * std::round(x / period);
      x = std::clamp(x, center - halfRPhi, center + halfRPhi);
    }
    consider(Vector2(x, z));
  }

  if (!coversFullAzimuth()) {
    // For a vertical edge, the distance minimized over z is convex in its
    // unrolled x coordinate. Its continuous minimum has z clamped to the
    // allowed interval. Check the two periodic copies bracketing that minimum.
    const double z = std::clamp(point.y(), -halfZ, halfZ);
    const double x = point.x() - metric(0, 1) / metric(0, 0) * (z - point.y());
    for (double edge : {-halfRPhi, halfRPhi}) {
      const double firstCopy = std::floor((x - edge) / period);
      for (double copy : {firstCopy, firstCopy + 1.}) {
        const double edgeX = edge + copy * period;
        const double edgeZ = std::clamp(
            point.y() - metric(0, 1) / metric(1, 1) * (edgeX - point.x()),
            -halfZ, halfZ);
        consider(Vector2(edgeX, edgeZ));
      }
    }
  }
  return lposition + (closest - point);
}

Vector2 CylinderBounds::center() const {
  return Vector2(get(eR) * get(eAveragePhi), 0.0);
}

std::ostream& CylinderBounds::toStream(std::ostream& sl) const {
  detail::OstreamStateGuard guard{sl};
  sl << std::fixed << std::setprecision(7);
  sl << "Acts::CylinderBounds: (radius, halfLengthZ, halfPhiSector, "
        "averagePhi) = ";
  sl << "(" << get(eR) << ", " << get(eHalfLengthZ) << ", ";
  sl << get(eHalfPhiSector) << ", " << get(eAveragePhi) << ")";
  return sl;
}

std::vector<Vector3> CylinderBounds::circleVertices(
    const Transform3 transform, unsigned int quarterSegments) const {
  std::vector<Vector3> vertices;

  double avgPhi = get(eAveragePhi);
  double halfPhi = get(eHalfPhiSector);

  std::vector<double> phiRef = {};
  if (bool fullCylinder = coversFullAzimuth(); fullCylinder) {
    phiRef = {avgPhi};
  }

  // Write the two bows/circles on either side
  std::vector<int> sides = {-1, 1};
  for (auto& side : sides) {
    // Helper method to create the segment
    auto svertices = detail::VerticesHelper::segmentVertices(
        {get(eR), get(eR)}, avgPhi - halfPhi, avgPhi + halfPhi, phiRef,
        quarterSegments, Vector3(0., 0., side * get(eHalfLengthZ)), transform);
    vertices.insert(vertices.end(), svertices.begin(), svertices.end());
  }

  return vertices;
}

void CylinderBounds::checkConsistency() noexcept(false) {
  if (get(eR) <= 0.) {
    throw std::invalid_argument(
        "CylinderBounds: invalid radial setup: radius is negative");
  }
  if (get(eHalfLengthZ) <= 0.) {
    throw std::invalid_argument(
        "CylinderBounds: invalid length setup: half length is negative");
  }
  if (get(eHalfPhiSector) <= 0. || get(eHalfPhiSector) > std::numbers::pi) {
    throw std::invalid_argument("CylinderBounds: invalid phi sector setup.");
  }
  if (get(eAveragePhi) != detail::radian_sym(get(eAveragePhi)) &&
      std::abs(std::abs(get(eAveragePhi)) - std::numbers::pi) > s_epsilon) {
    throw std::invalid_argument("CylinderBounds: invalid phi positioning.");
  }
}

}  // namespace Acts

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/SurfaceBoundsImpl.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"

namespace Acts {

inline SurfaceMultiIntersection DiscSurface::intersectImpl(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryTolerance& boundaryTolerance,
    ActsScalar tolerance) const {
  // Get the contextual transform
  auto gctxTransform = transform(gctx);
  // Use the intersection helper for planar surfaces
  auto intersection =
      PlanarHelper::intersect(gctxTransform, position, direction, tolerance);
  auto status = intersection.status();
  // Evaluate boundary check if requested (and reachable)
  if (intersection.status() != Intersection3D::Status::unreachable &&
      m_bounds != nullptr && !boundaryTolerance.isInfinite()) {
    // Built-in local to global for speed reasons
    const auto& tMatrix = gctxTransform.matrix();
    const Vector3 vecLocal(intersection.position() - tMatrix.block<3, 1>(0, 3));
    const Vector2 lcartesian = tMatrix.block<3, 2>(0, 0).transpose() * vecLocal;
    if (auto absoluteBound = boundaryTolerance.asAbsoluteBoundOpt();
        absoluteBound.has_value() && m_bounds->coversFullAzimuth()) {
      double modifiedTolerance = tolerance + absoluteBound->tolerance0;
      if (!m_bounds->insideRadialBounds(VectorHelpers::perp(lcartesian),
                                        modifiedTolerance)) {
        status = Intersection3D::Status::missed;
      }
    } else if (!insideBounds(localCartesianToPolar(lcartesian),
                             boundaryTolerance)) {
      status = Intersection3D::Status::missed;
    }
  }
  return {{Intersection3D(intersection.position(), intersection.pathLength(),
                          status),
           Intersection3D::invalid()},
          this};
}

}  // namespace Acts

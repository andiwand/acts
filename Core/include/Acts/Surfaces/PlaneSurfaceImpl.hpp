// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/SurfaceBoundsImpl.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"

namespace Acts {

inline SurfaceMultiIntersection PlaneSurface::intersectImpl(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryTolerance& boundaryTolerance,
    ActsScalar tolerance) const {
  // Get the contextual transform
  const auto& gctxTransform = transform(gctx);
  // Use the intersection helper for planar surfaces
  auto intersection =
      PlanarHelper::intersect(gctxTransform, position, direction, tolerance);
  auto status = intersection.status();
  // Evaluate boundary check if requested (and reachable)
  if (intersection.status() != Intersection3D::Status::unreachable) {
    // Built-in local to global for speed reasons
    const auto& tMatrix = gctxTransform.matrix();
    // Create the reference vector in local
    const Vector3 vecLocal(intersection.position() - tMatrix.block<3, 1>(0, 3));
    if (!bounds().insideImpl(tMatrix.block<3, 2>(0, 0).transpose() * vecLocal,
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

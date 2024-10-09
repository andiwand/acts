// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/DiscSurfaceImpl.hpp"
#include "Acts/Surfaces/PlaneSurfaceImpl.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

inline SurfaceMultiIntersection Surface::intersectImpl(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryTolerance& boundaryTolerance,
    ActsScalar tolerance) const {
  switch (type()) {
    case SurfaceType::Plane:
      return static_cast<const PlaneSurface*>(this)->intersectImpl(
          gctx, position, direction, boundaryTolerance, tolerance);
    case SurfaceType::Disc:
      return static_cast<const DiscSurface*>(this)->intersectImpl(
          gctx, position, direction, boundaryTolerance, tolerance);
    default:
      return intersect(gctx, position, direction, boundaryTolerance, tolerance);
  }
}

}  // namespace Acts

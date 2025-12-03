// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/RectangleSurface.hpp"

namespace Acts {

RectangleSurface::RectangleSurface(const RectangleSurface& other)
    : GeometryObject(), PlaneSurface(other), m_bounds(other.m_bounds) {}

RectangleSurface::RectangleSurface(const GeometryContext& gctx,
                                   const RectangleSurface& other,
                                   const Transform3& transform)
    : GeometryObject(),
      PlaneSurface(gctx, other, transform),
      m_bounds(other.m_bounds) {}

RectangleSurface::RectangleSurface(
    std::shared_ptr<const RectangleBounds> bounds,
    const DetectorElementBase& detelement)
    : PlaneSurface(bounds, detelement), m_bounds(*bounds) {}

RectangleSurface::RectangleSurface(
    const Transform3& transform, std::shared_ptr<const RectangleBounds> bounds)
    : PlaneSurface(transform, bounds), m_bounds(*bounds) {}

RectangleSurface& RectangleSurface::operator=(const RectangleSurface& other) {
  if (this != &other) {
    Surface::operator=(other);
    m_bounds = other.m_bounds;
  }
  return *this;
}

const SurfaceBounds& RectangleSurface::bounds() const {
  return m_bounds;
}

const std::shared_ptr<const RectangleBounds>& RectangleSurface::boundsPtr()
    const {
  static thread_local std::shared_ptr<const RectangleBounds> ptr =
      std::make_shared<const RectangleBounds>(m_bounds);
  return ptr;
}

void RectangleSurface::assignSurfaceBounds(
    std::shared_ptr<const RectangleBounds> newBounds) {
  m_bounds = *newBounds;
  PlaneSurface::assignSurfaceBounds(newBounds);
}

}  // namespace Acts

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/RadialBoundsImpl.hpp"
#include "Acts/Surfaces/RectangleBoundsImpl.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <iostream>

namespace Acts {

inline bool SurfaceBounds::insideImpl(
    const Vector2& localPosition,
    const BoundaryTolerance& boundaryTolerance) const {
  switch (type()) {
    case SurfaceBounds::eRectangle:
      return static_cast<const RectangleBounds*>(this)->insideImpl(
          localPosition, boundaryTolerance);
    case SurfaceBounds::eDisc:
      return static_cast<const RadialBounds*>(this)->insideImpl(
          localPosition, boundaryTolerance);
    default:
      return inside(localPosition, boundaryTolerance);
  }
}

}  // namespace Acts

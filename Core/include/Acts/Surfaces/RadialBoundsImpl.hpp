// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/detail/BoundaryCheckHelper.hpp"

namespace Acts {

inline bool RadialBounds::insideImpl(
    const Vector2& lposition,
    const BoundaryTolerance& boundaryTolerance) const {
  return detail::insideAlignedBox(Vector2(get(eMinR), -get(eHalfPhiSector)),
                                  Vector2(get(eMaxR), get(eHalfPhiSector)),
                                  boundaryTolerance, shifted(lposition),
                                  std::nullopt);
}

}  // namespace Acts

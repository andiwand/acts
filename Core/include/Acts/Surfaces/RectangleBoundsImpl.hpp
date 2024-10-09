// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/detail/BoundaryCheckHelper.hpp"

namespace Acts {

inline bool RectangleBounds::insideImpl(
    const Vector2& localPosition,
    const BoundaryTolerance& boundaryTolerance) const {
  return detail::insideAlignedBox(m_min, m_max, boundaryTolerance,
                                  localPosition, std::nullopt);
}

}  // namespace Acts

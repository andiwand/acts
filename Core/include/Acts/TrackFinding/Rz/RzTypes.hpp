// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstdint>
#include <limits>

namespace Acts::Experimental {

/// Index used when an RZ module, hit, layer, or stop is absent.
constexpr std::uint32_t kRzNone = std::numeric_limits<std::uint32_t>::max();

/// Components of the spatial free state; time is carried separately.
enum RzIndices : std::uint8_t {
  eRzPos0 = 0,
  eRzPos1 = 1,
  eRzPos2 = 2,
  eRzDir0 = 3,
  eRzDir1 = 4,
  eRzDir2 = 5,
  eRzQOverP = 6,
  eRzSize = 7,
};

using RzVector = Eigen::Matrix<double, eRzSize, 1>;
using RzMatrix = Eigen::Matrix<double, eRzSize, eRzSize>;

}  // namespace Acts::Experimental

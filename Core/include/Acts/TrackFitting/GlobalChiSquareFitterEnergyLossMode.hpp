// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>

namespace Acts::Experimental {

/// @addtogroup track_fitting
/// @{

/// Estimator for the deterministic energy loss correction.
/// Both estimators include radiation; the KF uses only Bethe ionisation loss.
/// Mode is an empirical approximation, not the Landau most probable value.
enum class Gx2fEnergyLossMode : std::uint8_t {
  /// Bethe (ionisation) + radiative, see @ref Acts::computeEnergyLossMean
  Mean,
  /// 0.9 * Bethe + 0.15 * radiative, see @ref Acts::computeEnergyLossMode
  Mode,
};

/// @}

}  // namespace Acts::Experimental

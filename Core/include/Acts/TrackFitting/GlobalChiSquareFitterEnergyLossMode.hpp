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

/// Selects the estimator used for the deterministic energy loss correction of
/// the GX2F
///
/// @c Mean is the default, so that the GX2F agrees with the KF and the CKF,
/// which both apply the mean loss. Note that they use
/// @ref Acts::computeEnergyLossBethe, i.e. the mean *ionisation* loss only,
/// whereas @c Mean here also carries the radiative term.
///
/// @note Neither of these is the Landau most probable value. @c Mode is the
/// empirical approximation of the most probable loss from
/// ATL-SOFT-PUB-2008-003. The true most probable ionisation loss is
/// @ref Acts::computeEnergyLossLandau.
///
/// @note This lives in its own header, so that it can be used without pulling
/// in the whole GX2F template header.
enum class Gx2fEnergyLossMode : std::uint8_t {
  /// Bethe (ionisation) + radiative, see @ref Acts::computeEnergyLossMean
  Mean,
  /// 0.9 * Bethe + 0.15 * radiative, see @ref Acts::computeEnergyLossMode
  Mode,
};

/// @}

}  // namespace Acts::Experimental

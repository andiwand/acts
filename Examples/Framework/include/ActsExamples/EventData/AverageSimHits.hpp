// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"

#include <tuple>

namespace ActsExamples {

/// Create (average) truth representation for selected simulated hits.
///
/// @param gCtx The geometry context for this
/// @param surface The reference surface of the measurement
/// @param simHits Selection of simulated hits from the container
/// @return a local position, a 4D global position, a direction
///
/// If more than one simulated hit is selected, the average truth information is
/// returned.
inline std::tuple<Acts::Vector2, Acts::Vector4, Acts::Vector3> averageSimHits(
    const Acts::GeometryContext& gCtx, const Acts::Surface& surface,
    const SimHitContainer::ConstSubset& simHits, const Acts::Logger& logger) {
  Acts::Vector2 avgLocal = Acts::Vector2::Zero();
  Acts::Vector4 avgPos4 = Acts::Vector4::Zero();
  Acts::Vector3 avgDir = Acts::Vector3::Zero();

  std::size_t n = 0u;
  for (const auto& simHit : simHits) {
    n += 1u;

    // We use the thickness of the surface as tolerance, because Geant4
    // treats the Surfaces as volumes and thus it is not ensured, that each hit
    // lies exactly on the Acts::Surface
    const double tolerance = surface.isSensitive() ? surface.thickness()
                                                   : Acts::s_onSurfaceTolerance;

    // transforming first to local positions and average that ensures that the
    // averaged position is still on the surface. the averaged global position
    // might not be on the surface anymore.
    const Acts::Result<Acts::Vector2> result = surface.globalToLocal(
        gCtx, simHit.position(), simHit.direction(), tolerance);
    if (result.ok()) {
      avgLocal += result.value();
    } else {
      ACTS_WARNING(
          "While averaging simhit, hit is not on the corresponding surface "
          << surface.geometryId() << "; use [0,0] as local position");
    }
    // global position should already be at the intersection. no need to perform
    // an additional intersection call.
    avgPos4 += simHit.fourPosition();
    avgDir += simHit.direction();
  }

  // only need to average if there are at least two inputs
  if (2u <= n) {
    double scale = 1.0 / n;
    avgLocal *= scale;
    avgPos4 *= scale;
    avgDir.normalize();
  }

  return {avgLocal, avgPos4, avgDir};
}

}  // namespace ActsExamples

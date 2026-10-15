// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Truth bookkeeping for a synthetic event, cheap enough to run alongside a
/// benchmark.

#include "Acts/EventData/SeedContainer.hpp"
#include "ActsFatras/Synthetic/SyntheticEvent.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace ActsFatras::Synthetic {

/// Truth content of an event, so that a caller can report what it is running
/// on.
struct EventSummary {
  /// Number of space points, both collections together
  std::size_t spacePoints{};
  /// How many of them the strip layers left
  std::size_t stripSpacePoints{};
  /// Number of generated primaries
  std::size_t primaries{};
  /// Number of generated secondaries
  std::size_t secondaries{};
  /// Primaries above the momentum threshold that leave enough of the selected
  /// space points to be seedable at all
  std::size_t seedablePrimaries{};
  /// Selected space points left by primaries
  std::size_t primaryHits{};
  /// Selected space points left by secondaries
  std::size_t secondaryHits{};
  /// Number of daughters of displaced vertices, which are counted among the
  /// secondaries as well
  std::size_t displaced{};
  /// Daughters of displaced vertices that are seedable by the same measure as
  /// `seedablePrimaries`
  std::size_t seedableDisplaced{};
  /// Selected space points per particle, indexed like `Event::particles`, so
  /// that a caller can bin the denominator of an efficiency itself
  std::vector<std::uint32_t> particleHits;
};

/// Summarise the truth content of an event.
///
/// @param event the generated event
/// @param ptThreshold the momentum threshold a primary must pass
/// @param selection which space points count towards a primary being seedable,
///        a primary crossing only the strips being no candidate for a
///        pixel-only pass
/// @return the summary
EventSummary summarize(
    const Event& event, float ptThreshold,
    SpacePointSelection selection = SpacePointSelection::Combined);

/// What the seeds of an event look like against the generator truth.
struct SeedingSummary {
  /// Number of seeds
  std::size_t seeds{};
  /// Seeds sharing enough space points with one primary above the momentum
  /// threshold
  std::size_t trueSeeds{};
  /// Primaries above the momentum threshold that have at least one true seed
  std::size_t matchedPrimaries{};
  /// Seeds sharing enough space points with one daughter of a displaced vertex
  /// above the momentum threshold
  std::size_t trueDisplacedSeeds{};
  /// Daughters of displaced vertices above the momentum threshold that have at
  /// least one seed of their own
  std::size_t matchedDisplaced{};
  /// Whether a particle has a seed of its own, indexed like `Event::particles`;
  /// only primaries and daughters of displaced vertices above the momentum
  /// threshold are ever set
  std::vector<bool> matched;
};

/// Match seeds against the generator truth.
/// @param event the generated event
/// @param seeds the seeds found on it, indexing the container they were made
///        from, which is the one they carry
/// @param ptThreshold the momentum threshold a primary, or a daughter of a
///        displaced vertex, must pass
/// @param minTrueSpacePoints space points a seed has to share with one particle
///        to be its seed; three constrain a helix
/// @return the summary
SeedingSummary evaluateSeeds(const Event& event,
                             const Acts::SeedContainer& seeds,
                             float ptThreshold,
                             std::size_t minTrueSpacePoints = 3);

}  // namespace ActsFatras::Synthetic

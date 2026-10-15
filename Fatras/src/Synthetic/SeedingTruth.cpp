// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/SeedingTruth.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace ActsFatras {

Synthetic::EventSummary Synthetic::summarize(
    const Event& event, const float ptThreshold,
    const SpacePointSelection selection) {
  // Counted off the collections rather than taken from `numHits`, so that a
  // pass seeded on one of them is scored against what it can actually see.
  std::vector<std::uint32_t> hits(event.particles.size(), 0);
  const auto count = [&](const Acts::SpacePointContainer& spacePoints) {
    const auto particleColumn = spacePoints.column<std::uint32_t>("particleId");
    for (const auto sp : spacePoints) {
      ++hits[sp.extra(particleColumn)];
    }
  };
  if (selection != SpacePointSelection::Strip) {
    count(event.spacePoints);
  }
  if (selection != SpacePointSelection::Pixel) {
    count(event.stripSpacePoints);
  }

  EventSummary summary;
  summary.stripSpacePoints = event.stripSpacePoints.size();
  summary.spacePoints = event.spacePoints.size() + summary.stripSpacePoints;
  for (std::size_t index = 0; index < event.particles.size(); ++index) {
    const GeneratedParticle& particle = event.particles[index];
    if (particle.primary()) {
      ++summary.primaries;
      summary.primaryHits += hits[index];
      if (particle.pt >= ptThreshold && hits[index] >= 3) {
        ++summary.seedablePrimaries;
      }
    } else {
      ++summary.secondaries;
      summary.secondaryHits += hits[index];
    }
  }
  return summary;
}

Synthetic::SeedingSummary Synthetic::evaluateSeeds(
    const Event& event, const Acts::SeedContainer& seeds,
    const float ptThreshold, const std::size_t minTrueSpacePoints) {
  // whichever selection the seeder was run on, see `selectSpacePoints`; every
  // one of them indexes the same particles
  const Acts::SpacePointContainer& spacePoints = seeds.spacePointContainer();
  const auto particleColumn = spacePoints.column<std::uint32_t>("particleId");

  SeedingSummary summary;
  summary.seeds = seeds.size();
  std::vector<bool> matched(event.particles.size(), false);
  for (const auto seed : seeds) {
    const auto indices = seed.spacePointIndices();
    if (indices.empty() || indices.size() < minTrueSpacePoints) {
      continue;
    }
    // The particle contributing the most space points. Scanning the seed
    // against itself beats a map, a seed holding a handful of them.
    std::uint32_t particle = 0;
    std::size_t best = 0;
    for (const Acts::SpacePointIndex index : indices) {
      const std::uint32_t candidate = spacePoints[index].extra(particleColumn);
      const auto count = static_cast<std::size_t>(
          std::ranges::count_if(indices, [&](Acts::SpacePointIndex other) {
            return spacePoints[other].extra(particleColumn) == candidate;
          }));
      if (count > best) {
        best = count;
        particle = candidate;
      }
    }
    if (best < minTrueSpacePoints) {
      continue;
    }
    if (const GeneratedParticle& info = event.particles[particle];
        !info.primary() || info.pt < ptThreshold) {
      continue;
    }
    ++summary.trueSeeds;
    matched[particle] = true;
  }
  summary.matchedPrimaries =
      static_cast<std::size_t>(std::ranges::count(matched, true));
  return summary;
}

}  // namespace ActsFatras

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>

#include "RzTrackFinderImpl.hpp"

namespace Acts::Experimental::detail::rz {

void Finder::beginWalk(const RzMeasurementAccessor& measurements,
                       const RzTrackStart& start, RzTrackCandidate& candidate,
                       Walk& walk) const {
  candidate.clear();
  walk = Walk{};
  walk.candidate = &candidate;
  State& state = walk.state;
  state.v = start.parameters;
  state.anchor = start.parameters;
  state.c = start.covariance;
  state.time = start.time;
  state.timeVariance = start.timeVariance;
  const float absCharge = m_cfg.particleHypothesis.absoluteCharge();
  state.massOverCharge =
      m_cfg.particleHypothesis.mass() / (absCharge > 0. ? absCharge : 1.);
  state.bz = m_bz;
  state.anchorBz = m_bz;

  const RzLayout& layout = *m_layout;
  for (const RzSeedMeasurement& entry : start.seedMeasurements) {
    if (entry.module == kRzNone || entry.index == kRzNone) {
      continue;
    }
    const std::uint32_t layer = layout.modules[entry.module].layer;
    if (layer != kRzNone) {
      walk.knownHits.emplace_back(layer, entry);
    }
  }
  if (start.module != kRzNone) {
    const std::uint32_t layer = layout.modules[start.module].layer;
    walk.propagation.startSurface = layout.layers[layer].surface;
    const RzSurface& surface = layout.surfaces[walk.propagation.startSurface];
    state.bz = bzAt(surface, alongCoordinate(surface, state.v), m_bz);
    state.br = surface.brAt(alongCoordinate(surface, state.v));
    state.anchorBz = state.bz;
    const RzSeedMeasurement known = walk.knownAt(layer);
    const bool took =
        known.module != kRzNone &&
        takeKnownHit(measurements, known, layer, kRzNone, state, candidate);
    ModuleList startModules;
    bool startOnModule = false;
    modulesAt(layer, state, startModules, startOnModule, candidate);
    if (!startModules.empty() &&
        searchLayer(measurements, layer, kRzNone, startModules, state,
                    candidate, took ? 1u : 0u,
                    took ? known.module : kRzNone) == 0 &&
        !took && startOnModule) {
      candidate.hits.push_back(
          {layer, kRzNone, kRzNone, kRzNone, startModules.front(), 0.});
    }
  }

  m_propagator.initialize(walk.propagation, state.v,
                          walk.propagation.startSurface);

  // Seed-layer measurements do not consume the hole budget.
  for (const RzTrackHit& hit : candidate.hits) {
    walk.holes += hit.isHole() ? 1 : 0;
  }
  walk.consecutiveHoles = walk.holes;
  walk.measurementsFound =
      static_cast<std::uint32_t>(candidate.hits.size()) - walk.holes;
  walk.lastHit = state;
}

void Finder::searchStop(const RzMeasurementAccessor& measurements,
                        const Crossing& crossing, Walk& walk) const {
  ++walk.layersCrossed;
  State& state = walk.state;
  RzTrackCandidate& candidate = *walk.candidate;
  // Try the known seed hit before searching for additional measurements.
  const RzSeedMeasurement known = walk.knownAt(crossing.layer);
  const bool took = known.module != kRzNone &&
                    takeKnownHit(measurements, known, crossing.layer,
                                 crossing.stop, state, candidate);
  // Use the same module search for measurement lookup and hole detection.
  bool onModule = false;
  modulesAt(crossing.layer, state, walk.crossedModules, onModule, candidate);
  if (walk.crossedModules.empty()) {
    if (took) {
      // Retain a known hit even if the search window contains no modules.
      ++walk.measurementsFound;
      walk.consecutiveHoles = 0;
      walk.lastHit = state;
    }
    return;
  }
  const std::uint32_t accepted =
      searchLayer(measurements, crossing.layer, crossing.stop,
                  walk.crossedModules, state, candidate, took ? 1u : 0u,
                  took ? known.module : kRzNone) +
      (took ? 1u : 0u);
  if (accepted > 0) {
    walk.measurementsFound += accepted;
    walk.consecutiveHoles = 0;
    walk.lastHit = state;
    // Apply the pT cut after a measurement updates the momentum estimate.
    if (m_cfg.ptMin > 0.) {
      const double qOverP = std::abs(state.v[eRzQOverP]);
      const double sinTheta = fastHypot(state.v[eRzDir0], state.v[eRzDir1]);
      if (qOverP > 0. && sinTheta / qOverP < m_cfg.ptMin) {
        walk.done = true;
        return;
      }
    }
  } else if (!onModule) {
    // passed between the modules: nothing was expected here
    return;
  } else {
    candidate.hits.push_back({crossing.layer, kRzNone, crossing.stop, kRzNone,
                              walk.crossedModules.front(), 0.});
    ++walk.holes;
    ++walk.consecutiveHoles;
    if (walk.holes > m_cfg.maxHoles ||
        walk.consecutiveHoles > m_cfg.maxConsecutiveHoles) {
      walk.done = true;
      return;
    }
  }
  // Stop candidates with too few measurements at the configured depth.
  if (m_cfg.layersForMinMeasurements > 0 &&
      walk.layersCrossed >= m_cfg.layersForMinMeasurements &&
      walk.measurementsFound < m_cfg.minMeasurementsAtLayer) {
    walk.done = true;
  }
}

bool Finder::finishWalk(const RzMeasurementAccessor& measurements,
                        Walk& walk) const {
  RzTrackCandidate& candidate = *walk.candidate;
  while (!candidate.hits.empty() && candidate.hits.back().isHole()) {
    candidate.hits.pop_back();
  }
  candidate.holes = 0;
  for (const RzTrackHit& hit : candidate.hits) {
    candidate.holes += hit.isHole();
  }
  candidate.time = walk.lastHit.time;
  candidate.timeVariance = walk.lastHit.timeVariance;
  candidate.parameters = walk.lastHit.v;
  candidate.covariance = walk.lastHit.c;
  if (candidate.measurements < m_cfg.minMeasurements) {
    return false;
  }
  if (m_cfg.backwardPass) {
    backwardPass(measurements, walk.lastHit, candidate);
  }
  return true;
}

void Finder::searchForward(const RzMeasurementAccessor& measurements,
                           Walk& walk) const {
  for (const Crossing& crossing : m_propagator.sensitiveCrossings(
           walk.state, walk.propagation, *walk.candidate)) {
    searchStop(measurements, crossing, walk);
    if (walk.done) {
      break;
    }
  }
}

bool Finder::findTrack(
    const RzMeasurementAccessor& measurements, const RzVector& start,
    const RzMatrix& startCovariance, std::uint32_t startModule,
    RzTrackCandidate& candidate,
    std::span<const RzSeedMeasurement> seedMeasurements) const {
  Walk walk;
  beginWalk(measurements,
            RzTrackStart{start, startCovariance, startModule, seedMeasurements},
            candidate, walk);
  searchForward(measurements, walk);
  return finishWalk(measurements, walk);
}

void Finder::findTracks(
    const RzMeasurementAccessor& measurements,
    std::span<const RzTrackStart> starts, std::size_t batch,
    const std::function<void(std::size_t, bool, RzTrackCandidate&)>& onTrack)
    const {
  const std::size_t width = std::max<std::size_t>(1, batch);
  std::vector<Walk> walks(std::min(width, starts.size()));
  std::vector<RzTrackCandidate> candidates(walks.size());
  for (std::size_t first = 0; first < starts.size(); first += width) {
    const std::size_t n = std::min(width, starts.size() - first);
    for (std::size_t i = 0; i < n; ++i) {
      beginWalk(measurements, starts[first + i], candidates[i], walks[i]);
    }
    if (n == 1) {
      searchForward(measurements, walks.front());
    } else {
      // Advance the batch to its next stops, then search those stops.
      bool any = true;
      while (any) {
        any = false;
        for (std::size_t i = 0; i < n; ++i) {
          Walk& walk = walks[i];
          walk.crossing =
              walk.done ? std::nullopt
                        : m_propagator.advance(walk.state, walk.propagation,
                                               *walk.candidate);
          walk.done = !walk.crossing.has_value();
          any = any || !walk.done;
        }
        for (std::size_t i = 0; i < n; ++i) {
          if (walks[i].crossing) {
            searchStop(measurements, *walks[i].crossing, walks[i]);
          }
        }
      }
    }
    for (std::size_t i = 0; i < n; ++i) {
      onTrack(first + i, finishWalk(measurements, walks[i]), candidates[i]);
    }
  }
}

}  // namespace Acts::Experimental::detail::rz

namespace Acts::Experimental {

RzTrackFinder::RzTrackFinder(const RzTrackFinderConfig& config,
                             const RzLayout& layout, double bz)
    : m_cfg(config), m_layout(&layout), m_bz(bz) {}

bool RzTrackFinder::findTrack(
    const RzMeasurementAccessor& measurements, const RzVector& start,
    const RzMatrix& startCovariance, std::uint32_t startModule,
    RzTrackCandidate& candidate,
    std::span<const RzSeedMeasurement> seedMeasurements) const {
  return detail::rz::Finder(m_cfg, *m_layout, m_bz)
      .findTrack(measurements, start, startCovariance, startModule, candidate,
                 seedMeasurements);
}

void RzTrackFinder::findTracks(
    const RzMeasurementAccessor& measurements,
    std::span<const RzTrackStart> starts, std::size_t batch,
    const std::function<void(std::size_t, bool, RzTrackCandidate&)>& onTrack)
    const {
  detail::rz::Finder(m_cfg, *m_layout, m_bz)
      .findTracks(measurements, starts, batch, onTrack);
}

}  // namespace Acts::Experimental

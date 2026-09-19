// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/container/small_vector.hpp>

#include "RzMeasurement.hpp"
#include "RzPropagation.hpp"

namespace Acts::Experimental::detail::rz {

using ModuleList = boost::container::small_vector<std::uint32_t, 8>;

/// State, navigation cursors and counters for one forward walk.
struct Walk {
  State state;
  /// Output state, excluding any transport beyond the last measurement.
  State lastHit;
  RzTrackCandidate* candidate{};
  /// Known seed measurements and their layers.
  boost::container::small_vector<std::pair<std::uint32_t, RzSeedMeasurement>, 8>
      knownHits;
  PropagationState propagation;
  std::uint32_t holes{};
  std::uint32_t consecutiveHoles{};
  std::uint32_t layersCrossed{};
  std::uint32_t measurementsFound{};
  /// the walk has ended, one way or another
  bool done{false};
  std::optional<Crossing> crossing;
  ModuleList crossedModules;

  RzSeedMeasurement knownAt(std::uint32_t layerIndex) const {
    for (const auto& [l, entry] : knownHits) {
      if (l == layerIndex) {
        return entry;
      }
    }
    return RzSeedMeasurement{};
  }
};

class Finder {
 public:
  Finder(const RzTrackFinderConfig& config, const RzLayout& layout, double bz)
      : m_cfg(config),
        m_layout(&layout),
        m_bz(bz),
        m_evaluator(config.maxModuleDistance,
                    config.gateFactor * config.chi2Cut, config.stripMargin),
        m_stepper(config.radialField, bz),
        m_propagator(config, layout, m_stepper) {}
  Finder(const Finder&) = delete;
  Finder& operator=(const Finder&) = delete;
  bool findTrack(
      const RzMeasurementAccessor& measurements, const RzVector& start,
      const RzMatrix& startCovariance, std::uint32_t startModule,
      RzTrackCandidate& candidate,
      std::span<const RzSeedMeasurement> seedMeasurements = {}) const;

  /// Follow starts in batches; `onTrack` receives each result in input order.
  /// A batch size of zero or one processes one track at a time.
  void findTracks(const RzMeasurementAccessor& measurements,
                  std::span<const RzTrackStart> starts, std::size_t batch,
                  const std::function<void(std::size_t, bool,
                                           RzTrackCandidate&)>& onTrack) const;

 private:
  // Forward finding.
  void searchForward(const RzMeasurementAccessor& measurements,
                     Walk& walk) const;
  void beginWalk(const RzMeasurementAccessor& measurements,
                 const RzTrackStart& start, RzTrackCandidate& candidate,
                 Walk& walk) const;
  void searchStop(const RzMeasurementAccessor& measurements,
                  const Crossing& crossing, Walk& walk) const;
  bool finishWalk(const RzMeasurementAccessor& measurements, Walk& walk) const;
  // Module and measurement search.
  void modulesAt(std::uint32_t layer, const State& state, ModuleList& modules,
                 bool& onModule, RzTrackCandidate& candidate) const;
  std::uint32_t searchLayer(const RzMeasurementAccessor& measurements,
                            std::uint32_t layer, std::uint32_t stop,
                            const ModuleList& modules, State& state,
                            RzTrackCandidate& candidate,
                            std::uint32_t skipRounds = 0,
                            std::uint32_t usedModule = kRzNone) const;
  Placed placeHit(const RzMeasurementAccessor& measurements,
                  const RzTrackHit& hit) const;
  bool takeKnownHit(const RzMeasurementAccessor& measurements,
                    const RzSeedMeasurement& seed, std::uint32_t layerIndex,
                    std::uint32_t stop, State& state,
                    RzTrackCandidate& candidate) const;

  // Filtering and backward refit.
  std::uint32_t saveForwardState(const State& state,
                                 RzTrackCandidate& candidate) const;
  void backwardPass(const RzMeasurementAccessor& measurements,
                    const State& forward, RzTrackCandidate& candidate) const;
  bool inwardSearch(const RzMeasurementAccessor& measurements, State& state,
                    RzTrackCandidate& candidate) const;

  const RzTrackFinderConfig& m_cfg;
  const RzLayout* m_layout;
  double m_bz;
  MeasurementEvaluator m_evaluator;
  Stepper m_stepper;
  Propagator m_propagator;
};

}  // namespace Acts::Experimental::detail::rz

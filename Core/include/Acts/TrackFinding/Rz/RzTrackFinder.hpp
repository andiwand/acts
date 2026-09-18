// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Kalman track finding on an RZ layout with closed-form helix transport.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/TrackFinding/Rz/RzLayout.hpp"
#include "Acts/TrackFinding/Rz/RzMeasurementGrid.hpp"
#include "Acts/TrackFinding/Rz/RzTransport.hpp"
#include "Acts/TrackFinding/Rz/RzTypes.hpp"

#include <cstdint>
#include <functional>
#include <numbers>
#include <optional>
#include <span>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts::Experimental {

struct RzTrackFinderConfig {
  /// Largest chi2 a measurement is accepted with
  double chi2Cut = 15.;
  /// Limit for the straight-line pre-gate, relative to `chi2Cut`.
  double gateFactor = 4.;
  /// Search window in units of the predicted position uncertainty
  double windowSigmas = 5.;
  /// Search window added on top, in length
  double windowMin = 1. * UnitConstants::mm;
  /// Furthest a module may sit from the RZ stop for a candidate to count
  double maxModuleDistance = 50. * UnitConstants::mm;
  /// Extra room along a strip on top of its half length
  double stripMargin = 2. * UnitConstants::mm;
  /// Edge tolerance for the hole decision.
  double moduleEdgeTolerance = 0.5 * UnitConstants::mm;
  std::uint32_t maxHoles = 3;
  std::uint32_t maxConsecutiveHoles = 2;
  std::uint32_t minMeasurements = 6;
  /// Measurements accepted on one layer, more than one for module overlaps
  std::uint32_t maxMeasurementsPerLayer = 2;
  /// Stop once the track has turned this far in the transverse plane
  double maxTurningAngle = std::numbers::pi;
  /// Minimum filtered transverse momentum after an update; zero disables it.
  double ptMin = 0.;
  /// Early hit-count stop; `layersForMinMeasurements == 0` disables it.
  std::uint32_t minMeasurementsAtLayer = 0;
  std::uint32_t layersForMinMeasurements = 0;
  bool applyMaterial = true;
  /// Search from the innermost hit to the beam line and closest approach.
  bool inwardSearch = true;
  /// Refilter inward so the inner state uses all found measurements.
  bool backwardPass = true;
  /// Apply the first-order radial-field correction at each step.
  bool radialField = true;
  /// Retain every filtered state for output. Otherwise retain only the
  /// checkpoint needed to start a partial backward pass.
  bool storeForwardStates = true;
  double backwardInflation = 100.;
  /// Number of inner measurements to refilter; zero refilters all.
  std::uint32_t backwardLayers = 6;
  /// Scale the forward q/p prior variance; zero freezes q/p during refitting.
  double backwardQOverPScale = 1.;
  ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();
};

/// One sensitive layer the track crossed: with a measurement, or a hole
struct RzTrackHit {
  std::uint32_t layer{kRzNone};
  /// Index of the measurement within its module, `kRzNone` for a hole
  std::uint32_t measurement{kRzNone};
  /// Index into `RzTrackCandidate::stopSurfaces` of the stop it was found
  /// at, `kRzNone` for the layer the track started on
  std::uint32_t stop{kRzNone};
  /// Index into `RzTrackCandidate::forwardStates` of the forward state after
  /// this measurement, `kRzNone` for a hole or when the state was not retained
  std::uint32_t forwardState{kRzNone};
  /// For a hole, the module the track crossed without leaving a measurement;
  /// for a measurement it is the one of the measurement's grid entry
  std::uint32_t module{kRzNone};
  double chi2{};

  bool isHole() const { return measurement == kRzNone; }
};

/// One of the measurements a seed is made of, as the finder addresses it.
struct RzSeedMeasurement {
  /// Index into `RzLayout::modules`
  std::uint32_t module{kRzNone};
  /// Index of the measurement within that module
  std::uint32_t index{kRzNone};
};

/// A filtered state, including the independent scalar time estimate.
struct RzTrackState {
  RzVector parameters;
  RzMatrix covariance;
  double time{};
  double timeVariance{};
  /// Accumulated time updates, used to transfer later timing information
  /// inward.
  double timeCorrection{};

  bool operator==(const RzTrackState&) const = default;
};

struct RzTrackCandidate {
  /// The state at the end of the forward pass, on the last measurement
  RzVector parameters{RzVector::Zero()};
  RzMatrix covariance{RzMatrix::Zero()};
  /// Backward-refitted inner state, at perigee after an inward search.
  RzVector innerParameters{RzVector::Zero()};
  RzMatrix innerCovariance{RzMatrix::Zero()};
  double time{};
  double timeVariance{};
  double innerTime{};
  double innerTimeVariance{};
  bool hasInner{false};
  /// Why the backward pass gave up, 0 if it did not
  std::uint32_t backwardFailure{};
  /// Measurements and holes in the order they were found, holes after the
  /// last measurement dropped
  std::vector<RzTrackHit> hits;
  std::uint32_t measurements{};
  std::uint32_t holes{};
  /// The RZ surfaces the forward pass stopped at, in order, for the backward
  /// pass to replay, and the path length from the stop before to each
  std::vector<std::uint32_t> stopSurfaces;
  std::vector<double> stopPaths;
  /// The coordinate along each stop's surface at the crossing, for its
  /// material band
  std::vector<double> stopAlong;
  /// Retained forward states: every measurement when requested for output,
  /// otherwise only the checkpoint for a partial backward pass
  std::vector<RzTrackState> forwardStates;
  double chi2{};
  /// Counters for the cost analysis
  std::uint32_t stops{};
  std::uint32_t candidatesTested{};
  /// Module and bin visits made by the layer lookup.
  std::uint32_t modulesTested{};
  std::uint32_t binsVisited{};
  /// Successful measurement evaluations, including known hits and the
  /// backward pass
  std::uint32_t exactEvaluated{};

  void clear() {
    hits.clear();
    stopSurfaces.clear();
    stopPaths.clear();
    stopAlong.clear();
    forwardStates.clear();
    measurements = 0;
    holes = 0;
    hasInner = false;
    backwardFailure = 0;
    chi2 = 0.;
    stops = 0;
    candidatesTested = 0;
    modulesTested = 0;
    binsVisited = 0;
    exactEvaluated = 0;
  }
};

/// Start state, module, and optional known seed measurements.
struct RzTrackStart {
  RzVector parameters{RzVector::Zero()};
  RzMatrix covariance{RzMatrix::Zero()};
  /// The module the start state sits on, or `kRzNone`; its layer is searched
  /// before any transport
  std::uint32_t module{kRzNone};
  /// The measurements the seed is made of. A layer that holds one is not
  /// searched: the measurement is taken.
  std::span<const RzSeedMeasurement> seedMeasurements{};
  double time{};
  double timeVariance{};
};

class RzTrackFinder {
 public:
  RzTrackFinder(const RzTrackFinderConfig& config, const RzLayout& layout,
                double bz);

  const RzTrackFinderConfig& config() const { return m_cfg; }
  /// The helix for a field value
  /// @param bz the field along z
  /// @return the helix
  RzHelix helixAt(double bz) const { return RzHelix{bz}; }
  /// The field the finder falls back on where a surface has no table
  double bz() const { return m_bz; }

  /// Follow one start state; known seed measurements are taken directly.
  /// @return true if the candidate has enough measurements.
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
  struct Pending;
  struct State;
  struct Evaluation;
  struct Prediction;
  struct Placed;
  struct Walk;
  using ModuleList = boost::container::small_vector<std::uint32_t, 8>;

  // Forward walk and navigation.
  void beginWalk(const RzMeasurementAccessor& measurements,
                 const RzTrackStart& start, RzTrackCandidate& candidate,
                 Walk& walk) const;
  bool advanceWalk(Walk& walk) const;
  void searchStop(const RzMeasurementAccessor& measurements, Walk& walk) const;
  bool finishWalk(const RzMeasurementAccessor& measurements, Walk& walk) const;
  std::optional<double> pathBackward(const RzHelix& helix, const RzVector& v,
                                     const RzSurface& surface,
                                     double guess) const;

  // Module and measurement search.
  void modulesAt(std::uint32_t layer, const State& state, ModuleList& modules,
                 bool& onModule, RzTrackCandidate& candidate) const;
  std::uint32_t searchLayer(const RzMeasurementAccessor& measurements,
                            std::uint32_t layer, std::uint32_t stop,
                            const ModuleList& modules, State& state,
                            RzTrackCandidate& candidate,
                            std::uint32_t skipRounds = 0,
                            std::uint32_t usedModule = kRzNone) const;
  Placed place(const RzModule& module, const RzMeasurement& measurement,
               const RzMeasurementFrame* frame) const;
  Placed placeHit(const RzMeasurementAccessor& measurements,
                  const RzTrackHit& hit) const;
  template <bool Cache = false>
  std::optional<Evaluation> evaluate(
      const State& state, const Placed& measurement, bool gate = true,
      bool useTime = true,
      std::optional<Prediction>* prediction = nullptr) const;
  bool takeKnownHit(const RzMeasurementAccessor& measurements,
                    const RzSeedMeasurement& seed, std::uint32_t layerIndex,
                    std::uint32_t stop, State& state,
                    RzTrackCandidate& candidate) const;

  // Filtering and material.
  void update(State& state, const Evaluation& evaluation) const;
  std::uint32_t saveForwardState(const State& state,
                                 RzTrackCandidate& candidate) const;
  void materialise(State& state, const Vector3& normal) const;
  bool applyMaterial(State& state, const MaterialSlab& slab,
                     const Vector3& normal, double direction = 1.) const;
  bool applyMaterial(State& state, const RzSurface& surface, std::int32_t band,
                     const Vector3& normal, double direction = 1.) const;
  void regainEnergy(State& state, const RzSurface& surface, std::int32_t band,
                    const Vector3& normal) const;
  void backwardPass(const RzMeasurementAccessor& measurements,
                    const State& forward, RzTrackCandidate& candidate) const;
  bool inwardSearch(const RzMeasurementAccessor& measurements, State& state,
                    RzTrackCandidate& candidate) const;

  RzTrackFinderConfig m_cfg;
  const RzLayout* m_layout{};
  double m_bz{};
};

}  // namespace Acts::Experimental

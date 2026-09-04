// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Track finding on an RZ layout: a free-frame Kalman filter walked from stop
/// to stop by closed-form helix transport, with the material of the stops
/// accumulated between updates.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/TrackFinding/Rz/RzLayout.hpp"
#include "Acts/TrackFinding/Rz/RzMeasurementGrid.hpp"
#include "Acts/TrackFinding/Rz/RzTransport.hpp"

#include <cstdint>
#include <numbers>
#include <optional>
#include <vector>

namespace Acts::Experimental {

struct RzTrackFinderConfig {
  /// Largest chi2 a measurement is accepted with
  double chi2Cut = 15.;
  /// A candidate whose chi2 along the straight line, against the covariance
  /// at the stop widened by the direction uncertainty over the module
  /// distance, exceeds this many times `chi2Cut` is dropped before the exact
  /// transport is built for it
  double gateFactor = 4.;
  /// Search window in units of the predicted position uncertainty
  double windowSigmas = 5.;
  /// Search window added on top, in length
  double windowMin = 1. * UnitConstants::mm;
  /// Furthest a module may sit from the RZ stop for a candidate to count
  double maxModuleDistance = 50. * UnitConstants::mm;
  /// Extra room along a strip on top of its half length
  double stripMargin = 2. * UnitConstants::mm;
  /// Room around a module's edge within which a crossing still counts as on
  /// the module, for the hole decision
  double moduleEdgeTolerance = 0.5 * UnitConstants::mm;
  std::uint32_t maxHoles = 3;
  std::uint32_t maxConsecutiveHoles = 2;
  std::uint32_t minMeasurements = 6;
  /// Measurements accepted on one layer, more than one for module overlaps
  std::uint32_t maxMeasurementsPerLayer = 2;
  /// Stop once the track has turned this far in the transverse plane
  double maxTurningAngle = std::numbers::pi;
  bool applyMaterial = true;
  /// Refilter the found measurements backwards from the forward result, so
  /// that the parameters at the inner end carry every hit's information: a
  /// filter run the other way, started from the diagonal of the forward
  /// covariance inflated by `backwardInflation`, i.e. from nothing.
  bool backwardPass = true;
  double backwardInflation = 100.;
  /// Run the backward pass over the innermost this many measurements only,
  /// from the forward state at the outermost of them, with the forward
  /// filter's final q/p and its variance as the prior: the impact parameters
  /// come from the inner hits, the momentum from the whole track. Zero runs
  /// it over every measurement.
  std::uint32_t backwardLayers = 6;
  /// In a partial backward pass, what the forward q/p variance is scaled by
  /// as the prior: the inner hits already went into it, so 1 double-counts
  /// them; 0 freezes q/p at the forward value and restores its variance after
  double backwardQOverPScale = 1.;
  ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();
};

/// One sensitive layer the track crossed: with a measurement, or a hole
struct RzTrackHit {
  std::uint32_t layer{kRzNone};
  /// Index into the measurement grid, `kRzNone` for a hole
  std::uint32_t measurement{kRzNone};
  /// Index into `RzTrackCandidate::stopSurfaces` of the stop it was found
  /// at, `kRzNone` for the layer the track started on
  std::uint32_t stop{kRzNone};
  /// Index into `RzTrackCandidate::forwardStates` of the forward state after
  /// this measurement, `kRzNone` for a hole
  std::uint32_t forwardState{kRzNone};
  double chi2{};

  bool isHole() const { return measurement == kRzNone; }
};

struct RzTrackCandidate {
  /// The state at the end of the forward pass, on the last measurement
  RzVector parameters{RzVector::Zero()};
  RzMatrix covariance{RzMatrix::Zero()};
  /// The state at the first measurement after the backward pass, if run
  RzVector innerParameters{RzVector::Zero()};
  RzMatrix innerCovariance{RzMatrix::Zero()};
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
  /// The forward state and covariance after each measurement, for a backward
  /// pass that starts part way in
  std::vector<std::pair<RzVector, RzMatrix>> forwardStates;
  double chi2{};
  double pathLength{};
  /// Counters for the cost analysis
  std::uint32_t stops{};
  std::uint32_t candidatesTested{};

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
    pathLength = 0.;
    stops = 0;
    candidatesTested = 0;
  }
};

class RzTrackFinder {
 public:
  RzTrackFinder(const RzTrackFinderConfig& config, const RzLayout& layout,
                double bz);

  const RzTrackFinderConfig& config() const { return m_cfg; }
  const RzHelix& helix() const { return m_helix; }

  /// Follow a track outward from a start state.
  /// @param grid the measurements of the event
  /// @param start the start state
  /// @param startCovariance its covariance
  /// @param startModule the module the start state sits on, or `kRzNone`; its
  ///        layer is searched before any transport
  /// @param candidate the result, cleared first
  /// @return true if the candidate has at least `minMeasurements` hits
  bool findTrack(const RzMeasurementGrid& grid, const RzVector& start,
                 const RzMatrix& startCovariance, std::uint32_t startModule,
                 RzTrackCandidate& candidate) const;

 private:
  /// The scalars multiple scattering and energy loss straggling accumulate in
  /// between two materialisations into the covariance
  struct Pending {
    double varAngle{};
    double varPosition{};
    double covAnglePosition{};
    double varQOverP{};

    bool empty() const { return varAngle == 0. && varQOverP == 0.; }
    void advance(double s) {
      varPosition += 2. * covAnglePosition * s + varAngle * s * s;
      covAnglePosition += varAngle * s;
    }
  };

  struct State {
    RzVector v;
    RzMatrix c;
    Pending pending;
    double turned{};
    /// Where the covariance sits: the state it was last moved to or updated
    /// at, and the path walked since. The covariance is moved once per
    /// sensitive stop by the Jacobian of that whole path, one helix step
    /// from the anchor: the passive stops in between only re-parametrise
    /// the same map, and the energy loss at them changes the curvature by
    /// too little to matter for a Jacobian.
    RzVector anchor;
    double pathSince{};

    /// Walk on without the covariance
    void travel(double s) { pathSince += s; }
    /// Bring the covariance to the state, on a surface with the given normal
    void moveCovariance(const RzHelix& helix, const Vector3& normal) {
      if (pathSince == 0.) {
        return;
      }
      c = helix.stepJacobianOnto(anchor, pathSince, v, normal).transport(c);
      anchor = v;
      pathSince = 0.;
    }
  };

  /// What evaluating a measurement against a state yields: the residual on
  /// the module, and the measurement pulled back to where the state is, as
  /// the rows of `H J` with `J` the transport to the module. The update then
  /// happens at the state's own stop, which for the linear model is the same
  /// as updating on the module and costs no covariance transport.
  struct Evaluation {
    double chi2{};
    /// `C (H J)^T`, one column per measured coordinate
    Eigen::Matrix<double, eRzSize, 2> ch;
    Eigen::Matrix<double, 2, 1> residual;
    Eigen::Matrix<double, 2, 2> sInv;
  };

  /// Take the residual of a measurement against the state brought to its
  /// module, and its chi2, with the state's covariance as the prediction's
  /// @param gate drop the measurement on the straight-line chi2 first; off
  ///        for a measurement the track is known to have
  /// @return nothing if the module cannot be reached or the strip is missed
  std::optional<Evaluation> evaluate(const State& state, const RzMeasurement& m,
                                     bool gate = true) const;

  /// Kalman update with an evaluated measurement, at the state's stop
  void update(State& state, const Evaluation& e) const;

  /// Apply the material of a stop: energy loss to the state now, scattering
  /// and straggling to `pending`
  /// @param direction +1 along the track, -1 against it (energy is regained)
  /// @return false if the track ranged out
  bool applyMaterial(State& state, const MaterialSlab& slab,
                     const Vector3& normal, double direction = 1.) const;

  /// The same from a band's table, the formulas as the fallback
  bool applyMaterial(State& state, const RzSurface& surface, int band,
                     const Vector3& normal, double direction = 1.) const;

  /// Refilter the candidate's measurements from the outer end inwards
  void backwardPass(const RzMeasurementGrid& grid, const State& forward,
                    RzTrackCandidate& candidate) const;

  /// Add `pending` to the covariance and project the position part onto the
  /// surface with the given normal
  void materialise(State& state, const Vector3& normal) const;

  /// Search a layer around the state and update with the best candidates
  /// @param stop the stop the layer is at, `kRzNone` for the start layer
  /// @return the number of measurements accepted
  std::uint32_t searchLayer(const RzMeasurementGrid& grid, std::uint32_t layer,
                            std::uint32_t stop, State& state,
                            RzTrackCandidate& candidate) const;

  /// Whether the state, brought to a module plane, lands on a module of the
  /// layer: a layer crossing without a measurement is a hole only then
  bool onModule(std::uint32_t layer, const State& state) const;

  /// Path length back to an RZ surface, negative, or nothing
  /// @param guess where to start looking, the forward path with its sign
  ///        flipped
  std::optional<double> pathBackward(const RzVector& v,
                                     const RzSurface& surface,
                                     double guess) const;

  RzTrackFinderConfig m_cfg;
  const RzLayout* m_layout{};
  RzHelix m_helix;
};

}  // namespace Acts::Experimental

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
#include <span>
#include <vector>

#include <boost/container/static_vector.hpp>

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
  /// Branch stopper: drop a candidate whose filtered transverse momentum has
  /// fallen below this, checked after each measurement. The seed's own
  /// estimate is not trusted for it — the first update has to have happened.
  /// Zero switches it off.
  double ptMin = 0.;
  /// A candidate that has not reached this many measurements by the time it
  /// has crossed this many sensitive layers is dropped. Zero switches it off.
  std::uint32_t minMeasurementsAtLayer = 0;
  std::uint32_t layersForMinMeasurements = 0;
  bool applyMaterial = true;
  /// Carry on inward past the innermost measurement, searching the layers
  /// between it and the beam line, and end at the closest approach. Without
  /// it a seed built from outer space points yields a track that starts where
  /// the seed did and never sees anything inside it, which both loses those
  /// hits and stops deduplication recognising the seeds that would have found
  /// them.
  bool inwardSearch = true;
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
  /// Index of the measurement within its module, `kRzNone` for a hole
  std::uint32_t measurement{kRzNone};
  /// Index into `RzTrackCandidate::stopSurfaces` of the stop it was found
  /// at, `kRzNone` for the layer the track started on
  std::uint32_t stop{kRzNone};
  /// Index into `RzTrackCandidate::forwardStates` of the forward state after
  /// this measurement, `kRzNone` for a hole
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

struct RzTrackCandidate {
  /// The state at the end of the forward pass, on the last measurement
  RzVector parameters{RzVector::Zero()};
  RzMatrix covariance{RzMatrix::Zero()};
  /// The state at the inner end after the backward pass, if run: at the
  /// closest approach to the beam axis when the inward search ran, otherwise
  /// at the first measurement
  RzVector innerParameters{RzVector::Zero()};
  RzMatrix innerCovariance{RzMatrix::Zero()};
  bool hasInner{false};
  /// Whether `innerParameters` are already at the closest approach, so that a
  /// caller has nothing left to extrapolate
  bool innerAtPerigee{false};
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
  /// Modules the layer lookups examined, and the bins they walked to find
  /// them: what the measurement binning costs per stop
  std::uint32_t modulesTested{};
  std::uint32_t binsVisited{};
  /// Candidates on polar modules, which the search cannot gate on the
  /// module's axes and so places and evaluates one by one
  std::uint32_t polarTested{};
  /// Exact transports made: the gate's pick per round, the known hits, the
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
    innerAtPerigee = false;
    backwardFailure = 0;
    chi2 = 0.;
    pathLength = 0.;
    stops = 0;
    candidatesTested = 0;
    modulesTested = 0;
    binsVisited = 0;
    polarTested = 0;
    exactEvaluated = 0;
  }
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

  /// Follow a track outward from a start state.
  /// @param measurements where to get a module's measurements
  /// @param start the start state
  /// @param startCovariance its covariance
  /// @param startModule the module the start state sits on, or `kRzNone`; its
  ///        layer is searched before any transport
  /// @param candidate the result, cleared first
  /// @return true if the candidate has at least `minMeasurements` hits
  /// @param seedMeasurements the measurements the seed is made of. A layer
  ///        that holds one is not searched: the measurement is taken.
  bool findTrack(
      const RzMeasurementAccessor& measurements, const RzVector& start,
      const RzMatrix& startCovariance, std::uint32_t startModule,
      RzTrackCandidate& candidate,
      std::span<const RzSeedMeasurement> seedMeasurements = {}) const;

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
    /// `Bz` the state moves in from here, and the one at the anchor
    double bz{};
    double anchorBz{};

    /// Walk on without the covariance
    void travel(double s) { pathSince += s; }
    /// Bring the covariance to the state, on a surface with the given normal
    void moveCovariance(const RzHelix& helix, const Vector3& normal) {
      if (pathSince == 0.) {
        return;
      }
      c = helix.stepJacobianOnto(anchor, pathSince, v, normal).transport(c);
      anchor = v;
      anchorBz = bz;
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

  /// A measurement placed in the global frame, as the exact transport needs
  /// it. The search holds measurements on their module's axes, which is all
  /// the gate reads; this is what the few that survive the gate are expanded
  /// into.
  struct Placed {
    Vector3 position{Vector3::Zero()};
    /// The direction the measured coordinate is taken along
    Vector3 u{Vector3::Zero()};
    /// The other one, which a strip does not measure
    Vector3 v{Vector3::Zero()};
    Vector3 normal{Vector3::Zero()};
    /// Variance along `u`
    double cov00{};
    double cov01{};
    /// Variance along `v`, unused by a strip
    double cov11{};
    double invLever{};
    /// Room along `v`, the coordinate a strip does not measure
    double halfV{};
    /// How far from the RZ stop the module may be met
    double maxDistance{};
    bool pixel{};
  };

  /// Place a hit's measurement in the global frame
  /// @param measurements where to get a module's measurements
  /// @param hit the hit, which names its module and the index within it
  /// @return the placed measurement
  Placed placeHit(const RzMeasurementAccessor& measurements,
                  const RzTrackHit& hit) const;

  /// Place a measurement in the global frame
  /// @param mod the module it sits on
  /// @param m the measurement
  /// @param frame its own axes, or nullptr for a cartesian module, whose
  ///        measurements all share the module's
  /// @return the placed measurement
  Placed place(const RzModule& mod, const RzMeasurement& m,
               const RzMeasurementFrame* frame) const;

  /// Take the residual of a measurement against the state brought to its
  /// module, and its chi2, with the state's covariance as the prediction's
  /// @param gate drop the measurement on the straight-line chi2 first; off
  ///        for a measurement the track is known to have
  /// @return nothing if the module cannot be reached or the strip is missed
  std::optional<Evaluation> evaluate(const State& state, const Placed& m,
                                     bool gate = true) const;

  /// Take a measurement the caller says the track is made of, without
  /// searching the layer for it. The seed's own measurements are known, and
  /// looking for them costs a module window opened against the seed's
  /// covariance - the widest the track ever has - and a full transport of
  /// every measurement the crossed modules carry.
  /// @return true if the measurement could be brought onto the track
  bool takeKnownHit(const RzMeasurementAccessor& measurements,
                    const RzSeedMeasurement& seed, std::uint32_t layerIndex,
                    std::uint32_t stop, State& state,
                    RzTrackCandidate& candidate) const;

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

  /// Give the state back the mean energy a stop took from it, and nothing
  /// else: the scattering and straggling of a stop the backward pass skips
  /// are not wanted, only the momentum the track had there
  void regainEnergy(State& state, const RzSurface& surface, int band,
                    const Vector3& normal) const;

  /// Refilter the candidate's measurements from the outer end inwards
  void backwardPass(const RzMeasurementAccessor& measurements,
                    const State& forward, RzTrackCandidate& candidate) const;

  /// Walk inward from the state, searching every sensitive layer between it
  /// and the beam line, and leave the state at the closest approach.
  /// @param measurements where to get a module's measurements
  /// @param state the state at the innermost measurement, moved to the
  ///        closest approach
  /// @param candidate the hits found are appended, outward to inward
  /// @return false if the closest approach could not be reached
  bool inwardSearch(const RzMeasurementAccessor& measurements, State& state,
                    RzTrackCandidate& candidate) const;

  /// Add `pending` to the covariance and project the position part onto the
  /// surface with the given normal
  void materialise(State& state, const Vector3& normal) const;

  /// How many modules of one layer a crossing can land on at once: a stereo
  /// pair, an overlap in phi or along, and room to spare
  static constexpr std::size_t kMaxModulesPerLayer = 8;
  using ModuleList =
      boost::container::static_vector<std::uint32_t, kMaxModulesPerLayer>;

  /// Search the modules the state crosses and update with the best candidates
  /// @param measurements where to get a module's measurements
  /// @param layer the layer the modules belong to
  /// @param stop the stop the layer is at, `kRzNone` for the start layer
  /// @param modules the modules the crossing landed on, from `modulesAt`
  /// @return the number of measurements accepted
  std::uint32_t searchLayer(const RzMeasurementAccessor& measurements,
                            std::uint32_t layer, std::uint32_t stop,
                            const ModuleList& modules, State& state,
                            RzTrackCandidate& candidate,
                            std::uint32_t skipRounds = 0,
                            std::uint32_t usedModule = kRzNone) const;

  /// The modules of the layer the state could have crossed, widened by where
  /// the state could be, which is what the search has to look at.
  /// @param layer the layer
  /// @param state the state at the stop
  /// @param modules filled with the modules, cleared first
  /// @param onModule set if the crossing lands on a module without the
  ///        widening — the hole decision, which must stay as tight as it was
  ///        or a track that merely passed near a module counts as having
  ///        missed one
  /// @param candidate its counters take what the lookup examined
  void modulesAt(std::uint32_t layer, const State& state, ModuleList& modules,
                 bool& onModule, RzTrackCandidate& candidate) const;

  /// Path length back to an RZ surface, negative, or nothing
  /// @param guess where to start looking, the forward path with its sign
  ///        flipped
  std::optional<double> pathBackward(const RzHelix& helix, const RzVector& v,
                                     const RzSurface& surface,
                                     double guess) const;

  RzTrackFinderConfig m_cfg;
  const RzLayout* m_layout{};
  double m_bz{};
};

}  // namespace Acts::Experimental

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFinding/Rz/RzTrackFinder.hpp"

#include <cmath>
#include <iterator>
#include <optional>

namespace Acts::Experimental::detail::rz {

/// Process noise accumulated between covariance materialisations.
struct Pending {
  double varAngle{};
  double varPosition{};
  double covAnglePosition{};
  double varQOverP{};

  bool empty() const { return varAngle == 0. && varQOverP == 0.; }
  // Signed path: position-direction correlations reverse when walking inward.
  void advance(double s) {
    varPosition += 2. * covAnglePosition * s + varAngle * s * s;
    covAnglePosition += varAngle * s;
  }
};

struct State {
  RzVector v;
  RzMatrix c;
  double time{};
  double timeVariance{};
  double timeCorrection{};
  double massOverCharge{};
  Pending pending;
  double turned{};
  /// Covariance anchor and path since its last transport.
  RzVector anchor;
  double pathSince{};
  /// `Bz` the state moves in from here, and the one at the anchor
  double bz{};
  double anchorBz{};
  /// The radial field where the state stands
  double br{};
  /// Accumulated radial-field correction to the Jacobian q/p column.
  Vector3 brQopPosition{Vector3::Zero()};
  Vector3 brQopDirection{Vector3::Zero()};

  /// Walk on without the covariance
  void travel(double s) {
    pathSince += s;
    const double mOverP = massOverCharge * v[eRzQOverP];
    time += s * std::sqrt(1. + mOverP * mOverP);
  }
  /// Bring the covariance to the state, on a surface with the given normal
  void moveCovariance(const RzHelix& helix, const Vector3& normal) {
    if (pathSince == 0.) {
      return;
    }
    c = helix
            .stepJacobianOnto(anchor, pathSince, v, normal,
                              detail::stepTrig(helix.kappa(anchor) * pathSince),
                              brQopPosition, brQopDirection)
            .transport(c);
    anchor = v;
    anchorBz = bz;
    pathSince = 0.;
    brQopPosition.setZero();
    brQopDirection.setZero();
  }
};

struct NavigationState {
  /// Next outward cylinder and next disc along z.
  std::size_t cyl{};
  std::int32_t disc{};
  std::int32_t discStep{1};
  bool cylindersLeft{true};
  bool discsLeft{true};
  /// Last stop type; probe the other type only when it could be nearer.
  bool inEndcap{false};
  /// Cylinder intersection cached until the state moves.
  std::uint32_t cylCached{kRzNone};
  std::optional<double> cylCachedPath;
};

enum class PropagationStatus : std::uint8_t {
  Active,
  NoTarget,
  TurningLimit,
  Escape,
  MaterialFailure,
};

struct PropagationState {
  NavigationState navigation;
  std::uint32_t startSurface{kRzNone};
  PropagationStatus status{PropagationStatus::Active};
};

struct Target {
  std::uint32_t surface;
  double path;
  bool cylinder;
};

struct Crossing {
  std::uint32_t layer;
  std::uint32_t stop;
};

/// Geometry traversal; caches remain valid only until the state moves.
class Navigator {
 public:
  explicit Navigator(const RzLayout& layout) : m_layout(layout) {}
  void initialize(NavigationState& navigation, const RzVector& start) const;
  std::optional<Target> next(NavigationState& navigation, const State& state,
                             double maxPath) const;
  static void moved(NavigationState& navigation, bool cylinder) {
    navigation.cylCached = kRzNone;
    navigation.inEndcap = !cylinder;
  }

 private:
  const RzLayout& m_layout;
};

/// Parameter transport, deferred covariance transport and material effects.
class Stepper {
 public:
  Stepper(const RzTrackFinderConfig& config, double bz)
      : m_cfg(config), m_bz(bz) {}
  Vector3 land(State& state, RzVector& landed, double path,
               const RzSurface& surface, double along) const;
  void materialise(State& state, const Vector3& normal) const;
  bool applyMaterial(State& state, const MaterialSlab& slab,
                     const Vector3& normal, double direction = 1.) const;
  bool applyMaterial(State& state, const RzSurface& surface, std::int32_t band,
                     const Vector3& normal, double direction = 1.) const;
  void regainEnergy(State& state, const RzSurface& surface, std::int32_t band,
                    const Vector3& normal) const;

 private:
  const RzTrackFinderConfig& m_cfg;
  double m_bz;
};

class SensitiveCrossings;

/// Advance through passive material to the next sensitive crossing.
class Propagator {
 public:
  Propagator(const RzTrackFinderConfig& config, const RzLayout& layout,
             const Stepper& stepper)
      : m_cfg(config),
        m_layout(layout),
        m_stepper(stepper),
        m_navigator(layout) {}
  void initialize(PropagationState& propagation, const RzVector& start,
                  std::uint32_t startSurface = kRzNone) const;
  std::optional<Crossing> advance(State& state, PropagationState& propagation,
                                  RzTrackCandidate& candidate) const;

  SensitiveCrossings sensitiveCrossings(State& state,
                                        PropagationState& propagation,
                                        RzTrackCandidate& candidate) const;

 private:
  const RzTrackFinderConfig& m_cfg;
  const RzLayout& m_layout;
  const Stepper& m_stepper;
  Navigator m_navigator;
};

/// Single-pass traversal. Each increment uses the current, possibly filtered
/// state.
class SensitiveCrossings {
 public:
  SensitiveCrossings(const Propagator& propagator, State& state,
                     PropagationState& propagation, RzTrackCandidate& candidate)
      : m_propagator(propagator),
        m_state(state),
        m_propagation(propagation),
        m_candidate(candidate) {}

  struct Iterator {
    using value_type = Crossing;
    using difference_type = std::ptrdiff_t;
    using iterator_concept = std::input_iterator_tag;

    SensitiveCrossings* range{};
    const Crossing& operator*() const { return *range->m_crossing; }
    Iterator& operator++() {
      range->advance();
      return *this;
    }
    void operator++(int) { ++*this; }
    bool operator==(std::default_sentinel_t) const {
      return !range->m_crossing;
    }
  };

  Iterator begin() {
    advance();
    return {this};
  }
  std::default_sentinel_t end() const { return {}; }

 private:
  void advance() {
    m_crossing = m_propagator.advance(m_state, m_propagation, m_candidate);
  }
  const Propagator& m_propagator;
  State& m_state;
  PropagationState& m_propagation;
  RzTrackCandidate& m_candidate;
  std::optional<Crossing> m_crossing;
};

inline SensitiveCrossings Propagator::sensitiveCrossings(
    State& state, PropagationState& propagation,
    RzTrackCandidate& candidate) const {
  return {*this, state, propagation, candidate};
}

Vector3 surfaceNormal(const RzSurface& surface, const RzVector& state);
double alongCoordinate(const RzSurface& surface, const RzVector& state);
double bzAt(const RzSurface& surface, double along, double fallback);

}  // namespace Acts::Experimental::detail::rz

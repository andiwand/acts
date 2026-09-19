// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFinding/Rz/RzTrackFinder.hpp"

#include <iterator>
#include <optional>

#include "RzState.hpp"

namespace Acts::Experimental::detail::rz {

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

/// Parameter transport and field updates.
class Stepper {
 public:
  Stepper(bool radialField, double bz) : m_radialField(radialField), m_bz(bz) {}
  Vector3 land(State& state, RzVector& landed, double path,
               const RzSurface& surface, double along) const;

 private:
  bool m_radialField;
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

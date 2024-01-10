// This file is part of the Acts project.
//
// Copyright (C) 2018-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/PropagatorOptions.hpp"

namespace Acts {

class Surface;

/// @brief The void navigator struct as a default navigator
///
/// It does not provide any navigation action, the compiler
/// should eventually optimise that the function call is not done
///
struct VoidNavigator {
  struct Options : public NavigatorPlainOptions {
    /// @brief Set the plain options
    ///
    /// @param pOptions The plain options
    void setPlainOptions(const NavigatorPlainOptions& pOptions) {
      // TODO is this safe?
      static_cast<NavigatorPlainOptions&>(*this) = pOptions;
    }
  };

  /// @brief Nested State struct, minimal requirement
  struct State {
    /// Indicator if the target is reached
    bool targetReached = false;

    /// Navigation state : a break has been detected
    bool navigationBreak = false;
  };

  /// Unique typedef to publish to the Propagator
  using state_type = State;

  State makeState(const Options& /*options*/) const {
    State result;
    return result;
  }

  const Surface* currentSurface(const State& /*state*/) const {
    return nullptr;
  }

  bool targetReached(const State& state) const { return state.targetReached; }

  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

  void currentSurface(State& /*state*/, const Surface* /*surface*/) const {}

  void targetReached(State& state, bool targetReached) const {
    state.targetReached = targetReached;
  }

  void navigationBreak(State& state, bool navigationBreak) const {
    state.navigationBreak = navigationBreak;
  }

  /// Navigation call - void
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t Type of the Stepper
  ///
  /// Empty call, compiler should optimise that
  template <typename propagator_state_t, typename stepper_t>
  void initialize(propagator_state_t& /*state*/,
                  const stepper_t& /*stepper*/) const {}

  /// Navigation call - void
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t Type of the Stepper
  ///
  /// Empty call, compiler should optimise that
  template <typename propagator_state_t, typename stepper_t>
  void preStep(propagator_state_t& /*state*/,
               const stepper_t& /*stepper*/) const {}

  /// Navigation call - void
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t Type of the Stepper
  ///
  /// Empty call, compiler should optimise that
  template <typename propagator_state_t, typename stepper_t>
  void postStep(propagator_state_t& /*state*/,
                const stepper_t& /*stepper*/) const {}
};

}  // namespace Acts

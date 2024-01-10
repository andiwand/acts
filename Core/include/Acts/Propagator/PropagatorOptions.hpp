// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"

namespace Acts {

class Surface;

/// @brief Class holding the trivial options for stepping
///
struct StepperPlainOptions {
  /// Tolerance for the error of the integration
  double stepTolerance = 1e-4;

  /// Cut-off value for the step size
  double stepSizeCutOff = 0.;

  /// Absolute maximum step size
  double maxStepSize = std::numeric_limits<double>::max();

  /// Maximum number of Runge-Kutta steps for the stepper step call
  unsigned int maxRungeKuttaStepTrials = 10000;
};

/// @brief Class holding the trivial options for navigation
///
struct NavigatorPlainOptions {
  /// Required tolerance to reach surface
  double surfaceTolerance = s_onSurfaceTolerance;

  /// The start surface
  const Surface* startSurface = nullptr;

  /// The target surface
  const Surface* targetSurface = nullptr;
};

/// @brief Class holding the trivial options for propagation
///
struct PropagatorPlainOptions {
  /// Propagation direction
  Direction direction = Direction::Forward;

  /// Maximum number of steps for one propagate call
  unsigned int maxSteps = 1000;

  /// Absolute maximum path length
  double pathLimit = std::numeric_limits<double>::max();

  /// Loop protection step, it adapts the pathLimit
  bool loopProtection = true;
  double loopFraction = 0.5;  ///< Allowed loop fraction, 1 is a full loop

  /// Configurations for Stepper
  StepperPlainOptions stepper;

  /// Configurations for Navigator
  NavigatorPlainOptions navigator;
};

/// @brief Options for propagate() call
///
/// @tparam action_list_t List of action types called after each
///    propagation step with the current propagation and stepper state
///
/// @tparam aborter_list_t List of abort conditions tested after each
///    propagation step using the current propagation and stepper state
///
template <typename propagator_t, typename action_list_t,
          typename aborter_list_t>
struct PropagatorOptions : public PropagatorPlainOptions {
  using action_list_type = action_list_t;
  using aborter_list_type = aborter_list_t;

  /// Delete default constructor
  PropagatorOptions() = delete;

  /// PropagatorOptions copy constructor
  PropagatorOptions(const PropagatorOptions& po) = default;

  /// PropagatorOptions copy except for aborters
  template <typename other_action_list_t, typename other_aborter_list_t>
  PropagatorOptions(const PropagatorOptions<propagator_t, other_action_list_t,
                                            other_aborter_list_t>& po)
      : PropagatorPlainOptions(po),
        geoContext{po.geoContext},
        magFieldContext{po.magFieldContext} {}

  /// PropagatorOptions with context
  PropagatorOptions(const GeometryContext& gctx,
                    const MagneticFieldContext& mctx)
      : geoContext(gctx), magFieldContext(mctx) {}

  /// @brief Expand the Options with extended aborters
  ///
  /// @tparam extended_aborter_list_t Type of the new aborter list
  ///
  /// @param aborters The new aborter list to be used (internally)
  template <typename extended_aborter_list_t>
  PropagatorOptions<propagator_t, action_list_t, extended_aborter_list_t>
  extend(extended_aborter_list_t aborters) const {
    PropagatorOptions<propagator_t, action_list_t, extended_aborter_list_t>
        eoptions(*this);
    eoptions.actionList = actionList;
    eoptions.abortList = std::move(aborters);
    return eoptions;
  }

  /// @brief Set the plain options
  ///
  /// @param pOptions The plain options
  void setPlainOptions(const PropagatorPlainOptions& pOptions) {
    // TODO is this safe?
    static_cast<PropagatorPlainOptions&>(*this) = pOptions;

    stepper.setPlainOptions(pOptions.stepper);
    navigator.setPlainOptions(pOptions.navigator);
  }

  /// The context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;

  /// The context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;

  /// Stepper options
  typename propagator_t::Stepper::Options stepper;

  /// Navigator options
  typename propagator_t::Navigator::Options navigator;

  /// List of actions
  action_list_t actionList;

  /// List of abort conditions
  aborter_list_t abortList;
};

}  // namespace Acts

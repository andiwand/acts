// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Propagator/PropagatorOptions.hpp"
#include "Acts/Propagator/detail/GenericDenseEnvironmentExtension.hpp"

namespace Acts {

/// @brief A typedef for the default GenericDenseEnvironmentExtension with
/// double.
using DenseEnvironmentExtension =
    detail::GenericDenseEnvironmentExtension<double>;

template <typename propagator_t, typename action_list_t,
          typename aborter_list_t>
struct DenseStepperPropagatorOptions
    : public PropagatorOptions<propagator_t, action_list_t, aborter_list_t> {
  /// Delete default constructor
  DenseStepperPropagatorOptions() = delete;

  /// PropagatorOptions copy except for aborters
  template <typename other_action_list_t, typename other_aborter_list_t>
  DenseStepperPropagatorOptions(
      const DenseStepperPropagatorOptions<propagator_t, other_action_list_t,
                                          other_aborter_list_t>& po)
      : PropagatorOptions<propagator_t, action_list_t, aborter_list_t>(po),
        meanEnergyLoss{po.meanEnergyLoss},
        includeGradient{po.includeGgradient},
        momentumCutOff{po.momentumCutOff} {}

  /// Copy Constructor
  DenseStepperPropagatorOptions(const DenseStepperPropagatorOptions& dspo) =
      default;

  /// Constructor with GeometryContext
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param mctx The current magnetic fielc context object
  DenseStepperPropagatorOptions(const GeometryContext& gctx,
                                const MagneticFieldContext& mctx)
      : PropagatorOptions<propagator_t, action_list_t, aborter_list_t>(gctx,
                                                                       mctx) {}

  /// @brief Expand the Options with extended aborters
  ///
  /// @tparam extended_aborter_list_t Type of the new aborter list
  ///
  /// @param aborters The new aborter list to be used (internally)
  template <typename extended_aborter_list_t>
  DenseStepperPropagatorOptions<propagator_t, action_list_t,
                                extended_aborter_list_t>
  extend(extended_aborter_list_t aborters) const {
    DenseStepperPropagatorOptions<propagator_t, action_list_t,
                                  extended_aborter_list_t>
        eoptions(*this);
    eoptions.actionList = actionList;
    eoptions.abortList = std::move(aborters);
    return eoptions;
  }

  using PropagatorOptions<propagator_t, action_list_t,
                          aborter_list_t>::setPlainOptions;

  using PropagatorOptions<propagator_t, action_list_t,
                          aborter_list_t>::actionList;

  /// Toggle between mean and mode evaluation of energy loss
  bool meanEnergyLoss = true;

  /// Boolean flag for inclusion of d(dEds)d(q/p) into energy loss
  bool includeGradient = true;

  /// Cut-off value for the momentum in SI units
  double momentumCutOff = 0.;
};

}  // namespace Acts

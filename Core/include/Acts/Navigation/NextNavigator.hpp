// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"

#include <iomanip>
#include <iterator>
#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/container/small_vector.hpp>

namespace Acts {

namespace Experimental {

class NextNavigator {
public:
  using Surfaces = std::vector<const Surface*>;
  using SurfaceIter = std::vector<const Surface*>::iterator;

  using NavigationSurfaces =
      boost::container::small_vector<SurfaceIntersection, 10>;
  using NavigationSurfaceIter = NavigationSurfaces::iterator;

  using NavigationLayers =
      boost::container::small_vector<LayerIntersection, 10>;
  using NavigationLayerIter = NavigationLayers::iterator;

  using NavigationBoundaries =
      boost::container::small_vector<BoundaryIntersection, 4>;
  using NavigationBoundaryIter = NavigationBoundaries::iterator;

  using ExternalSurfaces = std::multimap<uint64_t, GeometryIdentifier>;

  struct Config {
    /// Tracking Geometry for this Navigator
    std::shared_ptr<const TrackingGeometry> trackingGeometry;

    /// Configuration for this Navigator
    /// stop at every sensitive surface (whether it has material or not)
    bool resolveSensitive = true;
    /// stop at every material surface (whether it is passive or not)
    bool resolveMaterial = true;
    /// stop at every surface regardless what it is
    bool resolvePassive = false;

    /// The tolerance used to defined "reached"
    double tolerance = s_onSurfaceTolerance;
  };

  /// Nested State struct
  ///
  /// It acts as an internal state which is
  /// created for every propagation/extrapolation step
  /// and keep thread-local navigation information
  struct State : public NavigationState {
    /// Navigation state - external state: the start surface
    const Surface* startSurface = nullptr;
    /// Navigation state - external state: the current surface
    const Surface* currentSurface = nullptr;
    /// Navigation state - external state: the target surface
    const Surface* targetSurface = nullptr;
    /// Indicator if the target is reached
    bool targetReached = false;
    /// Navigation state : a break has been detected
    bool navigationBreak = false;

    /// Reset state
    ///
    /// @param geoContext is the geometry context
    /// @param pos is the global position
    /// @param dir is the momentum direction
    /// @param navDir is the navigation direction
    /// @param ssurface is the new starting surface
    /// @param tsurface is the target surface
    void reset(const GeometryContext& /*geoContext*/, const Vector3& /*pos*/,
               const Vector3& /*dir*/, NavigationDirection /*navDir*/,
               const Surface* /*ssurface*/, const Surface* /*tsurface*/) {
      // keep detector because we cannot get it back right now
      auto detector = currentDetector;

      // Reset everything first
      *this = State();

      currentDetector = detector;
    }
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param _logger a logger instance
  explicit NextNavigator(Config cfg,
                     std::shared_ptr<const Logger> _logger =
                         getDefaultLogger("NextNavigator", Logging::Level::INFO))
      : m_cfg{std::move(cfg)}, m_logger{std::move(_logger)} {}

  /// @brief Navigator status call
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void status(propagator_state_t& state, const stepper_t& stepper) const {
    ACTS_VERBOSE(volInfo(state) << "Entering navigator::status.");

    auto &nState = state.navigation;
    prepareNavigatorState(state, stepper);

    // Check if the navigator is inactive
    if (inactive(state, stepper)) {
      ACTS_VERBOSE(volInfo(state) << "navigator inactive");
      return;
    }

    // (a) Pre-stepping call from propgator
    if (state.navigation.currentVolume == nullptr) {
      // Initialize and return
      initialize(state, stepper);
      return;
    }

    if (nState.currentVolume == nullptr) {
      ACTS_ERROR("panic: no current volume");
      return;
    }

    nState.currentVolume->updateNavigationState(state.geoContext, nState);
  }

  /// @brief Navigator target call
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void target(propagator_state_t& state, const stepper_t& stepper) const {
    ACTS_VERBOSE(volInfo(state) << "Entering navigator::target.");

    auto &nState = state.navigation;

    // Check if the navigator is inactive
    if (inactive(state, stepper)) {
      ACTS_VERBOSE(volInfo(state) << "navigator inactive");
      return;
    }

    // Navigator target always resets the current surface
    state.navigation.currentSurface = nullptr;

    // Return to the propagator
  }
private:
  Config m_cfg;

  std::shared_ptr<const Logger> m_logger;

  template <typename propagator_state_t>
  std::string volInfo(const propagator_state_t& state) const {
    return (state.navigation.currentVolume
                ? state.navigation.currentVolume->name()
                : "No Volume") +
           " | ";
  }

  const Logger& logger() const { return *m_logger; }

  /// Initialize call - start of propagation
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  ///
  /// @return boolean return triggers exit to stepper
  template <typename propagator_state_t, typename stepper_t>
  void initialize(propagator_state_t& state, const stepper_t& stepper) const {
    // Call the navigation helper prior to actual navigation
    ACTS_VERBOSE(volInfo(state) << "Initialization.");

    auto &nState = state.navigation;

    if (state.navigation.currentDetector == nullptr) {
      // TODO how to get this from the tracking geometry?
      ACTS_ERROR(volInfo(state) << "detector not set");
    }

    if (state.navigation.currentVolume == nullptr) {
      state.navigation.currentVolume = state.navigation.currentDetector->findDetectorVolume(state.geoContext, nState.position);
    }
  }

  /// This checks if a navigation break had been triggered or navigator
  /// is misconfigured
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t, typename stepper_t>
  bool inactive(propagator_state_t& state, const stepper_t& stepper) const {
    // Void behavior in case no tracking geometry is present
    if (!m_cfg.trackingGeometry) {
      return true;
    }
    // turn the navigator into void when you are intructed to do nothing
    if (!m_cfg.resolveSensitive && !m_cfg.resolveMaterial &&
        !m_cfg.resolvePassive) {
      return true;
    }
    return false;
  }

  template <typename propagator_state_t, typename stepper_t>
  void prepareNavigatorState(propagator_state_t& state, const stepper_t& stepper) const {
    auto &nState = state.navigation;

    nState.position = stepper.position(state.stepping);
    nState.direction = stepper.direction(state.stepping);
    nState.absMomentum = stepper.momentum(state.stepping);
    nState.charge = stepper.charge(state.stepping);
    auto fieldResult = stepper.getField(state.stepping, nState.position);
    if (!fieldResult.ok()) {
      ACTS_ERROR("could not read from the magnetic field");
    }
    nState.magneticField = *fieldResult;
  }
};

}

}

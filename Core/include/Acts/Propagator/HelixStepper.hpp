// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/PropagatorTraits.hpp"
#include "Acts/Propagator/StepperOptions.hpp"
#include "Acts/Propagator/StepperStatistics.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"

#include <optional>
#include <string>
#include <tuple>

namespace Acts {

class IVolumeMaterial;

/// @brief Stepper that moves the track along a helix in closed form
///
/// Each step reads the magnetic field at the start position and moves the
/// track along the exact helix in that field. The transport jacobian is also
/// closed form. It is the exact derivative of the helix step, so it treats
/// the field as constant over the step and has no field gradient term.
///
/// The field at the end of a step is also the field at the start of the next
/// step. The stepper uses it to estimate the position error from the field
/// change along the step, and it adapts the step size to
/// @ref StepperPlainOptions::stepTolerance. An accepted step needs no extra
/// field lookup.
///
/// The navigation estimates the distance to a surface with a straight line.
/// A step of that length along the helix passes a surface perpendicular to the
/// track by about |omega x T|^2 h^3 / 3, with omega = (q/p) B, and the next
/// step goes back. A target aborter must accept an intersection that far behind
/// the track. The default near limit of @ref SurfaceReached is 100 um.
///
/// @note The stepper propagates in vacuum only. It ignores volume material.
class HelixStepper final {
 public:
  /// Type alias for bound track parameters
  using BoundParameters = BoundTrackParameters;
  /// Type alias for jacobian matrix
  using Jacobian = BoundMatrix;
  /// Type alias for covariance matrix
  using Covariance = BoundMatrix;
  /// Bound state tuple containing parameters, Jacobian, and path length
  using BoundState = std::tuple<BoundParameters, Jacobian, double>;

  /// Configuration for the helix stepper.
  struct Config {
    /// Magnetic field provider
    std::shared_ptr<const MagneticFieldProvider> bField;
  };

  /// Runtime options for helix propagation.
  struct Options : public StepperPlainOptions {
    /// Constructor
    /// @param gctx Geometry context
    /// @param mctx Magnetic field context
    Options(const GeometryContext& gctx, const MagneticFieldContext& mctx)
        : StepperPlainOptions(gctx, mctx) {}

    /// Set plain stepper options
    /// @param options Plain stepper options
    void setPlainOptions(const StepperPlainOptions& options) {
      static_cast<StepperPlainOptions&>(*this) = options;
    }
  };

  /// @brief State for track parameter propagation
  ///
  /// It contains the stepping information and is provided thread local
  /// by the propagator
  struct State {
    /// Constructor from the options and the field cache
    ///
    /// @param [in] optionsIn is the configuration of the stepper
    /// @param [in] fieldCacheIn is the cache object for the magnetic field
    State(const Options& optionsIn, MagneticFieldProvider::Cache fieldCacheIn)
        : options(optionsIn), fieldCache(std::move(fieldCacheIn)) {}

    /// Configuration options for the stepper
    Options options;

    /// Internal free vector parameters
    FreeVector pars = FreeVector::Zero();

    /// The propagation derivative
    FreeVector derivative = FreeVector::Zero();

    /// Jacobian from local to the global frame
    BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();

    /// Pure transport jacobian part from the helix steps
    FreeMatrix jacTransport = FreeMatrix::Identity();

    /// The full jacobian of the transport entire transport
    Jacobian jacobian = Jacobian::Identity();

    /// Covariance matrix for error propagation
    Covariance cov = Covariance::Zero();

    /// Covariance matrix (and indicator)
    /// associated with the initial error on track parameters
    bool covTransport = false;

    /// Particle hypothesis
    ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();

    /// Adaptive step size of the helix steps
    ConstrainedStep stepSize;

    /// Last performed step (for overstep limit calculation)
    double previousStepSize = 0.;

    /// Magnetic field at the current position, reused as the next step's
    /// field. Reset whenever the position is set from outside.
    std::optional<Vector3> field;

    /// Accumulated path length state
    double pathAccumulated = 0.;

    /// Total number of performed steps
    std::size_t nSteps = 0;

    /// Total number of attempted steps
    std::size_t nStepTrials = 0;

    /// Statistics of the stepper
    StepperStatistics statistics;

    /// This caches the current magnetic field cell and stays
    /// (and interpolates) within it as long as this is valid.
    MagneticFieldProvider::Cache fieldCache;
  };

  /// Constructor requires knowledge of the detector's magnetic field
  /// @param bField The magnetic field provider
  explicit HelixStepper(std::shared_ptr<const MagneticFieldProvider> bField);

  /// @brief Constructor with configuration
  /// @param config The configuration of the stepper
  explicit HelixStepper(const Config& config);

  /// Create a state object
  /// @param options Stepper options
  /// @return State object
  State makeState(const Options& options) const;

  /// Initialize the state from bound track parameters
  /// @param state The state to initialize
  /// @param par The bound track parameters
  void initialize(State& state, const BoundParameters& par) const;

  /// Initialize the state from bound parameters
  /// @param state The state to initialize
  /// @param boundParams Bound track parameters vector
  /// @param cov Covariance matrix
  /// @param particleHypothesis Particle hypothesis
  /// @param surface Reference surface
  void initialize(State& state, const BoundVector& boundParams,
                  const std::optional<BoundMatrix>& cov,
                  ParticleHypothesis particleHypothesis,
                  const Surface& surface) const;

  /// Get the field for the stepping, it checks first if the access is still
  /// within the Cell, and updates the cell if necessary.
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  /// @return Magnetic field vector
  Result<Vector3> getField(State& state, const Vector3& pos) const {
    return m_bField->getField(pos, state.fieldCache);
  }

  /// Global particle position accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Position vector
  Vector3 position(const State& state) const {
    return state.pars.template segment<3>(eFreePos0);
  }

  /// Momentum direction accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Direction vector
  Vector3 direction(const State& state) const {
    return state.pars.template segment<3>(eFreeDir0);
  }

  /// QoP direction accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Charge over momentum
  double qOverP(const State& state) const { return state.pars[eFreeQOverP]; }

  /// Absolute momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Absolute momentum
  double absoluteMomentum(const State& state) const {
    return particleHypothesis(state).extractMomentum(qOverP(state));
  }

  /// Momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Momentum vector
  Vector3 momentum(const State& state) const {
    return absoluteMomentum(state) * direction(state);
  }

  /// Charge access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Particle charge
  double charge(const State& state) const {
    return particleHypothesis(state).extractCharge(qOverP(state));
  }

  /// Particle hypothesis
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Particle hypothesis
  const ParticleHypothesis& particleHypothesis(const State& state) const {
    return state.particleHypothesis;
  }

  /// Time access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Time
  double time(const State& state) const { return state.pars[eFreeTime]; }

  /// Update surface status
  ///
  /// It checks the status to the reference surface & updates
  /// the step size accordingly
  ///
  /// @param [in,out] state The stepping state (thread-local cache)
  /// @param [in] surface The surface provided
  /// @param [in] index The surface intersection index
  /// @param [in] navDir The navigation direction
  /// @param [in] boundaryTolerance The boundary check for this status update
  /// @param [in] surfaceTolerance Surface tolerance used for intersection
  /// @param [in] stype The step size type to be set
  /// @param [in] logger A @c Logger instance
  /// @return Intersection status
  IntersectionStatus updateSurfaceStatus(
      State& state, const Surface& surface, std::uint8_t index,
      Direction navDir, const BoundaryTolerance& boundaryTolerance,
      double surfaceTolerance, ConstrainedStep::Type stype,
      const Logger& logger = getDummyLogger()) const {
    return detail::updateSingleSurfaceStatus<HelixStepper>(
        *this, state, surface, index, navDir, boundaryTolerance,
        surfaceTolerance, stype, logger);
  }

  /// Update step size
  ///
  /// This method intersects the provided surface and update the navigation
  /// step estimation accordingly (hence it changes the state). It also
  /// returns the status of the intersection to trigger onSurface in case
  /// the surface is reached.
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param target [in] The NavigationTarget
  /// @param direction [in] The propagation direction
  /// @param stype [in] The step size type to be set
  void updateStepSize(State& state, const NavigationTarget& target,
                      Direction direction, ConstrainedStep::Type stype) const {
    static_cast<void>(direction);
    double stepSize = target.pathLength();
    updateStepSize(state, stepSize, stype);
  }

  /// Update step size - explicitly with a double
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param stepSize [in] The step size value
  /// @param stype [in] The step size type to be set
  void updateStepSize(State& state, double stepSize,
                      ConstrainedStep::Type stype) const {
    state.previousStepSize = state.stepSize.value();
    state.stepSize.update(stepSize, stype);
  }

  /// Get the step size
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @param stype [in] The step size type to be returned
  /// @return Step size
  double getStepSize(const State& state, ConstrainedStep::Type stype) const {
    return state.stepSize.value(stype);
  }

  /// Release the Step size
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param [in] stype The step size type to be released
  void releaseStepSize(State& state, ConstrainedStep::Type stype) const {
    state.stepSize.release(stype);
  }

  /// Output the Step Size - single component
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @return String representation of step size
  std::string outputStepSize(const State& state) const {
    return state.stepSize.toString();
  }

  /// Create and return the bound state at the current position
  ///
  /// @brief This transports (if necessary) the covariance
  /// to the surface and creates a bound state. It does not check
  /// if the transported state is at the surface, this needs to
  /// be guaranteed by the propagator
  ///
  /// @param [in] state State that will be presented as @c BoundState
  /// @param [in] surface The surface to which we bind the state
  /// @param [in] transportCov Flag steering covariance transport
  /// @param [in] freeToBoundCorrection Correction for non-linearity effect during transform from free to bound
  ///
  /// @return A bound state:
  ///   - the parameters at the surface
  ///   - the stepwise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  Result<BoundState> boundState(
      State& state, const Surface& surface, bool transportCov = true,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const;

  /// @brief If necessary fill additional members needed for curvilinearState
  ///
  /// Compute path length derivatives in case they have not been computed
  /// yet, which is the case if no step has been executed yet.
  ///
  /// @param [in, out] state State of the stepper
  /// @return true if nothing is missing after this call, false otherwise.
  bool prepareCurvilinearState(State& state) const;

  /// Create and return a curvilinear state at the current position
  ///
  /// @brief This transports (if necessary) the covariance
  /// to the current position and creates a curvilinear state.
  ///
  /// @param [in] state State that will be presented as @c CurvilinearState
  /// @param [in] transportCov Flag steering covariance transport
  ///
  /// @return A curvilinear state:
  ///   - the curvilinear parameters at given position
  ///   - the stepweise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  BoundState curvilinearState(State& state, bool transportCov = true) const;

  /// Method to update a stepper state to the some parameters
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] freeParams Free parameters that will be written into @p state
  /// @param [in] boundParams Corresponding bound parameters used to update jacToGlobal in @p state
  /// @param [in] covariance The covariance that will be written into @p state
  /// @param [in] surface The surface used to update the jacToGlobal
  void update(State& state, const FreeVector& freeParams,
              const BoundVector& boundParams, const Covariance& covariance,
              const Surface& surface) const;

  /// Method to update the stepper state
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] uposition the updated position
  /// @param [in] udirection the updated direction
  /// @param [in] qop the updated qop value
  /// @param [in] time the updated time value
  void update(State& state, const Vector3& uposition, const Vector3& udirection,
              double qop, double time) const;

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state
  ///
  /// @param [in,out] state State of the stepper
  void transportCovarianceToCurvilinear(State& state) const;

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current position,
  /// or direction of the state
  ///
  /// @param [in,out] state State of the stepper
  /// @param [in] surface is the surface to which the covariance is forwarded to
  /// @param [in] freeToBoundCorrection Correction for non-linearity effect during transform from free to bound
  /// @note no check is done if the position is actually on the surface
  /// @return Failure if the parameters cannot be expressed on the surface
  Result<void> transportCovarianceToBound(
      State& state, const Surface& surface,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const;

  /// Perform a helix track parameter propagation step
  ///
  /// @param [in,out] state State of the stepper
  /// @param propDir is the direction of propagation
  /// @param material is ignored, the stepper propagates in vacuum only
  /// @return the result of the step
  ///
  /// @note The state contains the desired step size. It can be negative during
  ///       backwards track propagation, and since the step size is adaptive,
  ///       it can be modified by the stepper class during propagation.
  Result<double> step(State& state, Direction propDir,
                      const IVolumeMaterial* material) const;

 private:
  /// Magnetic field inside of the detector
  std::shared_ptr<const MagneticFieldProvider> m_bField;
};

template <>
struct SupportsBoundParameters<HelixStepper> : public std::true_type {};

}  // namespace Acts

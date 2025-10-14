// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiComponentTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/PropagatorState.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterError.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterExtensions.hpp"
#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"
#include "Acts/TrackFitting/detail/GsfComponentMerging.hpp"
#include "Acts/TrackFitting/detail/GsfUtils.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <limits>
#include <memory>
#include <type_traits>

namespace Acts {

/// Type alias for component reducer delegate function
using ComponentReducer =
    Delegate<void(std::vector<GsfComponent>&, std::size_t, const Surface&)>;

using BetheHeitlerApprox = AtlasBetheHeitlerApprox<6, 5>;

/// Combined options for the combinatorial Kalman filter.
///
/// @tparam source_link_iterator_t Type of the source link iterator
/// @tparam track_container_t Type of the track container
template <typename track_container_t>
struct CombinatorialKalmanFilter2Options {
  /// Type alias for track state container backend
  using TrackStateContainerBackend =
      typename track_container_t::TrackStateContainerBackend;
  /// Type alias for track state proxy from the container
  using TrackStateProxy = typename track_container_t::TrackStateProxy;

  /// PropagatorOptions with context
  ///
  /// @param gctx The geometry context for this track finding/fitting
  /// @param mctx The magnetic context for this track finding/fitting
  /// @param cctx The calibration context for this track finding/fitting
  /// @param extensions_ The extension struct
  /// @param pOptions The plain propagator options
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  CombinatorialKalmanFilter2Options(
      const GeometryContext& gctx, const MagneticFieldContext& mctx,
      std::reference_wrapper<const CalibrationContext> cctx,
      CombinatorialKalmanFilterExtensions<track_container_t> extensions_,
      const PropagatorPlainOptions& pOptions, bool mScattering = true,
      bool eLoss = true)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        extensions(extensions_),
        propagatorPlainOptions(pOptions),
        multipleScattering(mScattering),
        energyLoss(eLoss) {}

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  /// The filter extensions
  CombinatorialKalmanFilterExtensions<track_container_t> extensions;

  /// The trivial propagator options
  PropagatorPlainOptions propagatorPlainOptions;

  /// The target surface
  /// @note This is useful if the filtering should be terminated at a
  ///       certain surface
  const Surface* targetSurface = nullptr;

  /// Whether to consider multiple scattering.
  bool multipleScattering = true;

  /// Whether to consider energy loss.
  bool energyLoss = true;

  /// Skip the pre propagation call. This effectively skips the first surface
  /// @note This is useful if the first surface should not be considered in a second reverse pass
  bool skipPrePropagationUpdate = false;

  /// Maximum number of components which the GSF should handle
  std::size_t maxComponents = 16;

  /// When to discard components
  double weightCutoff = 1.0e-4;

  /// Takes a vector of components and reduces its number
  ComponentReducer mixtureReducer;

  /// How to reduce the states that are stored in the multi trajectory
  ComponentMergeMethod mergeMethod = ComponentMergeMethod::eMaxWeight;

  BetheHeitlerApprox betheHeitlerApprox = makeDefaultBetheHeitlerApprox();
};

/// Combinatorial Kalman filter to find tracks.
///
/// @tparam propagator_t Type of the propagator
///
/// The CombinatorialKalmanFilter contains an Actor and a Sequencer sub-class.
/// The Sequencer has to be part of the Navigator of the Propagator in order to
/// initialize and provide the measurement surfaces.
///
/// The Actor is part of the Propagation call and does the Kalman update.
/// Updater and Calibrator are given to the Actor for further use:
/// - The Updater is the implemented kalman updater formalism, it
///   runs via a visitor pattern through the measurements.
///
/// Measurements are not required to be ordered for the
/// CombinatorialKalmanFilter, measurement ordering needs to be figured out by
/// the navigation of the propagator.
///
/// The void components are provided mainly for unit testing.
///
template <typename propagator_t, typename track_container_t>
class CombinatorialKalmanFilter2 {
 public:
  /// Constructor with propagator and logging level
  /// @param pPropagator The propagator used for the track finding
  /// @param _logger The logger for messages
  explicit CombinatorialKalmanFilter2(
      propagator_t pPropagator,
      std::unique_ptr<const Logger> _logger = getDefaultLogger("CKF",
                                                               Logging::INFO))
      : m_propagator(std::move(pPropagator)),
        m_logger(std::move(_logger)),
        m_actorLogger{m_logger->cloneWithSuffix("Actor")},
        m_updaterLogger{m_logger->cloneWithSuffix("Updater")} {}

 private:
  using TrackStateContainerBackend =
      typename track_container_t::TrackStateContainerBackend;
  using TrackProxy = typename track_container_t::TrackProxy;
  using TrackStateProxy = typename track_container_t::TrackStateProxy;

  struct TemporaryStates {
    TrackStateContainerBackend traj;
    std::vector<MultiTrajectoryTraits::IndexType> tips;
    std::map<MultiTrajectoryTraits::IndexType, double> weights;
  };

  using FiltProjector =
      detail::MultiTrajectoryProjector<detail::StatesType::eFiltered,
                                       TrackStateContainerBackend>;

  using ComponentCache = GsfComponent;

  /// The propagator for the transport and material update
  propagator_t m_propagator;

  std::unique_ptr<const Logger> m_logger;
  std::shared_ptr<const Logger> m_actorLogger;
  std::shared_ptr<const Logger> m_updaterLogger;

  const Logger& logger() const { return *m_logger; }

  /// @brief Propagator Actor plugin for the CombinatorialKalmanFilter
  ///
  /// The CombinatorialKalmanFilter Actor does not rely on the measurements to
  /// be sorted along the track.
  class Actor {
   public:
    using SingleBoundState =
        std::tuple<BoundTrackParameters, BoundMatrix, double>;
    using MultieBoundState =
        std::tuple<MultiComponentBoundTrackParameters, BoundMatrix, double>;
    /// Broadcast the result_type
    using result_type = CombinatorialKalmanFilterResult<track_container_t>;

    using BranchStopperResult = CombinatorialKalmanFilterBranchStopperResult;

    /// The target surface aborter
    SurfaceReached targetReached{std::numeric_limits<double>::lowest()};

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// Whether to consider energy loss.
    bool energyLoss = true;

    /// Skip the pre propagation call. This effectively skips the first surface
    bool skipPrePropagationUpdate = false;

    /// Calibration context for the finding run
    const CalibrationContext* calibrationContextPtr{nullptr};

    /// Maximum number of components which the GSF should handle
    std::size_t maxComponents = 16;

    /// When to discard components
    double weightCutoff = 1.0e-4;

    /// Takes a vector of components and reduces its number
    ComponentReducer mixtureReducer;

    /// How to reduce the states that are stored in the multi trajectory
    ComponentMergeMethod mergeMethod = ComponentMergeMethod::eMaxWeight;

    BetheHeitlerApprox betheHeitlerApprox = makeDefaultBetheHeitlerApprox();

    /// @brief CombinatorialKalmanFilter actor operation
    ///
    /// @tparam propagator_state_t Type of the Propagator state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper is the stepper in use
    /// @param navigator is the navigator in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void act(propagator_state_t& state, const stepper_t& stepper,
             const navigator_t& navigator, result_type& result,
             const Logger& /*logger*/) const {
      assert(result.trackStates && "No MultiTrajectory set");

      if (result.finished) {
        return;
      }

      if (state.stage == PropagatorStage::prePropagation &&
          skipPrePropagationUpdate) {
        ACTS_VERBOSE("Skip pre-propagation update (first surface)");
        return;
      }

      ACTS_VERBOSE("CombinatorialKalmanFilter step");

      assert(!result.activeBranches.empty() && "No active branches");

      // Initialize path limit reached aborter
      if (result.pathLimitReached.internalLimit ==
          std::numeric_limits<double>::max()) {
        detail::setupLoopProtection(state, stepper, result.pathLimitReached,
                                    true, logger());
      }

      // Update:
      // - Waiting for a current surface
      if (auto surface = navigator.currentSurface(state.navigation);
          surface != nullptr) {
        // There are three scenarios:
        // 1) The surface is in the measurement map
        // -> Select source links
        // -> Perform the kalman update for selected non-outlier source links
        // -> Add track states in multitrajectory. Multiple states mean branch
        // splitting.
        // -> Call branch stopper to justify each branch
        // -> If there is non-outlier state, update stepper information
        // 2) The surface is not in the measurement map but with material or is
        // an active surface
        // -> Add a hole or passive material state in multitrajectory
        // -> Call branch stopper to justify the branch
        // 3) The surface is neither in the measurement map nor with material
        // -> Do nothing
        ACTS_VERBOSE("Perform filter step");
        auto res = filter(surface, state, stepper, navigator, result);
        if (!res.ok()) {
          ACTS_ERROR("Error in filter: " << res.error());
          result.lastError = res.error();
        }

        if (result.finished) {
          return;
        }
      }

      assert(!result.activeBranches.empty() && "No active branches");

      const bool isEndOfWorldReached =
          endOfWorldReached.checkAbort(state, stepper, navigator, logger());
      const bool isVolumeConstraintReached = volumeConstraintAborter.checkAbort(
          state, stepper, navigator, logger());
      const bool isPathLimitReached = result.pathLimitReached.checkAbort(
          state, stepper, navigator, logger());
      const bool isTargetReached =
          targetReached.checkAbort(state, stepper, navigator, logger());
      if (isEndOfWorldReached || isVolumeConstraintReached ||
          isPathLimitReached || isTargetReached) {
        if (isEndOfWorldReached) {
          ACTS_VERBOSE("End of world reached");
        } else if (isVolumeConstraintReached) {
          ACTS_VERBOSE("Volume constraint reached");
        } else if (isPathLimitReached) {
          ACTS_VERBOSE("Path limit reached");
        } else if (isTargetReached) {
          ACTS_VERBOSE("Target surface reached");

          // Bind the parameter to the target surface
          auto res = stepper.boundState(state.stepping, *targetReached.surface);
          if (!res.ok()) {
            ACTS_ERROR("Error while acquiring bound state for target surface: "
                       << res.error() << " " << res.error().message());
            result.lastError = res.error();
          } else {
            const auto& [boundParams, jacobian, pathLength] = *res;
            auto currentBranch = result.activeBranches.back();
            // Assign the fitted parameters
            currentBranch.parameters() = boundParams.parameters();
            currentBranch.covariance() = *boundParams.covariance();
            currentBranch.setReferenceSurface(
                boundParams.referenceSurface().getSharedPtr());
          }

          stepper.releaseStepSize(state.stepping,
                                  ConstrainedStep::Type::Navigator);
        }

        // Record the active branch and remove it from the list
        storeLastActiveBranch(result);
        result.activeBranches.pop_back();

        // Reset propagation state to track state at next active branch
        reset(state, stepper, navigator, result);
      }
    }

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    bool checkAbort(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, const result_type& result,
                    const Logger& /*logger*/) const {
      return !result.lastError.ok() || result.finished;
    }

    /// @brief CombinatorialKalmanFilter actor operation: reset propagation
    ///
    /// @tparam propagator_state_t Type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper is the stepper in use
    /// @param navigator is the navigator in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void reset(propagator_state_t& state, const stepper_t& stepper,
               const navigator_t& navigator, result_type& result) const {
      if (result.activeBranches.empty()) {
        ACTS_VERBOSE("Stop CKF with " << result.collectedTracks.size()
                                      << " found tracks");
        result.finished = true;

        return;
      }

      auto currentBranch = result.activeBranches.back();
      auto currentState = currentBranch.outermostTrackState();

      const auto& surface = currentState.referenceSurface();

      ACTS_VERBOSE("Propagation jumps to branch with tip = "
                   << currentBranch.tipIndex());

      MultiComponentBoundTrackParameters multiParams(
          surface.getSharedPtr(), currentState.filtered(),
          currentState.filteredCovariance(),
          stepper.particleHypothesis(state.stepping));

      // Reset the stepping state
      stepper.initialize(state.stepping, multiParams);

      // Reset the navigation state
      // Set targetSurface to nullptr for forward filtering
      state.navigation.options.startSurface = &surface;
      state.navigation.options.targetSurface = nullptr;
      auto navInitRes = navigator.initialize(
          state.navigation, stepper.position(state.stepping),
          stepper.direction(state.stepping), state.options.direction);
      if (!navInitRes.ok()) {
        ACTS_ERROR("Navigation initialization failed: " << navInitRes.error());
        result.lastError = navInitRes.error();
      }

      // No Kalman filtering for the starting surface, but still need
      // to consider the material effects here
      applyMultipleScattering(navigator.currentSurface(state.navigation), state,
                              stepper, navigator,
                              MaterialUpdateStage::PostUpdate);

      // Set path limit based on loop protection
      detail::setupLoopProtection(state, stepper, result.pathLimitReached, true,
                                  logger());

      // Set path limit based on target surface
      targetReached.checkAbort(state, stepper, navigator, logger());
    }

    /// @brief CombinatorialKalmanFilter actor operation:
    /// - filtering for all measurement(s) on surface
    /// - store selected track states in multiTrajectory
    /// - update propagator state to the (last) selected track state
    ///
    /// @tparam propagator_state_t Type of the Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param result The mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    Result<void> filter(const Surface* surface, propagator_state_t& state,
                        const stepper_t& stepper, const navigator_t& navigator,
                        result_type& result) const {
      using PM = TrackStatePropMask;

      bool isSensitive = surface->associatedDetectorElement() != nullptr;
      bool hasMaterial = surface->surfaceMaterial() != nullptr;
      bool isMaterialOnly = hasMaterial && !isSensitive;
      bool expectMeasurements = isSensitive;

      if (isSensitive) {
        ACTS_VERBOSE("Measurement surface " << surface->geometryId()
                                            << " detected.");
      } else if (isMaterialOnly) {
        ACTS_VERBOSE("Material surface " << surface->geometryId()
                                         << " detected.");
      } else {
        ACTS_VERBOSE("Passive surface " << surface->geometryId()
                                        << " detected.");
        return Result<void>::success();
      }

      // TODO GSF does an early exit here if no material and no measurement is
      // found

      // Transport the covariance to the surface
      for (auto cmp : stepper.componentIterable(state.stepping)) {
        if (isMaterialOnly) {
          cmp.singleStepper(stepper).transportCovarianceToCurvilinear(
              cmp.state());
        } else {
          cmp.singleStepper(stepper).transportCovarianceToBound(cmp.state(),
                                                                *surface);
        }
      }

      if (isMaterialOnly) {
        applyMultipleScattering(surface, state, stepper, navigator,
                                MaterialUpdateStage::FullUpdate);

        TemporaryStates tmpStates;

        std::vector<ComponentCache> componentCache;

        convoluteComponents(*surface, state, stepper, componentCache);

        if (componentCache.empty()) {
          ACTS_WARNING(
              "No components left after applying energy loss. "
              "Is the weight cutoff "
              << weightCutoff << " too high?");
          ACTS_WARNING("Return to propagator without applying energy loss");
          return Result<void>::success();
        }

        // reduce component number
        const auto finalCmpNumber = std::min(
            static_cast<std::size_t>(stepper.maxComponents), maxComponents);
        mixtureReducer(componentCache, finalCmpNumber, *surface);

        removeLowWeightComponents(componentCache);

        updateStepper(state, stepper, navigator, componentCache);

        return Result<void>::success();
      }

      applyMultipleScattering(surface, state, stepper, navigator,
                              MaterialUpdateStage::PreUpdate);

      // Bind the transported state to the current surface
      auto boundStateRes = stepper.boundState(state.stepping, *surface);
      if (!boundStateRes.ok()) {
        return boundStateRes.error();
      }
      auto& boundState = *boundStateRes;
      auto& [boundParams, jacobian, pathLength] = boundState;

      auto [singleParams, singleCov] = detail::gaussianMixtureMeanCov(
          boundParams.components(),
          [](const auto& cmp) -> std::tuple<double, BoundVector, BoundMatrix> {
            return {std::get<0>(cmp), std::get<1>(cmp), *std::get<2>(cmp)};
          });
      BoundTrackParameters singleBoundParams(
          boundParams.referenceSurface().getSharedPtr(), singleParams,
          singleCov, stepper.particleHypothesis(state.stepping));
      SingleBoundState singleState(singleBoundParams, jacobian, pathLength);

      auto currentBranch = result.activeBranches.back();
      TrackIndexType prevTip = currentBranch.tipIndex();

      using TrackStatesResult = Result<CkfTypes::BranchVector<TrackIndexType>>;
      TrackStatesResult tsRes = TrackStatesResult::success({});
      if (isSensitive) {
        // extend trajectory with measurements associated to the current surface
        // which may create extra trajectory branches if more than one
        // measurement is selected.
        tsRes = m_extensions.createTrackStates(
            state.geoContext, *calibrationContextPtr, *surface, singleState,
            prevTip, result.trackStateCandidates, *result.trackStates,
            logger());
      }

      TemporaryStates tmpStates;

      if (tsRes.ok() && !(*tsRes).empty()) {
        const CkfTypes::BranchVector<TrackIndexType>& newTrackStateList =
            *tsRes;
        Result<unsigned int> procRes = processNewTrackStates(
            state, stepper, newTrackStateList, result, tmpStates);
        if (!procRes.ok()) {
          ACTS_ERROR("Processing of selected track states failed: "
                     << procRes.error());
          return procRes.error();
        }
        unsigned int nBranchesOnSurface = *procRes;

        if (nBranchesOnSurface == 0) {
          ACTS_VERBOSE("All branches on surface " << surface->geometryId()
                                                  << " have been stopped");

          reset(state, stepper, navigator, result);

          return Result<void>::success();
        }

        // `currentBranch` is invalidated after `processNewTrackStates`
        currentBranch = result.activeBranches.back();
        prevTip = currentBranch.tipIndex();
      } else {
        if (!tsRes.ok()) {
          if (static_cast<CombinatorialKalmanFilterError>(
                  tsRes.error().value()) ==
              CombinatorialKalmanFilterError::NoMeasurementExpected) {
            // recoverable error returned by track state creator
            expectMeasurements = false;
          } else {
            ACTS_ERROR("Track state creation failed on surface "
                       << surface->geometryId() << ": " << tsRes.error());
            return tsRes.error();
          }
        }

        if (expectMeasurements) {
          ACTS_VERBOSE("Detected hole after measurement selection on surface "
                       << surface->geometryId());
        }

        auto stateMask = PM::Predicted | PM::Jacobian;

        // Add a hole or material track state to the multitrajectory
        TrackIndexType currentTip =
            addNonSourcelinkState(stateMask, singleState, result, isSensitive,
                                  expectMeasurements, prevTip);
        currentBranch.tipIndex() = currentTip;
        auto currentState = currentBranch.outermostTrackState();
        if (expectMeasurements) {
          currentBranch.nHoles()++;
        }

        BranchStopperResult branchStopperResult =
            m_extensions.branchStopper(currentBranch, currentState);

        // Check the branch
        if (branchStopperResult == BranchStopperResult::Continue) {
          // Remembered the active branch and its state
        } else {
          // No branch on this surface
          if (branchStopperResult == BranchStopperResult::StopAndKeep) {
            storeLastActiveBranch(result);
          }
          // Remove the branch from list
          result.activeBranches.pop_back();

          // Branch on the surface has been stopped - reset
          ACTS_VERBOSE("Branch on surface " << surface->geometryId()
                                            << " has been stopped");

          reset(state, stepper, navigator, result);

          return Result<void>::success();
        }
      }

      auto currentState = currentBranch.outermostTrackState();

      if (currentState.typeFlags().test(TrackStateFlag::OutlierFlag)) {
        // We don't need to update the stepper given an outlier state
        ACTS_VERBOSE("Outlier state detected on surface "
                     << surface->geometryId());
      } else if (currentState.typeFlags().test(
                     TrackStateFlag::MeasurementFlag)) {
        // If there are measurement track states on this surface
        // Update stepping state using filtered parameters of last track
        // state on this surface
        updateStepper(state, stepper, tmpStates);
        ACTS_VERBOSE("Stepping state is updated with filtered parameter:");
        ACTS_VERBOSE("-> " << currentState.filtered().transpose()
                           << " of track state with tip = "
                           << currentState.index());
      }

      applyMultipleScattering(surface, state, stepper, navigator,
                              MaterialUpdateStage::PostUpdate);

      return Result<void>::success();
    }

    /// Remove components with low weights and renormalize from the component
    /// cache
    /// TODO This function does not expect normalized components, but this
    /// could be redundant work...
    void removeLowWeightComponents(std::vector<ComponentCache>& cmps) const {
      auto proj = [](auto& cmp) -> double& { return cmp.weight; };

      detail::normalizeWeights(cmps, proj);

      auto new_end = std::remove_if(cmps.begin(), cmps.end(), [&](auto& cmp) {
        return proj(cmp) < weightCutoff;
      });

      // In case we would remove all components, keep only the largest
      if (std::distance(cmps.begin(), new_end) == 0) {
        cmps = {*std::max_element(
            cmps.begin(), cmps.end(),
            [&](auto& a, auto& b) { return proj(a) < proj(b); })};
        cmps.front().weight = 1.0;
      } else {
        cmps.erase(new_end, cmps.end());
        detail::normalizeWeights(cmps, proj);
      }
    }

    /// Function that updates the stepper from the MultiTrajectory
    template <typename propagator_state_t, typename stepper_t>
    void updateStepper(propagator_state_t& state, const stepper_t& stepper,
                       const TemporaryStates& tmpStates) const {
      auto cmps = stepper.componentIterable(state.stepping);

      for (auto [idx, cmp] : zip(tmpStates.tips, cmps)) {
        // we set ignored components to missed, so we can remove them after
        // the loop
        if (tmpStates.weights.at(idx) < weightCutoff) {
          cmp.status() = IntersectionStatus::unreachable;
          continue;
        }

        auto proxy = tmpStates.traj.getTrackState(idx);

        cmp.pars() =
            MultiTrajectoryHelpers::freeFiltered(state.geoContext, proxy);
        cmp.cov() = proxy.filteredCovariance();
        cmp.weight() = tmpStates.weights.at(idx);
      }

      stepper.removeMissedComponents(state.stepping);

      // TODO we have two normalization passes here now, this can probably be
      // optimized
      detail::normalizeWeights(
          cmps, [&](auto cmp) -> double& { return cmp.weight(); });
    }

    /// Function that updates the stepper from the ComponentCache
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void updateStepper(
        propagator_state_t& state, const stepper_t& stepper,
        const navigator_t& navigator,
        const std::vector<ComponentCache>& componentCache) const {
      const auto& surface = *navigator.currentSurface(state.navigation);

      // Clear components before adding new ones
      stepper.clearComponents(state.stepping);

      // Finally loop over components
      for (const auto& [weight, pars, cov] : componentCache) {
        // Add the component to the stepper
        BoundTrackParameters bound(surface.getSharedPtr(), pars, cov,
                                   stepper.particleHypothesis(state.stepping));

        auto res =
            stepper.addComponent(state.stepping, std::move(bound), weight);

        if (!res.ok()) {
          ACTS_ERROR("Error adding component to MultiStepper");
          continue;
        }

        auto& cmp = *res;
        auto freeParams = cmp.pars();
        cmp.jacToGlobal() = surface.boundToFreeJacobian(
            state.geoContext, freeParams.template segment<3>(eFreePos0),
            freeParams.template segment<3>(eFreeDir0));
        cmp.pathAccumulated() = state.stepping.pathAccumulated;
        cmp.jacobian() = BoundMatrix::Identity();
        cmp.derivative() = FreeVector::Zero();
        cmp.jacTransport() = FreeMatrix::Identity();
      }
    }

    template <typename propagator_state_t, typename stepper_t>
    void convoluteComponents(
        const Surface& surface, propagator_state_t& state,
        const stepper_t& stepper,
        std::vector<ComponentCache>& componentCache) const {
      for (auto cmp : stepper.componentIterable(state.stepping)) {
        auto& singleState = cmp.singleState(state).stepping;
        const auto& singleStepper = cmp.singleStepper(stepper);

        auto res = singleStepper.boundState(singleState, surface);
        if (!res.ok()) {
          ACTS_ERROR("Propagate to surface " << surface.geometryId()
                                             << " failed: " << res.error());
          return;
        }
        const auto& [boundParams, jacobian, pathLength] = *res;

        BoundTrackParameters bound(surface.getSharedPtr(),
                                   boundParams.parameters(),
                                   *boundParams.covariance(),
                                   stepper.particleHypothesis(state.stepping));

        applyBetheHeitler(surface, state, bound, cmp.weight(), componentCache);
      }
    }

    template <typename propagator_state_t>
    double applyBetheHeitler(
        const Surface& surface, const propagator_state_t& state,
        const BoundTrackParameters& old_bound, const double old_weight,
        std::vector<ComponentCache>& componentCaches) const {
      const auto p_prev = old_bound.absoluteMomentum();
      const auto& particleHypothesis = old_bound.particleHypothesis();

      // Evaluate material slab
      auto slab = surface.surfaceMaterial()->materialSlab(
          old_bound.position(state.geoContext), state.options.direction,
          MaterialUpdateStage::FullUpdate);

      const auto pathCorrection = surface.pathCorrection(
          state.geoContext, old_bound.position(state.geoContext),
          old_bound.direction());
      slab.scaleThickness(pathCorrection);

      const double pathXOverX0 = slab.thicknessInX0();

      // Emit a warning if the approximation is not valid for this x/x0
      if (!betheHeitlerApprox.validXOverX0(pathXOverX0)) {
        ACTS_DEBUG(
            "Bethe-Heitler approximation encountered invalid value for x/x0="
            << pathXOverX0 << " at surface " << surface.geometryId());
      }

      // Get the mixture
      const auto mixture = betheHeitlerApprox.mixture(pathXOverX0);

      // Create all possible new components
      for (const auto& gaussian : mixture) {
        // Here we combine the new child weight with the parent weight.
        // However, this must be later re-adjusted
        const auto new_weight = gaussian.weight * old_weight;

        if (new_weight < weightCutoff) {
          ACTS_VERBOSE("Skip component with weight " << new_weight);
          continue;
        }

        if (gaussian.mean < 1.e-8) {
          ACTS_WARNING("Skip component with gaussian "
                       << gaussian.mean << " +- " << gaussian.var);
          continue;
        }

        // compute delta p from mixture and update parameters
        auto new_pars = old_bound.parameters();

        const auto delta_p = [&]() {
          if (state.options.direction == Direction::Forward()) {
            return p_prev * (gaussian.mean - 1.);
          } else {
            return p_prev * (1. / gaussian.mean - 1.);
          }
        }();

        assert(p_prev + delta_p > 0. && "new momentum must be > 0");
        new_pars[eBoundQOverP] =
            particleHypothesis.qOverP(p_prev + delta_p, old_bound.charge());

        // compute inverse variance of p from mixture and update covariance
        auto new_cov = old_bound.covariance().value();

        const auto varInvP = [&]() {
          if (state.options.direction == Direction::Forward()) {
            const auto f = 1. / (p_prev * gaussian.mean);
            return f * f * gaussian.var;
          } else {
            return gaussian.var / (p_prev * p_prev);
          }
        }();

        new_cov(eBoundQOverP, eBoundQOverP) += varInvP;
        assert(std::isfinite(new_cov(eBoundQOverP, eBoundQOverP)) &&
               "new cov not finite");

        // Set the remaining things and push to vector
        componentCaches.push_back({new_weight, new_pars, new_cov});
      }

      return pathXOverX0;
    }

    /// Process new, incompomplete track states and set the filtered state
    ///
    /// @note will process the given list of new states, run the updater
    ///     or share the predicted state for states flagged as outliers
    ///     and add them to the list of active branches
    ///
    /// @param newTrackStateList index list of new track states
    /// @param result which contains among others the new states, and the list of active branches
    /// @return the number of newly added branches or an error
    template <typename propagator_state_t, typename stepper_t>
    Result<unsigned int> processNewTrackStates(
        propagator_state_t& state, const stepper_t& stepper,
        const CkfTypes::BranchVector<TrackIndexType>& newTrackStateList,
        result_type& result, TemporaryStates& tmpStates) const {
      using PM = TrackStatePropMask;

      unsigned int nBranchesOnSurface = 0;

      auto rootBranch = result.activeBranches.back();

      // Build the new branches by forking the root branch. Reverse the order
      // to process the best candidate first
      CkfTypes::BranchVector<TrackProxy> newBranches;
      for (auto it = newTrackStateList.rbegin(); it != newTrackStateList.rend();
           ++it) {
        // Keep the root branch as the first branch, make a copy for the
        // others
        auto shallowCopy = [&] {
          auto sc = rootBranch.container().makeTrack();
          sc.copyFromShallow(rootBranch);
          return sc;
        };
        auto newBranch =
            (it == newTrackStateList.rbegin()) ? rootBranch : shallowCopy();
        newBranch.tipIndex() = *it;
        newBranches.push_back(newBranch);
      }

      // Remove the root branch
      result.activeBranches.pop_back();

      // Update and select from the new branches
      for (TrackProxy newBranch : newBranches) {
        auto trackState = newBranch.outermostTrackState();
        TrackStateType typeFlags = trackState.typeFlags();

        if (typeFlags.test(TrackStateFlag::OutlierFlag)) {
          // No Kalman update for outlier
          // Set the filtered parameter index to be the same with predicted
          // parameter
          trackState.shareFrom(PM::Predicted, PM::Filtered);
          // Increment number of outliers
          newBranch.nOutliers()++;
        } else if (typeFlags.test(TrackStateFlag::MeasurementFlag)) {
          // Kalman update
          auto updateRes = kalmanUpdate(trackState, state, stepper);
          if (!updateRes.ok()) {
            ACTS_ERROR("Update step failed: " << updateRes.error());
            return updateRes.error();
          }
          if (nBranchesOnSurface == 0) {
            tmpStates = *updateRes;
          }
          ACTS_VERBOSE("Appended measurement track state with tip = "
                       << newBranch.tipIndex());
          // Set the measurement flag
          typeFlags.set(TrackStateFlag::MeasurementFlag);
          // Increment number of measurements
          newBranch.nMeasurements()++;
          newBranch.nDoF() += trackState.calibratedSize();
          newBranch.chi2() += trackState.chi2();
        } else {
          ACTS_WARNING("Cannot handle this track state flags");
          continue;
        }

        result.activeBranches.push_back(newBranch);

        BranchStopperResult branchStopperResult =
            m_extensions.branchStopper(newBranch, trackState);

        // Check if need to stop this branch
        if (branchStopperResult == BranchStopperResult::Continue) {
          // Record the number of branches on surface
          nBranchesOnSurface++;
        } else {
          // Record the number of stopped branches
          if (branchStopperResult == BranchStopperResult::StopAndKeep) {
            storeLastActiveBranch(result);
          }
          // Remove the branch from list
          result.activeBranches.pop_back();
        }
      }

      return nBranchesOnSurface;
    }

    /// @brief CombinatorialKalmanFilter actor operation: add a hole or material track state
    ///
    /// @param stateMask The bitmask that instructs which components to allocate
    /// @param boundState The bound state on current surface
    /// @param result is the mutable result state object and which to leave invalid
    /// @param isSensitive The surface is sensitive or passive
    /// @param expectMeasurements True if measurements where expected for this surface
    /// @param prevTip The index of the previous state
    ///
    /// @return The tip of added state
    TrackIndexType addNonSourcelinkState(TrackStatePropMask stateMask,
                                         const SingleBoundState& boundState,
                                         result_type& result, bool isSensitive,
                                         bool expectMeasurements,
                                         TrackIndexType prevTip) const {
      using PM = TrackStatePropMask;

      // Add a track state
      auto trackStateProxy =
          result.trackStates->makeTrackState(stateMask, prevTip);
      ACTS_VERBOSE("Create "
                   << (isSensitive
                           ? (expectMeasurements ? "Hole"
                                                 : "noMeasurementExpected")
                           : "Material")
                   << " output track state #" << trackStateProxy.index()
                   << " with mask: " << stateMask);

      const auto& [boundParams, jacobian, pathLength] = boundState;
      // Fill the track state
      trackStateProxy.predicted() = boundParams.parameters();
      trackStateProxy.predictedCovariance() = boundParams.covariance().value();
      trackStateProxy.jacobian() = jacobian;
      trackStateProxy.pathLength() = pathLength;
      // Set the surface
      trackStateProxy.setReferenceSurface(
          boundParams.referenceSurface().getSharedPtr());

      // Set the track state flags
      auto typeFlags = trackStateProxy.typeFlags();
      if (trackStateProxy.referenceSurface().surfaceMaterial() != nullptr) {
        typeFlags.set(TrackStateFlag::MaterialFlag);
      }
      typeFlags.set(TrackStateFlag::ParameterFlag);
      if (isSensitive) {
        typeFlags.set(expectMeasurements ? TrackStateFlag::HoleFlag
                                         : TrackStateFlag::NoExpectedHitFlag);
      }

      // Set the filtered parameter index to be the same with predicted
      // parameter
      trackStateProxy.shareFrom(PM::Predicted, PM::Filtered);

      return trackStateProxy.index();
    }

    /// This function performs the kalman update, computes the new posterior
    /// weights, renormalizes all components, and does some statistics.
    template <typename propagator_state_t, typename stepper_t>
    Result<TemporaryStates> kalmanUpdate(TrackStateProxy& trackState,
                                         propagator_state_t& state,
                                         const stepper_t& stepper) const {
      const auto& surface = trackState.referenceSurface();

      TemporaryStates tmpStates;

      auto cmps = stepper.componentIterable(state.stepping);
      for (auto cmp : cmps) {
        auto& singleState = cmp.singleState(state).stepping;
        const auto& singleStepper = cmp.singleStepper(stepper);

        TrackStatePropMask mask =
            TrackStatePropMask::Predicted | TrackStatePropMask::Filtered |
            TrackStatePropMask::Jacobian | TrackStatePropMask::Calibrated;
        TrackStateProxy trackStateProxy =
            tmpStates.traj.makeTrackState(mask, kTrackIndexInvalid);

        // TODO call calibrator again?

        trackStateProxy.setReferenceSurface(surface.getSharedPtr());
        // Bind the transported state to the current surface
        auto res = singleStepper.boundState(singleState, surface);
        if (!res.ok()) {
          ACTS_ERROR("Propagate to surface " << surface.geometryId()
                                             << " failed: " << res.error());
          return res.error();
        }
        const auto& [boundParams, jacobian, pathLength] = *res;

        // Fill the track state
        trackStateProxy.predicted() = boundParams.parameters();
        trackStateProxy.predictedCovariance() = *boundParams.covariance();
        trackStateProxy.allocateCalibrated(trackState.calibratedSize());
        trackStateProxy.setProjectorSubspaceIndices(
            trackState.projectorSubspaceIndices());
        trackStateProxy.effectiveCalibrated() =
            trackState.effectiveCalibrated();
        trackStateProxy.effectiveCalibratedCovariance() =
            trackState.effectiveCalibratedCovariance();

        auto updateRes = m_extensions.updater(state.geoContext, trackStateProxy,
                                              *updaterLogger);

        if (!updateRes.ok()) {
          return updateRes.error();
        }

        tmpStates.tips.push_back(trackStateProxy.index());
        tmpStates.weights[tmpStates.tips.back()] = cmp.weight();
      }

      detail::computePosteriorWeights(tmpStates.traj, tmpStates.tips,
                                      tmpStates.weights);

      detail::normalizeWeights(tmpStates.tips, [&](auto idx) -> double& {
        return tmpStates.weights.at(idx);
      });

      updateMultiTrajectory(trackState, tmpStates, surface);

      // Return success
      return Result<TemporaryStates>::success(tmpStates);
    }

    void updateMultiTrajectory(TrackStateProxy& state,
                               const TemporaryStates& tmpStates,
                               const Surface& surface) const {
      using PrtProjector =
          detail::MultiTrajectoryProjector<detail::StatesType::ePredicted,
                                           TrackStateContainerBackend>;
      using FltProjector =
          detail::MultiTrajectoryProjector<detail::StatesType::eFiltered,
                                           TrackStateContainerBackend>;

      const auto isMeasurement = state.typeFlags().test(MeasurementFlag);

      auto [prtMean, prtCov] =
          mergeGaussianMixture(tmpStates.tips, surface, mergeMethod,
                               PrtProjector{tmpStates.traj, tmpStates.weights});
      state.predicted() = prtMean;
      state.predictedCovariance() = prtCov;

      if (isMeasurement) {
        auto [fltMean, fltCov] = mergeGaussianMixture(
            tmpStates.tips, surface, mergeMethod,
            FltProjector{tmpStates.traj, tmpStates.weights});
        state.filtered() = fltMean;
        state.filteredCovariance() = fltCov;
      } else {
        state.shareFrom(TrackStatePropMask::Predicted,
                        TrackStatePropMask::Filtered);
      }
    }

    /// Apply the multiple scattering to the state
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void applyMultipleScattering(const Surface* surface,
                                 propagator_state_t& state,
                                 const stepper_t& stepper,
                                 const navigator_t& navigator,
                                 const MaterialUpdateStage& updateStage =
                                     MaterialUpdateStage::FullUpdate) const {
      if (!multipleScattering) {
        return;
      }
      if (surface->surfaceMaterial() == nullptr) {
        return;
      }

      for (auto cmp : stepper.componentIterable(state.stepping)) {
        auto singleState = cmp.singleState(state);
        const auto& singleStepper = cmp.singleStepper(stepper);

        detail::PointwiseMaterialInteraction interaction(surface, singleState,
                                                         singleStepper);
        if (interaction.evaluateMaterialSlab(singleState, navigator,
                                             updateStage)) {
          // In the Gsf we only need to handle the multiple scattering
          interaction.evaluatePointwiseMaterialInteraction(multipleScattering,
                                                           false);

          // Screen out material effects info
          ACTS_VERBOSE("Material effects on surface: "
                       << surface->geometryId()
                       << " at update stage: " << updateStage << " are :");
          ACTS_VERBOSE("eLoss = "
                       << interaction.Eloss << ", "
                       << "variancePhi = " << interaction.variancePhi << ", "
                       << "varianceTheta = " << interaction.varianceTheta
                       << ", "
                       << "varianceQoverP = " << interaction.varianceQoverP);

          // Update the state and stepper with material effects
          interaction.updateState(singleState, singleStepper, addNoise);

          assert(singleState.stepping.cov.array().isFinite().all() &&
                 "covariance not finite after multi scattering");
        }
      }
    }

    void storeLastActiveBranch(result_type& result) const {
      auto currentBranch = result.activeBranches.back();
      TrackIndexType currentTip = currentBranch.tipIndex();

      ACTS_VERBOSE("Storing track "
                   << currentBranch.index() << " with tip index " << currentTip
                   << ". nMeasurements = " << currentBranch.nMeasurements()
                   << ", nOutliers = " << currentBranch.nOutliers()
                   << ", nHoles = " << currentBranch.nHoles());

      result.collectedTracks.push_back(currentBranch);
    }

    CombinatorialKalmanFilterExtensions<track_container_t> m_extensions;

    /// End of world aborter
    EndOfWorldReached endOfWorldReached;

    /// Volume constraint aborter
    VolumeConstraintAborter volumeConstraintAborter;

    /// Actor logger instance
    const Logger* actorLogger{nullptr};
    /// Updater logger instance
    const Logger* updaterLogger{nullptr};

    const Logger& logger() const { return *actorLogger; }
  };

  /// Void path limit reached aborter to replace the default since the path
  /// limit is handled in the CKF actor internally.
  struct StubPathLimitReached {
    double internalLimit{};

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    bool checkAbort(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/,
                    const Logger& /*logger*/) const {
      return false;
    }
  };

 public:
  /// Combinatorial Kalman Filter implementation, calls the Kalman filter
  ///
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param initialParameters The initial track parameters
  /// @param tfOptions CombinatorialKalmanFilter2Options steering the track
  ///                  finding
  /// @param trackContainer Track container in which to store the results
  /// @param rootBranch The track to be used as the root branch
  ///
  /// @note The input measurements are given in the form of @c SourceLinks.
  ///       It's @c calibrator_t's job to turn them into calibrated measurements
  ///       used in the track finding.
  ///
  /// @return a container of track finding result for all the initial track
  /// parameters
  auto findTracks(
      const BoundTrackParameters& initialParameters,
      const CombinatorialKalmanFilter2Options<track_container_t>& tfOptions,
      track_container_t& trackContainer,
      typename track_container_t::TrackProxy rootBranch) const
      -> Result<std::vector<
          typename std::decay_t<decltype(trackContainer)>::TrackProxy>> {
    // Create the ActorList
    using Actors = ActorList<Actor>;

    // Create relevant options for the propagation options
    using PropagatorOptions = typename propagator_t::template Options<Actors>;
    PropagatorOptions propOptions(tfOptions.geoContext,
                                  tfOptions.magFieldContext);

    // Set the trivial propagator options
    propOptions.setPlainOptions(tfOptions.propagatorPlainOptions);

    // Catch the actor
    auto& combKalmanActor = propOptions.actorList.template get<Actor>();
    combKalmanActor.targetReached.surface = tfOptions.targetSurface;
    combKalmanActor.multipleScattering = tfOptions.multipleScattering;
    combKalmanActor.energyLoss = tfOptions.energyLoss;
    combKalmanActor.skipPrePropagationUpdate =
        tfOptions.skipPrePropagationUpdate;
    combKalmanActor.actorLogger = m_actorLogger.get();
    combKalmanActor.updaterLogger = m_updaterLogger.get();
    combKalmanActor.calibrationContextPtr = &tfOptions.calibrationContext.get();

    // copy delegates to calibrator, updater, branch stopper
    combKalmanActor.m_extensions = tfOptions.extensions;

    combKalmanActor.maxComponents = tfOptions.maxComponents;
    combKalmanActor.weightCutoff = tfOptions.weightCutoff;
    combKalmanActor.mixtureReducer = tfOptions.mixtureReducer;
    combKalmanActor.mergeMethod = tfOptions.mergeMethod;
    combKalmanActor.betheHeitlerApprox = tfOptions.betheHeitlerApprox;

    auto propState =
        m_propagator
            .template makeState<PropagatorOptions, StubPathLimitReached>(
                propOptions);

    MultiComponentBoundTrackParameters multiInitialParameters(
        initialParameters.referenceSurface().getSharedPtr(),
        initialParameters.parameters(), initialParameters.covariance(),
        initialParameters.particleHypothesis());

    auto initResult =
        m_propagator.template initialize<decltype(propState),
                                         MultiComponentBoundTrackParameters,
                                         StubPathLimitReached>(
            propState, multiInitialParameters);
    if (!initResult.ok()) {
      ACTS_ERROR("Propagation initialization failed: " << initResult.error());
      return initResult.error();
    }

    auto& r =
        propState
            .template get<CombinatorialKalmanFilterResult<track_container_t>>();
    r.tracks = &trackContainer;
    r.trackStates = &trackContainer.trackStateContainer();

    r.activeBranches.push_back(rootBranch);

    auto propagationResult = m_propagator.propagate(propState);

    auto result = m_propagator.makeResult(
        std::move(propState), propagationResult, propOptions, false);

    if (!result.ok()) {
      ACTS_ERROR("Propagation failed: " << result.error() << " "
                                        << result.error().message()
                                        << " with the initial parameters: \n"
                                        << initialParameters.parameters());
      return result.error();
    }

    auto& propRes = *result;

    // Get the result of the CombinatorialKalmanFilter
    auto combKalmanResult =
        std::move(propRes.template get<
                  CombinatorialKalmanFilterResult<track_container_t>>());

    Result<void> error = combKalmanResult.lastError;
    if (error.ok() && !combKalmanResult.finished) {
      error = Result<void>(
          CombinatorialKalmanFilterError::PropagationReachesMaxSteps);
    }
    if (!error.ok()) {
      ACTS_ERROR("CombinatorialKalmanFilter failed: "
                 << combKalmanResult.lastError.error() << " "
                 << combKalmanResult.lastError.error().message()
                 << " with the initial parameters: "
                 << initialParameters.parameters().transpose());
      return error.error();
    }

    return std::move(combKalmanResult.collectedTracks);
  }

  /// Combinatorial Kalman Filter implementation, calls the Kalman filter
  ///
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param initialParameters The initial track parameters
  /// @param tfOptions CombinatorialKalmanFilter2Options steering the track
  ///                  finding
  /// @param trackContainer Track container in which to store the results
  /// @note The input measurements are given in the form of @c SourceLinks.
  ///       It's @c calibrator_t's job to turn them into calibrated measurements
  ///       used in the track finding.
  ///
  /// @return a container of track finding result for all the initial track
  /// parameters
  auto findTracks(
      const BoundTrackParameters& initialParameters,
      const CombinatorialKalmanFilter2Options<track_container_t>& tfOptions,
      track_container_t& trackContainer) const
      -> Result<std::vector<
          typename std::decay_t<decltype(trackContainer)>::TrackProxy>> {
    auto rootBranch = trackContainer.makeTrack();
    return findTracks(initialParameters, tfOptions, trackContainer, rootBranch);
  }
};

}  // namespace Acts

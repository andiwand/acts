// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/Propagator/MultiStepperLoop.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter2.hpp"
#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"
#include "Acts/TrackFitting/GsfMixtureReduction.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"

#include <memory>
#include <utility>

namespace {

using Stepper = Acts::MultiStepperLoop<Acts::SympyStepper>;
using Navigator = Acts::Navigator;
using Propagator = Acts::Propagator<Stepper, Navigator>;
using CKF =
    Acts::CombinatorialKalmanFilter2<Propagator, ActsExamples::TrackContainer>;

struct TrackFinderFunction2Impl
    : public ActsExamples::TrackFindingAlgorithm::TrackFinderFunction {
  CKF trackFinder;

  explicit TrackFinderFunction2Impl(CKF&& f) : trackFinder(std::move(f)) {}

  ActsExamples::TrackFindingAlgorithm::TrackFinderResult operator()(
      const ActsExamples::TrackParameters& initialParameters,
      const ActsExamples::TrackFindingAlgorithm::TrackFinderOptions& options,
      ActsExamples::TrackContainer& tracks,
      ActsExamples::TrackProxy rootBranch) const override {
    Acts::CombinatorialKalmanFilter2Options<ActsExamples::TrackContainer>
        options2(options.geoContext, options.magFieldContext,
                 options.calibrationContext, options.extensions,
                 options.propagatorPlainOptions, options.multipleScattering,
                 options.energyLoss);

    options2.mixtureReducer.connect<&Acts::reduceMixtureWithKLDistance>();
    options2.betheHeitlerApprox = Acts::makeDefaultBetheHeitlerApprox();

    return trackFinder.findTracks(initialParameters, options2, tracks,
                                  rootBranch);
  }
};

}  // namespace

std::shared_ptr<ActsExamples::TrackFindingAlgorithm::TrackFinderFunction>
ActsExamples::TrackFindingAlgorithm::makeTrackFinderFunction2(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    const Acts::Logger& logger) {
  Stepper stepper(std::move(magneticField));
  Navigator::Config cfg{std::move(trackingGeometry)};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Navigator navigator(cfg, logger.cloneWithSuffix("Navigator2"));
  Propagator propagator(std::move(stepper), std::move(navigator),
                        logger.cloneWithSuffix("Propagator2"));
  CKF trackFinder(std::move(propagator), logger.cloneWithSuffix("Finder2"));

  // build the track finder functions. owns the track finder object.
  return std::make_shared<TrackFinderFunction2Impl>(std::move(trackFinder));
}

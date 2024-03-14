// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/MyTrackFindingAlgorithm.hpp"

#include "Acts/EventData/ProxyAccessor.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFinding/MyTrackFinding.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"

namespace ActsExamples {

MyTrackFindingAlgorithm::MyTrackFindingAlgorithm(Config config,
                                                 Acts::Logging::Level level)
    : IAlgorithm("MyTrackFindingAlgorithm", level), m_cfg(std::move(config)) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements input collection");
  }
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links input collection");
  }
  if (m_cfg.inputInitialTrackParameters.empty()) {
    throw std::invalid_argument(
        "Missing initial track parameters input collection");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing tracks output collection");
  }

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_inputInitialTrackParameters.initialize(m_cfg.inputInitialTrackParameters);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ProcessCode MyTrackFindingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Read input data
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& sourceLinks = m_inputSourceLinks(ctx);
  const auto& initialParameters = m_inputInitialTrackParameters(ctx);

  // Construct a perigee surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});

  Acts::PropagatorPlainOptions pOptions;
  pOptions.maxSteps = m_cfg.maxSteps;
  pOptions.direction = Acts::Direction::Forward;

  using Stepper = Acts::EigenStepper<>;
  using Navigator = Acts::Navigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;

  Stepper stepper(m_cfg.magneticField);
  Navigator::Config cfg{m_cfg.trackingGeometry};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Navigator navigator(cfg, logger().cloneWithSuffix("Navigator"));
  Propagator propagator(std::move(stepper), std::move(navigator),
                        logger().cloneWithSuffix("Propagator"));

  PassThroughCalibrator pcalibrator;
  MeasurementCalibratorAdapter calibrator(pcalibrator, measurements);

  IndexSourceLinkAccessor slAccessor;
  slAccessor.container = &sourceLinks;
  Acts::SourceLinkAccessorDelegate<IndexSourceLinkAccessor::Iterator>
      slAccessorDelegate;
  slAccessorDelegate.connect<&IndexSourceLinkAccessor::range>(&slAccessor);

  Acts::MyTrackFinding<Propagator, Acts::VectorMultiTrajectory> finder({});

  // Perform the track finding for all initial parameters
  ACTS_DEBUG("Invoke track finding with " << initialParameters.size()
                                          << " seeds.");

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();

  auto trackContainerTemp = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainerTemp =
      std::make_shared<Acts::VectorMultiTrajectory>();

  TrackContainer tracks(trackContainer, trackStateContainer);
  TrackContainer tracksTemp(trackContainerTemp, trackStateContainerTemp);

  tracks.addColumn<unsigned int>("trackGroup");
  tracksTemp.addColumn<unsigned int>("trackGroup");
  Acts::ProxyAccessor<unsigned int> seedNumber("trackGroup");

  unsigned int nSeed = 0;

  for (std::size_t iseed = 0; iseed < initialParameters.size(); ++iseed) {
    // Clear trackContainerTemp and trackStateContainerTemp
    tracksTemp.clear();

    auto result = finder.findTracks(initialParameters.at(iseed), tracksTemp);
    m_nTotalSeeds++;
    nSeed++;

    if (!result.ok()) {
      m_nFailedSeeds++;
      ACTS_WARNING("Track finding failed for seed " << iseed << " with error"
                                                    << result.error());
      continue;
    }

    auto& tracksForSeed = result.value();
    for (auto& track : tracksForSeed) {
      // Set the seed number, this number decrease by 1 since the seed number
      // has already been updated
      seedNumber(track) = nSeed - 1;
      auto destProxy = tracks.getTrack(tracks.addTrack());
      destProxy.copyFrom(track, true);  // make sure we copy track states!
    }
  }

  ACTS_DEBUG("Finalized track finding with " << tracks.size()
                                             << " track candidates.");

  m_memoryStatistics.local().hist +=
      tracks.trackStateContainer().statistics().hist;

  auto constTrackStateContainer =
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer));

  auto constTrackContainer = std::make_shared<Acts::ConstVectorTrackContainer>(
      std::move(*trackContainer));

  ConstTrackContainer constTracks{constTrackContainer,
                                  constTrackStateContainer};

  m_outputTracks(ctx, std::move(constTracks));
  return ProcessCode::SUCCESS;
}

ProcessCode MyTrackFindingAlgorithm::finalize() {
  ACTS_INFO("MyTrackFindingAlgorithm statistics:");
  ACTS_INFO("- total seeds: " << m_nTotalSeeds);
  ACTS_INFO("- failed seeds: " << m_nFailedSeeds);
  ACTS_INFO("- failure ratio: " << static_cast<double>(m_nFailedSeeds) /
                                       m_nTotalSeeds);

  auto memoryStatistics =
      m_memoryStatistics.combine([](const auto& a, const auto& b) {
        Acts::VectorMultiTrajectory::Statistics c;
        c.hist = a.hist + b.hist;
        return c;
      });
  std::stringstream ss;
  memoryStatistics.toStream(ss);
  ACTS_DEBUG("Track State memory statistics (averaged):\n" << ss.str());
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

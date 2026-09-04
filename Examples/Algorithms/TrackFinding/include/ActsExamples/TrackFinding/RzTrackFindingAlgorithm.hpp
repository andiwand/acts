// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/TrackFinding/Rz/RzLayout.hpp"
#include "Acts/TrackFinding/Rz/RzTrackFinder.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cstdint>
#include <memory>
#include <string>

namespace ActsExamples {

/// Track finding on the RZ skeleton of the tracking geometry, one track per
/// initial parameter set, forward only. Tracks come out with their parameters
/// at the perigee and their measurements as track states, so the truth
/// matcher and the performance writers apply unchanged.
class RzTrackFindingAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    std::string inputMeasurements;
    std::string inputInitialTrackParameters;
    std::string outputTracks;

    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;

    /// Measurement binning per layer, see `Acts::Experimental::RzLayoutOptions`
    std::uint32_t phiBins = 64;
    double alongBinWidth = 20 * Acts::UnitConstants::mm;

    /// The finder's cuts, see `Acts::Experimental::RzTrackFinderConfig`
    double chi2Cut = 15.;
    double windowSigmas = 5.;
    double windowMin = 1. * Acts::UnitConstants::mm;
    std::uint32_t maxHoles = 3;
    std::uint32_t maxConsecutiveHoles = 2;
    std::uint32_t minMeasurements = 6;
    std::uint32_t maxMeasurementsPerLayer = 2;
    bool applyMaterial = true;
    double backwardInflation = 100.;
    std::uint32_t backwardLayers = 6;
    double backwardQOverPScale = 1.;
  };

  RzTrackFindingAlgorithm(const Config& config,
                          std::unique_ptr<const Acts::Logger> log);

  ProcessCode execute(const AlgorithmContext& ctx) const override;
  ProcessCode finalize() override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  Acts::Experimental::RzLayout m_layout;
  std::shared_ptr<const Acts::PerigeeSurface> m_perigee;

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};
  ReadDataHandle<TrackParametersContainer> m_inputInitialTrackParameters{
      this, "InputInitialTrackParameters"};
  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};

  mutable std::atomic<std::size_t> m_nSeeds{0};
  mutable std::atomic<std::size_t> m_nTracks{0};
  mutable std::atomic<std::size_t> m_nStops{0};
  mutable std::atomic<std::size_t> m_nCandidates{0};
  mutable std::atomic<std::size_t> m_nMeasurementsOnTracks{0};
  mutable std::atomic<std::size_t> m_nHolesOnTracks{0};
  mutable std::atomic<std::size_t> m_nBackwardFailures{0};
};

}  // namespace ActsExamples

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"

#include "ActsExamples/EventData/ProtoTrack.hpp"

#include <ostream>
#include <stdexcept>
#include <utility>

namespace ActsExamples {

TruthTrackFinder::TruthTrackFinder(const Config& config,
                                   std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("TruthTrackFinder", std::move(logger)), m_cfg(config) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputParticleMeasurementsMap.empty()) {
    throw std::invalid_argument("Missing input hit-particles map collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing output proto tracks collection");
  }
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing input measurements collection");
  }
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing input simulated hits collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing input simulated hits measurements map");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputParticleMeasurementsMap.initialize(m_cfg.inputParticleMeasurementsMap);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);
}

ProcessCode TruthTrackFinder::execute(const AlgorithmContext& ctx) const {
  // prepare input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& particleMeasurementsMap = m_inputParticleMeasurementsMap(ctx);
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& simHits = m_inputSimHits(ctx);
  const auto& measurementSimHitsMap = m_inputMeasurementSimHitsMap(ctx);

  // prepare output collection
  ProtoTrackContainer tracks;
  tracks.reserve(particles.size());

  ACTS_VERBOSE("Create proto tracks for " << particles.size() << " particles");
  for (const auto& particle : particles) {
    // find the corresponding hits for this particle
    const auto& particleMeasurements =
        makeRange(particleMeasurementsMap.equal_range(particle.index()));

    ACTS_VERBOSE(" - From " << particleMeasurements.size()
                            << " measurements for particle " << particle);

    // fill hit indices to create the proto track
    ProtoTrack track;
    std::vector<const SimHit*> hits;
    track.reserve(particleMeasurements.size());
    hits.reserve(particleMeasurements.size());
    for (const auto& [_, measurementId] : particleMeasurements) {
      ConstVariableBoundMeasurementProxy measurement =
          measurements.getMeasurement(measurementId);

      ACTS_VERBOSE("   - Measurement " << measurementId << " with barcode "
                                       << particle.barcode() << " at "
                                       << measurement.geometryId());

      const auto measurementSimHits =
          measurementSimHitsMap.equal_range(measurementId);
      if (measurementSimHits.first == measurementSimHits.second) {
        ACTS_WARNING("No sim hit found for measurement index "
                     << measurementId);
        continue;
      }

      const SimHit& firstSimHit = simHits.at(measurementSimHits.first->second);

      track.emplace_back(measurementId);
      hits.emplace_back(&firstSimHit);
    }

    std::vector<std::size_t> indices(hits.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::ranges::sort(indices, [&hits](std::size_t a, std::size_t b) {
      return hits[a]->time() < hits[b]->time();
    });
    ProtoTrack sortedTrack;
    sortedTrack.reserve(hits.size());
    for (const std::size_t idx : indices) {
      sortedTrack.emplace_back(track[idx]);
    }

    // add proto track to the output collection
    tracks.emplace_back(std::move(sortedTrack));
  }

  m_outputProtoTracks(ctx, std::move(tracks));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

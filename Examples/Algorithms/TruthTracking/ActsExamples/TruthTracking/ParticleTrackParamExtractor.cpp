// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/ParticleTrackParamExtractor.hpp"

#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <stdexcept>
#include <utility>

namespace ActsExamples {

ParticleTrackParamExtractor::ParticleTrackParamExtractor(
    const Config& config, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("ParticleTrackParamExtractor", std::move(logger)),
      m_cfg(config) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing output track parameters collection");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
}

ProcessCode ParticleTrackParamExtractor::execute(
    const AlgorithmContext& ctx) const {
  const SimParticleContainer& particles = m_inputParticles(ctx);

  std::unordered_map<SimBarcode, std::shared_ptr<Acts::PerigeeSurface>>
      perigeeSurfaces;

  // create track parameters from the particles
  TrackParametersContainer trackParameters;

  for (const auto& particle : particles) {
    const auto vtxId = particle.barcode().vertexId();
    const auto particleHypothesis = particle.hypothesis();
    const double phi =
        Acts::VectorHelpers::phi(particle.generationState().direction());
    const double theta =
        Acts::VectorHelpers::theta(particle.generationState().direction());
    const double qOverP = particle.generationState().qOverP();
    const double time = particle.generationState().time();

    std::shared_ptr<Acts::PerigeeSurface>& perigee = perigeeSurfaces[vtxId];
    if (perigee == nullptr) {
      // assume that all particles with the same vertex ID originate from the
      // same position and use it to as the reference position for the perigee
      // frame.
      perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(
          particle.generationState().position());
    }

    trackParameters.emplace_back(
        perigee, Acts::BoundVector{0, 0, phi, theta, qOverP, time},
        std::nullopt, particleHypothesis);
  }

  m_outputTrackParameters(ctx, std::move(trackParameters));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/Csv/CsvInputOutput.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/GenerationProcessType.hpp"

#include <stdexcept>
#include <string>

#include "CsvOutputData.hpp"

namespace ActsExamples {

CsvParticleReader::CsvParticleReader(const Config& config,
                                     Acts::Logging::Level level)
    : m_cfg(config),
      m_eventsRange(
          determineEventFilesRange(m_cfg.inputDir, m_cfg.inputStem + ".csv")),
      m_logger(Acts::getDefaultLogger("CsvParticleReader", level)) {
  if (m_cfg.inputStem.empty()) {
    throw std::invalid_argument("Missing input filename stem");
  }
  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument("Missing output collection");
  }

  m_outputParticles.initialize(m_cfg.outputParticles);
}

std::string CsvParticleReader::name() const {
  return "CsvParticleReader";
}

std::pair<std::size_t, std::size_t> CsvParticleReader::availableEvents() const {
  return m_eventsRange;
}

ProcessCode CsvParticleReader::read(const AlgorithmContext& ctx) {
  SimParticleContainer particles(SimParticleColumns::Generated);

  const std::string path = perEventFilepath(
      m_cfg.inputDir, m_cfg.inputStem + ".csv", ctx.eventNumber);

  // vt and m are an optional columns
  BoostDescribeCsvReader<ParticleData> reader(path, {"vt", "m"});
  ParticleData data;
  while (reader.read(data)) {
    MutableSimParticle particle = particles.createParticle();
    particle.assignParentIndices({});
    particle.barcode() = ActsFatras::Barcode()
                             .withVertexPrimary(data.particle_id_pv)
                             .withVertexSecondary(data.particle_id_sv)
                             .withParticle(data.particle_id_part)
                             .withGeneration(data.particle_id_gen)
                             .withSubParticle(data.particle_id_subpart);
    particle.pdg() = Acts::PdgParticle{data.particle_type};
    particle.charge() = data.q * Acts::UnitConstants::e;
    particle.mass() = data.m * Acts::UnitConstants::GeV;
    particle.generationProcess() =
        static_cast<ActsFatras::GenerationProcessType>(data.process);
    particle.generationState().fourPosition() = Acts::Vector4(
        data.vx * Acts::UnitConstants::mm, data.vy * Acts::UnitConstants::mm,
        data.vz * Acts::UnitConstants::mm, data.vt * Acts::UnitConstants::mm);
    particle.generationState().absoluteMomentum() =
        Acts::fastHypot(data.px, data.py, data.pz);
    particle.generationState().direction() =
        Acts::Vector3(data.px, data.py, data.pz);
  }

  m_outputParticles(ctx, std::move(particles));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

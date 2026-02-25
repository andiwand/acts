// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SpacePointConverter.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"

namespace ActsExamples {

SpacePointConverter::SpacePointConverter(const Config& cfg,
                                         Acts::Logging::Level lvl)
    : IAlgorithm("SpacePointConverter", lvl), m_cfg(cfg) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point input collection");
  }
  if (m_cfg.outputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point output collection");
  }
}

ProcessCode SpacePointConverter::execute(const AlgorithmContext& ctx) const {
  const auto& inputSpacePoints = m_inputSpacePoints(ctx);

  Acts::SpacePointContainer2 outputSpacePoints(
      Acts::SpacePointColumns::SourceLinks | Acts::SpacePointColumns::X |
      Acts::SpacePointColumns::Y | Acts::SpacePointColumns::Z |
      Acts::SpacePointColumns::R | Acts::SpacePointColumns::Time |
      Acts::SpacePointColumns::VarianceZ | Acts::SpacePointColumns::VarianceR);
  outputSpacePoints.reserve(inputSpacePoints.size());
  for (const auto& sp : inputSpacePoints) {
    auto newSp = outputSpacePoints.createSpacePoint();
    newSp.assignSourceLinks(sp.sourceLinks());
    newSp.x() = sp.x();
    newSp.y() = sp.y();
    newSp.z() = sp.z();
    newSp.r() = sp.r();
    if (sp.t().has_value()) {
      newSp.time() = sp.t().value();
    }
    newSp.varianceZ() = sp.varianceZ();
    newSp.varianceR() = sp.varianceR();
  }

  m_outputSpacePoints(ctx, std::move(outputSpacePoints));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

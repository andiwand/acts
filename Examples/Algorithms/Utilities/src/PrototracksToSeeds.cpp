// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/PrototracksToSeeds.hpp"

#include "ActsExamples/EventData/Seed.hpp"
#include "ActsExamples/Utilities/EventDataTransforms.hpp"

#include <algorithm>

namespace ActsExamples {

PrototracksToSeeds::PrototracksToSeeds(Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("PrototracksToSeeds", lvl), m_cfg(std::move(cfg)) {
  m_outputSeeds.initialize(m_cfg.outputSeeds);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
}

ProcessCode PrototracksToSeeds::execute(const AlgorithmContext& ctx) const {
  auto prototracks = m_inputProtoTracks(ctx);
  const auto& spacePointContainer = m_inputSpacePoints(ctx);

  const auto nBefore = prototracks.size();
  prototracks.erase(std::remove_if(prototracks.begin(), prototracks.end(),
                                   [](const auto& t) { return t.size() < 3; }),
                    prototracks.end());
  ACTS_DEBUG("Discarded " << prototracks.size() - nBefore
                          << " prototracks with less then 3 hits");

  SeedContainer seeds;
  seeds.reserve(prototracks.size());

  for (const auto& pt : prototracks) {
    prototrackToSeed(pt, spacePointContainer, seeds);
  }

  m_outputSeeds(ctx, std::move(seeds));
  m_outputProtoTracks(ctx, std::move(prototracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

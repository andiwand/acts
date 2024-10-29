// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ActsExamples {
struct AlgorithmContext;

class SegmentSeedingAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    std::vector<std::string> inputSpacePoints;
    std::string outputSeeds;
  };

  SegmentSeedingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  std::vector<std::unique_ptr<ReadDataHandle<SimSpacePointContainer>>>
      m_inputSpacePoints{};

  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};
};

}  // namespace ActsExamples

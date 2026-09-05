// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/Logging.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/SequenceElement.hpp"

#include <string>

namespace ActsExamples {

/// Print track parameters.
class TrackParametersPrinter : public IAlgorithm {
 public:
  struct Config {
    /// Logger for this component. Unnamed by default, in which case it is
    /// named after the component. Assign a named logger to override.
    std::shared_ptr<const Acts::Logger> logger = makeDefaultLogger();

    /// Input tracks parameters collection.
    std::string inputTrackParameters;
  };

  explicit TrackParametersPrinter(const Config& cfg);

  ProcessCode execute(const AlgorithmContext& ctx) const override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<TrackParametersContainer> m_inputTrackParameters{
      this, "InputTrackParameters"};
};

}  // namespace ActsExamples

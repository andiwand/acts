// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/Logging.hpp"
#include "ActsExamples/Io/Podio/PodioInputConverter.hpp"

#include <string>

namespace ActsExamples {

/// Read in a track collection as EDM4hep from a @c podio::Frame.
class EDM4hepTrackInputConverter : public PodioInputConverter {
 public:
  struct Config {
    /// Logger for this component. Unnamed by default, in which case it is
    /// named after the component. Assign a named logger to override.
    std::shared_ptr<const Acts::Logger> logger = makeDefaultLogger();

    std::string inputFrame;
    /// Input track collection name in edm4hep
    std::string inputTracks = "ActsTracks";
    /// Output track collection
    std::string outputTracks;
    /// Magnetic field along the z axis (needed for the conversion of
    /// parameters)
    double Bz{};
  };

  /// constructor
  /// @param config is the configuration object
  explicit EDM4hepTrackInputConverter(const Config& config);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  ProcessCode convert(const AlgorithmContext& ctx,
                      const podio::Frame& frame) const final;

 private:
  Config m_cfg;

  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples

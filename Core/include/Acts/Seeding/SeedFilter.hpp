// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/EventData/SpacePointMutableData.hpp"
#include "Acts/Seeding/CandidatesForMiddleSp.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <vector>

namespace Acts {

/// Filter seeds at various stages with the currently
/// available information.
class SeedFilter final {
 public:
  using Config = SeedFilterConfig;

  struct State {
    // longitudinal impact parameter as defined by bottom and middle space point
    float zOrigin = 0;
    // number of minimum top SPs in seed confirmation
    std::size_t nTopSeedConf = 0;
    // radius of bottom component of seed that is used to define the number of
    // compatible top required
    float rMaxSeedConf = std::numeric_limits<float>::max();
  };

  explicit SeedFilter(const Config& config,
                      std::unique_ptr<const Logger> logger =
                          getDefaultLogger("SeedFilter", Logging::Level::INFO));

  /// Create Seeds for the all seeds with the same bottom and middle
  /// space point and discard all others.
  /// @param mutableData Container for mutable variables used in the seeding
  /// @param bottomSpIndex fixed bottom space point
  /// @param middleSpIndex fixed middle space point
  /// @param topSpIndices vector containing all space points that may be compatible
  ///                 with both bottom and middle space point
  /// @param invHelixDiameters vector containing 1/(2*r) values where r is the helix radius
  /// @param impactParameters vector containing the impact parameters
  /// @param seedFilterState holds quantities used in seed filter
  /// @param candidates_collector container for the seed candidates
  void filterSeeds_2SpFixed(const SpacePointContainer& spacePointContainer,
                            const SpacePointMutableData& mutableData,
                            SpacePointIndex bottomSpIndex,
                            SpacePointIndex middleSpIndex,
                            const std::vector<SpacePointIndex>& topSpIndices,
                            const std::vector<float>& invHelixDiameters,
                            const std::vector<float>& impactParameters,
                            State& seedFilterState,
                            CandidatesForMiddleSp& candidates_collector) const;

  /// Filter seeds once all seeds for one middle space point have been created
  /// @param mutableData Container for mutable variables used in the seeding
  /// @param candidates_collector collection of seed candidates
  /// @param outputCollection Output container for the seeds
  /// for all seeds with the same middle space point
  void filterSeeds_1SpFixed(const SpacePointContainer& spacePointContainer,
                            SpacePointMutableData& mutableData,
                            CandidatesForMiddleSp& candidatesCollector,
                            SeedContainer& outputCollection) const;

  /// Filter seeds once all seeds for one middle space point have been created
  /// @param mutableData Container for mutable variables used in the seeding
  /// @param candidates collection of seed candidates
  /// @param numQualitySeeds number of high quality seeds in seed confirmation
  /// @param outputCollection Output container for the seeds
  /// for all seeds with the same middle space point
  void filterSeeds_1SpFixed(SpacePointMutableData& mutableData,
                            std::vector<TripletCandidate>& candidates,
                            const std::size_t numQualitySeeds,
                            SeedContainer& outputCollection) const;

  const Config getConfig() const { return m_cfg; }

 private:
  const Config m_cfg;
  std::unique_ptr<const Logger> m_logger;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

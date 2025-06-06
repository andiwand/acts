// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TripletSeedingAlgorithm2.hpp"

#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Seeding2/TripletSeedFinder2.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <cmath>
#include <csignal>
#include <cstddef>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace ActsExamples {

namespace {

std::vector<std::array<Acts::GeometryIdentifier, 3>> readTriplets(
    const std::string& path) {
  std::vector<std::array<Acts::GeometryIdentifier, 3>> result;

  std::ifstream in(path);

  std::string line;
  std::getline(in, line);
  while (in) {
    std::getline(in, line);

    std::stringstream lineIn(line);
    std::vector<std::string> row;

    std::string cell;
    while (std::getline(lineIn, cell, ',')) {
      row.push_back(cell);
    }

    if (row.empty()) {
      continue;
    }
    if (row.size() != 3) {
      throw std::runtime_error(
          "TripletSeedingAlgorithm2: Invalid triplet line: " + line);
    }
    if (row[0].empty() || row[1].empty() || row[2].empty()) {
      throw std::runtime_error(
          "TripletSeedingAlgorithm2: Empty triplet line: " + line);
    }
    result.push_back({Acts::GeometryIdentifier(std::stoul(row[0])),
                      Acts::GeometryIdentifier(std::stoul(row[1])),
                      Acts::GeometryIdentifier(std::stoul(row[2]))});
  }

  return result;
}

}  // namespace

static std::vector<std::array<Acts::GeometryIdentifier, 3>> module_triplets;

TripletSeedingAlgorithm2::TripletSeedingAlgorithm2(const Config& cfg,
                                                   Acts::Logging::Level lvl)
    : IAlgorithm("SeedingAlgorithm2", lvl), m_cfg(cfg) {
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputSeeds.initialize(m_cfg.outputSeeds);

  static_assert(std::numeric_limits<
                    decltype(m_cfg.finderConfig.deltaRMaxTopSP)>::has_quiet_NaN,
                "Value of deltaRMaxTopSP must support NaN values");

  static_assert(std::numeric_limits<
                    decltype(m_cfg.finderConfig.deltaRMinTopSP)>::has_quiet_NaN,
                "Value of deltaRMinTopSP must support NaN values");

  static_assert(
      std::numeric_limits<
          decltype(m_cfg.finderConfig.deltaRMaxBottomSP)>::has_quiet_NaN,
      "Value of deltaRMaxBottomSP must support NaN values");

  static_assert(
      std::numeric_limits<
          decltype(m_cfg.finderConfig.deltaRMinBottomSP)>::has_quiet_NaN,
      "Value of deltaRMinBottomSP must support NaN values");

  if (m_cfg.gridConfig.cotThetaMax != m_cfg.finderConfig.cotThetaMax) {
    throw std::invalid_argument("Inconsistent config cotThetaMax");
  }

  if (m_cfg.gridConfig.minPt != m_cfg.finderConfig.minPt) {
    throw std::invalid_argument("Inconsistent config minPt");
  }

  if (m_cfg.gridConfig.bFieldInZ != m_cfg.finderOptions.bFieldInZ) {
    throw std::invalid_argument("Inconsistent config bFieldInZ");
  }

  if (m_cfg.gridConfig.zBinEdges.size() - 1 != m_cfg.zBinNeighborsTop.size() &&
      m_cfg.zBinNeighborsTop.empty() == false) {
    throw std::invalid_argument("Inconsistent config zBinNeighborsTop");
  }

  if (m_cfg.gridConfig.zBinEdges.size() - 1 !=
          m_cfg.zBinNeighborsBottom.size() &&
      m_cfg.zBinNeighborsBottom.empty() == false) {
    throw std::invalid_argument("Inconsistent config zBinNeighborsBottom");
  }

  if (!m_cfg.zBinsCustomLooping.empty()) {
    // check that the bins required in the custom bin looping
    // are contained in the bins defined by the total number of edges

    for (std::size_t i : m_cfg.zBinsCustomLooping) {
      if (i >= m_cfg.gridConfig.zBinEdges.size()) {
        throw std::invalid_argument(
            "Inconsistent config zBinsCustomLooping does not contain a subset "
            "of bins defined by zBinEdges");
      }
    }
  }

  if (m_cfg.useExtraCuts) {
    // This function will be applied to select space points during grid filling
    m_spacePointSelector.connect<itkFastTrackingSPselect>();

    // This function will be applied to the doublet compatibility selection
    m_cfg.finderConfig.experimentCuts.connect<itkFastTrackingCuts>();
  }

  m_cfg.gridConfig.bottomBinFinder.emplace(m_cfg.numPhiNeighbors,
                                           m_cfg.zBinNeighborsBottom, 0);
  m_cfg.gridConfig.topBinFinder.emplace(m_cfg.numPhiNeighbors,
                                        m_cfg.zBinNeighborsTop, 0);
  m_cfg.gridConfig.navigation[1ul] = m_cfg.zBinsCustomLooping;

  m_cfg.finderConfig.filter = std::make_unique<Acts::TripletSeedFilter2>(
      m_cfg.filterConfig, logger().cloneWithSuffix("Filter"));

  m_seedFinder = Acts::TripletSeedFinder2(m_cfg.finderConfig.derive(),
                                          logger().cloneWithSuffix("Finder"));

  module_triplets =
      readTriplets("/Users/andreas/cern/scripts/notebooks/module_triplets.csv");
  std::cout << module_triplets.size()
            << " triplets read from module_triplets.csv" << std::endl;
}

ProcessCode TripletSeedingAlgorithm2::execute(
    const AlgorithmContext& ctx) const {
  const SimSpacePointContainer& spacePoints = m_inputSpacePoints(ctx);

  Acts::SpacePointContainer2 coreSpacePoints;
  auto& rColumn = coreSpacePoints.createDenseExtraColumn<float>("r");
  auto& phiColumn = coreSpacePoints.createDenseExtraColumn<float>("phi");
  auto& varianceRColumn =
      coreSpacePoints.createDenseExtraColumn<float>("varianceR");
  auto& varianceZColumn =
      coreSpacePoints.createDenseExtraColumn<float>("varianceZ");
  coreSpacePoints.reserve(spacePoints.size());
  for (const auto& sp : spacePoints) {
    // check if the space point passes the selection
    if (m_spacePointSelector(sp)) {
      auto newSp = coreSpacePoints.createSpacePoint(
          std::array<Acts::SourceLink, 1>{Acts::SourceLink(&sp)}, sp.x(),
          sp.y(), sp.z());
      newSp.extra(rColumn) = sp.r();
      newSp.extra(phiColumn) = std::atan2(sp.y(), sp.x());
      newSp.extra(varianceRColumn) = sp.varianceR();
      newSp.extra(varianceZColumn) = sp.varianceZ();
    }
  }

  std::vector<Acts::SpacePointIndex2> sortedCoreSpacePoints;
  sortedCoreSpacePoints.resize(coreSpacePoints.size());
  std::iota(sortedCoreSpacePoints.begin(), sortedCoreSpacePoints.end(), 0);
  std::ranges::sort(
      sortedCoreSpacePoints,
      [&](const Acts::SpacePointIndex2& a, const Acts::SpacePointIndex2& b) {
        auto spA = coreSpacePoints.at(a);
        auto spB = coreSpacePoints.at(b);
        auto sourceLinkA = spA.sourceLinks()[0]
                               .get<const SimSpacePoint*>()
                               ->sourceLinks()[0]
                               .get<IndexSourceLink>();
        auto sourceLinkB = spB.sourceLinks()[0]
                               .get<const SimSpacePoint*>()
                               ->sourceLinks()[0]
                               .get<IndexSourceLink>();
        if (sourceLinkA.geometryId() == sourceLinkB.geometryId()) {
          // float cotThetaA = spA.z() / spA.extra(rColumn);
          // float cotThetaB = spB.z() / spB.extra(rColumn);
          // return cotThetaA < cotThetaB;
          return spA.extra(rColumn) < spB.extra(rColumn);
        }
        return sourceLinkA.geometryId() < sourceLinkB.geometryId();
      });

  Acts::CylindricalSpacePointGrid2 grid(m_cfg.gridConfig.derive(),
                                        logger().cloneWithSuffix("Grid"));

  grid.fill(coreSpacePoints, phiColumn, rColumn);

  // Compute radius Range
  // we rely on the fact the grid is storing the proxies
  // with a sorting in the radius
  const Acts::Range1D<float> rRange =
      grid.computeRadiusRange(coreSpacePoints, rColumn);

  Acts::TripletSeedFinder2::Options finderOptions = m_cfg.finderOptions;

  /// variable middle SP radial region of interest
  finderOptions.rMiddleSpRange = {
      std::floor(rRange.min() / 2) * 2 + m_cfg.deltaRMiddleMinSPRange,
      std::floor(rRange.max() / 2) * 2 - m_cfg.deltaRMiddleMaxSPRange};

  // run the seeding
  Acts::SeedContainer2 seeds;
  Acts::TripletSeedFinder2::State state;
  Acts::TripletSeedFinder2::Cache cache;

  auto derivedOptions = finderOptions.derive(m_seedFinder->config());
  m_seedFinder->initialize(state, derivedOptions);

  std::vector<Acts::SpacePointIndex2> bottomSp;
  std::vector<Acts::SpacePointIndex2> middleSp;
  std::vector<Acts::SpacePointIndex2> topSp;

  // for (const auto [bottom, middle, top] : grid.binnedGround()) {
  //   ACTS_VERBOSE("Process middle " << middle);

  //   bottomSp.clear();
  //   for (const auto b : bottom) {
  //     bottomSp.insert(bottomSp.end(), grid.at(b).begin(), grid.at(b).end());
  //   }
  //   middleSp.clear();
  //   middleSp = grid.at(middle);
  //   topSp.clear();
  //   for (const auto t : top) {
  //     topSp.insert(topSp.end(), grid.at(t).begin(), grid.at(t).end());
  //   }

  //   m_seedFinder->createSeeds(
  //       state, cache,
  //       Acts::TripletSeedFinder2::ContainerPointers(
  //           coreSpacePoints, rColumn, varianceRColumn, varianceZColumn),
  //       bottomSp, middleSp, topSp, seeds);
  // }

  for (const auto& triplet : module_triplets) {
    auto bottomSpRange = std::ranges::equal_range(
        sortedCoreSpacePoints, triplet[0].value(), {},
        [&](Acts::SpacePointIndex2 i) {
          auto sp = coreSpacePoints.at(i);
          auto sourceLink = sp.sourceLinks()[0]
                                .get<const SimSpacePoint*>()
                                ->sourceLinks()[0]
                                .get<IndexSourceLink>();
          return sourceLink.geometryId().value();
        });
    auto middleSpRange = std::ranges::equal_range(
        sortedCoreSpacePoints, triplet[1].value(), {},
        [&](Acts::SpacePointIndex2 i) {
          auto sp = coreSpacePoints.at(i);
          auto sourceLink = sp.sourceLinks()[0]
                                .get<const SimSpacePoint*>()
                                ->sourceLinks()[0]
                                .get<IndexSourceLink>();
          return sourceLink.geometryId().value();
        });
    auto topSpRange = std::ranges::equal_range(
        sortedCoreSpacePoints, triplet[2].value(), {},
        [&](Acts::SpacePointIndex2 i) {
          auto sp = coreSpacePoints.at(i);
          auto sourceLink = sp.sourceLinks()[0]
                                .get<const SimSpacePoint*>()
                                ->sourceLinks()[0]
                                .get<IndexSourceLink>();
          return sourceLink.geometryId().value();
        });

    if (bottomSpRange.empty() || middleSpRange.empty() || topSpRange.empty()) {
      continue;
    }

    bottomSp.clear();
    middleSp.clear();
    topSp.clear();

    std::ranges::transform(
        bottomSpRange.begin(), bottomSpRange.end(),
        std::back_inserter(bottomSp),
        [&](Acts::SpacePointIndex2 i) { return sortedCoreSpacePoints[i]; });
    std::ranges::transform(
        middleSpRange.begin(), middleSpRange.end(),
        std::back_inserter(middleSp),
        [&](Acts::SpacePointIndex2 i) { return sortedCoreSpacePoints[i]; });
    std::ranges::transform(
        topSpRange.begin(), topSpRange.end(), std::back_inserter(topSp),
        [&](Acts::SpacePointIndex2 i) { return sortedCoreSpacePoints[i]; });

    m_seedFinder->createSeeds(
        state, cache,
        Acts::TripletSeedFinder2::ContainerPointers(
            coreSpacePoints, rColumn, varianceRColumn, varianceZColumn),
        bottomSp, middleSp, topSp, seeds);
  }

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePoints.size() << " space points");

  // we have seeds of proxies
  // convert them to seed of external space points
  SimSeedContainer seedContainerForStorage;
  seedContainerForStorage.reserve(seeds.size());
  for (const auto& seed : seeds) {
    auto sps = seed.spacePointIndices();
    seedContainerForStorage.emplace_back(*coreSpacePoints.at(sps[0])
                                              .sourceLinks()[0]
                                              .get<const SimSpacePoint*>(),
                                         *coreSpacePoints.at(sps[1])
                                              .sourceLinks()[0]
                                              .get<const SimSpacePoint*>(),
                                         *coreSpacePoints.at(sps[2])
                                              .sourceLinks()[0]
                                              .get<const SimSpacePoint*>());
    seedContainerForStorage.back().setVertexZ(seed.vertexZ());
    seedContainerForStorage.back().setQuality(seed.quality());
  }

  m_outputSeeds(ctx, std::move(seedContainerForStorage));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

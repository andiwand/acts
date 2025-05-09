// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"

#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"
#include "ActsExamples/EventData/Seed.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"

#include <cmath>
#include <csignal>
#include <cstddef>
#include <limits>
#include <stdexcept>

using namespace Acts::HashedStringLiteral;

namespace ActsExamples {

SeedingAlgorithm::SeedingAlgorithm(SeedingAlgorithm::Config cfg,
                                   Acts::Logging::Level lvl)
    : IAlgorithm("SeedingAlgorithm", lvl), m_cfg(std::move(cfg)) {
  // Seed Finder config requires Seed Filter object before conversion to
  // internal units
  m_cfg.seedFilterConfig = m_cfg.seedFilterConfig.toInternalUnits();
  m_cfg.seedFinderConfig.seedFilter = std::make_unique<Acts::SeedFilter>(
      m_cfg.seedFilterConfig, logger().cloneWithSuffix("SeedFilter"));
  m_cfg.seedFinderConfig =
      m_cfg.seedFinderConfig.toInternalUnits().calculateDerivedQuantities();
  m_cfg.seedFinderOptions =
      m_cfg.seedFinderOptions.toInternalUnits().calculateDerivedQuantities(
          m_cfg.seedFinderConfig);
  m_cfg.gridConfig = m_cfg.gridConfig.toInternalUnits();
  m_cfg.gridOptions = m_cfg.gridOptions.toInternalUnits();
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point input collections");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds output collection");
  }

  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputSeeds.initialize(m_cfg.outputSeeds);

  if (m_cfg.gridConfig.rMax != m_cfg.seedFinderConfig.rMax &&
      m_cfg.allowSeparateRMax == false) {
    throw std::invalid_argument(
        "Inconsistent config rMax: using different values in gridConfig and "
        "seedFinderConfig. If values are intentional set allowSeparateRMax to "
        "true");
  }

  if (m_cfg.seedFilterConfig.deltaRMin != m_cfg.seedFinderConfig.deltaRMin) {
    throw std::invalid_argument("Inconsistent config deltaRMin");
  }

  if (m_cfg.gridConfig.deltaRMax != m_cfg.seedFinderConfig.deltaRMax) {
    throw std::invalid_argument("Inconsistent config deltaRMax");
  }

  static_assert(
      std::numeric_limits<
          decltype(m_cfg.seedFinderConfig.deltaRMaxTopSP)>::has_quiet_NaN,
      "Value of deltaRMaxTopSP must support NaN values");

  static_assert(
      std::numeric_limits<
          decltype(m_cfg.seedFinderConfig.deltaRMinTopSP)>::has_quiet_NaN,
      "Value of deltaRMinTopSP must support NaN values");

  static_assert(
      std::numeric_limits<
          decltype(m_cfg.seedFinderConfig.deltaRMaxBottomSP)>::has_quiet_NaN,
      "Value of deltaRMaxBottomSP must support NaN values");

  static_assert(
      std::numeric_limits<
          decltype(m_cfg.seedFinderConfig.deltaRMinBottomSP)>::has_quiet_NaN,
      "Value of deltaRMinBottomSP must support NaN values");

  if (std::isnan(m_cfg.seedFinderConfig.deltaRMaxTopSP)) {
    m_cfg.seedFinderConfig.deltaRMaxTopSP = m_cfg.seedFinderConfig.deltaRMax;
  }

  if (std::isnan(m_cfg.seedFinderConfig.deltaRMinTopSP)) {
    m_cfg.seedFinderConfig.deltaRMinTopSP = m_cfg.seedFinderConfig.deltaRMin;
  }

  if (std::isnan(m_cfg.seedFinderConfig.deltaRMaxBottomSP)) {
    m_cfg.seedFinderConfig.deltaRMaxBottomSP = m_cfg.seedFinderConfig.deltaRMax;
  }

  if (std::isnan(m_cfg.seedFinderConfig.deltaRMinBottomSP)) {
    m_cfg.seedFinderConfig.deltaRMinBottomSP = m_cfg.seedFinderConfig.deltaRMin;
  }

  if (m_cfg.gridConfig.zMin != m_cfg.seedFinderConfig.zMin) {
    throw std::invalid_argument("Inconsistent config zMin");
  }

  if (m_cfg.gridConfig.zMax != m_cfg.seedFinderConfig.zMax) {
    throw std::invalid_argument("Inconsistent config zMax");
  }

  if (m_cfg.seedFilterConfig.maxSeedsPerSpM !=
      m_cfg.seedFinderConfig.maxSeedsPerSpM) {
    throw std::invalid_argument("Inconsistent config maxSeedsPerSpM");
  }

  if (m_cfg.gridConfig.cotThetaMax != m_cfg.seedFinderConfig.cotThetaMax) {
    throw std::invalid_argument("Inconsistent config cotThetaMax");
  }

  if (m_cfg.gridConfig.minPt != m_cfg.seedFinderConfig.minPt) {
    throw std::invalid_argument("Inconsistent config minPt");
  }

  if (m_cfg.gridOptions.bFieldInZ != m_cfg.seedFinderOptions.bFieldInZ) {
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

  if (!m_cfg.seedFinderConfig.zBinsCustomLooping.empty()) {
    // check that the bins required in the custom bin looping
    // are contained in the bins defined by the total number of edges

    for (std::size_t i : m_cfg.seedFinderConfig.zBinsCustomLooping) {
      if (i >= m_cfg.gridConfig.zBinEdges.size()) {
        throw std::invalid_argument(
            "Inconsistent config zBinsCustomLooping does not contain a subset "
            "of bins defined by zBinEdges");
      }
    }
  }

  if (m_cfg.useExtraCuts) {
    // This function will be applied to select space points during grid filling
    m_cfg.seedFinderConfig.spacePointSelector
        .connect<itkFastTrackingSPselect>();

    // This function will be applied to the doublet compatibility selection
    m_cfg.seedFinderConfig.experimentCuts.connect<itkFastTrackingCuts>();
  }

  m_bottomBinFinder = std::make_unique<const Acts::GridBinFinder<3ul>>(
      m_cfg.numPhiNeighbors, cfg.zBinNeighborsBottom, 0);
  m_topBinFinder = std::make_unique<const Acts::GridBinFinder<3ul>>(
      m_cfg.numPhiNeighbors, m_cfg.zBinNeighborsTop, 0);

  m_seedFinder = std::make_unique<SeedFinder>(
      m_cfg.seedFinderConfig, logger().cloneWithSuffix("SeedFinder"));
}

ProcessCode SeedingAlgorithm::execute(const AlgorithmContext& ctx) const {
  const SpacePointContainer& spacePointContainer = m_inputSpacePoints(ctx);

  Acts::CylindricalSpacePointGrid grid =
      Acts::CylindricalSpacePointGridCreator::createGrid(
          m_cfg.gridConfig, m_cfg.gridOptions, logger());

  Acts::CylindricalSpacePointGridCreator::fillGrid(
      m_cfg.seedFinderConfig, grid, spacePointContainer, logger());

  // Compute radius Range
  // we rely on the fact the grid is storing the proxies
  // with a sorting in the radius
  float minRange = std::numeric_limits<float>::max();
  float maxRange = std::numeric_limits<float>::lowest();
  for (const auto& coll : grid) {
    if (coll.empty()) {
      continue;
    }
    auto firstEl = spacePointContainer.at(coll.front());
    auto lastEl = spacePointContainer.at(coll.back());
    minRange = std::min(firstEl.radius(), minRange);
    maxRange = std::max(lastEl.radius(), maxRange);
  }

  std::array<std::vector<std::size_t>, 3ul> navigation;
  navigation[1ul] = m_cfg.seedFinderConfig.zBinsCustomLooping;

  auto spacePointsGrouping =
      Acts::CylindricalBinnedGroup(std::move(grid), *m_bottomBinFinder,
                                   *m_topBinFinder, std::move(navigation));

  /// variable middle SP radial region of interest
  const Acts::Range1D<float> rMiddleSPRange(
      std::floor(minRange / 2) * 2 +
          m_cfg.seedFinderConfig.deltaRMiddleMinSPRange,
      std::floor(maxRange / 2) * 2 -
          m_cfg.seedFinderConfig.deltaRMiddleMaxSPRange);

  // run the seeding
  static thread_local SeedContainer seeds;
  seeds.clear();
  static thread_local SeedFinder::State state;
  state.spacePointMutableData.resize(spacePointContainer.size());

  for (const auto [bottom, middle, top] : spacePointsGrouping) {
    m_seedFinder->createSeedsForGroup(
        spacePointContainer, m_cfg.seedFinderOptions, state,
        spacePointsGrouping.grid(), seeds, bottom, middle, top, rMiddleSPRange);
  }

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePointContainer.size() << " space points");

  m_outputSeeds(ctx, SeedContainer(seeds));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

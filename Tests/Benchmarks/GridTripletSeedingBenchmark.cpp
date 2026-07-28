// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Throughput benchmark for the grid triplet seeder on an ITk-like synthetic
/// event.
///
/// The detector and the event come from SyntheticItkEvent.hpp, which is shared
/// with the other seeding benchmarks. What is added here is the configuration
/// and the driving loop, both taken from how ATLAS runs this seeder on the ITk
/// pixel detector: `ActsTrk::GridTripletSeedingTool` together with the
/// overrides that `ActsPixelSeedingToolCfg` applies to it, and the ITk main
/// pass values of the tracking flags it reads (a seed momentum threshold of
/// 900 MeV and a maximum primary impact parameter of 5 mm).
///
/// Everything the experiment does per event is inside the benchmarked region:
/// filling the grid, sorting each bin in radius, copying into the packed
/// container the seeder consumes, and the seeding itself.

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding/CylindricalSpacePointGrid.hpp"
#include "Acts/Seeding/DoubletSeedFinder.hpp"
#include "Acts/Seeding/TripletSeedFinder.hpp"
#include "Acts/Seeding/TripletSeeder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"
#include "ActsTests/CommonHelpers/BenchmarkTools.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <vector>

#include <boost/program_options.hpp>

#include "SyntheticItkEvent.hpp"

namespace po = boost::program_options;

using namespace Acts::UnitLiterals;
using namespace ActsTests;
using namespace ActsTests::SyntheticItk;

namespace {

/// The ATLAS ITk pixel seeding configuration.
///
/// The values are those of `ActsTrk::GridTripletSeedingTool` in Athena after
/// `ActsPixelSeedingToolCfg` has applied its pixel overrides. Only the
/// properties that reach ACTS are kept; the Athena-side extras (the Hough
/// vertex collision region, the strip measurement details) are left out.
struct ItkPixelConfig {
  // -- shared between the grid, the doublet finders and the filter
  float minPt = 900_MeV;
  float cotThetaMax = 27.2899f;  // eta = 4
  float impactMax = 5_mm;        // ITk main pass maxPrimaryImpactSeed
  float bFieldInZ = 2_T;

  // -- grid
  float gridRMax = 320_mm;
  float zMin = -3000_mm;
  float zMax = 3000_mm;
  float deltaRMax = 280_mm;
  float phiMin = -std::numbers::pi_v<float>;
  float phiMax = std::numbers::pi_v<float>;
  int phiBinDeflectionCoverage = 3;
  int maxPhiBins = 200;
  int numPhiNeighbors = 1;
  std::vector<float> zBinEdges{-3000., -2700., -2500., -1400., -925.,
                               -500.,  -250.,  250.,   500.,   925.,
                               1400.,  2500.,  2700.,  3000.};
  std::vector<float> rBinEdges{0., 320.};
  std::vector<std::size_t> zBinsCustomLooping{2,  3, 4, 5, 12, 11,
                                              10, 9, 7, 6, 8};
  std::vector<std::size_t> rBinsCustomLooping{1};
  std::vector<std::pair<int, int>> zBinNeighborsTop{
      {0, 0}, {-1, 0}, {-2, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {-1, 1},
      {0, 1}, {0, 1},  {0, 1},  {0, 2},  {0, 1},  {0, 0}};
  std::vector<std::pair<int, int>> zBinNeighborsBottom{
      {0, 0},  {0, 1},  {0, 1},  {0, 1},  {0, 1},  {0, 1}, {0, 0},
      {-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {0, 0}};
  std::vector<std::pair<int, int>> rBinNeighborsTop{{0, 0}};
  std::vector<std::pair<int, int>> rBinNeighborsBottom{{0, 0}};

  // -- doublet finders
  float deltaRMinTopSP = 6_mm;
  float deltaRMaxTopSP = 280_mm;
  float deltaRMinBottomSP = 6_mm;
  float deltaRMaxBottomSP = 150_mm;
  // the pixel configuration switches the longitudinal doublet cut off
  float deltaZMax = std::numeric_limits<float>::infinity();
  float collisionRegionMin = -200_mm;
  float collisionRegionMax = 200_mm;
  bool interactionPointCut = true;
  float helixCutTolerance = 1.f;

  // -- triplet finder
  float sigmaScattering = 2.f;
  float radLengthPerSeed = 0.098045f;
  float toleranceParam = 1.1_mm;

  // -- middle space point range
  bool useVariableMiddleSPRange = false;
  float deltaRMiddleMinSPRange = 10.f;
  float deltaRMiddleMaxSPRange = 10.f;
  std::vector<std::pair<float, float>> rRangeMiddleSP{
      {0, 0},    {140, 260}, {40, 260}, {40, 260}, {40, 260},
      {40, 260}, {70, 260},  {40, 260}, {40, 260}, {40, 260},
      {40, 260}, {140, 260}, {0, 0}};

  // -- filter
  float deltaRMin = 20_mm;
  float deltaInvHelixDiameter = 0.00003f;
  float compatSeedWeight = 100.f;
  float impactWeightFactor = 100.f;
  float zOriginWeightFactor = 1.f;
  std::size_t maxSeedsPerSpM = 4;
  std::size_t compatSeedLimit = 3;
  float seedWeightIncrement = 0.f;
  // the pixel configuration disables the increment by putting it out of reach
  float numSeedIncrement = std::numeric_limits<float>::infinity();
  bool seedConfirmation = true;
  std::size_t maxSeedsPerSpMConf = 5;
  std::size_t maxQualitySeedsPerSpMConf = 5;
  bool useDeltaRinsteadOfTopRadius = true;
  float absDeltaEtaWeightFactor = 0.f;
  float absDeltaEtaMinImpact = 2_mm;
  /// Drop seeds that are the best one for none of their three space points
  bool seedQualitySelection = true;
};

Acts::CylindricalSpacePointGrid::Config makeGridConfig(
    const ItkPixelConfig& cfg) {
  Acts::CylindricalSpacePointGrid::Config gridCfg;
  gridCfg.minPt = cfg.minPt;
  gridCfg.rMin = 0;
  gridCfg.rMax = cfg.gridRMax;
  gridCfg.zMin = cfg.zMin;
  gridCfg.zMax = cfg.zMax;
  gridCfg.deltaRMax = cfg.deltaRMax;
  gridCfg.cotThetaMax = cfg.cotThetaMax;
  gridCfg.impactMax = cfg.impactMax;
  gridCfg.phiMin = cfg.phiMin;
  gridCfg.phiMax = cfg.phiMax;
  gridCfg.phiBinDeflectionCoverage = cfg.phiBinDeflectionCoverage;
  gridCfg.maxPhiBins = cfg.maxPhiBins;
  gridCfg.zBinEdges = cfg.zBinEdges;
  gridCfg.rBinEdges = cfg.rBinEdges;
  gridCfg.bFieldInZ = cfg.bFieldInZ;
  gridCfg.bottomBinFinder.emplace(cfg.numPhiNeighbors, cfg.zBinNeighborsBottom,
                                  cfg.rBinNeighborsBottom);
  gridCfg.topBinFinder.emplace(cfg.numPhiNeighbors, cfg.zBinNeighborsTop,
                               cfg.rBinNeighborsTop);
  gridCfg.navigation[0ul] = {};
  gridCfg.navigation[1ul] = cfg.zBinsCustomLooping;
  gridCfg.navigation[2ul] = cfg.rBinsCustomLooping;
  return gridCfg;
}

Acts::DoubletSeedFinder::Config makeBottomDoubletConfig(
    const ItkPixelConfig& cfg) {
  Acts::DoubletSeedFinder::Config doubletCfg;
  doubletCfg.spacePointsSortedByRadius = true;
  doubletCfg.candidateDirection = Acts::Direction::Backward();
  doubletCfg.deltaRMin = cfg.deltaRMinBottomSP;
  doubletCfg.deltaRMax = cfg.deltaRMaxBottomSP;
  doubletCfg.deltaZMin = -cfg.deltaZMax;
  doubletCfg.deltaZMax = cfg.deltaZMax;
  doubletCfg.impactMax = cfg.impactMax;
  doubletCfg.interactionPointCut = cfg.interactionPointCut;
  doubletCfg.collisionRegionMin = cfg.collisionRegionMin;
  doubletCfg.collisionRegionMax = cfg.collisionRegionMax;
  doubletCfg.cotThetaMax = cfg.cotThetaMax;
  doubletCfg.minPt = cfg.minPt;
  doubletCfg.helixCutTolerance = cfg.helixCutTolerance;
  return doubletCfg;
}

Acts::TripletSeedFinder::Config makeTripletConfig(const ItkPixelConfig& cfg) {
  Acts::TripletSeedFinder::Config tripletCfg;
  tripletCfg.useStripInfo = false;
  tripletCfg.sortedByCotTheta = true;
  tripletCfg.minPt = cfg.minPt;
  tripletCfg.sigmaScattering = cfg.sigmaScattering;
  tripletCfg.radLengthPerSeed = cfg.radLengthPerSeed;
  tripletCfg.impactMax = cfg.impactMax;
  tripletCfg.helixCutTolerance = cfg.helixCutTolerance;
  tripletCfg.toleranceParam = cfg.toleranceParam;
  return tripletCfg;
}

Acts::BroadTripletSeedFilter::Config makeFilterConfig(
    const ItkPixelConfig& cfg) {
  Acts::BroadTripletSeedFilter::Config filterCfg;
  filterCfg.deltaInvHelixDiameter = cfg.deltaInvHelixDiameter;
  filterCfg.deltaRMin = cfg.deltaRMin;
  filterCfg.compatSeedWeight = cfg.compatSeedWeight;
  filterCfg.impactWeightFactor = cfg.impactWeightFactor;
  filterCfg.zOriginWeightFactor = cfg.zOriginWeightFactor;
  filterCfg.maxSeedsPerSpM = cfg.maxSeedsPerSpM;
  filterCfg.compatSeedLimit = cfg.compatSeedLimit;
  filterCfg.seedWeightIncrement = cfg.seedWeightIncrement;
  filterCfg.numSeedIncrement = cfg.numSeedIncrement;
  filterCfg.seedConfirmation = cfg.seedConfirmation;
  filterCfg.maxSeedsPerSpMConf = cfg.maxSeedsPerSpMConf;
  filterCfg.maxQualitySeedsPerSpMConf = cfg.maxQualitySeedsPerSpMConf;
  filterCfg.useDeltaRinsteadOfTopRadius = cfg.useDeltaRinsteadOfTopRadius;
  filterCfg.absDeltaEtaWeightFactor = cfg.absDeltaEtaWeightFactor;
  filterCfg.absDeltaEtaMinImpact = cfg.absDeltaEtaMinImpact;

  // central seed confirmation
  filterCfg.centralSeedConfirmationRange.zMinSeedConf = -250_mm;
  filterCfg.centralSeedConfirmationRange.zMaxSeedConf = 250_mm;
  filterCfg.centralSeedConfirmationRange.rMaxSeedConf = 140_mm;
  filterCfg.centralSeedConfirmationRange.nTopForLargeR = 1;
  filterCfg.centralSeedConfirmationRange.nTopForSmallR = 2;
  filterCfg.centralSeedConfirmationRange.seedConfMinBottomRadius = 60_mm;
  filterCfg.centralSeedConfirmationRange.seedConfMaxZOrigin = 150_mm;
  filterCfg.centralSeedConfirmationRange.minImpactSeedConf = 1_mm;
  // forward seed confirmation, identical but without the z restriction
  filterCfg.forwardSeedConfirmationRange =
      filterCfg.centralSeedConfirmationRange;
  filterCfg.forwardSeedConfirmationRange.zMinSeedConf = -3000_mm;
  filterCfg.forwardSeedConfirmationRange.zMaxSeedConf = 3000_mm;

  return filterCfg;
}

/// The radial range middle space point candidates are taken from, looked up in
/// the z bin the candidates sit in.
/// @param cfg the seeding configuration
/// @param z the longitudinal position of the middle space points in mm
/// @param variableRange the range derived from the extent of the event
/// @return the radial range in mm
std::pair<float, float> radiusRangeForMiddle(
    const ItkPixelConfig& cfg, float z,
    const Acts::Range1D<float>& variableRange) {
  if (cfg.useVariableMiddleSPRange) {
    return {variableRange.min(), variableRange.max()};
  }
  auto edge = std::ranges::lower_bound(cfg.zBinEdges, z);
  std::size_t zBin = std::distance(cfg.zBinEdges.begin(), edge);
  if (zBin > 0) {
    --zBin;
  }
  const auto& range = cfg.rRangeMiddleSP.at(zBin);
  return {range.first, range.second};
}

/// Everything that is reused between benchmark iterations, so that the timed
/// region allocates as little as the experiment's own does.
struct SeedingCache {
  Acts::CylindricalSpacePointGrid::Config gridConfig;
  Acts::BroadTripletSeedFilter::Config filterConfig;
  std::shared_ptr<const Acts::DoubletSeedFinder> bottomDoubletFinder;
  std::shared_ptr<const Acts::DoubletSeedFinder> topDoubletFinder;
  std::shared_ptr<const Acts::TripletSeedFinder> tripletFinder;
  Acts::TripletSeeder seeder;
  Acts::TripletSeeder::Cache seederCache;
  std::unique_ptr<const Acts::Logger> logger;
};

/// Run one full seeding pass over an event, the way the ATLAS tool does.
/// @param cfg the seeding configuration
/// @param cache the reusable finders and their state
/// @param spacePoints the space points of the event
/// @param seeds receives the seeds
void createSeeds(const ItkPixelConfig& cfg, SeedingCache& cache,
                 const Acts::SpacePointContainer& spacePoints,
                 Acts::SeedContainer& seeds) {
  Acts::CylindricalSpacePointGrid grid(cache.gridConfig, cache.logger->clone());

  for (std::size_t i = 0; i < spacePoints.size(); ++i) {
    const auto& sp = spacePoints[i];
    grid.insert(i, sp.phi(), sp.z(), sp.r());
  }
  for (std::size_t i = 0; i < grid.numberOfBins(); ++i) {
    std::ranges::sort(grid.at(i),
                      [&](Acts::SpacePointIndex a, Acts::SpacePointIndex b) {
                        return spacePoints[a].r() < spacePoints[b].r();
                      });
  }

  // the seeder works on a packed container laid out in grid bin order
  Acts::SpacePointContainer gridSpacePoints(
      Acts::SpacePointColumns::CopiedFromIndex |
      Acts::SpacePointColumns::PackedXY | Acts::SpacePointColumns::PackedZR |
      Acts::SpacePointColumns::VarianceZ | Acts::SpacePointColumns::VarianceR);
  gridSpacePoints.reserve(grid.numberOfSpacePoints());
  std::vector<Acts::SpacePointIndexRange> gridSpacePointRanges;
  gridSpacePointRanges.reserve(grid.numberOfBins());
  for (std::size_t i = 0; i < grid.numberOfBins(); ++i) {
    const std::uint32_t begin = gridSpacePoints.size();
    for (const Acts::SpacePointIndex spIndex : grid.at(i)) {
      const auto& sp = spacePoints[spIndex];
      auto newSp = gridSpacePoints.createSpacePoint();
      newSp.copiedFromIndex() = sp.index();
      newSp.xy() = std::array<float, 2>{sp.x(), sp.y()};
      newSp.zr() = std::array<float, 2>{sp.z(), sp.r()};
      newSp.varianceZ() = sp.varianceZ();
      newSp.varianceR() = sp.varianceR();
    }
    gridSpacePointRanges.emplace_back(begin, gridSpacePoints.size());
  }

  // radial extent of the event, exploiting that every bin is sorted in radius
  const Acts::Range1D<float> rRange = [&]() -> Acts::Range1D<float> {
    float minRange = std::numeric_limits<float>::max();
    float maxRange = std::numeric_limits<float>::lowest();
    for (const Acts::SpacePointIndexRange& range : gridSpacePointRanges) {
      if (range.first == range.second) {
        continue;
      }
      minRange = std::min(gridSpacePoints[range.first].zr()[1], minRange);
      maxRange = std::max(gridSpacePoints[range.second - 1].zr()[1], maxRange);
    }
    return {minRange, maxRange};
  }();
  const Acts::Range1D<float> variableMiddleRange(
      std::floor(rRange.min() / 2) * 2 + cfg.deltaRMiddleMinSPRange,
      std::floor(rRange.max() / 2) * 2 - cfg.deltaRMiddleMaxSPRange);

  Acts::BroadTripletSeedFilter::State filterState;
  Acts::BroadTripletSeedFilter::Cache filterCache;
  Acts::BroadTripletSeedFilter filter(cache.filterConfig, filterState,
                                      filterCache, *cache.logger);

  std::vector<Acts::SpacePointContainer::ConstRange> bottomSpRanges;
  std::vector<Acts::SpacePointContainer::ConstRange> topSpRanges;

  Acts::SeedContainer candidates;

  for (const auto [bottom, middle, top] : grid.binnedGroup()) {
    const Acts::SpacePointContainer::ConstRange middleSpRange =
        gridSpacePoints.range(gridSpacePointRanges[middle]).asConst();
    if (middleSpRange.empty()) {
      continue;
    }

    bottomSpRanges.clear();
    for (const std::size_t b : bottom) {
      bottomSpRanges.push_back(
          gridSpacePoints.range(gridSpacePointRanges[b]).asConst());
    }
    topSpRanges.clear();
    for (const std::size_t t : top) {
      topSpRanges.push_back(
          gridSpacePoints.range(gridSpacePointRanges[t]).asConst());
    }

    // all middle candidates of a group share a z bin, so this is a per-group
    // lookup
    const std::pair<float, float> middleRange = radiusRangeForMiddle(
        cfg, middleSpRange.front().zr()[0], variableMiddleRange);

    cache.seeder.createSeedsFromGroups(
        cache.seederCache, *cache.bottomDoubletFinder, *cache.topDoubletFinder,
        *cache.tripletFinder, filter, gridSpacePoints, bottomSpRanges,
        middleSpRange, topSpRanges, middleRange, candidates);
  }

  // Keep a seed only if it is the best one found for at least one of its three
  // space points. ATLAS applies this on top of the ACTS filter.
  seeds.assignSpacePointContainer(spacePoints);
  seeds.reserve(candidates.size());
  for (Acts::MutableSeedProxy candidate : candidates) {
    const auto& indices = candidate.spacePointIndices();
    if (cfg.seedQualitySelection) {
      const float quality = candidate.quality();
      const bool best =
          std::ranges::any_of(indices, [&](Acts::SpacePointIndex index) {
            return filterState.bestSeedQualityMap.at(index) <= quality;
          });
      if (!best) {
        continue;
      }
    }
    const std::array<Acts::SpacePointIndex, 3> original{
        gridSpacePoints.at(indices[0]).copiedFromIndex(),
        gridSpacePoints.at(indices[1]).copiedFromIndex(),
        gridSpacePoints.at(indices[2]).copiedFromIndex()};
    auto seed = seeds.createSeed();
    seed.assignSpacePointIndices(original);
    seed.quality() = candidate.quality();
    seed.vertexZ() = candidate.vertexZ();
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  std::size_t pileup = EventConfig{}.pileup;
  float secondaryRate = EventConfig{}.secondaryRate;
  std::size_t numRuns = 10;
  std::string dumpPrefix;
  bool verbose = false;

  try {
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")(
        "pileup", po::value<std::size_t>(&pileup)->default_value(pileup),
        "number of overlaid minimum-bias interactions")(
        "secondary-rate",
        po::value<float>(&secondaryRate)->default_value(secondaryRate),
        "mean number of secondaries produced per primary crossing")(
        "runs", po::value<std::size_t>(&numRuns)->default_value(numRuns),
        "number of benchmark runs")(
        "dump", po::value<std::string>(&dumpPrefix),
        "write the generated event to <arg>_spacepoints.csv and "
        "<arg>_particles.csv")("verbose", po::bool_switch(&verbose),
                               "log the seeder's own statistics");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help") > 0) {
      std::cout << desc << std::endl;
      return 0;
    }
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }

  const DetectorLayout layout = makePixelLayout();

  EventConfig eventConfig;
  eventConfig.pileup = pileup;
  eventConfig.secondaryRate = secondaryRate;
  const Event event = generateEvent(layout, eventConfig);
  if (!dumpPrefix.empty()) {
    writeEventCsv(event, layout, dumpPrefix);
  }

  ItkPixelConfig cfg;
  cfg.bFieldInZ = eventConfig.bFieldZ * 1_T;

  SeedingCache cache;
  cache.logger = Acts::getDefaultLogger(
      "GridTriplet",
      verbose ? Acts::Logging::Level::DEBUG : Acts::Logging::Level::FATAL);
  cache.gridConfig = makeGridConfig(cfg);
  cache.filterConfig = makeFilterConfig(cfg);
  const Acts::DoubletSeedFinder::Config bottomConfig =
      makeBottomDoubletConfig(cfg);
  Acts::DoubletSeedFinder::Config topConfig = bottomConfig;
  topConfig.candidateDirection = Acts::Direction::Forward();
  topConfig.deltaRMin = cfg.deltaRMinTopSP;
  topConfig.deltaRMax = cfg.deltaRMaxTopSP;
  cache.bottomDoubletFinder = Acts::DoubletSeedFinder::create(
      Acts::DoubletSeedFinder::DerivedConfig(bottomConfig, cfg.bFieldInZ));
  cache.topDoubletFinder = Acts::DoubletSeedFinder::create(
      Acts::DoubletSeedFinder::DerivedConfig(topConfig, cfg.bFieldInZ));
  cache.tripletFinder =
      Acts::TripletSeedFinder::create(Acts::TripletSeedFinder::DerivedConfig(
          makeTripletConfig(cfg), cfg.bFieldInZ));
  cache.seeder = Acts::TripletSeeder(cache.logger->clone());

  const EventSummary summary = summarize(event, cfg.minPt / 1_GeV);
  std::cout << "spacePoints=" << summary.spacePoints
            << " primaryHits=" << summary.primaryHits
            << " secondaryHits=" << summary.secondaryHits
            << " primaries=" << summary.primaries
            << " secondaries=" << summary.secondaries
            << " seedable=" << summary.seedablePrimaries << std::endl;

  // matched outside the timed region, so that the truth lookup does not show
  // up in the measurement
  Acts::SeedContainer reference;
  createSeeds(cfg, cache, event.spacePoints, reference);
  const SeedingSummary seedSummary =
      evaluateSeeds(event, reference, cfg.minPt / 1_GeV);

  const auto result = microBenchmark(
      [&] {
        Acts::SeedContainer seeds;
        createSeeds(cfg, cache, event.spacePoints, seeds);
        assumeRead(seeds);
      },
      1, numRuns);

  std::cout << "seeds=" << seedSummary.seeds
            << " trueSeeds=" << seedSummary.trueSeeds << " efficiency="
            << static_cast<float>(seedSummary.matchedPrimaries) /
                   static_cast<float>(
                       std::max<std::size_t>(1, summary.seedablePrimaries))
            << "\n"
            << result << std::endl;

  return 0;
}

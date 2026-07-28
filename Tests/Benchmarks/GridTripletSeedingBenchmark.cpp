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
/// The same loop runs on either of the two space point grids ACTS offers, so
/// that they can be compared on equal terms:
///
///  - `--grid cylindrical` bins in (phi, z, r), which is what ATLAS runs.
///  - `--grid spherical` bins in (phi, eta, r) instead. Since the middle axis
///    now measures a direction rather than a position, its bins group space
///    points the way the doublet cuts do, and far fewer candidate pairs have to
///    be formed and rejected. What it costs is that the eta a space point is
///    binned at, cot(theta) = z / r, is the eta seen from the origin and not
///    the one of the track, so a vertex displaced by z0 moves a space point at
///    radius r by z0 / r. The neighbour windows are sized from that offset,
///    see `makeNeighbors`.
///
/// On the default event the spherical grid finds 89.8% of the seedable
/// primaries in 0.71 s, against 86.0% in 0.92 s for the cylindrical one, and it
/// returns a comparable number of seeds while doing so, 15.7k against 15.2k.
/// Two things had to be right for that. The bins have to be wide enough that a
/// window is one neighbour either side, because the doublet finder is entered
/// once per middle space point and neighbour bin whatever the bin holds, and a
/// finer binning loses more to those entries than it gains in rejected pairs.
/// And the longitudinal doublet cut, which ATLAS switches off because its z
/// bins already impose one, has to be switched back on, since the eta bins do
/// not. Giving the cylindrical grid the same cut does not help it: it costs
/// efficiency there and saves 4% of the time.
///
/// Everything the experiment does per event is inside the benchmarked region:
/// filling the grid, sorting each bin in radius, copying into the packed
/// container the seeder consumes, and the seeding itself. What each grid offers
/// the seeder to work with is reported alongside, so that a tuning run can be
/// read without a profiler: `bottomPairs` and `topPairs` count the space point
/// pairs the binning lets through to the doublet finder, and `ranges` counts
/// how often it is entered to do so.

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding/CylindricalSpacePointGrid.hpp"
#include "Acts/Seeding/DoubletSeedFinder.hpp"
#include "Acts/Seeding/SphericalSpacePointGrid.hpp"
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
#include <string>
#include <type_traits>
#include <vector>

#include <boost/program_options.hpp>

#include "SyntheticItkEvent.hpp"

namespace po = boost::program_options;
namespace Exp = Acts::Experimental;

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

/// The (phi, eta, r) grid, which has no ATLAS counterpart to copy from. The phi
/// axis, the doublet finders, the triplet finder and the filter are shared with
/// the cylindrical setup above; what differs is the middle axis, and that is
/// what the numbers here describe. Every one of them was scanned on the
/// synthetic event, and each is at the point where efficiency stops paying for
/// the time it costs.
struct SphericalGridConfig {
  /// Pseudorapidity acceptance, matching `cotThetaMax`
  float etaMin = -4.f;
  float etaMax = 4.f;
  /// Bin width on the middle axis, in cot(theta) rather than in eta. What sets
  /// it is the offset below and not a resolution: the window of a bin has to
  /// reach that far regardless, so binning any finer only splits one range into
  /// several and enters the doublet finder again for each without rejecting
  /// anything more.
  float cotBinWidth = 0.6f;
  /// Widening of the bins proportional to |cot(theta)|, which would give the
  /// forward region a relative rather than an absolute resolution. It is where
  /// most of the space points are, but every value tried cost more in extra
  /// candidate pairs than it saved in bins, so this is off.
  float cotBinGrowth = 0.f;
  /// Radial bin edges in mm. Splitting the radial range does cut the number of
  /// candidate pairs, by a third for the four bin version, but a bin the
  /// doublet finder is entered for costs the same whether it holds one space
  /// point or a hundred, and the extra entries came out more expensive than the
  /// pairs they saved. Hence a single bin.
  std::vector<float> rBinEdges{0.f, 320.f};
  /// How far a bottom or top space point of the same track may sit from the
  /// middle one on the cot(theta) axis. A track from z0 crosses radius r at
  /// cot(theta) = cot(theta)_track + z0 / r, so two space points of one track
  /// are at most |z0| * (1 / r_inner - 1 / r_outer) apart on this axis, which
  /// `makeNeighbors` turns into per-bin windows. The geometric worst case is
  /// several units for both, but the beam spot is 45 mm long and the bound is
  /// only reached by a track that is both far off centre and seeded at the
  /// innermost radius. What the scan finds is that reaching further than this
  /// towards the tops buys nothing at all, and that the bottoms, which sit at
  /// the small radii where the offset is largest, want about twice the reach.
  float cotSpreadBottom = 1.5f;
  float cotSpreadTop = 0.6f;
  /// Radial range middle candidates are taken from, as in `rRangeMiddleSP`
  float middleRMin = 40.f;
  float middleRMax = 260.f;
  /// Longitudinal reach of middle candidates in mm. The z bins ATLAS loops over
  /// stop short of the outermost ones, which is how it avoids seeding from the
  /// last disk, where nothing can sit beyond a middle space point. A bin of
  /// this grid holds a range of cot(theta) rather than of z, but since z = r *
  /// cot(theta) the same statement is a cap on the radius of the middle
  /// candidates of a bin, which is what `makeSphericalMiddleAxis` turns it
  /// into.
  float middleZMax = 2700.f;
  /// Longitudinal reach of a doublet in mm. The pixel configuration switches
  /// this cut off because with a cylindrical grid the z bin windows already
  /// bound how far apart in z the two space points of a doublet can be. A
  /// spherical grid has no z bins to do that, so the cut has to be made
  /// explicitly again, and it is what the comparison turns on: without it the
  /// seeder returns 28k seeds instead of 15k for the same 7.3k true ones, and
  /// spends longer than the cylindrical grid does producing them. The ATLAS
  /// tool itself defaults to 600 mm, which costs efficiency here; 900 mm is
  /// where it stops doing so.
  float deltaZMax = 900.f;
};

/// Turn a reach along an axis into the per-bin neighbour windows the grid bin
/// finder wants: the window of a bin is the smallest one covering the bin
/// itself grown by `below` and `above`. The bins are not uniform on either
/// axis of the spherical grid, so stating the reach once in the coordinate of
/// the axis is what keeps the windows honest everywhere.
/// @param edges the bin edges of the axis
/// @param below how far below the lower edge of a bin to reach
/// @param above how far above the upper edge of a bin to reach
/// @return one window per bin
std::vector<std::pair<int, int>> makeNeighbors(const std::vector<float>& edges,
                                               float below, float above) {
  const int numBins = static_cast<int>(edges.size()) - 1;
  std::vector<std::pair<int, int>> neighbors(numBins);
  for (int bin = 0; bin < numBins; ++bin) {
    int down = 0;
    while (bin - down > 0 && edges[bin] - below < edges[bin - down]) {
      ++down;
    }
    int up = 0;
    while (bin + up + 1 < numBins &&
           edges[bin + 1] + above > edges[bin + up + 1]) {
      ++up;
    }
    neighbors[bin] = {-down, up};
  }
  return neighbors;
}

/// The binning of the middle axis, in whichever coordinate the grid bins it.
/// The driving loop uses it to look up the radial range middle candidates are
/// taken from without having to know which grid it is running on.
struct MiddleAxis {
  /// Bin edges, in mm for the cylindrical grid and in cot(theta) for the
  /// spherical one
  std::vector<float> edges;
  /// Radial range of middle candidates, one entry per bin
  std::vector<std::pair<float, float>> rRanges;
};

Acts::CylindricalSpacePointGrid::Config makeCylindricalGridConfig(
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

/// The middle axis of the cylindrical grid: its z bins, with the radial range
/// ATLAS takes middle candidates from in each of them.
/// @param cfg the seeding configuration
/// @return the middle axis binning
MiddleAxis makeCylindricalMiddleAxis(const ItkPixelConfig& cfg) {
  return {cfg.zBinEdges, cfg.rRangeMiddleSP};
}

/// The bin edges of the middle axis of the spherical grid, in cot(theta), which
/// is the coordinate space points are binned in. Bins may grow with
/// |cot(theta)|, so they are built outwards from the centre in both directions.
/// @param sph the spherical grid configuration
/// @return the bin edges in cot(theta)
std::vector<float> makeCotThetaBinEdges(const SphericalGridConfig& sph) {
  std::vector<float> positive{0.f};
  const float cotMax = std::sinh(std::max(-sph.etaMin, sph.etaMax));
  while (positive.back() < cotMax) {
    positive.push_back(positive.back() + sph.cotBinWidth +
                       sph.cotBinGrowth * positive.back());
  }

  std::vector<float> edges;
  edges.reserve(2 * positive.size() - 1);
  for (auto it = positive.rbegin(); it != positive.rend(); ++it) {
    edges.push_back(-*it);
  }
  edges.insert(edges.end(), std::next(positive.begin()), positive.end());
  return edges;
}

/// The same edges in pseudorapidity, which is how the grid takes them.
/// @param sph the spherical grid configuration
/// @return the bin edges in eta
std::vector<float> makeEtaBinEdges(const SphericalGridConfig& sph) {
  std::vector<float> edges = makeCotThetaBinEdges(sph);
  for (float& edge : edges) {
    edge = std::asinh(edge);
  }
  return edges;
}

Exp::SphericalSpacePointGrid::Config makeSphericalGridConfig(
    const ItkPixelConfig& cfg, const SphericalGridConfig& sph) {
  Exp::SphericalSpacePointGrid::Config gridCfg;
  gridCfg.minPt = cfg.minPt;
  gridCfg.rMin = 0;
  gridCfg.rMax = cfg.gridRMax;
  gridCfg.etaMin = sph.etaMin;
  gridCfg.etaMax = sph.etaMax;
  gridCfg.etaBinEdges = makeEtaBinEdges(sph);
  gridCfg.deltaRMax = cfg.deltaRMax;
  gridCfg.impactMax = cfg.impactMax;
  gridCfg.phiMin = cfg.phiMin;
  gridCfg.phiMax = cfg.phiMax;
  gridCfg.phiBinDeflectionCoverage = cfg.phiBinDeflectionCoverage;
  gridCfg.maxPhiBins = cfg.maxPhiBins;
  gridCfg.rBinEdges = sph.rBinEdges;
  gridCfg.bFieldInZ = cfg.bFieldInZ;

  // both axes take their windows from the reach of what looks for them: the
  // radial one from how far the doublet finders go, the middle one from how far
  // the beam spot can move a space point along it
  const std::vector<float> cotThetaEdges = makeCotThetaBinEdges(sph);
  gridCfg.bottomBinFinder.emplace(
      cfg.numPhiNeighbors,
      makeNeighbors(cotThetaEdges, sph.cotSpreadBottom, sph.cotSpreadBottom),
      makeNeighbors(sph.rBinEdges, cfg.deltaRMaxBottomSP, 0.f));
  gridCfg.topBinFinder.emplace(
      cfg.numPhiNeighbors,
      makeNeighbors(cotThetaEdges, sph.cotSpreadTop, sph.cotSpreadTop),
      makeNeighbors(sph.rBinEdges, 0.f, cfg.deltaRMaxTopSP));
  gridCfg.navigation[0ul] = {};
  gridCfg.navigation[1ul] = {};
  gridCfg.navigation[2ul] = {};
  return gridCfg;
}

/// The middle axis of the spherical grid. Every bin offers the same range,
/// except that the more forward ones have it cut back so that no middle
/// candidate sits beyond `middleZMax`.
/// @param sph the spherical grid configuration
/// @return the middle axis binning
MiddleAxis makeSphericalMiddleAxis(const SphericalGridConfig& sph) {
  const std::vector<float> edges = makeCotThetaBinEdges(sph);
  std::vector<std::pair<float, float>> ranges;
  ranges.reserve(edges.size() - 1);
  for (std::size_t bin = 0; bin + 1 < edges.size(); ++bin) {
    const float cotMax =
        std::max(std::abs(edges[bin]), std::abs(edges[bin + 1]));
    ranges.emplace_back(sph.middleRMin,
                        std::min(sph.middleRMax, sph.middleZMax / cotMax));
  }
  return {edges, ranges};
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
/// the middle axis bin the candidates sit in.
/// @param cfg the seeding configuration
/// @param axis the binning of the middle axis
/// @param coordinate the middle axis coordinate of the middle space points
/// @param variableRange the range derived from the extent of the event
/// @return the radial range in mm
std::pair<float, float> radiusRangeForMiddle(
    const ItkPixelConfig& cfg, const MiddleAxis& axis, float coordinate,
    const Acts::Range1D<float>& variableRange) {
  if (cfg.useVariableMiddleSPRange) {
    return {variableRange.min(), variableRange.max()};
  }
  auto edge = std::ranges::lower_bound(axis.edges, coordinate);
  std::size_t bin = std::distance(axis.edges.begin(), edge);
  if (bin > 0) {
    --bin;
  }
  return axis.rRanges.at(std::min(bin, axis.rRanges.size() - 1));
}

/// Everything that is reused between benchmark iterations, so that the timed
/// region allocates as little as the experiment's own does.
struct SeedingCache {
  Acts::BroadTripletSeedFilter::Config filterConfig;
  std::shared_ptr<const Acts::DoubletSeedFinder> bottomDoubletFinder;
  std::shared_ptr<const Acts::DoubletSeedFinder> topDoubletFinder;
  std::shared_ptr<const Acts::TripletSeedFinder> tripletFinder;
  Acts::TripletSeeder seeder;
  Acts::TripletSeeder::Cache seederCache;
  std::unique_ptr<const Acts::Logger> logger;
};

/// Run one full seeding pass over an event, the way the ATLAS tool does.
/// @tparam GridType the space point grid to bin the event in
/// @param cfg the seeding configuration
/// @param gridConfig the configuration of the grid
/// @param axis the binning of the middle axis of the grid
/// @param cache the reusable finders and their state
/// @param spacePoints the space points of the event
/// @param seeds receives the seeds
template <typename GridType>
void createSeeds(const ItkPixelConfig& cfg,
                 const typename GridType::Config& gridConfig,
                 const MiddleAxis& axis, SeedingCache& cache,
                 const Acts::SpacePointContainer& spacePoints,
                 Acts::SeedContainer& seeds) {
  constexpr bool spherical =
      std::is_same_v<GridType, Exp::SphericalSpacePointGrid>;

  GridType grid(gridConfig, cache.logger->clone());

  for (std::size_t i = 0; i < spacePoints.size(); ++i) {
    const auto& sp = spacePoints[i];
    if constexpr (spherical) {
      grid.insert(i, sp.phi(), sp.z() / sp.r(), sp.r());
    } else {
      grid.insert(i, sp.phi(), sp.z(), sp.r());
    }
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

    // all middle candidates of a group share a middle axis bin, so this is a
    // per-group lookup
    const std::array<float, 2> front = middleSpRange.front().zr();
    const std::pair<float, float> middleRange = radiusRangeForMiddle(
        cfg, axis, spherical ? front[0] / front[1] : front[0],
        variableMiddleRange);

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

/// Count the space point pairs the grid offers to the doublet finders. This is
/// the combinatorial load a grid configuration imposes, and it is what the two
/// grids have to be compared on: everything downstream only ever sees pairs
/// that were offered here.
/// @tparam GridType the space point grid to bin the event in
/// @param cfg the seeding configuration
/// @param gridConfig the configuration of the grid
/// @param axis the binning of the middle axis of the grid
/// @param spacePoints the space points of the event
/// @param logger the logger to build the grid with
template <typename GridType>
void printGridStatistics(const ItkPixelConfig& cfg,
                         const typename GridType::Config& gridConfig,
                         const MiddleAxis& axis,
                         const Acts::SpacePointContainer& spacePoints,
                         const Acts::Logger& logger) {
  constexpr bool spherical =
      std::is_same_v<GridType, Exp::SphericalSpacePointGrid>;

  GridType grid(gridConfig, logger.clone());
  for (std::size_t i = 0; i < spacePoints.size(); ++i) {
    const auto& sp = spacePoints[i];
    if constexpr (spherical) {
      grid.insert(i, sp.phi(), sp.z() / sp.r(), sp.r());
    } else {
      grid.insert(i, sp.phi(), sp.z(), sp.r());
    }
  }

  std::size_t occupied = 0;
  for (std::size_t i = 0; i < grid.numberOfBins(); ++i) {
    occupied += grid.at(i).empty() ? 0 : 1;
  }

  std::size_t middles = 0;
  std::size_t bottomPairs = 0;
  std::size_t topPairs = 0;
  // the doublet finder is entered once per middle space point and neighbour
  // bin, so this counts the fixed cost the binning imposes independently of how
  // many space points the bins hold
  std::size_t ranges = 0;
  for (const auto [bottom, middle, top] : grid.binnedGroup()) {
    std::size_t numMiddle = 0;
    for (const Acts::SpacePointIndex index : grid.at(middle)) {
      const auto& sp = spacePoints[index];
      const std::pair<float, float> range = radiusRangeForMiddle(
          cfg, axis, spherical ? sp.z() / sp.r() : sp.z(), {0, 0});
      numMiddle += (sp.r() >= range.first && sp.r() <= range.second) ? 1 : 0;
    }
    middles += numMiddle;
    ranges += numMiddle * (bottom.size() + top.size());
    for (const std::size_t b : bottom) {
      bottomPairs += numMiddle * grid.at(b).size();
    }
    for (const std::size_t t : top) {
      topPairs += numMiddle * grid.at(t).size();
    }
  }

  std::cout << "bins=" << grid.numberOfBins() << " occupied=" << occupied
            << " middles=" << middles << " bottomPairs=" << bottomPairs
            << " topPairs=" << topPairs << " ranges=" << ranges << std::endl;
}

}  // namespace

int main(int argc, char* argv[]) {
  std::size_t pileup = EventConfig{}.pileup;
  float secondaryRate = EventConfig{}.secondaryRate;
  std::size_t numRuns = 10;
  float minPt = ItkPixelConfig{}.minPt / 1_MeV;
  // Efficiency is counted over a harder threshold than the seeder is cut at,
  // which keeps the turn-on out of it: a seeder run at 900 MeV loses tracks
  // just above its own threshold for reasons that say nothing about the seeder,
  // and where exactly it loses them moves with the momentum resolution.
  float truthPt = 1000.f;
  std::string gridName = "cylindrical";
  std::string dumpPrefix;
  std::vector<float> middleRange;
  bool allMiddleBins = false;
  bool forceDeltaZMax = false;
  bool verbose = false;

  SphericalGridConfig sph;

  try {
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")(
        "grid", po::value<std::string>(&gridName)->default_value(gridName),
        "space point grid to bin the event in, cylindrical or spherical")(
        "pileup", po::value<std::size_t>(&pileup)->default_value(pileup),
        "number of overlaid minimum-bias interactions")(
        "secondary-rate",
        po::value<float>(&secondaryRate)->default_value(secondaryRate),
        "mean number of secondaries produced per primary crossing")(
        "runs", po::value<std::size_t>(&numRuns)->default_value(numRuns),
        "number of benchmark runs")(
        "min-pt", po::value<float>(&minPt)->default_value(minPt),
        "seed momentum threshold in MeV")(
        "truth-pt", po::value<float>(&truthPt)->default_value(truthPt),
        "momentum in MeV a primary has to reach to be counted in the "
        "efficiency")(
        "cot-bin-width",
        po::value<float>(&sph.cotBinWidth)->default_value(sph.cotBinWidth),
        "spherical grid: width of the middle axis bins in cot(theta)")(
        "cot-bin-growth",
        po::value<float>(&sph.cotBinGrowth)->default_value(sph.cotBinGrowth),
        "spherical grid: widening of the bins per unit of |cot(theta)|")(
        "cot-spread-bottom",
        po::value<float>(&sph.cotSpreadBottom)
            ->default_value(sph.cotSpreadBottom),
        "spherical grid: cot(theta) reach towards bottom candidates")(
        "cot-spread-top",
        po::value<float>(&sph.cotSpreadTop)->default_value(sph.cotSpreadTop),
        "spherical grid: cot(theta) reach towards top candidates")(
        "delta-z-max",
        po::value<float>(&sph.deltaZMax)->default_value(sph.deltaZMax),
        "spherical grid: longitudinal reach of a doublet in mm")(
        "r-bins", po::value<std::vector<float>>(&sph.rBinEdges)->multitoken(),
        "spherical grid: radial bin edges in mm")(
        "middle-range",
        po::value<std::vector<float>>(&middleRange)->multitoken(),
        "override the radial range of middle candidates in every bin of either "
        "grid, so that both see the same middle space points")(
        "all-middle-bins", po::bool_switch(&allMiddleBins),
        "cylindrical grid: drop the custom looping order, which otherwise "
        "leaves the two outermost z bins out of the middle bins entirely")(
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
    // the cut belongs to the spherical setup, but asking for it explicitly
    // applies it to either grid, so that the two can be told apart from it
    forceDeltaZMax = !vm["delta-z-max"].defaulted();
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }

  if (gridName != "cylindrical" && gridName != "spherical") {
    std::cerr << "error: unknown grid '" << gridName << "'" << std::endl;
    return 1;
  }
  const bool spherical = gridName == "spherical";

  // the middle axis is walked outwards one bin at a time, so a width that does
  // not advance would never terminate
  if (sph.cotBinWidth <= 0.f || sph.cotBinGrowth < 0.f) {
    std::cerr << "error: the middle axis bins have to grow" << std::endl;
    return 1;
  }
  if (sph.rBinEdges.size() < 2 || !std::ranges::is_sorted(sph.rBinEdges)) {
    std::cerr << "error: --r-bins takes at least two increasing edges"
              << std::endl;
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
  // reaches the phi binning of the grid, both doublet finders and the triplet
  // finder, so it has to be set before any of them is built
  cfg.minPt = minPt * 1_MeV;
  if (spherical || forceDeltaZMax) {
    cfg.deltaZMax = sph.deltaZMax;
  }

  SeedingCache cache;
  cache.logger = Acts::getDefaultLogger(
      "GridTriplet",
      verbose ? Acts::Logging::Level::DEBUG : Acts::Logging::Level::FATAL);
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

  if (allMiddleBins) {
    cfg.zBinsCustomLooping.clear();
  }
  const Acts::CylindricalSpacePointGrid::Config cylindricalConfig =
      makeCylindricalGridConfig(cfg);
  const Exp::SphericalSpacePointGrid::Config sphericalConfig =
      makeSphericalGridConfig(cfg, sph);
  MiddleAxis axis =
      spherical ? makeSphericalMiddleAxis(sph) : makeCylindricalMiddleAxis(cfg);
  if (middleRange.size() == 2) {
    std::ranges::fill(axis.rRanges,
                      std::pair<float, float>{middleRange[0], middleRange[1]});
  } else if (!middleRange.empty()) {
    std::cerr << "error: --middle-range takes two values" << std::endl;
    return 1;
  }

  // dispatch once, so that the grid type is a compile time property of the
  // seeding pass and not a branch inside it
  const auto seed = [&](Acts::SeedContainer& seeds) {
    if (spherical) {
      createSeeds<Exp::SphericalSpacePointGrid>(
          cfg, sphericalConfig, axis, cache, event.spacePoints, seeds);
    } else {
      createSeeds<Acts::CylindricalSpacePointGrid>(
          cfg, cylindricalConfig, axis, cache, event.spacePoints, seeds);
    }
  };

  if (spherical) {
    printGridStatistics<Exp::SphericalSpacePointGrid>(
        cfg, sphericalConfig, axis, event.spacePoints, *cache.logger);
  } else {
    printGridStatistics<Acts::CylindricalSpacePointGrid>(
        cfg, cylindricalConfig, axis, event.spacePoints, *cache.logger);
  }

  const EventSummary summary = summarize(event, truthPt * 1_MeV / 1_GeV);
  std::cout << "grid=" << gridName << "\nspacePoints=" << summary.spacePoints
            << " primaryHits=" << summary.primaryHits
            << " secondaryHits=" << summary.secondaryHits
            << " primaries=" << summary.primaries
            << " secondaries=" << summary.secondaries
            << " seedable=" << summary.seedablePrimaries << std::endl;

  // matched outside the timed region, so that the truth lookup does not show
  // up in the measurement
  Acts::SeedContainer reference;
  seed(reference);
  const SeedingSummary seedSummary =
      evaluateSeeds(event, reference, truthPt * 1_MeV / 1_GeV);

  const auto result = microBenchmark(
      [&] {
        Acts::SeedContainer seeds;
        seed(seeds);
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

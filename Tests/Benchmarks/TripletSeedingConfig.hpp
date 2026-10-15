// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// The triplet seeding cuts the seeding benchmarks share. The grid and the
/// orthogonal seeder build the same doublet finders, triplet finder and filter
/// and differ only in which space points they offer each other, so sharing the
/// cuts makes a comparison between them one of the binning alone.

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding/DoubletSeedFinder.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Seeding/TripletSeedFinder.hpp"

#include <cstddef>
#include <limits>
#include <numbers>
#include <utility>
#include <vector>

namespace ActsTests::TripletSeedingConfig {

using namespace Acts::UnitLiterals;

/// The ATLAS ITk seeding configuration: the values of
/// `ActsTrk::GridTripletSeedingTool` in Athena, at the ITk main pass values of
/// the tracking flags it reads (900 MeV, 5 mm, 1100 mm). Only the properties
/// that reach ACTS are kept.
///
/// The defaults below are what `ActsPixelSeedingToolCfg` leaves, i.e. the pixel
/// pass; `makeItkStripConfig` applies `ActsStripSeedingToolCfg` instead. Both
/// configure one and the same tool, which is why they are one struct.
struct ItkSeedingConfig {
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
  /// Whether a space point is a pair of strips, in which case the triplet is
  /// fitted with each of its three points free to move along its own strip
  bool useStripInfo = false;
  /// Cheap cut on the doublets of a strip triplet, ahead of that fit
  float cotThetaDiffMax = std::numeric_limits<float>::infinity();

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
  /// Seed confirmation in the central region
  Acts::SeedConfirmationRangeConfig centralSeedConfirmationRange{
      .zMinSeedConf = -250_mm,
      .zMaxSeedConf = 250_mm,
      .rMaxSeedConf = 140_mm,
      .nTopForLargeR = 1,
      .nTopForSmallR = 2,
      .seedConfMinBottomRadius = 60_mm,
      .seedConfMaxZOrigin = 150_mm,
      .minImpactSeedConf = 1_mm};
  /// Seed confirmation in the forward region, identical but without the z
  /// restriction
  Acts::SeedConfirmationRangeConfig forwardSeedConfirmationRange{
      .zMinSeedConf = -3000_mm,
      .zMaxSeedConf = 3000_mm,
      .rMaxSeedConf = 140_mm,
      .nTopForLargeR = 1,
      .nTopForSmallR = 2,
      .seedConfMinBottomRadius = 60_mm,
      .seedConfMaxZOrigin = 150_mm,
      .minImpactSeedConf = 1_mm};

  // -- doublet azimuth cut, which ATLAS adds through the experiment cuts of
  //    the doublet finders; see `DoubletDPhiCut` in the grid benchmark
  bool doubletDPhiCut = false;
  /// Azimuth difference allowed whatever the radii
  float doubletDPhiConst = 0.015f;
  /// Azimuth difference allowed per unit of radial separation, for curvature
  float doubletDPhiSlope = 2e-4f / 1_mm;
  /// Cap on the impact parameter term
  float doubletDPhiCap = 0.10f;
  /// Impact parameter the swing term is computed for; negative takes
  /// `impactMax`
  float doubletDPhiD0Max = -1.f;
};

/// The ATLAS ITk strip seeding configuration, i.e. the same tool after
/// `ActsStripSeedingToolCfg`. The strips reach three times as far out as the
/// pixels, sit four layers deep rather than nine, and each of their space
/// points is a crossing of two strips rather than a point, so what changes is
/// the extent, how far a doublet may reach, and that the triplet is fitted
/// against the strips instead of through the crossings.
/// @return the configuration
inline ItkSeedingConfig makeItkStripConfig() {
  ItkSeedingConfig cfg;

  cfg.impactMax = 20_mm;

  cfg.gridRMax = 1000_mm;
  cfg.deltaRMax = 600_mm;
  cfg.rBinEdges = {0.f, 1100.f};
  // The strips are seeded from the middle of the detector outwards, which the
  // pixel pass does the other way round. Both orders leave the two outermost z
  // bins out of the middle bins entirely.
  cfg.zBinsCustomLooping = {7, 8, 6, 9, 5, 10, 4, 11, 3, 12, 2};
  cfg.zBinNeighborsBottom = {{0, 0},  {0, 1},  {0, 1},  {0, 1},  {0, 2},
                             {0, 1},  {0, 0},  {-1, 0}, {-2, 0}, {-1, 0},
                             {-1, 0}, {-1, 0}, {0, 0}};

  cfg.deltaRMinTopSP = 20_mm;
  cfg.deltaRMaxTopSP = 300_mm;
  cfg.deltaRMinBottomSP = 20_mm;
  cfg.deltaRMaxBottomSP = 300_mm;
  cfg.deltaZMax = 900_mm;
  // A strip space point of a primary can sit far off the line to the beam spot,
  // the crossing moving along the strips, so the cut is left to the triplet
  cfg.interactionPointCut = false;

  cfg.useStripInfo = true;

  // Middle candidates are taken from a range derived from the extent of the
  // event rather than from a table, four layers being too few to write one for
  cfg.useVariableMiddleSPRange = true;
  cfg.deltaRMiddleMinSPRange = 30_mm;
  cfg.deltaRMiddleMaxSPRange = 150_mm;

  cfg.impactWeightFactor = 1.f;
  cfg.compatSeedLimit = 4;
  cfg.seedWeightIncrement = 10100.f;
  cfg.numSeedIncrement = 1.f;
  cfg.useDeltaRinsteadOfTopRadius = false;
  // Seed confirmation is a pixel notion: its ranges are radii inside the pixel
  // detector, and the strips have too few layers to ask for a fourth point.
  cfg.seedConfirmation = false;
  cfg.maxSeedsPerSpMConf = 100;
  cfg.maxQualitySeedsPerSpMConf = 100;
  cfg.seedQualitySelection = false;

  return cfg;
}

/// The ATLAS ITk large radius tracking (LRT) configuration: the same tool after
/// `ActsLargeRadiusStripSeedingToolCfg`, which applies its LRT cuts on top of
/// `ActsStripSeedingToolCfg`, at the values of the LRT pass flags it reads
/// (1000 MeV, 300 mm, 500 mm). LRT seeds on the strips alone, so this is a
/// strip configuration opened up for tracks that do not come from the beam
/// spot.
///
/// The ATLAS pass runs second and on the space points the primary pass left
/// unused; nothing here removes any, so it sees every one of them.
/// @return the configuration
inline ItkSeedingConfig makeItkLrtStripConfig() {
  ItkSeedingConfig cfg = makeItkStripConfig();

  cfg.minPt = 1000_MeV;
  cfg.impactMax = 300_mm;
  cfg.collisionRegionMin = -500_mm;
  cfg.collisionRegionMax = 500_mm;
  // cut on the impact parameter of a doublet, at the opened-up `impactMax`
  cfg.interactionPointCut = true;
  // what is left of the strip acceptance once a track no longer has to
  // point back at the beam spot
  cfg.cotThetaMax = 5.f;

  // only sizes the phi bins of the grid, which the impact parameter widens
  cfg.deltaRMax = 400_mm;
  cfg.zBinEdges = {-3000., -2500., -1400., -910., -500., -250.,
                   250.,   500.,   910.,   1400., 2500., 3000.};
  cfg.zBinsCustomLooping = {2, 10, 3, 9, 6, 4, 8, 5, 7};
  cfg.zBinNeighborsTop = {{0, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {-1, 1},
                          {0, 1}, {0, 1},  {0, 1},  {0, 1},  {0, 0}};
  cfg.zBinNeighborsBottom = {{0, 0},  {0, 1},  {0, 1},  {0, 1},  {0, 1}, {0, 0},
                             {-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {0, 0}};

  cfg.deltaRMinTopSP = 50_mm;
  cfg.deltaRMaxTopSP = 250_mm;
  cfg.deltaRMinBottomSP = 50_mm;
  cfg.deltaRMaxBottomSP = 250_mm;
  cfg.deltaZMax = 850_mm;

  cfg.cotThetaDiffMax = 0.3f;

  // a table again, rather than the extent of the event
  cfg.useVariableMiddleSPRange = false;
  cfg.deltaRMiddleMaxSPRange = 75_mm;
  cfg.rRangeMiddleSP = {{0, 0},     {400, 850}, {500, 800}, {500, 800},
                        {500, 800}, {500, 800}, {500, 800}, {500, 800},
                        {500, 800}, {400, 850}, {0, 0}};

  cfg.doubletDPhiCut = true;
  cfg.doubletDPhiCap = 0.10f;

  cfg.maxSeedsPerSpM = 1;
  cfg.absDeltaEtaWeightFactor = 10.f;
  cfg.absDeltaEtaMinImpact = 2_mm;
  cfg.zOriginWeightFactor = 1.f;
  cfg.seedConfirmation = true;
  cfg.maxSeedsPerSpMConf = 1;
  cfg.maxQualitySeedsPerSpMConf = 1;
  cfg.centralSeedConfirmationRange = {.zMinSeedConf = -1400_mm,
                                      .zMaxSeedConf = 1400_mm,
                                      .rMaxSeedConf = 140_mm,
                                      .nTopForLargeR = 1,
                                      .nTopForSmallR = 0,
                                      .seedConfMinBottomRadius = 0_mm,
                                      .seedConfMaxZOrigin = 1500_mm,
                                      .minImpactSeedConf = 200_mm};
  cfg.forwardSeedConfirmationRange = {.zMinSeedConf = -3000_mm,
                                      .zMaxSeedConf = 3000_mm,
                                      .rMaxSeedConf = 140_mm,
                                      .nTopForLargeR = 1,
                                      .nTopForSmallR = 0,
                                      .seedConfMinBottomRadius = 350_mm,
                                      .seedConfMaxZOrigin = 1500_mm,
                                      .minImpactSeedConf = 200_mm};

  return cfg;
}

inline Acts::DoubletSeedFinder::Config makeBottomDoubletConfig(
    const ItkSeedingConfig& cfg) {
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

inline Acts::TripletSeedFinder::Config makeTripletConfig(
    const ItkSeedingConfig& cfg) {
  Acts::TripletSeedFinder::Config tripletCfg;
  tripletCfg.useStripInfo = cfg.useStripInfo;
  tripletCfg.cotThetaDiffMax = cfg.cotThetaDiffMax;
  tripletCfg.sortedByCotTheta = true;
  tripletCfg.minPt = cfg.minPt;
  tripletCfg.sigmaScattering = cfg.sigmaScattering;
  tripletCfg.radLengthPerSeed = cfg.radLengthPerSeed;
  tripletCfg.impactMax = cfg.impactMax;
  tripletCfg.helixCutTolerance = cfg.helixCutTolerance;
  tripletCfg.toleranceParam = cfg.toleranceParam;
  return tripletCfg;
}

inline Acts::BroadTripletSeedFilter::Config makeFilterConfig(
    const ItkSeedingConfig& cfg) {
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

  filterCfg.centralSeedConfirmationRange = cfg.centralSeedConfirmationRange;
  filterCfg.forwardSeedConfirmationRange = cfg.forwardSeedConfirmationRange;

  return filterCfg;
}

}  // namespace ActsTests::TripletSeedingConfig

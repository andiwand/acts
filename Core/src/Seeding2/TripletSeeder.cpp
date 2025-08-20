// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/TripletSeeder.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding2/DoubletSeedFinder.hpp"
#include "Acts/Seeding2/TripletSeedFinder.hpp"

#include <Eigen/Dense>

#include <chrono>

namespace Acts::Experimental {

namespace {

template <typename DoubletCollections>
void createAndFilterTriplets(TripletSeeder::Cache& cache,
                             const TripletSeedFinder& tripletFinder,
                             const ITripletSeedFilter& filter,
                             const SpacePointContainer2& spacePoints,
                             DoubletCollections bottomDoublets,
                             const ConstSpacePointProxy2& spM,
                             DoubletCollections topDoublets) {
  for (auto bottomDoublet : bottomDoublets) {
    if (topDoublets.empty()) {
      break;
    }

    cache.tripletTopCandidates.clear();
    tripletFinder.createTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                             topDoublets,
                                             cache.tripletTopCandidates);

    cache.counterTriplet0 += cache.tripletTopCandidates.counterTriplet0;
    cache.counterTriplet1 += cache.tripletTopCandidates.counterTriplet1;
    cache.counterTriplet2 += cache.tripletTopCandidates.counterTriplet2;
    cache.counterTriplet3 += cache.tripletTopCandidates.counterTriplet3;
    cache.counterTriplet4 += cache.tripletTopCandidates.counterTriplet4;
    cache.counterTriplet5 += cache.tripletTopCandidates.counterTriplet5;
    cache.counterTriplet6 += cache.tripletTopCandidates.counterTriplet6;
    cache.counterTriplet7 += cache.tripletTopCandidates.counterTriplet7;
    cache.counterTriplet8 += cache.tripletTopCandidates.counterTriplet8;
    cache.counterTriplet9 += cache.tripletTopCandidates.counterTriplet9;
    cache.counterTriplet10 += cache.tripletTopCandidates.counterTriplet10;

    auto startTime1 = std::chrono::high_resolution_clock::now();

    filter.filterTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                      cache.tripletTopCandidates);

    auto endTime1 = std::chrono::high_resolution_clock::now();
    cache.totalFilter1Time += std::chrono::duration<double, std::nano>(endTime1 - startTime1).count() * 1e-9;
  }
}

template <typename SpacePointCollections>
void createSeedsFromGroupsImpl(
    const Logger& logger, TripletSeeder::Cache& cache,
    const DoubletSeedFinder& bottomFinder, const DoubletSeedFinder& topFinder,
    const TripletSeedFinder& tripletFinder, const ITripletSeedFilter& filter,
    const SpacePointContainer2& spacePoints,
    SpacePointCollections& bottomSpGroups,
    const ConstSpacePointProxy2& middleSp, SpacePointCollections& topSpGroups,
    SeedContainer2& outputSeeds) {
  MiddleSpInfo middleSpInfo = DoubletSeedFinder::computeMiddleSpInfo(middleSp);

  auto startTime1 = std::chrono::high_resolution_clock::now();

  // create middle-top doublets
  cache.topDoublets.clear();
  for (auto& topSpGroup : topSpGroups) {
    topFinder.createDoublets(middleSp, middleSpInfo, topSpGroup,
                             cache.topDoublets);
  }

  auto endTime1 = std::chrono::high_resolution_clock::now();
  cache.totalTopDoubletTime += std::chrono::duration<double, std::nano>(endTime1 - startTime1).count() * 1e-9;

  // no top SP found -> cannot form any triplet
  if (cache.topDoublets.empty()) {
    ACTS_VERBOSE("No compatible Tops, returning");
    return;
  }

  if (!filter.sufficientTopDoublets(spacePoints, middleSp, cache.topDoublets)) {
    return;
  }

  auto startTime2 = std::chrono::high_resolution_clock::now();

  // create middle-bottom doublets
  cache.bottomDoublets.clear();
  for (auto& bottomSpGroup : bottomSpGroups) {
    bottomFinder.createDoublets(middleSp, middleSpInfo, bottomSpGroup,
                                cache.bottomDoublets);
  }

  auto endTime2 = std::chrono::high_resolution_clock::now();
  cache.totalBottomDoubletTime += std::chrono::duration<double, std::nano>(endTime2 - startTime2).count() * 1e-9;

  // no bottom SP found -> cannot form any triplet
  if (cache.bottomDoublets.empty()) {
    ACTS_VERBOSE("No compatible Bottoms, returning");
    return;
  }

  ACTS_VERBOSE("Candidates: " << cache.bottomDoublets.size() << " bottoms and "
                              << cache.topDoublets.size()
                              << " tops for middle candidate indexed "
                              << middleSp.index());

  auto startTime3 = std::chrono::high_resolution_clock::now();

  // combine doublets to triplets
  if (tripletFinder.config().sortedByCotTheta) {
    cache.bottomDoublets.sortByCotTheta({0, cache.bottomDoublets.size()},
                                        cache.sortedBottoms);
    cache.topDoublets.sortByCotTheta({0, cache.topDoublets.size()},
                                     cache.sortedTops);

    createAndFilterTriplets(cache, tripletFinder, filter, spacePoints,
                            cache.bottomDoublets.subset(cache.sortedBottoms),
                            middleSp,
                            cache.topDoublets.subset(cache.sortedTops));
  } else {
    createAndFilterTriplets(cache, tripletFinder, filter, spacePoints,
                            cache.bottomDoublets.range(), middleSp,
                            cache.topDoublets.range());
  }

  auto endTime3 = std::chrono::high_resolution_clock::now();
  cache.totalTripletTime += std::chrono::duration<double, std::nano>(endTime3 - startTime3).count() * 1e-9;

  auto startTime4 = std::chrono::high_resolution_clock::now();

  filter.filterTripletsMiddleFixed(spacePoints, outputSeeds);

  auto endTime4 = std::chrono::high_resolution_clock::now();
  cache.totalFilter2Time += std::chrono::duration<double, std::nano>(endTime4 - startTime4).count() * 1e-9;
}

}  // namespace

TripletSeeder::TripletSeeder(std::unique_ptr<const Logger> logger_)
    : m_logger(std::move(logger_)) {
  if (m_logger == nullptr) {
    throw std::invalid_argument("TripletSeeder: logger cannot be null");
  }
}

void TripletSeeder::createSeedsFromGroup(
    Cache& cache, const DoubletSeedFinder& bottomFinder,
    const DoubletSeedFinder& topFinder, const TripletSeedFinder& tripletFinder,
    const ITripletSeedFilter& filter, const SpacePointContainer2& spacePoints,
    SpacePointContainer2::ConstSubset& bottomSps,
    const ConstSpacePointProxy2& middleSp,
    SpacePointContainer2::ConstSubset& topSps,
    SeedContainer2& outputSeeds) const {
  assert((bottomFinder.config().spacePointsSortedByRadius ==
          topFinder.config().spacePointsSortedByRadius) &&
         "Inconsistent space point sorting");

  std::array<SpacePointContainer2::ConstSubset, 1> bottomSpGroups{bottomSps};
  std::array<SpacePointContainer2::ConstSubset, 1> topSpGroups{topSps};

  createSeedsFromGroupsImpl(*m_logger, cache, bottomFinder, topFinder,
                            tripletFinder, filter, spacePoints, bottomSpGroups,
                            middleSp, topSpGroups, outputSeeds);
}

void TripletSeeder::createSeedsFromGroups(
    Cache& cache, const DoubletSeedFinder& bottomFinder,
    const DoubletSeedFinder& topFinder, const TripletSeedFinder& tripletFinder,
    const ITripletSeedFilter& filter, const SpacePointContainer2& spacePoints,
    const std::span<SpacePointContainer2::ConstRange>& bottomSpGroups,
    const SpacePointContainer2::ConstRange& middleSpGroup,
    const std::span<SpacePointContainer2::ConstRange>& topSpGroups,
    const std::pair<float, float>& radiusRangeForMiddle,
    SeedContainer2& outputSeeds) const {
  assert((bottomFinder.config().spacePointsSortedByRadius ==
          topFinder.config().spacePointsSortedByRadius) &&
         "Inconsistent space point sorting");
  const bool spacePointsSortedByRadius =
      bottomFinder.config().spacePointsSortedByRadius;

  if (middleSpGroup.empty()) {
    return;
  }

  if (spacePointsSortedByRadius) {
    // Initialize initial offsets for bottom and top space points with binary
    // search. This requires at least one middle space point to be present which
    // is already checked above.
    const ConstSpacePointProxy2 firstMiddleSp = middleSpGroup.front();
    const float firstMiddleSpR = firstMiddleSp.zr()[1];

    for (auto& bottomSpGroup : bottomSpGroups) {
      // Find the first bottom space point that is within the deltaRMax of the
      // first middle space point.
      auto low = std::ranges::lower_bound(
          bottomSpGroup, firstMiddleSpR - bottomFinder.config().deltaRMax, {},
          [&](const ConstSpacePointProxy2& sp) { return sp.zr()[1]; });
      bottomSpGroup = bottomSpGroup.subrange(low - bottomSpGroup.begin());
    }

    for (auto& topSpGroup : topSpGroups) {
      // Find the first top space point that is within the deltaRMin of the
      // first middle space point.
      auto low = std::ranges::lower_bound(
          topSpGroup, firstMiddleSpR + topFinder.config().deltaRMin, {},
          [&](const ConstSpacePointProxy2& sp) { return sp.zr()[1]; });
      topSpGroup = topSpGroup.subrange(low - topSpGroup.begin());
    }
  }

  for (ConstSpacePointProxy2 spM : middleSpGroup) {
    const float rM = spM.zr()[1];

    if (spacePointsSortedByRadius) {
      // check if spM is outside our radial region of interest
      if (rM < radiusRangeForMiddle.first) {
        continue;
      }
      if (rM > radiusRangeForMiddle.second) {
        // break because SPs are sorted in r
        break;
      }
    }

    createSeedsFromGroupsImpl(*m_logger, cache, bottomFinder, topFinder,
                              tripletFinder, filter, spacePoints,
                              bottomSpGroups, spM, topSpGroups, outputSeeds);
  }
}

}  // namespace Acts::Experimental

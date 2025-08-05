// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/SeedFilter.hpp"

#include "Acts/EventData/Seed.hpp"
#include "Acts/Seeding/detail/UtilityFunctions.hpp"

#include <algorithm>
#include <utility>

namespace Acts {

template <typename external_spacepoint_t>
SeedFilter<external_spacepoint_t>::SeedFilter(
    const SeedFilterConfig& config,
    IExperimentCuts<external_spacepoint_t>* expCuts /* = 0*/)
    : m_cfg(config), m_experimentCuts(expCuts) {}

template <typename external_spacepoint_t>
SeedFilter<external_spacepoint_t>::SeedFilter(
    const SeedFilterConfig& config, std::unique_ptr<const Logger> logger,
    IExperimentCuts<external_spacepoint_t>* expCuts /* = 0*/)
    : m_cfg(config), m_logger(std::move(logger)), m_experimentCuts(expCuts) {}

template <typename external_spacepoint_t>
void SeedFilter<external_spacepoint_t>::filterSeeds_2SpFixed(
    const SpacePointMutableData& mutableData,
    const external_spacepoint_t& bottomSp,
    const external_spacepoint_t& middleSp,
    const std::vector<const external_spacepoint_t*>& topSpVec,
    const std::vector<float>& invHelixDiameterVec,
    const std::vector<float>& impactParametersVec,
    SeedFilterState& seedFilterState,
    CandidatesForMiddleSp<const external_spacepoint_t>& candidatesCollector)
    const {
  // seed confirmation
  SeedConfirmationRangeConfig seedConfRange;
  if (m_cfg.seedConfirmation) {
    // check if bottom SP is in the central or forward region
    seedConfRange =
        (bottomSp.z() > m_cfg.centralSeedConfirmationRange.zMaxSeedConf ||
         bottomSp.z() < m_cfg.centralSeedConfirmationRange.zMinSeedConf)
            ? m_cfg.forwardSeedConfirmationRange
            : m_cfg.centralSeedConfirmationRange;
    // set the minimum number of top SP depending on whether the bottom SP is
    // in the central or forward region
    seedFilterState.nTopSeedConf =
        bottomSp.radius() > seedConfRange.rMaxSeedConf
            ? seedConfRange.nTopForLargeR
            : seedConfRange.nTopForSmallR;
  }

  std::size_t maxWeightSeedIndex = 0;
  bool maxWeightSeed = false;
  float weightMax = std::numeric_limits<float>::lowest();
  float zOrigin = seedFilterState.zOrigin;

  // initialize original index locations
  std::vector<std::size_t> topSpIndexVec(topSpVec.size());
  for (std::size_t i(0); i < topSpIndexVec.size(); ++i) {
    topSpIndexVec[i] = i;
  }

  if (topSpVec.size() > 2) {
    // sort indexes based on comparing values in invHelixDiameterVec
    std::ranges::sort(topSpIndexVec, {},
                      [&invHelixDiameterVec](const std::size_t t) {
                        return invHelixDiameterVec[t];
                      });
  }

  std::size_t beginCompTopIndex = 0;
  // loop over top SPs and other compatible top SP candidates
  for (const std::size_t topSpIndex : topSpIndexVec) {
    float invHelixDiameter = invHelixDiameterVec[topSpIndex];
    float lowerLimitCurv = invHelixDiameter - m_cfg.deltaInvHelixDiameter;
    float upperLimitCurv = invHelixDiameter + m_cfg.deltaInvHelixDiameter;
    // use deltaR instead of top radius
    float currentTopR = m_cfg.useDeltaRorTopRadius
                            ? mutableData.deltaR(topSpVec[topSpIndex]->index())
                            : topSpVec[topSpIndex]->radius();
    float impact = impactParametersVec[topSpIndex];

    float weight = -(impact * m_cfg.impactWeightFactor);

    // loop over compatible top SP candidates
    for (std::size_t variableCompTopIndex = beginCompTopIndex;
         variableCompTopIndex < topSpIndexVec.size(); variableCompTopIndex++) {
      std::size_t compatibletopSpIndex = topSpIndexVec[variableCompTopIndex];
      if (compatibletopSpIndex == topSpIndex) {
        continue;
      }

      float otherTopR =
          m_cfg.useDeltaRorTopRadius
              ? mutableData.deltaR(topSpVec[compatibletopSpIndex]->index())
              : topSpVec[compatibletopSpIndex]->radius();

      // curvature difference within limits?
      if (invHelixDiameterVec[compatibletopSpIndex] < lowerLimitCurv) {
        // the SPs are sorted in curvature so we skip unnecessary iterations
        beginCompTopIndex = variableCompTopIndex + 1;
        continue;
      }
      if (invHelixDiameterVec[compatibletopSpIndex] > upperLimitCurv) {
        // the SPs are sorted in curvature so we skip unnecessary iterations
        break;
      }
      // compared top SP should have at least deltaRMin distance
      float deltaR = currentTopR - otherTopR;
      if (std::abs(deltaR) < m_cfg.deltaRMin) {
        continue;
      }
      weight += m_cfg.compatSeedWeight;
      break;
    }

    if (m_experimentCuts != nullptr) {
      // add detector specific considerations on the seed weight
      weight += m_experimentCuts->seedWeight(bottomSp, middleSp,
                                             *topSpVec[topSpIndex]);
      // discard seeds according to detector specific cuts (e.g.: weight)
      if (!m_experimentCuts->singleSeedCut(weight, bottomSp, middleSp,
                                           *topSpVec[topSpIndex])) {
        continue;
      }
    }

    // keep the normal behavior without seed quality confirmation
    // if we have not yet reached our max number of seeds we add the new seed
    // to outCont

    candidatesCollector.push(bottomSp, middleSp, *topSpVec[topSpIndex],
                              weight, zOrigin, false);
  }  // loop on tops
  // if no high quality seed was found for a certain middle+bottom SP pair,
  // lower quality seeds can be accepted
  if (m_cfg.seedConfirmation && maxWeightSeed &&
      candidatesCollector.nHighQualityCandidates() == 0) {
    // if we have not yet reached our max number of seeds we add the new seed to
    // outCont

    candidatesCollector.push(bottomSp, middleSp, *topSpVec[maxWeightSeedIndex],
                             weightMax, zOrigin, false);
  }
}

// after creating all seeds with a common middle space point, filter again

template <typename external_spacepoint_t>
template <typename collection_t>
void SeedFilter<external_spacepoint_t>::filterSeeds_1SpFixed(
    SpacePointMutableData& mutableData,
    CandidatesForMiddleSp<const external_spacepoint_t>& candidatesCollector,
    collection_t& outputCollection) const {
  // retrieve all candidates
  // this collection is already sorted
  // higher weights first
  std::size_t numQualitySeeds = candidatesCollector.nHighQualityCandidates();
  auto extended_collection = candidatesCollector.storage();
  filterSeeds_1SpFixed(mutableData, extended_collection, numQualitySeeds,
                       outputCollection);
}

template <typename external_spacepoint_t>
template <typename collection_t>
void SeedFilter<external_spacepoint_t>::filterSeeds_1SpFixed(
    SpacePointMutableData& mutableData,
    std::vector<typename CandidatesForMiddleSp<
        const external_spacepoint_t>::value_type>& candidates,
    const std::size_t numQualitySeeds, collection_t& outputCollection) const {
  if (m_experimentCuts != nullptr) {
    candidates = m_experimentCuts->cutPerMiddleSP(std::move(candidates));
  }

  unsigned int maxSeeds = candidates.size();

  if (maxSeeds > m_cfg.maxSeedsPerSpM) {
    maxSeeds = m_cfg.maxSeedsPerSpM + 1;
  }

  // default filter removes the last seeds if maximum amount exceeded
  // ordering by weight by filterSeeds_2SpFixed means these are the lowest
  // weight seeds
  unsigned int numTotalSeeds = 0;
  for (const auto& [bottom, medium, top, bestSeedQuality, zOrigin,
                    qualitySeed] : candidates) {
    // stop if we reach the maximum number of seeds
    if (numTotalSeeds >= maxSeeds) {
      break;
    }

    if (m_cfg.seedConfirmation) {
      // continue if higher-quality seeds were found
      if (numQualitySeeds > 0 && !qualitySeed) {
        continue;
      }
      if (bestSeedQuality < mutableData.quality(bottom->index()) &&
          bestSeedQuality < mutableData.quality(medium->index()) &&
          bestSeedQuality < mutableData.quality(top->index())) {
        continue;
      }
    }

    // set quality of seed components
    mutableData.setQuality(bottom->index(), bestSeedQuality);
    mutableData.setQuality(medium->index(), bestSeedQuality);
    mutableData.setQuality(top->index(), bestSeedQuality);

    Seed<external_spacepoint_t> seed{*bottom, *medium, *top};
    seed.setVertexZ(zOrigin);
    seed.setQuality(bestSeedQuality);

    ACTS_VERBOSE("Adding seed: [b=" << bottom->index() << ", m="
                                    << medium->index() << ", t=" << top->index()
                                    << "], quality=" << bestSeedQuality
                                    << ", vertexZ=" << zOrigin);
    detail::pushBackOrInsertAtEnd(outputCollection, std::move(seed));
    ++numTotalSeeds;
  }
  ACTS_VERBOSE("Identified " << numTotalSeeds << " seeds");
}

}  // namespace Acts

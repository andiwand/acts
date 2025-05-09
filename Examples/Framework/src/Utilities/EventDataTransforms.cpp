// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/EventDataTransforms.hpp"

#include "Acts/EventData/SpacePointContainer.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"

#include <algorithm>
#include <optional>

ActsExamples::ProtoTrack ActsExamples::seedToPrototrack(
    const ActsExamples::SeedProxy& seed,
    const SpacePointContainer& spacePointContainer) {
  ProtoTrack track;
  track.reserve(seed.size());
  for (auto spacePoint : seed.spacePoints(spacePointContainer)) {
    for (const auto& slink : spacePoint.sourceLinks()) {
      const auto& islink = slink.get<IndexSourceLink>();
      track.emplace_back(islink.index());
    }
  }
  return track;
}

std::optional<ActsExamples::SpacePointProxy>
ActsExamples::findSpacePointByMeasurementIndex(
    ActsExamples::Index measurementIndex,
    const SpacePointContainer& spacePointContainer) {
  auto match = [&](const SpacePointProxy& sp) {
    auto sourceLinks = sp.sourceLinks();
    return std::ranges::any_of(sourceLinks, [&](const Acts::SourceLink& sl) {
      return sl.template get<IndexSourceLink>().index() == measurementIndex;
    });
  };

  auto found = std::ranges::find_if(spacePointContainer, match);

  if (found == spacePointContainer.end()) {
    return std::nullopt;
  }

  return *found;
}

ActsExamples::SeedProxy ActsExamples::prototrackToSeed(
    const ActsExamples::ProtoTrack& track,
    const ActsExamples::SpacePointContainer& spacePointContainer,
    SeedContainer& seedContainer) {
  auto findSpacePointIndex =
      [&](ActsExamples::Index measurementIndex) -> Acts::SpacePointIndex {
    auto found =
        findSpacePointByMeasurementIndex(measurementIndex, spacePointContainer);
    if (!found.has_value()) {
      throw std::runtime_error("No spacepoint found for measurement index " +
                               std::to_string(measurementIndex));
    }
    return found->index();
  };

  const auto s = track.size();
  if (s < 3) {
    throw std::runtime_error(
        "Cannot convert track with less then 3 spacePoints to seed");
  }

  std::vector<Acts::SpacePointIndex> spacePointIndices;
  spacePointIndices.reserve(track.size());

  std::transform(track.begin(), track.end(),
                 std::back_inserter(spacePointIndices), findSpacePointIndex);
  std::ranges::sort(spacePointIndices, {}, [&](const auto& spacePointIndex) {
    return spacePointContainer.at(spacePointIndex).radius();
  });

  // Simply use r = m*z + t and solve for r=0 to find z vertex position...
  // Probably not the textbook way to do
  auto front = spacePointContainer.at(spacePointIndices.front());
  auto back = spacePointContainer.at(spacePointIndices.back());
  const auto m = (back.radius() - front.radius()) / (back.z() - front.z());
  const auto t = front.radius() - m * front.z();
  const auto z_vertex = -t / m;

  auto seed = seedContainer.makeSeed(
      spacePointIndices[0], spacePointIndices[s / 2], spacePointIndices[s - 1]);
  seed.vertexZ() = z_vertex;
  return SeedProxy(seed);
}

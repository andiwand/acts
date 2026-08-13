// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/SyntheticEvent.hpp"

#include <cstdint>

namespace ActsFatras::Synthetic {

namespace {

/// Copy one collection onto the end of another. The columns the generator adds
/// beyond the standard ones are carried across by hand, `copyFrom` knowing only
/// the standard ones.
/// @param out the collection to append to, with the columns already created
/// @param in the collection to take
void append(Acts::SpacePointContainer& out,
            const Acts::SpacePointContainer& in) {
  const Acts::SpacePointColumns columns =
      kSpacePointColumns & ~Acts::SpacePointColumns::CopiedFromIndex;
  auto outLayer = out.column<std::uint32_t>("layerId");
  auto outParticle = out.column<std::uint32_t>("particleId");
  auto outWidth = out.column<float>("clusterWidth");
  auto outLocalY = out.column<float>("localPositionY");
  const auto inLayer = in.column<std::uint32_t>("layerId");
  const auto inParticle = in.column<std::uint32_t>("particleId");
  const auto inWidth = in.column<float>("clusterWidth");
  const auto inLocalY = in.column<float>("localPositionY");
  const bool strips =
      in.hasColumns(Acts::SpacePointColumns::StripCalibrationDetails);

  for (const auto sp : in) {
    auto copy = out.createSpacePoint();
    out.copyFrom(copy.index(), in, sp.index(), columns);
    // not the index in `out`: what a caller wants of it is the point this one
    // was made from, which sits in the collection it came from
    copy.copiedFromIndex() = sp.index();
    copy.extra(outLayer) = sp.extra(inLayer);
    copy.extra(outParticle) = sp.extra(inParticle);
    copy.extra(outWidth) = sp.extra(inWidth);
    copy.extra(outLocalY) = sp.extra(inLocalY);
    if (strips) {
      copy.outerStripCalibrationDetails() = sp.outerStripCalibrationDetails();
    }
  }
}

}  // namespace

Acts::SpacePointContainer selectSpacePoints(
    const Event& event, const SpacePointSelection selection) {
  const bool pixel = selection != SpacePointSelection::Strip;
  const bool strip = selection != SpacePointSelection::Pixel;

  Acts::SpacePointContainer out;
  out.createColumns(strip ? kSpacePointColumns |
                                Acts::SpacePointColumns::StripCalibrationDetails
                          : kSpacePointColumns);
  out.createColumn<std::uint32_t>("layerId");
  out.createColumn<std::uint32_t>("particleId");
  out.createColumn<float>("clusterWidth");
  out.createColumn<float>("localPositionY");
  out.reserve(
      static_cast<std::uint32_t>((pixel ? event.spacePoints.size() : 0) +
                                 (strip ? event.stripSpacePoints.size() : 0)));

  if (pixel) {
    append(out, event.spacePoints);
  }
  if (strip) {
    append(out, event.stripSpacePoints);
  }
  return out;
}

}  // namespace ActsFatras::Synthetic

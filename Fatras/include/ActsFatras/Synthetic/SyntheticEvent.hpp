// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"

#include <cstdint>
#include <vector>

namespace ActsFatras::Synthetic {

/// A generated particle, kept so that a caller can report the truth content of
/// an event without a truth-matching stage.
struct GeneratedParticle {
  /// Transverse momentum
  float pt{};
  /// Pseudorapidity
  float eta{};
  /// Azimuth of the momentum at the perigee, in [-pi, pi]
  float phi{};
  /// Transverse impact parameter with respect to the beam axis
  float d0{};
  /// Longitudinal position at the perigee
  float z0{};
  /// Charge, in units of the elementary one
  float charge{};
  /// Radius at which the particle was produced, zero for a primary
  float productionRadius{};
  /// Longitudinal position at which the particle was produced, zero for a
  /// primary
  float productionZ{};
  /// Index into `DetectorLayout::surfaces` of the surface it was produced on,
  /// `kNoSurface` for one from the luminous region or from a decay, both of
  /// which happen away from any surface
  std::uint32_t productionSurface{kNoSurface};
  /// Number of space points the particle left in the detector
  std::uint32_t numHits{};
  /// How many interactions back the luminous region is: zero for a particle
  /// out of it, one for what such a particle made, and so on. See
  /// `SecondaryConfig::maxGenerations`.
  std::uint8_t generation{};

  /// Whether the particle comes from the luminous region rather than from a
  /// material interaction or a decay.
  /// @return whether it is a primary
  bool primary() const { return generation == 0; }
};

/// The standard columns both collections carry; the strip one adds
/// `Acts::SpacePointColumns::StripCalibrationDetails` on top of them.
constexpr Acts::SpacePointColumns kSpacePointColumns =
    Acts::SpacePointColumns::CopiedFromIndex | Acts::SpacePointColumns::X |
    Acts::SpacePointColumns::Y | Acts::SpacePointColumns::Z |
    Acts::SpacePointColumns::R | Acts::SpacePointColumns::Phi |
    Acts::SpacePointColumns::VarianceZ | Acts::SpacePointColumns::VarianceR;

/// A generated event.
///
/// Beyond the standard columns and the variances a triplet seeder needs, the
/// space points carry `layerId` into `DetectorLayout::layers`, `particleId`
/// into `Event::particles`, and the `clusterWidth` and `localPositionY` the
/// GBTS seeder expects.
///
/// The two collections are kept apart rather than flagged within one, because
/// what separates them is a *column*: only a strip space point has a stereo
/// pair to carry, and a consumer that wants the pair should not have to test
/// every point for whether its entry means anything. Both index the same
/// `particles`, and their `layerId`s are the same index space. A seeder takes
/// one container, so see `selectSpacePoints`.
struct Event {
  /// Space points of the pixel layers, i.e. of every layer the description
  /// gives no `StripSensor`
  Acts::SpacePointContainer spacePoints;
  /// Space points of the strip layers, which carry
  /// `Acts::SpacePointColumns::StripCalibrationDetails` as well
  Acts::SpacePointContainer stripSpacePoints;
  /// The generating particles, indexed by the `particleId` column of either
  std::vector<GeneratedParticle> particles;
};

/// Which of an event's space points a seeder is to run on.
enum class SpacePointSelection {
  /// The pixel layers only, which is what a seeder tuned on pixels takes
  Pixel,
  /// The strip layers only
  Strip,
  /// Both of them, which is what a combined pass takes
  Combined,
};

/// Gather a selection of an event's space points into the single container a
/// seeder takes.
///
/// Columns are the union of what the selection holds, so a combined container
/// carries a stereo pair on every point and the entry means nothing on a pixel
/// one -- which is the price of putting both in one place, and why the event
/// itself does not. `copiedFromIndex` keeps pointing into the collection the
/// point came from, so a caller holding the event can find what it started as.
///
/// @param event the generated event
/// @param selection which of its collections to gather
/// @return the container
Acts::SpacePointContainer selectSpacePoints(const Event& event,
                                            SpacePointSelection selection);

}  // namespace ActsFatras::Synthetic

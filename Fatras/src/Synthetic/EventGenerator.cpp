// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/EventGenerator.hpp"

#include "ActsFatras/Synthetic/detail/StripReadout.hpp"

#include <cmath>
#include <cstdint>
#include <numbers>
#include <string>
#include <vector>

namespace ActsFatras::Synthetic {

namespace {

/// Create the dynamic columns the generator fills, unless they are already
/// there: `SpacePointContainer::clear` keeps them, and `createColumn` throws.
struct Columns {
  Acts::MutableSpacePointColumnProxy<std::uint32_t> layer;
  Acts::MutableSpacePointColumnProxy<float> clusterWidth;
  Acts::MutableSpacePointColumnProxy<float> localPositionY;
  Acts::MutableSpacePointColumnProxy<std::uint32_t> particle;
};

template <typename T>
Acts::MutableSpacePointColumnProxy<T> ensureColumn(
    Acts::SpacePointContainer& spacePoints, const std::string& name) {
  if (spacePoints.hasColumn(name)) {
    return spacePoints.column<T>(name);
  }
  return spacePoints.createColumn<T>(name);
}

Columns ensureColumns(Acts::SpacePointContainer& spacePoints) {
  return Columns{ensureColumn<std::uint32_t>(spacePoints, "layerId"),
                 ensureColumn<float>(spacePoints, "clusterWidth"),
                 ensureColumn<float>(spacePoints, "localPositionY"),
                 ensureColumn<std::uint32_t>(spacePoints, "particleId")};
}

/// What both collections carry; the strip one adds the pair on top.
constexpr Acts::SpacePointColumns kStandardColumns =
    Acts::SpacePointColumns::CopiedFromIndex | Acts::SpacePointColumns::X |
    Acts::SpacePointColumns::Y | Acts::SpacePointColumns::Z |
    Acts::SpacePointColumns::R | Acts::SpacePointColumns::Phi |
    Acts::SpacePointColumns::VarianceZ | Acts::SpacePointColumns::VarianceR;

/// Append one hit, everything but the variances it carries and whatever the
/// collection holds beyond them.
/// @param spacePoints the collection to append to
/// @param columns its dynamic columns
/// @param hit the hit
/// @return the appended space point
Acts::MutableSpacePointProxy fill(Acts::SpacePointContainer& spacePoints,
                                  Columns& columns,
                                  const detail::SmearedHit& hit) {
  // the helix azimuth is not wrapped while it is propagated
  const float phi = std::remainder(hit.phi, 2.f * std::numbers::pi_v<float>);
  auto sp = spacePoints.createSpacePoint();
  sp.x() = hit.r * std::cos(phi);
  sp.y() = hit.r * std::sin(phi);
  sp.z() = hit.z;
  sp.r() = hit.r;
  sp.phi() = phi;
  sp.copiedFromIndex() = sp.index();
  sp.extra(columns.layer) = hit.layer;
  // GBTS reads both; nothing here resolves a cluster, so they stay at zero
  sp.extra(columns.clusterWidth) = 0.f;
  sp.extra(columns.localPositionY) = 0.f;
  sp.extra(columns.particle) = hit.particle;
  return sp;
}

}  // namespace

EventGenerator::EventGenerator(const DetectorLayout& layout,
                               const EventConfig& config)
    : m_layout(&layout),
      m_cfg(config),
      m_primaries(config),
      m_propagator(layout, config),
      m_rng(config.seed) {}

void EventGenerator::reset(const std::uint32_t seed) {
  m_rng.seed(seed);
}

Event EventGenerator::generate() {
  Event event;
  generate(event);
  return event;
}

void EventGenerator::generate(Event& event) {
  const DetectorLayout& layout = *m_layout;
  const MeasurementConfig& cfg = m_cfg.simulation.measurement;

  m_scratch.tracks.clear();
  event.particles.clear();
  m_primaries.generate(m_rng, m_scratch.tracks, event.particles);
  m_propagator.propagate(m_rng, m_scratch, event.particles);

  Acts::SpacePointContainer& spacePoints = event.spacePoints;
  spacePoints.clear();
  spacePoints.createColumns(kStandardColumns);
  // not const: `SpacePointProxy::extra` takes a mutable column proxy by
  // non-const reference
  Columns columns = ensureColumns(spacePoints);
  spacePoints.reserve(static_cast<std::uint32_t>(m_scratch.hits.size()));

  const float measuredVariance = cfg.positionSmearing * cfg.positionSmearing;
  for (const detail::SmearedHit& hit : m_scratch.hits) {
    const bool cylinder =
        layout.layers[hit.layer].shape == SurfaceShape::Cylinder;
    auto sp = fill(spacePoints, columns, hit);
    // the hit sits at the nominal position along the normal, so that direction
    // carries no error: a cylinder measures z, a disc r
    sp.varianceZ() = cylinder ? measuredVariance : 0.f;
    sp.varianceR() = cylinder ? 0.f : measuredVariance;
  }

  Acts::SpacePointContainer& stripSpacePoints = event.stripSpacePoints;
  stripSpacePoints.clear();
  stripSpacePoints.createColumns(
      kStandardColumns | Acts::SpacePointColumns::StripCalibrationDetails);
  Columns stripColumns = ensureColumns(stripSpacePoints);
  stripSpacePoints.reserve(
      static_cast<std::uint32_t>(m_scratch.stripHits.size()));

  const std::vector<detail::StripLayer> strips = detail::stripLayers(layout);
  for (const detail::StripHit& hit : m_scratch.stripHits) {
    const detail::StripLayer& strip = strips[hit.hit.layer];
    const bool cylinder =
        layout.layers[hit.hit.layer].shape == SurfaceShape::Cylinder;
    auto sp = fill(stripSpacePoints, stripColumns, hit.hit);
    // The projection error lands along the strip, which runs along z on a
    // barrel and radially on an endcap; across it two strips measure one
    // coordinate between them and do better than either.
    sp.varianceZ() = cylinder ? strip.varianceAlong : 0.f;
    sp.varianceR() = cylinder ? 0.f : strip.varianceAlong;
    sp.outerStripCalibrationDetails() = hit.strips;
  }
}

}  // namespace ActsFatras::Synthetic

namespace ActsFatras {

Synthetic::Event Synthetic::generateEvent(const DetectorLayout& layout,
                                          const EventConfig& config) {
  return EventGenerator(layout, config).generate();
}

}  // namespace ActsFatras

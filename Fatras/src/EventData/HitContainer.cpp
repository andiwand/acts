// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/EventData/HitContainer.hpp"

namespace ActsFatras {

void HitContainer::reserve(std::uint32_t size) noexcept {
  m_hits.reserve(size);
}

void HitContainer::clear() noexcept {
  m_hits.clear();
  m_hitsByParticle.clear();
  m_hitsBySurface.clear();
}

Hit& HitContainer::emplace_back(const Acts::GeometryIdentifier geometryId,
                                const ParticleIndex particleId,
                                const Acts::Vector4& pos4,
                                const Acts::Vector4& before4,
                                const Acts::Vector4& after4,
                                const std::int32_t index) {
  m_hits.emplace_back(geometryId, particleId, pos4, before4, after4, index);
  const auto hitIndex = static_cast<HitIndex>(m_hits.size() - 1);
  m_hitsByParticle[particleId].push_back(hitIndex);
  m_hitsBySurface[geometryId].push_back(hitIndex);
  return m_hits.back();
}

std::span<const HitIndex> HitContainer::hitIndicesByParticle(
    ParticleIndex particleId) const {
  auto it = m_hitsByParticle.find(particleId);
  if (it != m_hitsByParticle.end()) {
    return it->second;
  }
  return {};
}

std::span<const HitIndex> HitContainer::hitIndicesBySurface(
    Acts::GeometryIdentifier geometryId) const {
  auto it = m_hitsBySurface.find(geometryId);
  if (it != m_hitsBySurface.end()) {
    return it->second;
  }
  return {};
}

}  // namespace ActsFatras

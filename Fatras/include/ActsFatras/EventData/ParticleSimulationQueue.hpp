// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/EventData/ParticleContainer.hpp"

#include <cstddef>

namespace ActsFatras {

// TODO split into frontend backend?
class ParticleSimulationQueue {
 public:
  explicit ParticleSimulationQueue(ParticleContainer &particles)
      : m_particles(&particles) {}

  std::size_t size() { return m_committed.size(); }
  bool empty() { return m_committed.empty(); }

  MutableParticleProxy createParticle() const {
    return m_particles->createParticle();
  }

  void selectParticle(MutableParticleProxy particle) const {
    m_selected.push_back(particle.index());
  }

  ConstParticleSubset selectedParticles() const {
    return m_particles->subset(m_selected).asConst();
  }

  void commitParticle(MutableParticleProxy particle) {
    m_committed.push_back(particle.index());
  }
  void commitParticle(ConstParticleProxy particle) {
    m_committed.push_back(particle.index());
  }

  MutableParticleProxy popParticle() {
    if (m_committed.empty()) {
      throw std::runtime_error(
          "ParticleSimulationQueue: no more particles to pop");
    }
    const auto index = m_committed.back();
    m_committed.pop_back();
    return m_particles->at(index);
  }

 private:
  ParticleContainer *m_particles{};
  mutable std::vector<ParticleIndex> m_selected;
  std::vector<ParticleIndex> m_committed;
};

}  // namespace ActsFatras

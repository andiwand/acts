// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/detail/PrimaryGenerator.hpp"

#include "ActsFatras/Synthetic/detail/Helix.hpp"
#include "ActsFatras/Synthetic/detail/Sampling.hpp"

#include <cmath>
#include <numbers>

namespace ActsFatras::Synthetic::detail {

PrimaryGenerator::PrimaryGenerator(const EventConfig& config)
    : m_cfg(config), m_mass(config.particleMass()) {}

void PrimaryGenerator::generate(
    std::mt19937& rng, std::vector<Helix>& tracks,
    std::vector<GeneratedParticle>& particles) const {
  constexpr float pi = std::numbers::pi_v<float>;

  const GenerationConfig& cfg = m_cfg.generation;
  const std::size_t numPrimaries = cfg.numPrimaries();
  tracks.reserve(tracks.size() + numPrimaries);
  particles.reserve(particles.size() + numPrimaries);

  for (std::size_t t = 0; t < numPrimaries; ++t) {
    const float pt =
        samplePt(sampleUniform01(rng), cfg.minPt, cfg.ptScale, cfg.ptExponent);
    const float rapidity =
        sampleRapidity(rng, cfg.minRapidity, cfg.maxRapidity, cfg.rapidityEdge,
                       cfg.rapidityEdgeWidth);
    // cot(theta) = pz / pT = (mT / pT) sinh(y); the mass is the whole of the
    // difference between a rapidity and a pseudorapidity
    const float cotTheta = std::hypot(1.f, m_mass / pt) * std::sinh(rapidity);
    const float d0 = cfg.d0Sigma * sampleNormal(rng);
    const float phi = sampleUniform(rng, -pi, pi);
    const float z0 = cfg.beamspotSigmaZ * sampleNormal(rng);
    // as many of one sign as of the other, there being no charge asymmetry
    // worth describing in a minimum-bias event
    const float charge = sampleCharge(rng);
    const Helix helix =
        makeHelix(d0, phi, z0, cotTheta, pt, charge, m_cfg.bFieldZ);
    tracks.push_back(helix);

    // out of the luminous region, so generation zero and no production vertex
    // to speak of
    particles.push_back(makeParticle(helix, pt));
  }

  generateDisplaced(rng, tracks, particles);
}

void PrimaryGenerator::generateDisplaced(
    std::mt19937& rng, std::vector<Helix>& tracks,
    std::vector<GeneratedParticle>& particles) const {
  constexpr float pi = std::numbers::pi_v<float>;

  const DisplacedConfig& cfg = m_cfg.displaced;
  const float logMinPt = std::log(cfg.minPt);
  const float logMaxPt = std::log(cfg.maxPt);

  for (std::size_t v = 0; v < cfg.vertices; ++v) {
    // the parent is neutral, so it flies in a straight line out of the beam
    // spot and its decay radius alone places the vertex
    const float parentPhi = sampleUniform(rng, -pi, pi);
    const float parentEta = sampleUniform(rng, -cfg.maxAbsEta, cfg.maxAbsEta);
    const float radius = sampleUniform(rng, cfg.minRadius, cfg.maxRadius);
    const float z = m_cfg.generation.beamspotSigmaZ * sampleNormal(rng) +
                    radius * std::sinh(parentEta);

    for (std::size_t t = 0; t < cfg.tracksPerVertex; ++t) {
      // A daughter pointing back towards the beam line would be an incoming
      // track, which nothing downstream expects of a decay. Only an opening
      // angle above pi/2 can make one, so the loop is bounded in practice.
      float phi = parentPhi;
      for (int attempt = 0; attempt < 100; ++attempt) {
        phi = parentPhi +
              sampleUniform(rng, -cfg.maxOpeningAngle, cfg.maxOpeningAngle);
        if (std::cos(phi - parentPhi) > 0.f) {
          break;
        }
      }
      const float eta = parentEta + sampleUniform(rng, -cfg.maxOpeningAngle,
                                                  cfg.maxOpeningAngle);
      const float pt = std::exp(sampleUniform(rng, logMinPt, logMaxPt));
      // a vertex is neutral, so its daughters alternate in charge; the odd
      // one out of an odd number is drawn
      const bool oddOneOut =
          cfg.tracksPerVertex % 2 == 1 && t + 1 == cfg.tracksPerVertex;
      const float charge =
          oddOneOut ? sampleCharge(rng) : (t % 2 == 0 ? 1.f : -1.f);
      const Helix helix =
          makeHelixFromPoint(std::max(radius, 1e-3f), parentPhi, z, phi,
                             std::sinh(eta), pt, charge, m_cfg.bFieldZ);
      tracks.push_back(helix);

      // a decay product, so one generation down: it is not counted as a
      // primary, and it makes no secondaries of its own at the default
      // `SecondaryConfig::maxGenerations`
      GeneratedParticle particle = makeParticle(helix, pt);
      particle.productionRadius = radius;
      particle.productionZ = z;
      particle.generation = 1;
      particle.displaced = true;
      particles.push_back(particle);
    }
  }
}

}  // namespace ActsFatras::Synthetic::detail

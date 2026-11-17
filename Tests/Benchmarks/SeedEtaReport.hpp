// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SeedContainer.hpp"
#include "ActsFatras/Synthetic/SyntheticEvent.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <numbers>
#include <random>
#include <set>
#include <string>
#include <vector>

/// Seeding efficiency and seed sizes per |eta| bin of the matched primary.
/// A seed belongs to the particle owning at least 3 of its space points.
inline void printSeedEtaReport(const ActsFatras::Synthetic::Event& event,
                               const Acts::SeedContainer& seeds,
                               const float ptThreshold, const std::string& tag) {
  const std::vector<float> edges = {0.0, 0.3, 0.6, 1.0, 1.5,
                                    2.0, 2.5, 3.0, 4.0};
  const std::size_t nBins = edges.size() - 1;
  const auto binOf = [&](float eta) -> int {
    const float a = std::abs(eta);
    for (std::size_t i = 0; i < nBins; ++i) {
      if (a < edges[i + 1]) {
        return static_cast<int>(i);
      }
    }
    return -1;
  };

  const auto particleColumn =
      event.spacePoints.column<std::uint32_t>("particleId");
  const auto layerColumn = event.spacePoints.column<std::uint32_t>("layerId");

  const std::size_t nP = event.particles.size();
  // distinct logical layers per particle
  std::vector<std::set<std::uint32_t>> layers(nP);
  for (std::size_t i = 0; i < event.spacePoints.size(); ++i) {
    const auto sp = event.spacePoints[i];
    const std::uint32_t p = sp.extra(particleColumn);
    if (p < nP) {
      layers[p].insert(sp.extra(layerColumn));
    }
  }

  // per particle: best (largest) true seed size and number of true seeds
  std::vector<std::uint32_t> best(nP, 0);
  std::vector<std::uint32_t> nSeeds(nP, 0);
  std::array<std::array<std::size_t, 8>, 16> sizeHist{};  // [bin][size]
  for (const auto seed : seeds) {
    const auto idx = seed.spacePointIndices();
    std::uint32_t owner = 0;
    std::size_t nBest = 0;
    for (const auto i : idx) {
      const std::uint32_t c = event.spacePoints[i].extra(particleColumn);
      const auto n = static_cast<std::size_t>(std::ranges::count_if(
          idx, [&](auto o) {
            return event.spacePoints[o].extra(particleColumn) == c;
          }));
      if (n > nBest) {
        nBest = n;
        owner = c;
      }
    }
    if (nBest < 3 || owner >= nP) {
      continue;
    }
    const auto& info = event.particles[owner];
    if (!info.primary() || info.pt < ptThreshold) {
      continue;
    }
    best[owner] = std::max<std::uint32_t>(best[owner],
                                          static_cast<std::uint32_t>(nBest));
    ++nSeeds[owner];
    if (const int b = binOf(info.eta); b >= 0) {
      ++sizeHist[b][std::min<std::size_t>(idx.size(), 7)];
    }
  }

  std::printf(
      "[%s] |eta| bin   seedable  eff(>=3hits)  seedable4L  eff(>=4layers)  "
      "only3SP  dup(>1 seed)  true seeds by size 3/4/5/6+\n",
      tag.c_str());
  for (std::size_t b = 0; b < nBins; ++b) {
    std::size_t s3 = 0, m3 = 0, s4 = 0, m4 = 0, only3 = 0, dup = 0;
    for (std::size_t p = 0; p < nP; ++p) {
      const auto& info = event.particles[p];
      if (!info.primary() || info.pt < ptThreshold || info.numHits < 3 ||
          binOf(info.eta) != static_cast<int>(b)) {
        continue;
      }
      ++s3;
      const bool matched = best[p] >= 3;
      m3 += matched;
      if (layers[p].size() >= 4) {
        ++s4;
        m4 += matched;
      }
      only3 += (best[p] == 3);
      dup += (nSeeds[p] > 1);
    }
    const auto& h = sizeHist[b];
    std::printf(
        "[%s] %.1f-%.1f  %8zu  %11.4f  %10zu  %14.4f  %7.4f  %12.4f  "
        "%zu/%zu/%zu/%zu\n",
        tag.c_str(), edges[b], edges[b + 1], s3,
        s3 ? double(m3) / s3 : 0., s4, s4 ? double(m4) / s4 : 0.,
        m3 ? double(only3) / m3 : 0., m3 ? double(dup) / m3 : 0., h[3], h[4],
        h[5], h[6] + h[7]);
  }
}

/// Emulate hit inefficiency: move a fraction `GBTS_HIT_INEFF` of the space
/// points to a random azimuth on their own layer and hand them to a dummy
/// secondary, so they are noise and no longer on their particle.
inline void applyHitInefficiency(ActsFatras::Synthetic::Event& event,
                                 std::uint32_t seed) {
  const char* env = std::getenv("GBTS_HIT_INEFF");
  const float frac = env != nullptr ? std::strtof(env, nullptr) : 0.f;
  if (frac <= 0.f) {
    return;
  }
  ActsFatras::Synthetic::GeneratedParticle dummy;
  dummy.generation = 1;
  dummy.productionRadius = 1.f;
  const auto dummyId = static_cast<std::uint32_t>(event.particles.size());
  event.particles.push_back(dummy);
  auto particleColumn = event.spacePoints.column<std::uint32_t>("particleId");
  std::mt19937 rng(seed * 7919u + 13u);
  std::uniform_real_distribution<float> u(0.f, 1.f);
  std::size_t moved = 0;
  for (std::size_t i = 0; i < event.spacePoints.size(); ++i) {
    if (u(rng) >= frac) {
      continue;
    }
    auto sp = event.spacePoints.at(static_cast<Acts::SpacePointIndex>(i));
    const float phi = std::numbers::pi_v<float> * (2.f * u(rng) - 1.f);
    const float r = sp.r();
    sp.x() = r * std::cos(phi);
    sp.y() = r * std::sin(phi);
    sp.phi() = phi;
    sp.extra(particleColumn) = dummyId;
    ++moved;
  }
  std::printf("hit inefficiency %.3f: moved %zu space points\n", frac, moved);
}

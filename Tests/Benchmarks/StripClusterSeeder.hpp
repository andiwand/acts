// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Exploration: a seeder for single strip clusters, transverse cuts first.
///
/// A seeding point is one cluster of the inner sensor. It knows its module
/// row exactly, its position along the strip only to the strip length, and it
/// carries a list of candidate partners on the outer sensor, each of which
/// resolves it to a precise point once a direction is given.
///
///  1. Modules (surface, strip row) have exact r and z intervals, so which
///     modules can form a doublet is decided once, from geometry, by interval
///     arithmetic: the connection table.
///  2. A doublet runs through the connected modules, the azimuth window and the
///     transverse cuts first. Only then are both clusters resolved against the
///     chord with their partners, and the r-z cuts run on resolved points.
///  3. A triplet requires both doublets to resolve the middle with the same
///     partner, passes the circle cuts, and is resolved again with the circle
///     tangent at each point for a precise r-z compatibility cut.

#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/SpacePointFormation/detail/StripSpacePointCalibrationImpl.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <unordered_map>
#include <vector>

#include "StripClusterSeeding.hpp"

namespace ActsTests::StripClusterSeeding {

struct SeederConfig {
  /// Minimum transverse momentum, as a helix radius in mm
  float minHelixRadius = 1667.f;
  float impactMax = 300.f;
  float zOriginMax = 500.f;
  float cotThetaMax = 5.f;
  float deltaRMin = 50.f;
  float deltaRMax = 250.f;
  float deltaZMax = 850.f;
  /// Athena `doubletDPhiCut`
  float dPhiConst = 0.015f;
  float dPhiSlope = 2e-4f;
  float dPhiCap = 0.10f;
  /// Resolution tolerance in strip half lengths, doublet and triplet stage
  float doubletTolerance = 1.75f;
  float tripletTolerance = 1.1f;
  /// cot(theta) difference of the two calibrated doublets, pre-filter
  float cotThetaDiffMax = 0.1f;
  /// cot(theta) difference of the two halves of the triplet, resolved with
  /// the circle tangent
  float tripletCotThetaDiffMax = 0.03f;
  std::size_t maxSeedsPerMiddle = 2;
  /// Margin on the module intervals in mm
  float margin = 2.f;
};

struct SeederStats {
  std::size_t middles{};
  std::size_t candidates{};
  std::size_t afterDPhi{};
  std::size_t afterImpact{};
  std::size_t doublets{};
  std::size_t tripletCandidates{};
  std::size_t afterPartner{};
  std::size_t afterCircle{};
  std::size_t seeds{};
};

class StripClusterSeeder {
 public:
  explicit StripClusterSeeder(SeederConfig cfg) : m_cfg(cfg) {}

  /// Register the modules the clusters of an event sit in. Geometry: called
  /// for every event up front, outside the timed region.
  void addModules(const DetectorLayout& layout,
                  const std::vector<StripLayer>& strips,
                  const std::vector<StripCluster>& clusters) {
    for (const StripCluster& c : clusters) {
      if (c.side != 0) {
        continue;
      }
      const std::uint64_t key = moduleKey(c);
      if (m_moduleOf.contains(key)) {
        continue;
      }
      const StripLayer& layer = strips[c.layer];
      const bool barrel =
          layout.layers[c.layer].shape == SurfaceShape::Cylinder;
      const float along = layer.halfLength * m_cfg.doubletTolerance;
      Module m;
      m.barrel = barrel;
      const float r = clusterR(c);
      if (barrel) {
        m.rMin = r - m_cfg.margin;
        m.rMax = r + m_cfg.margin;
        m.zMin = c.center[2] - along - m_cfg.margin;
        m.zMax = c.center[2] + along + m_cfg.margin;
      } else {
        m.rMin = r - along - m_cfg.margin;
        m.rMax = r + along + m_cfg.margin;
        m.zMin = c.center[2] - m_cfg.margin;
        m.zMax = c.center[2] + m_cfg.margin;
      }
      m_moduleOf.emplace(key, static_cast<std::uint32_t>(m_modules.size()));
      m_modules.push_back(m);
    }
  }

  /// Build the connection table from the registered modules.
  void connect() {
    const std::size_t n = m_modules.size();
    m_bottoms.assign(n, {});
    m_tops.assign(n, {});
    m_numConnections = 0;
    for (std::size_t mi = 0; mi < n; ++mi) {
      const Module& m = m_modules[mi];
      for (std::size_t oi = 0; oi < n; ++oi) {
        if (oi == mi) {
          continue;
        }
        const Module& o = m_modules[oi];
        // the other one as the bottom, then as the top
        for (const bool asTop : {false, true}) {
          const Module& inner = asTop ? m : o;
          const Module& outer = asTop ? o : m;
          float window = 0.f;
          if (!compatible(inner, outer, window)) {
            continue;
          }
          (asTop ? m_tops : m_bottoms)[mi].push_back(
              {static_cast<std::uint32_t>(oi), window});
          ++m_numConnections;
        }
      }
    }
  }

  std::size_t numModules() const { return m_modules.size(); }
  std::size_t numConnections() const { return m_numConnections; }

  /// Seed one event.
  /// @param points the seeding points, inner sensor clusters
  /// @param clusters the clusters of the event
  /// @param clusterOf the cluster of each point
  /// @param partners the partner table, indexed like the points
  /// @param seeds receives the seeds, on `points`
  void seed(const Acts::SpacePointContainer& points,
            const std::vector<StripCluster>& clusters,
            const std::vector<std::uint32_t>& clusterOf,
            const PartnerTable& partners, Acts::SeedContainer& seeds,
            SeederStats& stats) {
    sortPoints(points, clusters, clusterOf, partners);
    seeds.assignSpacePointContainer(points);

    const float pi = std::numbers::pi_v<float>;
    for (std::uint32_t mod = 0; mod < m_modules.size(); ++mod) {
      for (std::uint32_t im = m_modBegin[mod]; im < m_modBegin[mod + 1]; ++im) {
        const float zM = m_z[im];
        const float rM = m_r[im];
        const float absZ = std::abs(zM);
        // Athena LRT rRangeMiddleSP
        const float rLo = absZ < 1400.f ? 500.f : 400.f;
        const float rHi = absZ < 1400.f ? 800.f : 850.f;
        if (absZ > 2500.f || rM < rLo || rM > rHi) {
          continue;
        }
        ++stats.middles;
        m_bottom.clear();
        m_top.clear();
        findDoublets(im, m_bottoms[mod], true, partners, stats);
        findDoublets(im, m_tops[mod], false, partners, stats);
        if (m_bottom.empty() || m_top.empty()) {
          continue;
        }
        std::ranges::sort(m_top, {}, &Doublet::cot);
        triplets(im, partners, stats);
        (void)pi;
      }
    }
    // write out
    for (const Seed& s : m_seeds) {
      auto seed = seeds.createSeed();
      seed.assignSpacePointIndices(std::array<Acts::SpacePointIndex, 3>{
          m_original[s.b], m_original[s.m], m_original[s.t]});
      seed.quality() = -s.score;
      seed.vertexZ() = s.zOrigin;
    }
    stats.seeds += m_seeds.size();
    m_seeds.clear();
  }

 private:
  struct Module {
    bool barrel{};
    float rMin{}, rMax{}, zMin{}, zMax{};
  };
  struct Connection {
    std::uint32_t module;
    float window;
  };
  struct Doublet {
    std::uint32_t other;
    std::uint32_t partnerM;
    std::uint32_t partnerO;
    float cot;
    float zOrigin;
    std::array<float, 3> resolvedO;
    std::array<float, 3> resolvedM;
  };
  struct Seed {
    std::uint32_t b, m, t;
    float score;
    float zOrigin;
  };

  /// Interval check whether an inner and an outer module can hold a doublet.
  bool compatible(const Module& in, const Module& out, float& window) const {
    const float dRMin = out.rMin - in.rMax;
    const float dRMax = out.rMax - in.rMin;
    if (dRMax < m_cfg.deltaRMin || dRMin > m_cfg.deltaRMax) {
      return false;
    }
    const float dZMin = out.zMin - in.zMax;
    const float dZMax = out.zMax - in.zMin;
    if (dZMin > m_cfg.deltaZMax || dZMax < -m_cfg.deltaZMax) {
      return false;
    }
    // smallest |dz| over the largest usable dr
    const float absDZ = (dZMin <= 0.f && dZMax >= 0.f)
                            ? 0.f
                            : std::min(std::abs(dZMin), std::abs(dZMax));
    const float usableDR = std::min(dRMax, m_cfg.deltaRMax);
    if (absDZ > m_cfg.cotThetaMax * usableDR) {
      return false;
    }
    // z origin over the corners of both boxes, with dr kept positive
    float z0Min = std::numeric_limits<float>::max();
    float z0Max = std::numeric_limits<float>::lowest();
    for (float rI : {in.rMin, in.rMax}) {
      for (float rO : {out.rMin, out.rMax}) {
        const float dR = rO - rI;
        if (dR < 1.f) {
          continue;
        }
        for (float zI : {in.zMin, in.zMax}) {
          for (float zO : {out.zMin, out.zMax}) {
            const float z0 = zI - rI * (zO - zI) / dR;
            z0Min = std::min(z0Min, z0);
            z0Max = std::max(z0Max, z0);
          }
        }
      }
    }
    if (z0Min > m_cfg.zOriginMax || z0Max < -m_cfg.zOriginMax) {
      return false;
    }
    // azimuth reach, as the Athena doublet cut at the widest corners
    const auto swing = [&](float r) {
      return std::asin(std::min(1.f, m_cfg.impactMax / std::max(r, 1.f)));
    };
    const float swingMax = std::max(std::abs(swing(in.rMin) - swing(out.rMax)),
                                    std::abs(swing(in.rMax) - swing(out.rMin)));
    window = m_cfg.dPhiConst + m_cfg.dPhiSlope * std::min(dRMax, 1000.f) +
             std::min(m_cfg.dPhiCap, swingMax);
    return true;
  }

  void sortPoints(const Acts::SpacePointContainer& points,
                  const std::vector<StripCluster>& clusters,
                  const std::vector<std::uint32_t>& clusterOf,
                  const PartnerTable& partners) {
    const std::size_t n = points.size();
    m_order.resize(n);
    m_moduleOfPoint.resize(n);
    for (std::uint32_t i = 0; i < n; ++i) {
      m_order[i] = i;
      const auto it = m_moduleOf.find(moduleKey(clusters[clusterOf[i]]));
      m_moduleOfPoint[i] = it == m_moduleOf.end()
                               ? std::numeric_limits<std::uint32_t>::max()
                               : it->second;
    }
    std::ranges::sort(m_order, [&](std::uint32_t a, std::uint32_t b) {
      return m_moduleOfPoint[a] != m_moduleOfPoint[b]
                 ? m_moduleOfPoint[a] < m_moduleOfPoint[b]
                 : points[a].phi() < points[b].phi();
    });
    m_x.resize(n);
    m_y.resize(n);
    m_z.resize(n);
    m_r.resize(n);
    m_phi.resize(n);
    m_asin.resize(n);
    m_original.resize(n);
    m_modBegin.assign(m_modules.size() + 1, 0);
    for (std::uint32_t k = 0; k < n; ++k) {
      const std::uint32_t i = m_order[k];
      const auto sp = points[i];
      m_x[k] = sp.x();
      m_y[k] = sp.y();
      m_z[k] = sp.z();
      m_r[k] = sp.r();
      m_phi[k] = sp.phi();
      m_asin[k] =
          std::asin(std::min(1.f, m_cfg.impactMax / std::max(sp.r(), 1.f)));
      m_original[k] = i;
      const std::uint32_t mod = m_moduleOfPoint[i];
      if (mod != std::numeric_limits<std::uint32_t>::max()) {
        ++m_modBegin[mod + 1];
      }
    }
    for (std::size_t mod = 0; mod < m_modules.size(); ++mod) {
      m_modBegin[mod + 1] += m_modBegin[mod];
    }
    // the partners in the same order, so that a module's points and their
    // partners sit next to each other in memory
    m_partnerBegin.resize(n + 1);
    m_partners.clear();
    m_partners.reserve(partners.size());
    for (std::uint32_t k = 0; k < n; ++k) {
      const std::uint32_t i = m_order[k];
      m_partnerBegin[k] = static_cast<std::uint32_t>(m_partners.size());
      for (std::uint32_t j = partners.begin(i); j < partners.end(i); ++j) {
        m_partners.push_back(partners.details(j));
      }
    }
    m_partnerBegin[n] = static_cast<std::uint32_t>(m_partners.size());
  }

  void findDoublets(std::uint32_t im, const std::vector<Connection>& conns,
                    bool bottom, const PartnerTable& partners,
                    SeederStats& stats) {
    const float pi = std::numbers::pi_v<float>;
    const float rM = m_r[im];
    const float phiM = m_phi[im];
    const float zM = m_z[im];
    const std::uint32_t pM = m_original[im];

    for (const Connection& conn : conns) {
      const std::uint32_t begin = m_modBegin[conn.module];
      const std::uint32_t end = m_modBegin[conn.module + 1];
      if (begin == end) {
        continue;
      }
      const auto scan = [&](float lo, float hi) {
        auto first =
            std::lower_bound(m_phi.begin() + begin, m_phi.begin() + end, lo);
        for (auto it = first; it != m_phi.begin() + end && *it <= hi; ++it) {
          const auto io = static_cast<std::uint32_t>(it - m_phi.begin());
          ++stats.candidates;
          const float rO = m_r[io];
          const float dR = bottom ? rM - rO : rO - rM;
          // exact Athena azimuth cut
          float dPhi = m_phi[io] - phiM;
          dPhi = dPhi > pi ? dPhi - 2.f * pi
                           : (dPhi < -pi ? dPhi + 2.f * pi : dPhi);
          const float bound =
              m_cfg.dPhiConst + m_cfg.dPhiSlope * std::abs(dR) +
              std::min(m_cfg.dPhiCap, std::abs(m_asin[io] - m_asin[im]));
          if (std::abs(dPhi) > bound) {
            continue;
          }
          ++stats.afterDPhi;
          // a loose transverse impact cut (line plus sagitta) rejected 0.2%
          // here: dropped, the circle cut of the triplet does the work
          ++stats.afterImpact;
          resolveDoublet(im, io, bottom, partners);
        }
      };
      const float lo = phiM - conn.window;
      const float hi = phiM + conn.window;
      scan(lo, hi);
      if (lo < -pi) {
        scan(lo + 2.f * pi, pi);
      }
      if (hi > pi) {
        scan(-pi, hi - 2.f * pi);
      }
    }
    (void)zM;
    (void)pM;
  }

  void resolveDoublet(std::uint32_t im, std::uint32_t io, bool bottom,
                      const PartnerTable& partners) {
    // chord, pointing outwards
    std::array<float, 3> chord{m_x[io] - m_x[im], m_y[io] - m_y[im],
                               m_z[io] - m_z[im]};
    if (bottom) {
      chord = {-chord[0], -chord[1], -chord[2]};
    }
    // the middle first: for a given middle most chords fail to resolve it,
    // and then the other point is never touched
    constexpr std::size_t kMax = 6;
    std::array<std::array<float, 3>, kMax> resolvedM{};
    std::array<float, kMax> radiusM{};
    std::array<std::uint32_t, kMax> partnerM{};
    std::size_t numM = 0;
    for (std::uint32_t kM = m_partnerBegin[im];
         kM < m_partnerBegin[im + 1] && numM < kMax; ++kM) {
      if (Acts::detail::calibrateOuterStripSpacePoint(
              chord, m_partners[kM], resolvedM[numM], m_cfg.doubletTolerance)) {
        radiusM[numM] = std::sqrt(resolvedM[numM][0] * resolvedM[numM][0] +
                                  resolvedM[numM][1] * resolvedM[numM][1]);
        partnerM[numM] = kM;
        ++numM;
      }
    }
    if (numM == 0) {
      return;
    }
    ++m_middleResolved;
    for (std::uint32_t kO = m_partnerBegin[io]; kO < m_partnerBegin[io + 1];
         ++kO) {
      std::array<float, 3> ro{};
      if (!Acts::detail::calibrateOuterStripSpacePoint(
              chord, m_partners[kO], ro, m_cfg.doubletTolerance)) {
        continue;
      }
      const float rOr = std::sqrt(ro[0] * ro[0] + ro[1] * ro[1]);
      for (std::size_t i = 0; i < numM; ++i) {
        const std::array<float, 3>& rm = resolvedM[i];
        const float rMr = radiusM[i];
        const std::uint32_t kM = partnerM[i];
        const float dR = bottom ? rMr - rOr : rOr - rMr;
        if (dR < m_cfg.deltaRMin || dR > m_cfg.deltaRMax) {
          continue;
        }
        const float dZ = bottom ? rm[2] - ro[2] : ro[2] - rm[2];
        if (std::abs(dZ) > m_cfg.deltaZMax ||
            std::abs(dZ) > m_cfg.cotThetaMax * dR) {
          continue;
        }
        const float cot = dZ / dR;
        const float zOrigin = rm[2] - rMr * cot;
        if (std::abs(zOrigin) > m_cfg.zOriginMax) {
          continue;
        }
        (bottom ? m_bottom : m_top)
            .push_back({io, kM, kO, cot, zOrigin, ro, rm});
      }
    }
  }

  void triplets(std::uint32_t im, const PartnerTable& partners,
                SeederStats& stats) {
    stats.doublets += m_bottom.size() + m_top.size();
    m_candidates.clear();
    const std::uint32_t pM = m_original[im];
    for (const Doublet& b : m_bottom) {
      const auto lo = std::ranges::lower_bound(
          m_top, b.cot - m_cfg.cotThetaDiffMax, {}, &Doublet::cot);
      for (auto it = lo;
           it != m_top.end() && it->cot <= b.cot + m_cfg.cotThetaDiffMax;
           ++it) {
        const Doublet& t = *it;
        ++stats.tripletCandidates;
        // both doublets have to resolve the middle with the same partner
        if (t.partnerM != b.partnerM) {
          continue;
        }
        ++stats.afterPartner;
        // circle through the resolved transverse positions
        const std::array<float, 2> pB{b.resolvedO[0], b.resolvedO[1]};
        const std::array<float, 2> pMid{
            0.5f * (b.resolvedM[0] + t.resolvedM[0]),
            0.5f * (b.resolvedM[1] + t.resolvedM[1])};
        const std::array<float, 2> pT{t.resolvedO[0], t.resolvedO[1]};
        const float ax = pMid[0] - pB[0];
        const float ay = pMid[1] - pB[1];
        const float bx = pT[0] - pB[0];
        const float by = pT[1] - pB[1];
        const float d = 2.f * (ax * by - ay * bx);
        if (std::abs(d) < 1e-6f) {
          continue;
        }
        const float a2 = ax * ax + ay * ay;
        const float b2 = bx * bx + by * by;
        const float cx = pB[0] + (by * a2 - ay * b2) / d;
        const float cy = pB[1] + (ax * b2 - bx * a2) / d;
        const float radius = std::hypot(pB[0] - cx, pB[1] - cy);
        if (radius < m_cfg.minHelixRadius) {
          continue;
        }
        const float impact = std::abs(std::hypot(cx, cy) - radius);
        if (impact > m_cfg.impactMax) {
          continue;
        }
        ++stats.afterCircle;
        // resolve all three again with the circle tangent
        const float cotAvg = 0.5f * (b.cot + t.cot);
        const auto tangent = [&](const std::array<float, 2>& p) {
          float ux = -(p[1] - cy);
          float uy = p[0] - cx;
          const float n = std::hypot(ux, uy);
          ux /= n;
          uy /= n;
          if (ux * bx + uy * by < 0.f) {
            ux = -ux;
            uy = -uy;
          }
          return std::array<float, 3>{ux, uy, cotAvg};
        };
        std::array<float, 3> qB{};
        std::array<float, 3> qM{};
        std::array<float, 3> qT{};
        if (!Acts::detail::calibrateOuterStripSpacePoint(
                tangent(pB), m_partners[b.partnerO], qB,
                m_cfg.tripletTolerance) ||
            !Acts::detail::calibrateOuterStripSpacePoint(
                tangent(pMid), m_partners[b.partnerM], qM,
                m_cfg.tripletTolerance) ||
            !Acts::detail::calibrateOuterStripSpacePoint(
                tangent(pT), m_partners[t.partnerO], qT,
                m_cfg.tripletTolerance)) {
          continue;
        }
        const float rB = std::hypot(qB[0], qB[1]);
        const float rMq = std::hypot(qM[0], qM[1]);
        const float rT = std::hypot(qT[0], qT[1]);
        const float cotB = (qM[2] - qB[2]) / (rMq - rB);
        const float cotT = (qT[2] - qM[2]) / (rT - rMq);
        const float diff = std::abs(cotB - cotT);
        if (diff > m_cfg.tripletCotThetaDiffMax) {
          continue;
        }
        const float zOrigin = qB[2] - rB * cotB;
        if (std::abs(zOrigin) > m_cfg.zOriginMax) {
          continue;
        }
        m_candidates.push_back({b.other, im, t.other, diff, zOrigin});
      }
    }
    // keep the best few per middle
    const std::size_t keep =
        std::min(m_cfg.maxSeedsPerMiddle, m_candidates.size());
    std::partial_sort(
        m_candidates.begin(), m_candidates.begin() + keep, m_candidates.end(),
        [](const Seed& a, const Seed& c) { return a.score < c.score; });
    m_seeds.insert(m_seeds.end(), m_candidates.begin(),
                   m_candidates.begin() + keep);
    (void)pM;
  }

  SeederConfig m_cfg;
  std::vector<Module> m_modules;
  std::unordered_map<std::uint64_t, std::uint32_t> m_moduleOf;
  std::vector<std::vector<Connection>> m_bottoms;
  std::vector<std::vector<Connection>> m_tops;
  std::size_t m_numConnections{};

  // per event, sorted by module and azimuth
  std::vector<std::uint32_t> m_order;
  std::vector<std::uint32_t> m_moduleOfPoint;
  std::vector<std::uint32_t> m_modBegin;
  std::vector<float> m_x, m_y, m_z, m_r, m_phi, m_asin;
  std::vector<std::uint32_t> m_original;
  std::vector<std::uint32_t> m_partnerBegin;
  std::vector<Acts::OuterStripSpacePointCalibrationDetailsDerived> m_partners;

  std::size_t m_middleResolved{};

 public:
  std::size_t middleResolved() const { return m_middleResolved; }
  void resetMiddleResolved() { m_middleResolved = 0; }

 private:
  std::vector<Doublet> m_bottom;
  std::vector<Doublet> m_top;
  std::vector<Seed> m_candidates;
  std::vector<Seed> m_seeds;
};

}  // namespace ActsTests::StripClusterSeeding

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Exploration: seed the strips from single clusters instead of space points.
///
/// Two inputs for the grid triplet seeder, both built from the strip clusters
/// of the synthetic event:
///
///  - `formSpacePoints` pairs the clusters of the two sensors of a module into
///    space points under the beam spot constraint, every compatible pair, so
///    that ghosts appear as they do in a real formation step. The gap
///    tolerance is the knob that recovers displaced crossings.
///  - `makeClusterSpacePoints` takes the clusters of one sensor only, as
///    points at the centre of their strip. The coordinate along the strip
///    (z on a barrel, r on an endcap) is only known to the strip length, and
///    the variance says so.
///
/// `refineSeeds` then turns a seed of single clusters into one of space points:
/// it finds the partner cluster on the other sensor of each module, resolves
/// the pair with the direction of the seed instead of the beam spot, and cuts
/// on the r-z consistency of the three resolved points.

#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/SpacePointFormation/StripSpacePointBuilder.hpp"
#include "Acts/SpacePointFormation/detail/StripSpacePointCalibrationImpl.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/SyntheticEvent.hpp"
#include "ActsFatras/Synthetic/detail/StripReadout.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <vector>

namespace ActsTests::StripClusterSeeding {

using ActsFatras::Synthetic::DetectorLayout;
using ActsFatras::Synthetic::StripCluster;
using ActsFatras::Synthetic::SurfaceShape;
using ActsFatras::Synthetic::detail::StripLayer;

/// @return the key that groups the clusters one module row can pair
inline std::uint64_t moduleKey(const StripCluster& c) {
  return (static_cast<std::uint64_t>(c.surface) << 32) |
         static_cast<std::uint32_t>(c.row);
}

/// @return the azimuth of a cluster centre
inline float clusterPhi(const StripCluster& c) {
  return std::atan2(c.center[1], c.center[0]);
}

/// @return the radius of a cluster centre
inline float clusterR(const StripCluster& c) {
  return std::hypot(c.center[0], c.center[1]);
}

/// The clusters of one sensor side, sorted by module row and azimuth, for
/// partner lookups in an azimuth window.
class ClusterIndex {
 public:
  /// @param clusters every cluster of the event
  /// @param side the sensor side to index
  void build(const std::vector<StripCluster>& clusters, std::uint8_t side) {
    m_entries.clear();
    m_entries.reserve(clusters.size() / 2);
    for (std::uint32_t i = 0; i < clusters.size(); ++i) {
      const StripCluster& c = clusters[i];
      if (c.side != side) {
        continue;
      }
      m_entries.push_back({moduleKey(c), clusterPhi(c), i});
    }
    std::ranges::sort(m_entries, [](const Entry& a, const Entry& b) {
      return a.key != b.key ? a.key < b.key : a.phi < b.phi;
    });
  }

  /// Call `f(index)` for every indexed cluster in the module row `key` with
  /// an azimuth within `halfWidth` of `phi`.
  template <typename F>
  void forEachInWindow(std::uint64_t key, float phi, float halfWidth,
                       F&& f) const {
    // a stereo strip near the end of its row meets the other sensor in the
    // neighbouring row
    for (std::int64_t dr = -1; dr <= 1; ++dr) {
      const std::uint64_t k =
          (key & 0xffffffff00000000ull) |
          static_cast<std::uint32_t>(static_cast<std::int32_t>(key) + dr);
      forEachInRow(k, phi, halfWidth, f);
    }
  }

  template <typename F>
  void forEachInRow(std::uint64_t key, float phi, float halfWidth, F& f) const {
    constexpr float pi = std::numbers::pi_v<float>;
    const auto scan = [&](float lo, float hi) {
      auto it = std::ranges::lower_bound(
          m_entries, Entry{key, lo, 0}, [](const Entry& a, const Entry& b) {
            return a.key != b.key ? a.key < b.key : a.phi < b.phi;
          });
      for (; it != m_entries.end() && it->key == key && it->phi <= hi; ++it) {
        f(it->index);
      }
    };
    scan(phi - halfWidth, phi + halfWidth);
    if (phi - halfWidth < -pi) {
      scan(phi - halfWidth + 2 * pi, pi);
    }
    if (phi + halfWidth > pi) {
      scan(-pi, phi + halfWidth - 2 * pi);
    }
  }

  std::size_t size() const { return m_entries.size(); }

 private:
  struct Entry {
    std::uint64_t key;
    float phi;
    std::uint32_t index;
  };
  std::vector<Entry> m_entries;
};

/// The columns every container built here carries.
constexpr Acts::SpacePointColumns kColumns =
    Acts::SpacePointColumns::CopiedFromIndex | Acts::SpacePointColumns::X |
    Acts::SpacePointColumns::Y | Acts::SpacePointColumns::Z |
    Acts::SpacePointColumns::R | Acts::SpacePointColumns::Phi |
    Acts::SpacePointColumns::VarianceZ | Acts::SpacePointColumns::VarianceR;

/// Fill the position columns of a space point.
inline void setPosition(Acts::MutableSpacePointProxy& sp, float x, float y,
                        float z) {
  sp.x() = x;
  sp.y() = y;
  sp.z() = z;
  sp.r() = std::hypot(x, y);
  sp.phi() = std::atan2(y, x);
  sp.copiedFromIndex() = sp.index();
}

/// What `formSpacePoints` did.
struct FormationStats {
  std::size_t pairsTried{};
  std::size_t formed{};
  std::size_t ghosts{};
  /// Pairs of the two clusters of one crossing
  std::size_t trueCrossings{};
};

/// Pair the clusters of both sensors into space points under the beam spot
/// constraint, every compatible pair.
///
/// @param layout the detector layout
/// @param strips the strip layer descriptors
/// @param clusters the clusters of the event
/// @param outerIndex the index of the outer sensor clusters
/// @param gapParameter half the curvature of the softest track a pair should
///        still resolve, as `EventConfig::stripGapParameter`; this is the
///        tolerance that is relaxed to recover displaced crossings
/// @param ghostParticle the particle id a ghost gets
/// @param out receives the space points, with the `particleId` column
/// @return statistics
inline FormationStats formSpacePoints(
    const DetectorLayout& layout, const std::vector<StripLayer>& strips,
    const std::vector<StripCluster>& clusters, const ClusterIndex& outerIndex,
    float gapParameter, std::uint32_t ghostParticle, bool dropGhosts,
    Acts::SpacePointContainer& out, bool bestPartnerOnly = false) {
  FormationStats stats;
  out.clear();
  if (!out.hasColumn("particleId")) {
    out.createColumns(kColumns |
                      Acts::SpacePointColumns::StripCalibrationDetails);
    out.createColumn<std::uint32_t>("particleId");
  }
  auto particleColumn = out.column<std::uint32_t>("particleId");

  Acts::StripSpacePointBuilder::ConstrainedOptions options;
  options.vertex = Acts::Vector3::Zero();
  const auto ends = [](const StripCluster& c) {
    Acts::StripSpacePointBuilder::StripEnds e;
    e.top = Acts::Vector3(c.center[0] + c.halfVector[0],
                          c.center[1] + c.halfVector[1],
                          c.center[2] + c.halfVector[2]);
    e.bottom = Acts::Vector3(c.center[0] - c.halfVector[0],
                             c.center[1] - c.halfVector[1],
                             c.center[2] - c.halfVector[2]);
    return e;
  };

  for (std::uint32_t i = 0; i < clusters.size(); ++i) {
    const StripCluster& inner = clusters[i];
    if (inner.side != 0) {
      continue;
    }
    const StripLayer& layer = strips[inner.layer];
    const bool cylinder =
        layout.layers[inner.layer].shape == SurfaceShape::Cylinder;
    const float r = clusterR(inner);
    const float slope = cylinder ? 1.f : r / std::abs(inner.center[2]);
    options.stripLengthGapTolerance =
        gapParameter * r * layer.gapOverStereo * slope;
    // how far along the strip the pair may resolve, in half lengths, and so
    // how far apart across the strips the two centres may be
    const float reach = 1.f + static_cast<float>(options.stripLengthTolerance) +
                        static_cast<float>(options.stripLengthGapTolerance) /
                            (2.f * layer.halfLength);
    const float sinStereo = 2.f * layer.sinHalfStereo * layer.cosHalfStereo;
    const float halfWidth =
        (layer.halfLength * reach * sinStereo + 2.f * layer.halfGap) / r;

    const Acts::StripSpacePointBuilder::ConstrainedStripCache innerCache =
        Acts::StripSpacePointBuilder::makeConstrainedStripCache(ends(inner),
                                                                options);
    // with one point per inner cluster, the pair resolving closest to the
    // middle of the inner strip wins
    float bestOffset = std::numeric_limits<float>::max();
    std::uint32_t best = 0;
    Acts::Vector3 bestPoint = Acts::Vector3::Zero();
    if (bestPartnerOnly) {
      outerIndex.forEachInWindow(
          moduleKey(inner), clusterPhi(inner), halfWidth, [&](std::uint32_t o) {
            ++stats.pairsTried;
            const auto outerCache =
                Acts::StripSpacePointBuilder::makeConstrainedStripCache(
                    ends(clusters[o]), options);
            Acts::Vector3 point;
            if (Acts::StripSpacePointBuilder::computeConstrainedSpacePoint(
                    innerCache, outerCache, options, point)
                    .has_value()) {
              return;
            }
            const float offset = static_cast<float>(
                (point - Acts::Vector3(inner.center[0], inner.center[1],
                                       inner.center[2]))
                    .norm());
            if (offset < bestOffset) {
              bestOffset = offset;
              best = o;
              bestPoint = point;
            }
          });
      if (bestOffset == std::numeric_limits<float>::max()) {
        continue;
      }
      const StripCluster& outer = clusters[best];
      const bool ghost = best != i + 1;
      ++stats.formed;
      stats.ghosts += ghost ? 1 : 0;
      stats.trueCrossings += ghost ? 0 : 1;
      auto sp = out.createSpacePoint();
      setPosition(sp, static_cast<float>(bestPoint.x()),
                  static_cast<float>(bestPoint.y()),
                  static_cast<float>(bestPoint.z()));
      sp.varianceZ() = cylinder ? layer.varianceAlong : 0.f;
      sp.varianceR() = cylinder ? 0.f : layer.varianceAlong;
      Acts::OuterStripSpacePointCalibrationDetails details{};
      details.outerCenter = outer.center;
      for (int k = 0; k < 3; ++k) {
        details.innerToOuterSeparation[k] = outer.center[k] - inner.center[k];
      }
      details.outerHalfVector = outer.halfVector;
      details.innerHalfVector = inner.halfVector;
      sp.outerStripCalibrationDetails() = details;
      // a ghost still carries the inner cluster's particle: the point is the
      // inner cluster, only its partner is wrong
      sp.extra(particleColumn) = inner.particle;
      continue;
    }
    outerIndex.forEachInWindow(
        moduleKey(inner), clusterPhi(inner), halfWidth, [&](std::uint32_t o) {
          const StripCluster& outer = clusters[o];
          ++stats.pairsTried;
          const Acts::StripSpacePointBuilder::ConstrainedStripCache outerCache =
              Acts::StripSpacePointBuilder::makeConstrainedStripCache(
                  ends(outer), options);
          Acts::Vector3 point;
          if (Acts::StripSpacePointBuilder::computeConstrainedSpacePoint(
                  innerCache, outerCache, options, point)
                  .has_value()) {
            return;
          }
          // the generator writes the two clusters of a crossing next to
          // each other; anything else is a ghost, even of the same particle
          const bool ghost = o != i + 1;
          if (ghost && dropGhosts) {
            return;
          }
          ++stats.formed;
          stats.ghosts += ghost ? 1 : 0;
          stats.trueCrossings += ghost ? 0 : 1;
          auto sp = out.createSpacePoint();
          setPosition(sp, static_cast<float>(point.x()),
                      static_cast<float>(point.y()),
                      static_cast<float>(point.z()));
          sp.varianceZ() = cylinder ? layer.varianceAlong : 0.f;
          sp.varianceR() = cylinder ? 0.f : layer.varianceAlong;
          Acts::OuterStripSpacePointCalibrationDetails details{};
          details.outerCenter = outer.center;
          for (int k = 0; k < 3; ++k) {
            details.innerToOuterSeparation[k] =
                outer.center[k] - inner.center[k];
          }
          details.outerHalfVector = outer.halfVector;
          details.innerHalfVector = inner.halfVector;
          sp.outerStripCalibrationDetails() = details;
          sp.extra(particleColumn) = ghost ? ghostParticle : inner.particle;
        });
  }
  return stats;
}

/// Formation used as a filter: for every inner cluster, the outer clusters it
/// resolves with under the beam spot constraint at the given gap tolerance.
/// A cluster with none is not seeded on; one with several stays one point.
/// @return per inner cluster index, the begin offset into `partners`, sized
///         clusters.size() + 1
inline std::vector<std::uint32_t> formationPartners(
    const DetectorLayout& layout, const std::vector<StripLayer>& strips,
    const std::vector<StripCluster>& clusters, const ClusterIndex& outerIndex,
    float gapParameter, std::vector<std::uint32_t>& partners) {
  std::vector<std::uint32_t> begin(clusters.size() + 1, 0);
  partners.clear();
  Acts::StripSpacePointBuilder::ConstrainedOptions options;
  const auto ends = [](const StripCluster& c) {
    Acts::StripSpacePointBuilder::StripEnds e;
    e.top = Acts::Vector3(c.center[0] + c.halfVector[0],
                          c.center[1] + c.halfVector[1],
                          c.center[2] + c.halfVector[2]);
    e.bottom = Acts::Vector3(c.center[0] - c.halfVector[0],
                             c.center[1] - c.halfVector[1],
                             c.center[2] - c.halfVector[2]);
    return e;
  };
  for (std::uint32_t i = 0; i < clusters.size(); ++i) {
    begin[i] = static_cast<std::uint32_t>(partners.size());
    const StripCluster& inner = clusters[i];
    if (inner.side != 0) {
      continue;
    }
    const StripLayer& layer = strips[inner.layer];
    const bool cylinder =
        layout.layers[inner.layer].shape == SurfaceShape::Cylinder;
    const float r = clusterR(inner);
    const float slope = cylinder ? 1.f : r / std::abs(inner.center[2]);
    options.stripLengthGapTolerance =
        gapParameter * r * layer.gapOverStereo * slope;
    const float reach = 1.f + static_cast<float>(options.stripLengthTolerance) +
                        static_cast<float>(options.stripLengthGapTolerance) /
                            (2.f * layer.halfLength);
    const float sinStereo = 2.f * layer.sinHalfStereo * layer.cosHalfStereo;
    const float halfWidth =
        (layer.halfLength * reach * sinStereo + 2.f * layer.halfGap) / r;
    const auto innerCache =
        Acts::StripSpacePointBuilder::makeConstrainedStripCache(ends(inner),
                                                                options);
    outerIndex.forEachInWindow(
        moduleKey(inner), clusterPhi(inner), halfWidth, [&](std::uint32_t o) {
          const auto outerCache =
              Acts::StripSpacePointBuilder::makeConstrainedStripCache(
                  ends(clusters[o]), options);
          Acts::Vector3 point;
          if (!Acts::StripSpacePointBuilder::computeConstrainedSpacePoint(
                   innerCache, outerCache, options, point)
                   .has_value()) {
            partners.push_back(o);
          }
        });
  }
  begin[clusters.size()] = static_cast<std::uint32_t>(partners.size());
  return begin;
}

/// Points at the centre of the strip of every cluster of one sensor side.
///
/// @param layout the detector layout
/// @param strips the strip layer descriptors
/// @param clusters the clusters of the event
/// @param side the sensor side to take
/// @param out receives the points, with the `particleId` column
/// @param clusterOf receives the cluster index of each point
inline void makeClusterSpacePoints(
    const DetectorLayout& layout, const std::vector<StripLayer>& strips,
    const std::vector<StripCluster>& clusters, std::uint8_t side,
    Acts::SpacePointContainer& out, std::vector<std::uint32_t>& clusterOf,
    int truthPosition = 0,
    const std::vector<std::uint32_t>* keepBegin = nullptr) {
  out.clear();
  clusterOf.clear();
  if (!out.hasColumn("particleId")) {
    out.createColumns(kColumns);
    out.createColumn<std::uint32_t>("particleId");
  }
  auto particleColumn = out.column<std::uint32_t>("particleId");
  out.reserve(static_cast<std::uint32_t>(clusters.size() / 2));
  for (std::uint32_t i = 0; i < clusters.size(); ++i) {
    const StripCluster& c = clusters[i];
    if (c.side != side) {
      continue;
    }
    if (keepBegin != nullptr && (*keepBegin)[i] == (*keepBegin)[i + 1]) {
      continue;
    }
    const StripLayer& layer = strips[c.layer];
    const bool cylinder =
        layout.layers[c.layer].shape == SurfaceShape::Cylinder;
    // flat over the strip length
    const float lengthVariance =
        (2.f * layer.halfLength) * (2.f * layer.halfLength) / 12.f;
    std::array<float, 3> position = c.center;
    if (truthPosition != 0) {
      // diagnostic: the point of the true pair, resolved with a loose
      // tolerance; 1 takes all of it, 2 only the coordinate along the strip,
      // 3 only the transverse position
      const StripCluster& inner = side == 0 ? c : clusters[i - 1];
      const StripCluster& outer = side == 0 ? clusters[i + 1] : c;
      std::array<float, 3> separation{};
      for (int j = 0; j < 3; ++j) {
        separation[j] = outer.center[j] - inner.center[j];
      }
      const auto derived =
          Acts::detail::deriveOuterStripSpacePointCalibrationDetails(
              inner.halfVector, outer.halfVector, separation, outer.center);
      const float r0 = clusterR(c);
      std::array<float, 3> direction{c.center[0] / r0, c.center[1] / r0,
                                     c.center[2] / r0};
      std::array<float, 3> point{};
      if (Acts::detail::calibrateOuterStripSpacePoint(direction, derived, point,
                                                      100.f)) {
        if (truthPosition == 1) {
          position = point;
        } else if (truthPosition == 2) {
          if (cylinder) {
            position[2] = point[2];
          } else {
            const float scale = std::hypot(point[0], point[1]) / r0;
            position[0] *= scale;
            position[1] *= scale;
          }
        } else if (truthPosition == 3) {
          const float rp = std::hypot(point[0], point[1]);
          const float phiP = std::atan2(point[1], point[0]);
          const float rc = cylinder ? r0 : rp;
          (void)rc;
          if (cylinder) {
            position[0] = r0 * std::cos(phiP);
            position[1] = r0 * std::sin(phiP);
          } else {
            position[0] = r0 * std::cos(phiP);
            position[1] = r0 * std::sin(phiP);
          }
        }
      }
    }
    auto sp = out.createSpacePoint();
    setPosition(sp, position[0], position[1], position[2]);
    sp.varianceZ() = cylinder ? lengthVariance : 0.f;
    sp.varianceR() = cylinder ? 0.f : lengthVariance;
    sp.extra(particleColumn) = c.particle;
    clusterOf.push_back(i);
  }
}

/// Knobs of `refineSeeds`.
struct RefineConfig {
  /// How far outside its strip a resolved point may land, in half lengths
  float tolerance = 1.05f;
  /// Maximum difference of cot(theta) between the two resolved doublets
  float cotThetaDiffMax = 0.15f;
  /// z origin window of the bottom doublet, on resolved positions
  float collisionRegionMin = -500.f;
  float collisionRegionMax = 500.f;
  /// Distance back along the seed direction to put the formation vertex at
  float lever = 1000.f;
};

/// What `refineSeeds` did.
struct RefineStats {
  std::size_t seedsIn{};
  std::size_t seedsOut{};
  std::size_t partnerCandidates{};
  std::size_t clustersWithoutPartner{};
  std::size_t combinations{};
  std::size_t truePartners{};
  std::size_t chosenPartners{};
  /// Diagnostics on the true seeds (three clusters of one particle)
  std::size_t trueSeeds{};
  std::size_t trueClusterPartnerInWindow{};
  std::size_t trueClusterPartnerResolved{};
  std::size_t trueClusters{};
  std::size_t trueSeedsOut{};
};

/// Resolve the three clusters of each seed with its own direction and keep the
/// seeds whose resolved points agree in r-z.
///
/// @param layout the detector layout
/// @param strips the strip layer descriptors
/// @param clusters the clusters of the event
/// @param clusterOf the cluster index of each point of `seedsIn`'s container
/// @param otherIndex the index of the clusters of the other sensor side
/// @param side the sensor side the seeds were built from
/// @param cfg the cuts
/// @param seedsIn the seeds of single clusters
/// @param seedsOut receives the kept seeds, on the same container
/// @param stats accumulates statistics
inline void refineSeeds(const DetectorLayout& /*layout*/,
                        const std::vector<StripLayer>& strips,
                        const std::vector<StripCluster>& clusters,
                        const std::vector<std::uint32_t>& clusterOf,
                        const ClusterIndex& otherIndex, std::uint8_t side,
                        const RefineConfig& cfg,
                        const Acts::SeedContainer& seedsIn,
                        Acts::SeedContainer& seedsOut, RefineStats& stats) {
  seedsOut.assignSpacePointContainer(seedsIn.spacePointContainer());
  seedsOut.reserve(seedsIn.size());
  const Acts::SpacePointContainer& points = seedsIn.spacePointContainer();

  constexpr std::size_t kMaxCandidates = 4;
  std::array<std::array<std::array<float, 3>, kMaxCandidates>, 3> resolved{};
  std::array<std::array<bool, kMaxCandidates>, 3> truePartner{};
  std::array<std::size_t, 3> numResolved{};

  for (const auto seed : seedsIn) {
    ++stats.seedsIn;
    const auto indices = seed.spacePointIndices();
    std::array<std::array<float, 3>, 3> p{};
    for (int k = 0; k < 3; ++k) {
      const auto sp = points[indices[k]];
      p[k] = {sp.x(), sp.y(), sp.z()};
    }

    // circle through the three transverse positions, which single clusters
    // measure well; the tangent at each point is the transverse direction
    const float ax = p[1][0] - p[0][0];
    const float ay = p[1][1] - p[0][1];
    const float bx = p[2][0] - p[0][0];
    const float by = p[2][1] - p[0][1];
    const float d = 2.f * (ax * by - ay * bx);
    const float chordX = p[2][0] - p[0][0];
    const float chordY = p[2][1] - p[0][1];
    // r-z slope from the coarse coordinates; it barely matters to the
    // resolution, which reads the direction across the strips
    const float chordXY = std::hypot(chordX, chordY);
    const float cot = (p[2][2] - p[0][2]) / chordXY;

    const std::uint32_t c0 = clusterOf[indices[0]];
    const bool isTrue =
        clusters[c0].particle == clusters[clusterOf[indices[1]]].particle &&
        clusters[c0].particle == clusters[clusterOf[indices[2]]].particle;
    stats.trueSeeds += isTrue ? 1 : 0;
    bool ok = true;
    for (int k = 0; k < 3 && ok; ++k) {
      float tx = chordX / chordXY;
      float ty = chordY / chordXY;
      if (std::abs(d) > 1e-6f) {
        const float a2 = ax * ax + ay * ay;
        const float b2 = bx * bx + by * by;
        const float cx = p[0][0] + (by * a2 - ay * b2) / d;
        const float cy = p[0][1] + (ax * b2 - bx * a2) / d;
        // perpendicular to the radius of the circle, oriented outwards
        float ux = -(p[k][1] - cy);
        float uy = p[k][0] - cx;
        const float norm = std::hypot(ux, uy);
        ux /= norm;
        uy /= norm;
        if (ux * chordX + uy * chordY < 0.f) {
          ux = -ux;
          uy = -uy;
        }
        tx = ux;
        ty = uy;
      }
      const float dn = std::sqrt(1.f + cot * cot);
      const std::array<float, 3> direction{tx / dn, ty / dn, cot / dn};

      const std::uint32_t ci = clusterOf[indices[k]];
      const StripCluster& c = clusters[ci];
      // the generator writes the two clusters of a crossing next to each other
      const std::uint32_t truePartnerIndex = side == 0 ? ci + 1 : ci - 1;
      stats.trueClusters += isTrue ? 1 : 0;
      const StripLayer& layer = strips[c.layer];
      const float r = clusterR(c);
      const float sinStereo = 2.f * layer.sinHalfStereo * layer.cosHalfStereo;
      const float halfWidth =
          (layer.halfLength * cfg.tolerance * sinStereo + 2.f * layer.halfGap) /
          r;
      numResolved[k] = 0;
      otherIndex.forEachInWindow(
          moduleKey(c), clusterPhi(c), halfWidth, [&](std::uint32_t o) {
            ++stats.partnerCandidates;
            if (numResolved[k] == kMaxCandidates) {
              return;
            }
            const StripCluster& other = clusters[o];
            const bool isTruePartner = isTrue && o == truePartnerIndex;
            stats.trueClusterPartnerInWindow += isTruePartner ? 1 : 0;
            const StripCluster& inner = side == 0 ? c : other;
            const StripCluster& outer = side == 0 ? other : c;
            std::array<float, 3> separation{};
            for (int j = 0; j < 3; ++j) {
              separation[j] = outer.center[j] - inner.center[j];
            }
            const auto derived =
                Acts::detail::deriveOuterStripSpacePointCalibrationDetails(
                    inner.halfVector, outer.halfVector, separation,
                    outer.center);
            std::array<float, 3> point{};
            if (!Acts::detail::calibrateOuterStripSpacePoint(
                    direction, derived, point, cfg.tolerance)) {
              return;
            }
            stats.trueClusterPartnerResolved += isTruePartner ? 1 : 0;
            resolved[k][numResolved[k]] = point;
            truePartner[k][numResolved[k]] = other.particle == c.particle;
            ++numResolved[k];
          });
      if (numResolved[k] == 0) {
        ++stats.clustersWithoutPartner;
        ok = false;
      }
    }
    if (!ok) {
      continue;
    }

    // the combination of partners whose two doublets agree best in r-z
    float best = std::numeric_limits<float>::max();
    std::array<std::size_t, 3> chosen{};
    for (std::size_t i = 0; i < numResolved[0]; ++i) {
      const auto& b = resolved[0][i];
      const float rb = std::hypot(b[0], b[1]);
      for (std::size_t j = 0; j < numResolved[1]; ++j) {
        const auto& m = resolved[1][j];
        const float rm = std::hypot(m[0], m[1]);
        const float cotB = (m[2] - b[2]) / (rm - rb);
        const float zOrigin = m[2] - rm * cotB;
        if (zOrigin < cfg.collisionRegionMin ||
            zOrigin > cfg.collisionRegionMax) {
          continue;
        }
        for (std::size_t l = 0; l < numResolved[2]; ++l) {
          ++stats.combinations;
          const auto& t = resolved[2][l];
          const float rt = std::hypot(t[0], t[1]);
          const float cotT = (t[2] - m[2]) / (rt - rm);
          const float diff = std::abs(cotT - cotB);
          if (diff < best) {
            best = diff;
            chosen = {i, j, l};
          }
        }
      }
    }
    if (best > cfg.cotThetaDiffMax) {
      continue;
    }
    for (int k = 0; k < 3; ++k) {
      ++stats.chosenPartners;
      stats.truePartners += truePartner[k][chosen[k]] ? 1 : 0;
    }
    ++stats.seedsOut;
    stats.trueSeedsOut += isTrue ? 1 : 0;
    auto out = seedsOut.createSeed();
    out.assignSpacePointIndices(std::array<Acts::SpacePointIndex, 3>{
        indices[0], indices[1], indices[2]});
    out.quality() = seed.quality();
    out.vertexZ() = seed.vertexZ();
  }
}

}  // namespace ActsTests::StripClusterSeeding

namespace ActsTests::StripClusterSeeding {

/// Every candidate partner of each cluster of one sensor side, with the
/// calibration terms that depend on the pair only. No direction is chosen
/// yet, so a cluster stays one point however many partners it has.
class PartnerTable {
 public:
  /// @param strips the strip layer descriptors
  /// @param clusters the clusters of the event
  /// @param clusterOf the cluster of each seeding point
  /// @param otherIndex the index of the other sensor side
  /// @param side the side the seeding points come from
  /// @param reach how far along the strip a pair may resolve, in half lengths
  void build(const std::vector<StripLayer>& strips,
             const std::vector<StripCluster>& clusters,
             const std::vector<std::uint32_t>& clusterOf,
             const ClusterIndex& otherIndex, std::uint8_t side, float reach) {
    m_begin.assign(clusterOf.size() + 1, 0);
    m_details.clear();
    m_true.clear();
    for (std::size_t p = 0; p < clusterOf.size(); ++p) {
      m_begin[p] = static_cast<std::uint32_t>(m_details.size());
      const std::uint32_t ci = clusterOf[p];
      const StripCluster& c = clusters[ci];
      const StripLayer& layer = strips[c.layer];
      const float sinStereo = 2.f * layer.sinHalfStereo * layer.cosHalfStereo;
      const float halfWidth =
          (layer.halfLength * reach * sinStereo + 2.f * layer.halfGap) /
          clusterR(c);
      otherIndex.forEachInWindow(
          moduleKey(c), clusterPhi(c), halfWidth, [&](std::uint32_t o) {
            const StripCluster& other = clusters[o];
            const StripCluster& inner = side == 0 ? c : other;
            const StripCluster& outer = side == 0 ? other : c;
            std::array<float, 3> separation{};
            for (int j = 0; j < 3; ++j) {
              separation[j] = outer.center[j] - inner.center[j];
            }
            m_details.push_back(
                Acts::detail::deriveOuterStripSpacePointCalibrationDetails(
                    inner.halfVector, outer.halfVector, separation,
                    outer.center));
            m_true.push_back(o == (side == 0 ? ci + 1 : ci - 1));
          });
    }
    m_begin[clusterOf.size()] = static_cast<std::uint32_t>(m_details.size());
  }

  /// The same from the partners formation resolved, inner side only.
  void build(const std::vector<StripCluster>& clusters,
             const std::vector<std::uint32_t>& clusterOf,
             const std::vector<std::uint32_t>& partnerBegin,
             const std::vector<std::uint32_t>& partnerList) {
    m_begin.assign(clusterOf.size() + 1, 0);
    m_details.clear();
    m_true.clear();
    for (std::size_t p = 0; p < clusterOf.size(); ++p) {
      m_begin[p] = static_cast<std::uint32_t>(m_details.size());
      const std::uint32_t ci = clusterOf[p];
      const StripCluster& inner = clusters[ci];
      for (std::uint32_t k = partnerBegin[ci]; k < partnerBegin[ci + 1]; ++k) {
        const std::uint32_t o = partnerList[k];
        const StripCluster& outer = clusters[o];
        std::array<float, 3> separation{};
        for (int j = 0; j < 3; ++j) {
          separation[j] = outer.center[j] - inner.center[j];
        }
        m_details.push_back(
            Acts::detail::deriveOuterStripSpacePointCalibrationDetails(
                inner.halfVector, outer.halfVector, separation, outer.center));
        m_true.push_back(o == ci + 1);
      }
    }
    m_begin[clusterOf.size()] = static_cast<std::uint32_t>(m_details.size());
  }

  std::uint32_t begin(std::size_t point) const { return m_begin[point]; }
  std::uint32_t end(std::size_t point) const { return m_begin[point + 1]; }
  const Acts::OuterStripSpacePointCalibrationDetailsDerived& details(
      std::uint32_t k) const {
    return m_details[k];
  }
  std::size_t size() const { return m_details.size(); }
  std::size_t numTrue() const {
    return static_cast<std::size_t>(std::ranges::count(m_true, true));
  }

 private:
  std::vector<std::uint32_t> m_begin;
  std::vector<Acts::OuterStripSpacePointCalibrationDetailsDerived> m_details;
  std::vector<bool> m_true;
};

/// The doublet cut of the cluster seeding: resolve both clusters against the
/// chord of the doublet with each of their candidate partners, and keep the
/// doublet if one combination passes the r-z cuts on resolved positions. Runs
/// after the transverse cuts of the doublet finder, so only doublets that
/// survive those pay for it.
struct PartnerCalibrationCut {
  const PartnerTable* partners{};
  float tolerance = 1.75f;
  float collisionRegionMin = -500.f;
  float collisionRegionMax = 500.f;
  float cotThetaMax = 5.f;
  mutable std::size_t calls{};
  mutable std::size_t passed{};

  /// @param middleIndex index of the middle seeding point
  /// @param otherIndex index of the other seeding point
  /// @param middle middle position (coarse)
  /// @param other other position (coarse)
  /// @param isBottom whether the other point is the bottom one
  bool operator()(std::uint32_t middleIndex, std::uint32_t otherIndex,
                  const std::array<float, 3>& middle,
                  const std::array<float, 3>& other, bool isBottom) const {
    ++calls;
    // the chord, pointing outwards along the track
    std::array<float, 3> chord{other[0] - middle[0], other[1] - middle[1],
                               other[2] - middle[2]};
    if (isBottom) {
      chord = {-chord[0], -chord[1], -chord[2]};
    }
    constexpr std::size_t kMax = 4;
    std::array<std::array<float, 3>, kMax> resolvedM{};
    std::array<std::array<float, 3>, kMax> resolvedO{};
    std::size_t numM = 0;
    std::size_t numO = 0;
    for (std::uint32_t k = partners->begin(middleIndex);
         k < partners->end(middleIndex) && numM < kMax; ++k) {
      if (Acts::detail::calibrateOuterStripSpacePoint(
              chord, partners->details(k), resolvedM[numM], tolerance)) {
        ++numM;
      }
    }
    if (numM == 0) {
      return false;
    }
    for (std::uint32_t k = partners->begin(otherIndex);
         k < partners->end(otherIndex) && numO < kMax; ++k) {
      if (Acts::detail::calibrateOuterStripSpacePoint(
              chord, partners->details(k), resolvedO[numO], tolerance)) {
        ++numO;
      }
    }
    for (std::size_t i = 0; i < numM; ++i) {
      const float rM = std::hypot(resolvedM[i][0], resolvedM[i][1]);
      for (std::size_t j = 0; j < numO; ++j) {
        const float rO = std::hypot(resolvedO[j][0], resolvedO[j][1]);
        const float deltaR = rM - rO;
        const float deltaZ = resolvedM[i][2] - resolvedO[j][2];
        if (std::abs(deltaZ) > cotThetaMax * std::abs(deltaR)) {
          continue;
        }
        const float zOrigin = resolvedM[i][2] - rM * deltaZ / deltaR;
        if (zOrigin < collisionRegionMin || zOrigin > collisionRegionMax) {
          continue;
        }
        ++passed;
        return true;
      }
    }
    return false;
  }
};

}  // namespace ActsTests::StripClusterSeeding

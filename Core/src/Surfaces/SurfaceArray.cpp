// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/SurfaceArray.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/Ranges.hpp"
#include "Acts/Utilities/detail/MultiAxisHelper.hpp"
#include "Acts/Utilities/detail/OstreamStateGuard.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <map>
#include <numbers>
#include <ranges>
#include <string_view>
#include <utility>

namespace Acts {

const bool s_gridPreFilterEnabled = [] {
  const char* env = std::getenv("ACTS_GRID_PREFILTER");
  return env == nullptr || std::string_view(env) != "0";
}();

/// Base interface for all surface lookups.
struct SurfaceArray::ISurfaceGridLookup {
  virtual ~ISurfaceGridLookup() = default;

  /// Fill provided surfaces into the contained @c Grid.
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surfaces Input surface pointers
  virtual void fill(const GeometryContext& gctx,
                    std::span<const Surface* const> surfaces) = 0;

  /// Get all surfaces in bin given by the global bin index
  /// @param bin the global bin index
  /// @return span of surface pointers of the bin at that position
  virtual std::span<const Surface* const> at(std::size_t bin) const = 0;

  /// Performs lookup at @c pos and returns bin content as const reference
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position Lookup position
  /// @param direction Lookup direction
  /// @return A span of surface pointers
  virtual std::span<const Surface* const> at(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction) const = 0;

  /// Get all surfaces in bin given by local grid indices and neighbor
  /// distance.
  /// @param gridIndices the local grid indices
  /// @param neighborDistance the neighbor distance to include in the lookup
  /// @return span of surface pointers of the bin at that position and its neighbors
  virtual std::span<const Surface* const> neighbors(
      std::array<std::size_t, 2> gridIndices,
      std::uint8_t neighborDistance) const = 0;

  /// Performs a lookup at @c pos, but returns neighbors as well
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position Lookup position
  /// @param direction Lookup direction
  /// @return A span of surface pointers
  virtual std::span<const Surface* const> neighbors(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction) const = 0;

  /// Neighbor lookup carrying the per-surface pre-rejection data
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position Lookup position
  /// @param direction Lookup direction
  /// @return the pack with its extents and the query point
  virtual SurfaceArray::NeighborQuery neighborQuery(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction) const = 0;

  /// Returns the total size of the grid (including under/overflow bins)
  /// @return Size of the grid data structure
  virtual std::size_t size() const = 0;

  /// Gets the center position of bin @c bin in global coordinates
  /// @param bin the global bin index
  /// @return The bin center
  virtual Vector3 getBinCenter(std::size_t bin) const = 0;

  /// Returns copies of the axes used in the grid as @c AnyAxis
  /// @return The axes
  /// @note This returns copies. Use for introspection and querying.
  virtual std::vector<const IAxis*> getAxes() const = 0;

  /// Get the representative surface used for this lookup
  /// @return Surface pointer
  virtual const Surface* surfaceRepresentation() const = 0;

  /// Checks if global bin is valid
  /// @param bin the global bin index
  /// @return bool if the bin is valid
  /// @note Valid means that the index points to a bin which is not a under
  ///       or overflow bin or out of range in any axis.
  virtual bool isValidBin(std::size_t bin) const = 0;

  /// The binning values described by this surface grid lookup. They are in
  /// order of the axes (optional) and empty for eingle lookups
  /// @return Vector of axis directions for binning
  virtual std::vector<AxisDirection> binningValues() const { return {}; }

  /// Get the number of local bins in each dimension. This is used to
  /// determine the size of the grid for neighbor lookups.
  /// @return Array of number of local bins in each dimension
  virtual std::array<std::size_t, 2> numLocalBins() const = 0;

  /// Get the maximum neighbor distance that is supported by this lookup. This
  /// is used to determine how many neighbors to include in neighbor lookups.
  /// @return Maximum neighbor distance
  virtual std::uint8_t maxNeighborDistance() const = 0;
};

namespace {

struct SingleElementLookupImpl final : SurfaceArray::ISurfaceGridLookup {
  explicit SingleElementLookupImpl(const Surface* element)
      : m_element({element}) {}

  std::span<const Surface* const> at(std::size_t bin) const override {
    if (bin != 0) {
      throw std::out_of_range(
          "SingleElementLookupImpl only contains one bin with index 0");
    }
    return m_element;
  }

  std::span<const Surface* const> at(
      const GeometryContext& /*gctx*/, const Vector3& /*position*/,
      const Vector3& /*direction*/) const override {
    return m_element;
  }

  std::span<const Surface* const> neighbors(
      std::array<std::size_t, 2> gridIndices,
      std::uint8_t neighborDistance) const override {
    if (gridIndices != std::array<std::size_t, 2>{0, 0} ||
        neighborDistance != 0) {
      throw std::out_of_range(
          "SingleElementLookupImpl only contains one bin with zero neighbor "
          "distance");
    }
    return m_element;
  }

  std::span<const Surface* const> neighbors(
      const GeometryContext& /*gctx*/, const Vector3& /*position*/,
      const Vector3& /*direction*/) const override {
    return m_element;
  }

  SurfaceArray::NeighborQuery neighborQuery(
      const GeometryContext& /*gctx*/, const Vector3& /*position*/,
      const Vector3& /*direction*/) const override {
    return {.surfaces = m_element, .extents = {}};
  }

  std::size_t size() const override { return 1; }

  Vector3 getBinCenter(std::size_t /*bin*/) const override {
    return Vector3(0, 0, 0);
  }

  std::vector<const IAxis*> getAxes() const override { return {}; }

  const Surface* surfaceRepresentation() const override { return nullptr; }

  void fill(const GeometryContext& /*gctx*/,
            std::span<const Surface* const> /*surfaces*/) override {}

  bool isValidBin(std::size_t bin) const override { return bin == 0; }

  std::array<std::size_t, 2> numLocalBins() const override { return {1, 1}; }

  std::uint8_t maxNeighborDistance() const override { return 0; }

 private:
  std::vector<const Surface*> m_element;
};

template <class Axis1, class Axis2>
struct SurfaceGridLookupImpl final : SurfaceArray::ISurfaceGridLookup {
  SurfaceGridLookupImpl(std::shared_ptr<RegularSurface> representative,
                        double tolerance, std::tuple<Axis1, Axis2> axes,
                        std::vector<AxisDirection> binValues = {},
                        std::uint8_t maxNeighborDistance = 1)
      : m_representative(std::move(representative)),
        m_tolerance(tolerance),
        m_axes(std::move(axes)),
        m_binValues(std::move(binValues)),
        m_maxNeighborDistance(maxNeighborDistance) {
    // The representative surface is fixed for the lifetime of the lookup, so
    // resolve the cylinder case once instead of dynamic_cast-ing on every
    // grid lookup.
    if (const auto* cylinderBounds =
            dynamic_cast<const CylinderBounds*>(&m_representative->bounds());
        cylinderBounds != nullptr) {
      m_cylinderRadius = cylinderBounds->get(CylinderBounds::eR);
    }
    m_fillingGrid.resize(size());
  }

  void fill(const GeometryContext& gctx,
            std::span<const Surface* const> surfaces) override {
    for (const Surface* surface : surfaces) {
      const std::optional<std::size_t> globalBin =
          fillSurfaceToBinMapping(gctx, *surface);
      if (!globalBin.has_value()) {
        continue;
      }

      fillBinToSurfaceMapping(gctx, *surface, *globalBin);
    }

    for (std::vector<const Surface*>& binSurfaces : m_fillingGrid) {
      std::ranges::sort(binSurfaces);
      const auto last = std::ranges::unique(binSurfaces);
      binSurfaces.erase(last.begin(), last.end());
      binSurfaces.shrink_to_fit();
    }

    checkGrid(surfaces);

    computeSurfaceExtents(gctx, surfaces);
    populateNeighborCache();
  }

  std::span<const Surface* const> at(std::size_t globalBin) const override {
    return m_fillingGrid.at(globalBin);
  }

  std::span<const Surface* const> at(const GeometryContext& gctx,
                                     const Vector3& position,
                                     const Vector3& direction) const override {
    const std::optional<GridIndex> localBins = findLocalBin2D(
        gctx, position, direction, std::numeric_limits<double>::infinity());
    if (!localBins.has_value()) {
      return {};
    }
    const std::size_t globalBin = globalBinFromLocalBins3D(*localBins, 0);
    return m_neighborSurfacePacks.at(globalBin);
  }

  std::span<const Surface* const> neighbors(
      std::array<std::size_t, 2> gridIndices,
      std::uint8_t neighborDistance) const override {
    return m_neighborSurfacePacks.at(
        globalBinFromLocalBins3D(gridIndices, neighborDistance));
  }

  std::span<const Surface* const> neighbors(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction) const override {
    std::size_t globalBin = 0;
    if (!lookupBin(gctx, position, direction, globalBin, nullptr, nullptr)) {
      return {};
    }
    return m_neighborSurfacePacks.at(globalBin);
  }

  SurfaceArray::NeighborQuery neighborQuery(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction) const override {
    std::size_t globalBin = 0;
    GridPoint gridLocal{};
    double cosIncidence = 1;
    if (!lookupBin(gctx, position, direction, globalBin, &gridLocal,
                   &cosIncidence)) {
      return {};
    }

    // Lateral distance a trajectory covers while travelling from the
    // representative surface to a surface offset from it by at most
    // m_maxNormalOffset. That is how far the pierce point on the surface can
    // sit from the pierce point on the representative.
    // A point of a surface sits at most m_maxNormalOffset off the
    // representative. The trajectory covers that perpendicular distance in
    // m_maxNormalOffset / cos(theta) of path, which displaces it by
    // m_maxNormalOffset * tan(theta) within the representative.
    const double c = std::max(std::abs(cosIncidence), 1e-3);
    // The safety factor covers what the flat-surface argument above ignores:
    // the representative surface curves away over the lateral reach, and the
    // trajectory is a helix rather than the straight line assumed here. An
    // audit over 4.7e6 rejections found the bare bound short by 0.06 mm once.
    constexpr double safety = 1.1;
    const double lateral =
        safety * m_maxNormalOffset * std::sqrt(std::max(0.0, 1 - c * c)) / c;

    return {.surfaces = m_neighborSurfacePacks.at(globalBin),
            .extents = m_neighborExtentPacks.at(globalBin),
            .gridLocal = gridLocal,
            .margin = {(lateral + s_onSurfaceTolerance) * m_marginScale[0],
                       (lateral + s_onSurfaceTolerance) * m_marginScale[1]},
            .wrap = {m_cylinderRadius > 0, m_cylinderRadius <= 0}};
  }

  std::size_t size() const override {
    const GridIndex nBins = numLocalBins2D();
    return (nBins[0] + 2) * (nBins[1] + 2);
  }

  std::vector<AxisDirection> binningValues() const override {
    return m_binValues;
  }

  Vector3 getBinCenter(std::size_t bin) const override {
    const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
    const GridPoint gridLocal = binCenter(localBinsFromGlobalBin2D(bin));
    const Vector2 surfaceLocal = gridToSurfaceLocal(gridLocal);
    return m_representative->localToGlobal(gctx, surfaceLocal);
  }

  std::vector<const IAxis*> getAxes() const override {
    return {&std::get<0>(m_axes), &std::get<1>(m_axes)};
  }

  const Surface* surfaceRepresentation() const override {
    return m_representative.get();
  }

  bool isValidBin(std::size_t globalBin) const override {
    const GridIndex indices = localBinsFromGlobalBin2D(globalBin);
    return isValidBin(indices);
  }

  std::array<std::size_t, 2> numLocalBins() const override {
    return numLocalBins2D();
  }

  std::uint8_t maxNeighborDistance() const override {
    return m_maxNeighborDistance;
  }

 private:
  using GridIndex = std::array<std::size_t, 2>;
  using GridPoint = std::array<double, 2>;

  std::shared_ptr<RegularSurface> m_representative;
  double m_tolerance{};
  // needs to be a tuple for the grid_helper functions
  std::tuple<Axis1, Axis2> m_axes;
  std::vector<AxisDirection> m_binValues;
  std::uint8_t m_maxNeighborDistance{};
  // radius of the representative surface if it is a cylinder, 0 otherwise
  double m_cylinderRadius{0};

  // legacy grid for filling and for deprecated lookup methods.
  // TODO: remove this once deprecated lookup methods are removed and filling is
  // done directly into the neighbor cache
  std::vector<std::vector<const Surface*>> m_fillingGrid;

  // containers to store the surfaces in the custom grid
  std::vector<const Surface*> m_surfacePacks;
  std::vector<std::span<const Surface* const>> m_neighborSurfacePacks;
  // grid-local extents of the surfaces, parallel to the two above
  std::vector<SurfaceArray::SurfaceExtent> m_extentPacks;
  std::vector<std::span<const SurfaceArray::SurfaceExtent>>
      m_neighborExtentPacks;
  // scratch, only used while filling
  std::map<const Surface*, SurfaceArray::SurfaceExtent> m_surfaceExtents;
  // largest distance of any surface from the representative along its normal
  double m_maxNormalOffset{0};
  // converts a lateral distance into grid units on each axis
  std::array<double, 2> m_marginScale{1, 1};

  bool isValidBin(const GridIndex& indices) const {
    const GridIndex nBins = numLocalBins2D();
    for (std::size_t i = 0; i < indices.size(); ++i) {
      const std::size_t idx = indices.at(i);
      if (idx <= 0 || idx >= nBins.at(i) + 1) {
        return false;
      }
    }
    return true;
  }

  GridIndex numLocalBins2D() const {
    return {std::get<0>(m_axes).getNBins(), std::get<1>(m_axes).getNBins()};
  }

  GridIndex localBinsFromPosition2D(const GridPoint& point) const {
    return detail::MultiAxisHelper::getLocalBinsFromPoint(point, m_axes);
  }

  GridIndex localBinsFromGlobalBin2D(std::size_t globalBin) const {
    return detail::MultiAxisHelper::getLocalBinsFromGlobalBin(globalBin,
                                                              m_axes);
  }

  std::size_t globalBinFromLocalBins2D(const GridIndex& localBins) const {
    return detail::MultiAxisHelper::getGlobalBinFromLocalBins(localBins,
                                                              m_axes);
  }

  std::size_t globalBinFromLocalBins3D(const GridIndex& localBins,
                                       std::uint8_t neighborDistance) const {
    const std::size_t globalGridBin =
        detail::MultiAxisHelper::getGlobalBinFromLocalBins(localBins, m_axes);
    return globalGridBin * (m_maxNeighborDistance + 1) + neighborDistance;
  }

  GridPoint binCenter(const GridIndex& localBins) const {
    return detail::MultiAxisHelper::getBinCenter(localBins, m_axes);
  }

  /// map surface center to grid
  std::optional<std::size_t> fillSurfaceToBinMapping(
      const GeometryContext& gctx, const Surface& surface) {
    const Vector3 position =
        surface.referencePosition(gctx, AxisDirection::AxisR);
    const Vector3 normal = m_representative->normal(gctx, position);
    const std::optional<Vector2> surfaceLocal =
        findSurfaceLocal(gctx, position, normal, m_tolerance);
    if (!surfaceLocal.has_value()) {
      return std::nullopt;
    }
    const GridPoint gridLocal = surfaceToGridLocal(*surfaceLocal);
    const GridIndex localBins = localBinsFromPosition2D(gridLocal);
    const std::size_t globalBin = globalBinFromLocalBins2D(localBins);
    m_fillingGrid.at(globalBin).push_back(&surface);
    return globalBin;
  }

  /// flood fill neighboring bins given a starting bin
  void fillBinToSurfaceMapping(const GeometryContext& gctx,
                               const Surface& surface, std::size_t startBin) {
    const GridIndex startIndices = localBinsFromGlobalBin2D(startBin);
    const auto startNeighborIndices =
        detail::MultiAxisHelper::neighborHoodIndices(startIndices, 1u, m_axes);

    std::set<std::size_t> visited({startBin});
    std::vector<std::size_t> queue(startNeighborIndices.begin(),
                                   startNeighborIndices.end());

    while (!queue.empty()) {
      const std::size_t current = queue.back();
      queue.pop_back();

      // Skip overflow bins as they do not produce a valid bin center
      if (!isValidBin(current)) {
        continue;
      }
      if (visited.contains(current)) {
        continue;
      }

      const GridIndex currentIndices = localBinsFromGlobalBin2D(current);
      visited.insert(current);

      const GridPoint gridLocal = binCenter(currentIndices);
      const Vector2 surfaceLocal = gridToSurfaceLocal(gridLocal);
      const Vector3 normal = m_representative->normal(gctx, surfaceLocal);
      const Vector3 global =
          m_representative->localToGlobal(gctx, surfaceLocal, normal);

      const Intersection3D intersection =
          surface.intersect(gctx, global, normal, BoundaryTolerance::None())
              .closest();
      if (!intersection.isValid() ||
          std::abs(intersection.pathLength()) > m_tolerance) {
        continue;
      }
      m_fillingGrid.at(current).push_back(&surface);

      const auto neighborIndices = detail::MultiAxisHelper::neighborHoodIndices(
          currentIndices, 1u, m_axes);
      queue.insert(queue.end(), neighborIndices.begin(), neighborIndices.end());
    }
  }

  /// calculate neighbors for every bin and store in map
  /// Precompute, per surface, its axis-aligned extent in grid-local
  /// coordinates and how far it sits off the representative surface. Both feed
  /// the conservative pre-rejection in NeighborQuery.
  void computeSurfaceExtents(const GeometryContext& gctx,
                             std::span<const Surface* const> surfaces) {
    // On a phi axis the extent of a surface straddling the seam is meaningless
    // as a min/max pair; such surfaces get an infinite extent so they are never
    // rejected.
    const std::size_t angleAxis = m_cylinderRadius > 0 ? 0 : 1;
    constexpr double inf = std::numeric_limits<double>::infinity();

    double rMin = inf;
    m_maxNormalOffset = 0;
    for (const Surface* surface : surfaces) {
      std::array<double, 4> extent = {inf, -inf, inf, -inf};
      const std::vector<Vector3> corners =
          surface->polyhedronRepresentation(gctx, 4u).vertices;

      // Both the radial grid coordinate of a disc and the radial offset from a
      // cylinder have their extremum in the interior of an edge, not at a
      // corner, whenever the edge passes closest to the beam axis in the
      // middle. Sample those points too or the extent comes out too small and
      // the filter rejects surfaces the trajectory really crosses.
      std::vector<Vector3> samples = corners;
      samples.reserve(2 * corners.size());
      for (std::size_t i = 0; i < corners.size(); ++i) {
        const Vector3& a = corners[i];
        const Vector3& b = corners[(i + 1) % corners.size()];
        const Vector2 d = (b - a).head<2>();
        const double dd = d.squaredNorm();
        if (dd <= 0) {
          continue;
        }
        const double t = std::clamp(-a.head<2>().dot(d) / dd, 0.0, 1.0);
        samples.push_back(a + t * (b - a));
      }

      for (const Vector3& sample : samples) {
        const GridPoint gridLocal = projectToGridLocal(gctx, sample);
        extent[0] = std::min(extent[0], gridLocal[0]);
        extent[1] = std::max(extent[1], gridLocal[0]);
        extent[2] = std::min(extent[2], gridLocal[1]);
        extent[3] = std::max(extent[3], gridLocal[1]);
        m_maxNormalOffset =
            std::max(m_maxNormalOffset, normalOffset(gctx, sample));
        rMin = std::min(rMin, sample.head<2>().norm());
      }
      if (extent[2 * angleAxis + 1] - extent[2 * angleAxis] >
          std::numbers::pi) {
        extent[2 * angleAxis] = -inf;
        extent[2 * angleAxis + 1] = inf;
      }
      // Round the half-width up when narrowing to float so the test can only
      // ever become more permissive, never less.
      const auto toCentre = [](double lo, double hi) {
        return static_cast<float>(0.5 * (lo + hi));
      };
      const auto toHalf = [](double lo, double hi) {
        const double half = 0.5 * (hi - lo);
        if (!std::isfinite(half)) {
          return std::numeric_limits<float>::max();
        }
        return std::nextafter(static_cast<float>(half),
                              std::numeric_limits<float>::max());
      };
      SurfaceArray::SurfaceExtent packed;
      packed.u0 = std::isfinite(extent[0]) ? toCentre(extent[0], extent[1]) : 0;
      packed.hu = toHalf(extent[0], extent[1]);
      packed.v0 = std::isfinite(extent[2]) ? toCentre(extent[2], extent[3]) : 0;
      packed.hv = toHalf(extent[2], extent[3]);
      m_surfaceExtents[surface] = packed;
    }

    // A lateral distance in mm becomes an angle on the phi axis. Divide by the
    // smallest radius present so the margin is never underestimated.
    const double rSafe = std::isfinite(rMin) && rMin > 0 ? rMin : 1;
    m_marginScale = {1, 1};
    m_marginScale[angleAxis] = 1 / rSafe;
  }

  void populateNeighborCache() {
    m_surfacePacks.clear();
    m_neighborSurfacePacks.clear();
    m_extentPacks.clear();
    m_neighborExtentPacks.clear();

    using SurfacePackRange = std::pair<std::size_t, std::size_t>;
    std::vector<SurfacePackRange> neighborSurfacePacks;
    neighborSurfacePacks.resize(size() * (m_maxNeighborDistance + 1));

    std::vector<const Surface*> surfacePack;
    std::map<std::vector<const Surface*>, SurfacePackRange> surfacesToPackRange;
    for (std::size_t inputGlobalBin = 0; inputGlobalBin < m_fillingGrid.size();
         ++inputGlobalBin) {
      const GridIndex indices = localBinsFromGlobalBin2D(inputGlobalBin);

      if (!isValidBin(indices)) {
        continue;
      }

      for (std::uint8_t neighborDistance = 0;
           neighborDistance <= m_maxNeighborDistance; ++neighborDistance) {
        surfacePack.clear();

        for (const std::size_t idx :
             detail::MultiAxisHelper::neighborHoodIndices(
                 indices, neighborDistance, m_axes)) {
          const std::vector<const Surface*>& binContent = m_fillingGrid.at(idx);
          std::copy(binContent.begin(), binContent.end(),
                    std::back_inserter(surfacePack));
        }

        std::ranges::sort(surfacePack);
        const auto last = std::ranges::unique(surfacePack);
        surfacePack.erase(last.begin(), last.end());

        const std::size_t outputGlobalBin =
            globalBinFromLocalBins3D(indices, neighborDistance);

        if (const auto it = surfacesToPackRange.find(surfacePack);
            it != surfacesToPackRange.end()) {
          neighborSurfacePacks[outputGlobalBin] = it->second;
        } else {
          const SurfacePackRange surfacePackRange = {
              m_surfacePacks.size(),
              m_surfacePacks.size() + surfacePack.size()};
          m_surfacePacks.insert(m_surfacePacks.end(), surfacePack.begin(),
                                surfacePack.end());
          for (const Surface* surface : surfacePack) {
            m_extentPacks.push_back(m_surfaceExtents.at(surface));
          }
          surfacesToPackRange[surfacePack] = surfacePackRange;
          neighborSurfacePacks[outputGlobalBin] = surfacePackRange;
        }
      }
    }

    m_surfacePacks.shrink_to_fit();
    m_extentPacks.shrink_to_fit();
    m_surfaceExtents.clear();

    m_neighborSurfacePacks.reserve(neighborSurfacePacks.size());
    std::ranges::transform(neighborSurfacePacks,
                           std::back_inserter(m_neighborSurfacePacks),
                           [this](const SurfacePackRange& range) {
                             return std::span<const Surface* const>(
                                 m_surfacePacks.data() + range.first,
                                 m_surfacePacks.data() + range.second);
                           });
    m_neighborExtentPacks.reserve(neighborSurfacePacks.size());
    std::ranges::transform(
        neighborSurfacePacks, std::back_inserter(m_neighborExtentPacks),
        [this](const SurfacePackRange& range) {
          return std::span<const SurfaceArray::SurfaceExtent>(
              m_extentPacks.data() + range.first,
              m_extentPacks.data() + range.second);
        });
  }

  void checkGrid(std::span<const Surface* const> surfaces) {
    const std::set<const Surface*> allSurfaces(surfaces.begin(),
                                               surfaces.end());

    std::set<const Surface*> seenSurface;
    for (std::size_t globalBin = 0; globalBin < m_fillingGrid.size();
         ++globalBin) {
      for (const Surface* surface : m_fillingGrid.at(globalBin)) {
        seenSurface.insert(surface);
      }
    }

    if (allSurfaces != seenSurface) {
      std::set<const Surface*> diff;
      std::ranges::set_difference(allSurfaces, seenSurface,
                                  std::inserter(diff, diff.begin()));

      throw std::logic_error(std::format(
          "SurfaceArray grid does not contain all surfaces provided! "
          "{} surfaces not seen",
          diff.size()));
    }
  }

  Vector2 gridToSurfaceLocal(const GridPoint& gridLocal) const {
    Vector2 surfaceLocal = {gridLocal[0], gridLocal[1]};
    if (m_cylinderRadius > 0) {
      surfaceLocal[0] *= m_cylinderRadius;
    }
    return surfaceLocal;
  }

  GridPoint surfaceToGridLocal(const Vector2& local) const {
    GridPoint gridLocal = {local[0], local[1]};
    if (m_cylinderRadius > 0) {
      gridLocal[0] /= m_cylinderRadius;
    }
    return gridLocal;
  }

  /// Locate the neighbor pack for a trajectory. Optionally reports the
  /// grid-local pierce point and the cosine of the incidence angle.
  bool lookupBin(const GeometryContext& gctx, const Vector3& position,
                 const Vector3& direction, std::size_t& globalBin,
                 GridPoint* gridLocalOut, double* cosIncidenceOut) const {
    const std::optional<Vector2> surfaceLocal = findSurfaceLocal(
        gctx, position, direction, std::numeric_limits<double>::infinity());
    if (!surfaceLocal.has_value()) {
      return false;
    }

    const GridPoint gridLocal = surfaceToGridLocal(*surfaceLocal);
    const GridIndex localBins = localBinsFromPosition2D(gridLocal);

    const Vector3 normal = m_representative->normal(gctx, *surfaceLocal);
    const double cosIncidence = normal.dot(direction);
    // using 1e-6 to avoid division by zero, the actual value does not matter as
    // long as it is small compared to the angles we want to distinguish
    const double neighborDistanceReal = std::min<double>(
        m_maxNeighborDistance,
        std::max<double>(1, 1 / (1e-6 + std::abs(cosIncidence))));
    // clamp value to range before converting to std::uint8_t to avoid overflow
    const std::uint8_t neighborDistance =
        clampValue<std::uint8_t>(neighborDistanceReal);

    globalBin = globalBinFromLocalBins3D(localBins, neighborDistance);
    if (gridLocalOut != nullptr) {
      *gridLocalOut = gridLocal;
    }
    if (cosIncidenceOut != nullptr) {
      *cosIncidenceOut = cosIncidence;
    }
    return true;
  }

  /// Move a global point onto the representative surface along its normal, so
  /// that globalToLocal accepts it.
  Vector3 projectOntoRepresentative(const GeometryContext& gctx,
                                    const Vector3& global) const {
    const Transform3 inv =
        m_representative->localToGlobalTransform(gctx).inverse(Eigen::Isometry);
    Vector3 local = inv * global;
    if (m_cylinderRadius > 0) {
      const double r = local.head<2>().norm();
      if (r > 0) {
        local.head<2>() *= m_cylinderRadius / r;
      }
    } else {
      local.z() = 0;
    }
    return m_representative->localToGlobalTransform(gctx) * local;
  }

  /// Project a global point into grid-local coordinates of the representative
  /// surface, ignoring how far off the surface the point sits.
  GridPoint projectToGridLocal(const GeometryContext& gctx,
                               const Vector3& global) const {
    const Vector3 local =
        m_representative->localToGlobalTransform(gctx).inverse(
            Eigen::Isometry) *
        global;
    if (m_cylinderRadius > 0) {
      return {std::atan2(local.y(), local.x()), local.z()};
    }
    return {local.head<2>().norm(), std::atan2(local.y(), local.x())};
  }

  /// Distance of a global point from the representative surface along its
  /// normal.
  double normalOffset(const GeometryContext& gctx,
                      const Vector3& global) const {
    const Vector3 local =
        m_representative->localToGlobalTransform(gctx).inverse(
            Eigen::Isometry) *
        global;
    if (m_cylinderRadius > 0) {
      return std::abs(local.head<2>().norm() - m_cylinderRadius);
    }
    return std::abs(local.z());
  }

  std::optional<Vector2> findSurfaceLocal(const GeometryContext& gctx,
                                          const Vector3& position,
                                          const Vector3& direction,
                                          double tolerance) const {
    const Intersection3D intersection =
        m_representative
            ->intersect(gctx, position, direction,
                        BoundaryTolerance::Infinite())
            .closest();
    if (!intersection.isValid() ||
        std::abs(intersection.pathLength()) > tolerance) {
      return std::nullopt;
    }
    const Vector2 surfaceLocal =
        m_representative
            ->globalToLocal(gctx,
                            position + intersection.pathLength() * direction,
                            direction)
            .value();
    return surfaceLocal;
  }

  std::optional<GridIndex> findLocalBin2D(const GeometryContext& gctx,
                                          const Vector3& position,
                                          const Vector3& direction,
                                          double tolerance) const {
    const std::optional<Vector2> surfaceLocal =
        findSurfaceLocal(gctx, position, direction, tolerance);
    if (!surfaceLocal.has_value()) {
      return std::nullopt;
    }
    const GridPoint gridLocal = surfaceToGridLocal(*surfaceLocal);
    return localBinsFromPosition2D(gridLocal);
  }
};

std::unique_ptr<SurfaceArray::ISurfaceGridLookup> makeSurfaceGridLookup(
    std::shared_ptr<RegularSurface> representative, double tolerance,
    std::tuple<const IAxis&, const IAxis&> axes,
    std::uint8_t maxNeighborDistance) {
  const auto& [iAxisA, iAxisB] = axes;

  return iAxisA.visit([&]<typename axis_a_t>(const axis_a_t& axisA) {
    return iAxisB.visit(
        [&]<typename axis_b_t>(const axis_b_t& axisB)
            -> std::unique_ptr<SurfaceArray::ISurfaceGridLookup> {
          return std::make_unique<SurfaceGridLookupImpl<axis_a_t, axis_b_t>>(
              std::move(representative), tolerance,
              std::tuple<axis_a_t, axis_b_t>{axisA, axisB},
              std::vector<AxisDirection>(), maxNeighborDistance);
        });
  });
}

}  // namespace

SurfaceArray::SurfaceArray(std::shared_ptr<const Surface> srf)
    : m_gridLookup(std::make_unique<SingleElementLookupImpl>(srf.get())),
      m_surfaces({std::move(srf)}) {
  m_surfacesRawPointers.push_back(m_surfaces.at(0).get());
}

SurfaceArray::SurfaceArray(const GeometryContext& gctx,
                           std::vector<std::shared_ptr<const Surface>> surfaces,
                           std::shared_ptr<RegularSurface> representative,
                           double tolerance,
                           std::tuple<const IAxis&, const IAxis&> axes,
                           std::uint8_t maxNeighborDistance) {
  m_gridLookup = makeSurfaceGridLookup(std::move(representative), tolerance,
                                       axes, maxNeighborDistance);
  m_surfaces = std::move(surfaces);
  m_surfacesRawPointers =
      m_surfaces |
      std::views::transform(
          [](const std::shared_ptr<const Surface>& sp) { return sp.get(); }) |
      Ranges::to<std::vector>;
  m_gridLookup->fill(gctx, m_surfacesRawPointers);
}

SurfaceArray::SurfaceArray(SurfaceArray&& other) noexcept = default;

SurfaceArray& SurfaceArray::operator=(SurfaceArray&& other) noexcept = default;

SurfaceArray::~SurfaceArray() = default;

std::span<const Surface* const> SurfaceArray::at(std::size_t bin) const {
  return m_gridLookup->at(bin);
}

std::span<const Surface* const> SurfaceArray::at(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  return m_gridLookup->at(gctx, position, direction);
}

std::span<const Surface* const> SurfaceArray::neighbors(
    std::array<std::size_t, 2> gridIndices,
    std::uint8_t neighborDistance) const {
  return m_gridLookup->neighbors(gridIndices, neighborDistance);
}

std::span<const Surface* const> SurfaceArray::neighbors(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  return m_gridLookup->neighbors(gctx, position, direction);
}

SurfaceArray::NeighborQuery SurfaceArray::neighborQuery(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  return m_gridLookup->neighborQuery(gctx, position, direction);
}

std::size_t SurfaceArray::size() const {
  return m_gridLookup->size();
}

Vector3 SurfaceArray::getBinCenter(std::size_t bin) const {
  return m_gridLookup->getBinCenter(bin);
}

std::vector<const IAxis*> SurfaceArray::getAxes() const {
  return m_gridLookup->getAxes();
}

bool SurfaceArray::isValidBin(std::size_t bin) const {
  return m_gridLookup->isValidBin(bin);
}

std::vector<AxisDirection> SurfaceArray::binningValues() const {
  return m_gridLookup->binningValues();
}

std::ostream& SurfaceArray::toStream(const GeometryContext& /*gctx*/,
                                     std::ostream& sl) const {
  detail::OstreamStateGuard guard{sl};
  sl << std::fixed << std::setprecision(4);
  sl << "SurfaceArray:" << std::endl;
  sl << " - no surfaces: " << m_surfaces.size() << std::endl;

  const std::vector<const IAxis*> axes = m_gridLookup->getAxes();

  for (const auto [j, axis] : enumerate(axes)) {
    const AxisBoundaryType bdt = axis->getBoundaryType();
    sl << " - axis " << (j + 1) << std::endl;
    sl << "   - boundary type: ";
    if (bdt == AxisBoundaryType::Open) {
      sl << "open";
    }
    if (bdt == AxisBoundaryType::Bound) {
      sl << "bound";
    }
    if (bdt == AxisBoundaryType::Closed) {
      sl << "closed";
    }
    sl << std::endl;
    sl << "   - type: " << (axis->isEquidistant() ? "equidistant" : "variable")
       << std::endl;
    sl << "   - n bins: " << axis->getNBins() << std::endl;
    sl << "   - bin edges: [ ";
    const std::vector<double> binEdges = axis->getBinEdges();
    for (const auto [i, binEdge] : enumerate(binEdges)) {
      if (i > 0) {
        sl << ", ";
      }
      // Do not display negative zeroes
      sl << ((std::abs(binEdge) >= 5e-4) ? binEdge : 0.0);
    }
    sl << " ]" << std::endl;
  }
  return sl;
}

const Surface* SurfaceArray::surfaceRepresentation() const {
  return m_gridLookup->surfaceRepresentation();
}

std::array<std::size_t, 2> SurfaceArray::numLocalBins() const {
  return m_gridLookup->numLocalBins();
}

std::uint8_t SurfaceArray::maxNeighborDistance() const {
  return m_gridLookup->maxNeighborDistance();
}

}  // namespace Acts

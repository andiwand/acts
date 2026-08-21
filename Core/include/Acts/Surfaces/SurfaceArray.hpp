// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/IAxis.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <span>
#include <vector>

// Forward declared friend class for private access in tests.
namespace ActsTests {
struct SurfaceArrayCreatorFixture;
}

namespace Acts {

/// TEMPORARY measurement switch for the grid pre-filter -- not for commit.
/// Set ACTS_GRID_PREFILTER=0 to fall back to testing every surface of the pack.
extern const bool s_gridPreFilterEnabled;

/// @brief Provides Surface binning in 2 dimensions
///
/// Uses @c Grid under the hood to implement the storage and lookup
/// Contains a lookup struct which talks to the @c Grid
/// and performs utility actions. This struct needs to be initialised
/// externally and passed to @c SurfaceArray on construction.
class SurfaceArray {
 public:
  /// Constructor with a single surface
  /// @param srf The one and only surface
  explicit SurfaceArray(std::shared_ptr<const Surface> srf);

  /// Constructor to create a surface grid lookup for a given representative
  /// surface, tolerance, and axes.
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surfaces The input vector of surfaces that will be accessible
  ///                 through this @ref SurfaceArray.
  /// @param representative The surface which is used as representative
  /// @param tolerance The tolerance used for intersection checks
  /// @param axes The axes used for the grid
  /// @param maxNeighborDistance Maximum next neighbor distance to be included in neighbor lookups
  SurfaceArray(const GeometryContext& gctx,
               std::vector<std::shared_ptr<const Surface>> surfaces,
               std::shared_ptr<RegularSurface> representative, double tolerance,
               std::tuple<const IAxis&, const IAxis&> axes,
               std::uint8_t maxNeighborDistance = 1);

  // non-copyable but movable due to unique_ptr member. deferred implementation
  // to source since the pimpl is not fully defined in the header.

  /// @param other the other SurfaceArray to copy from
  SurfaceArray(const SurfaceArray& other) = delete;

  /// @param other the other SurfaceArray to move from
  SurfaceArray(SurfaceArray&& other) noexcept;

  /// @param other the other SurfaceArray to copy-assign from
  /// @return reference to this SurfaceArray
  SurfaceArray& operator=(const SurfaceArray& other) = delete;

  /// @param other the other SurfaceArray to move-assign from
  /// @return reference to this SurfaceArray
  SurfaceArray& operator=(SurfaceArray&& other) noexcept;

  ~SurfaceArray();

  /// Get all surfaces in bin given by the global bin index
  /// @param bin the global bin index
  /// @return span of surface pointers of the bin at that position
  std::span<const Surface* const> at(std::size_t bin) const;

  /// Get all surfaces in bin given by position @p pos.
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position the lookup position
  /// @param direction the lookup direction
  /// @return span of surface pointers of the bin at that position
  std::span<const Surface* const> at(const GeometryContext& gctx,
                                     const Vector3& position,
                                     const Vector3& direction) const;

  /// Get all surfaces in bin given by local grid indices and neighbor
  /// distance.
  /// @param gridIndices the local grid indices
  /// @param neighborDistance the neighbor distance to include in the lookup
  /// @return span of surface pointers of the bin at that position and its neighbors
  std::span<const Surface* const> neighbors(
      std::array<std::size_t, 2> gridIndices,
      std::uint8_t neighborDistance) const;

  /// Get all surfaces in bin at @p pos and its neighbors
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position to lookup
  /// @param direction The direction to lookup
  /// @return span of surface pointers of neighbors and nominal
  std::span<const Surface* const> neighbors(const GeometryContext& gctx,
                                            const Vector3& position,
                                            const Vector3& direction) const;

  /// A neighbor lookup that also carries a cheap conservative test for whether
  /// the trajectory can cross a given surface of the pack at all.
  ///
  /// The pack of a bin is the union over everything reachable from anywhere in
  /// that bin, so it is much wider than what a single trajectory crosses. The
  /// lookup already knows where the trajectory pierces the representative
  /// surface; comparing that against each surface's precomputed extent in the
  /// same frame rejects most of the pack without a 3D intersection.
  /// Centre and half-width of a surface on both grid axes. Floats and a
  /// centre/half form rather than doubles and a min/max pair: this is streamed
  /// once per candidate in the rejection loop, so halving it to 16 bytes halves
  /// that loop's memory traffic, and the comparison needs no arithmetic. The
  /// precision is irrelevant against margins of millimetre scale.
  struct SurfaceExtent {
    /// Centre on the first grid axis
    float u0{};
    /// Half-width on the first grid axis
    float hu{};
    /// Centre on the second grid axis
    float v0{};
    /// Half-width on the second grid axis
    float hv{};
  };

  struct NeighborQuery {
    /// The surfaces of the pack
    std::span<const Surface* const> surfaces;
    /// Parallel to @c surfaces. Empty if no extents are available, in which
    /// case @c mayCross always accepts.
    std::span<const SurfaceExtent> extents;
    /// Where the trajectory pierces the representative surface, grid-local
    std::array<double, 2> gridLocal{};
    /// Per-axis slack covering how far a surface sits off the representative
    /// at this incidence angle
    std::array<double, 2> margin{};
    /// Which axis is an angle and therefore has to be compared modulo 2*pi
    std::array<bool, 2> wrap{};

    /// Fold a difference of two angles, each in (-pi, pi], back into that
    /// range. Both inputs are already in range so the difference is in
    /// (-2pi, 2pi) and a pair of comparisons is enough -- this runs once per
    /// candidate of every pack, where a divide and a round were showing up.
    /// @param d the difference
    /// @return the folded difference
    constexpr static double foldToPi(double d) noexcept {
      constexpr double pi = std::numbers::pi_v<double>;
      constexpr double twoPi = 2 * pi;
      if (d > pi) {
        return d - twoPi;
      }
      if (d < -pi) {
        return d + twoPi;
      }
      return d;
    }

    /// Conservative pre-rejection.
    /// @param i index into @c surfaces
    /// @return false only if the trajectory provably cannot cross the surface
    bool mayCross(std::size_t i) const {
      if (i >= extents.size() || !s_gridPreFilterEnabled) {
        return true;
      }
      const SurfaceExtent& e = extents[i];
      double du = gridLocal[0] - e.u0;
      if (wrap[0]) {
        du = foldToPi(du);
      }
      if (std::abs(du) > e.hu + margin[0]) {
        return false;
      }
      double dv = gridLocal[1] - e.v0;
      if (wrap[1]) {
        dv = foldToPi(dv);
      }
      return std::abs(dv) <= e.hv + margin[1];
    }
  };

  /// Neighbor lookup with the per-surface pre-rejection data
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position to lookup
  /// @param direction The direction to lookup
  /// @return the pack together with its extents and the query point
  NeighborQuery neighborQuery(const GeometryContext& gctx,
                              const Vector3& position,
                              const Vector3& direction) const;

  /// Get the size of the underlying grid structure including under/overflow
  /// bins
  /// @return the size
  std::size_t size() const;

  /// Get the center of the bin identified by global bin index @p bin
  /// @param bin the global bin index
  /// @return Center position of the bin in global coordinates
  Vector3 getBinCenter(std::size_t bin) const;

  /// Get all surfaces attached to this @c SurfaceArray
  /// @return Reference to vector of all surfaces
  /// @note This does not reflect the actual state of the grid. It only
  ///       returns what was given in the constructor, without any checks
  ///       if that is actually what's in the grid.
  const std::vector<const Surface*>& surfaces() const {
    return m_surfacesRawPointers;
  }

  /// Get vector of axes spanning the grid as @c AnyAxis
  /// @return vector of @c AnyAxis
  /// @note The axes in the vector are copies. Only use for introspection and
  ///       querying.
  std::vector<const IAxis*> getAxes() const;

  /// Checks if global bin is valid
  /// @param bin the global bin index
  /// @return bool if the bin is valid
  /// @note Valid means that the index points to a bin which is not a under
  ///       or overflow bin or out of range in any axis.
  bool isValidBin(std::size_t bin) const;

  /// The binning values described by this surface grid lookup. They are in
  /// order of the axes
  /// @return Vector of axis directions for binning
  std::vector<AxisDirection> binningValues() const;

  /// Get string representation of this @c SurfaceArray
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param sl Output stream to write to
  /// @return the output stream given as @p sl
  std::ostream& toStream(const GeometryContext& gctx, std::ostream& sl) const;

  /// Get the representative surface used for this surface array
  /// @return Surface pointer
  const Surface* surfaceRepresentation() const;

  /// Get the number of local bins in each dimension. This is used to
  /// determine the size of the grid for neighbor lookups.
  /// @return Array of number of local bins in each dimension
  std::array<std::size_t, 2> numLocalBins() const;

  /// Get the maximum neighbor distance that is supported by this lookup. This
  /// is used to determine how many neighbors to include in neighbor lookups.
  /// @return Maximum neighbor distance
  std::uint8_t maxNeighborDistance() const;

  /// Forward declaration of the internal lookup struct. The actual definition
  /// is in the source file.
  struct ISurfaceGridLookup;

 private:
  /// The actual grid lookup implementation
  std::unique_ptr<ISurfaceGridLookup> m_gridLookup;
  /// this vector makes sure we have shared ownership over the surfaces
  std::vector<std::shared_ptr<const Surface>> m_surfaces;
  /// this vector is returned, so that (expensive) copying of the shared_ptr
  /// vector does not happen by default
  std::vector<const Surface*> m_surfacesRawPointers;
};

}  // namespace Acts

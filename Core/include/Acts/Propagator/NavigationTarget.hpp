// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <variant>

namespace Acts {

/// @brief The navigation target
///
/// This struct represents a navigation target which is communicated from the
/// navigator to the stepper through the propagator.
///
/// @note This incorporates `std::optional` semantics as the next target might
///       not exist.
class NavigationTarget {
 public:
  /// Type alias for the intersection object
  using Intersection = Intersection3D;

  /// Create a surface intersection from a 3D intersection, intersection index,
  /// and a surface
  ///
  /// @param intersection is the intersection
  /// @param intersectionIndex is the intersection index
  /// @param target is the intersected target
  /// @param boundaryTolerance is the boundary tolerance used for this
  /// intersection
  constexpr NavigationTarget(
      const Intersection3D& intersection, IntersectionIndex intersectionIndex,
      const Surface& target,
      const BoundaryTolerance& boundaryTolerance) noexcept
      : m_pathLength(intersection.pathLength()),
        m_status(intersection.status()),
        m_intersectionIndex(intersectionIndex),
        m_target(&target),
        m_surfaceRepresentation(&target),
        m_boundaryTolerance(boundaryTolerance) {}

  /// Create a layer intersection from a 3D intersection, intersection index,
  /// and a layer
  ///
  /// @param intersection is the intersection
  /// @param intersectionIndex is the intersection index
  /// @param target is the intersected target
  /// @param surfaceRepresentation is the surface representation of the layer
  /// @param boundaryTolerance is the boundary tolerance used for this
  /// intersection
  constexpr NavigationTarget(
      const Intersection3D& intersection, IntersectionIndex intersectionIndex,
      const Layer& target, const Surface& surfaceRepresentation,
      const BoundaryTolerance& boundaryTolerance) noexcept
      : m_pathLength(intersection.pathLength()),
        m_status(intersection.status()),
        m_intersectionIndex(intersectionIndex),
        m_target(&target),
        m_surfaceRepresentation(&surfaceRepresentation),
        m_boundaryTolerance(boundaryTolerance) {}

  /// Create a boundary surface intersection from a 3D intersection,
  /// intersection index, and a boundary surface
  ///
  /// @param intersection is the intersection
  /// @param intersectionIndex is the intersection index
  /// @param target is the intersected target
  /// @param boundaryTolerance is the boundary tolerance used for this
  /// intersection
  constexpr NavigationTarget(
      const Intersection3D& intersection, IntersectionIndex intersectionIndex,
      const BoundarySurface& target,
      const BoundaryTolerance& boundaryTolerance) noexcept
      : m_pathLength(intersection.pathLength()),
        m_status(intersection.status()),
        m_intersectionIndex(intersectionIndex),
        m_target(&target),
        m_surfaceRepresentation(&target.surfaceRepresentation()),
        m_boundaryTolerance(boundaryTolerance) {}

  /// Create a portal intersection from a 3D intersection, intersection index,
  /// and a portal
  ///
  /// @param intersection is the intersection
  /// @param intersectionIndex is the intersection index
  /// @param target is the intersected target
  /// @param boundaryTolerance is the boundary tolerance used for this
  /// intersection
  NavigationTarget(const Intersection3D& intersection,
                   IntersectionIndex intersectionIndex, const Portal& target,
                   const BoundaryTolerance& boundaryTolerance) noexcept
      : m_pathLength(intersection.pathLength()),
        m_status(intersection.status()),
        m_intersectionIndex(intersectionIndex),
        m_target(&target),
        m_surfaceRepresentation(&target.surface()),
        m_boundaryTolerance(boundaryTolerance) {}

  /// Copy constructor
  constexpr NavigationTarget(const NavigationTarget&) noexcept = default;

  /// Move constructor
  constexpr NavigationTarget(NavigationTarget&&) noexcept = default;

  /// Copy assignment operator
  /// @return Reference to this object
  constexpr NavigationTarget& operator=(const NavigationTarget&) noexcept =
      default;

  /// Move assignment operator
  /// @return Reference to this object
  constexpr NavigationTarget& operator=(NavigationTarget&&) noexcept = default;

  /// Overwrite path length and status from an intersection
  /// @param intersection the intersection to take them from
  constexpr void setIntersection(const Intersection3D& intersection) noexcept {
    m_pathLength = intersection.pathLength();
    m_status = intersection.status();
  }

  /// Returns the intersection index
  /// @return the intersection index
  constexpr IntersectionIndex intersectionIndex() const noexcept {
    return m_intersectionIndex;
  }

  /// Mutable access to the intersection index
  /// @return the intersection index
  constexpr IntersectionIndex& intersectionIndex() noexcept {
    return m_intersectionIndex;
  }

  /// Returns the surface that has been intersected
  /// @return the surface
  constexpr const Surface& surface() const noexcept {
    return *m_surfaceRepresentation;
  }

  /// Returns the layer that has been intersected
  /// @return the layer
  constexpr const Layer& layer() const {
    return *std::get<const Layer*>(m_target);
  }

  /// Returns the boundary surface that has been intersected
  /// @return the boundary surface
  constexpr const BoundarySurface& boundarySurface() const {
    return *std::get<const BoundarySurface*>(m_target);
  }

  /// Returns the portal that has been intersected
  /// @return the portal
  constexpr const Portal& portal() const {
    return *std::get<const Portal*>(m_target);
  }

  /// Returns whether the target is a surface
  /// @return true if the target is a surface
  constexpr bool isSurfaceTarget() const noexcept {
    return std::holds_alternative<const Surface*>(m_target);
  }

  /// Returns whether the target is a layer
  /// @return true if the target is a layer
  constexpr bool isLayerTarget() const noexcept {
    return std::holds_alternative<const Layer*>(m_target);
  }

  /// Returns whether the target is a portal
  /// @return true if the target is a portal
  constexpr bool isPortalTarget() const noexcept {
    return std::holds_alternative<const BoundarySurface*>(m_target) ||
           std::holds_alternative<const Portal*>(m_target);
  }

  /// Returns the boundary tolerance used for this intersection
  /// @return the boundary tolerance
  constexpr const BoundaryTolerance& boundaryTolerance() const noexcept {
    return m_boundaryTolerance;
  }

  /// Returns whether the intersection was successful or not
  /// @return true if the intersection is valid
  constexpr bool isValid() const noexcept {
    return m_status != IntersectionStatus::unreachable;
  }

  /// Returns the path length to the interseciton
  /// @return the path length
  constexpr double pathLength() const noexcept { return m_pathLength; }

  /// Returns the status of the interseciton
  /// @return the status
  constexpr IntersectionStatus status() const noexcept { return m_status; }

  /// Returns whether this is a none target
  /// @return true if this is a none target
  constexpr bool isNone() const noexcept {
    return std::holds_alternative<std::monostate>(m_target);
  }

  /// Factory method to create a none target
  /// @return a none target
  constexpr static NavigationTarget None() noexcept {
    return NavigationTarget();
  }

  /// Comparison operator by path length
  /// @param aIntersection First navigation target
  /// @param bIntersection Second navigation target
  /// @return true if aIntersection is before bIntersection
  constexpr static bool pathLengthOrder(
      const NavigationTarget& aIntersection,
      const NavigationTarget& bIntersection) noexcept {
    return aIntersection.pathLength() < bIntersection.pathLength();
  }

  /// Comparison operator by closest distance to the reference point
  /// @param aIntersection First navigation target
  /// @param bIntersection Second navigation target
  /// @return true if aIntersection is closer than bIntersection
  constexpr static bool closestOrder(
      const NavigationTarget& aIntersection,
      const NavigationTarget& bIntersection) noexcept {
    using enum IntersectionStatus;
    if (aIntersection.status() == unreachable) {
      return false;
    }
    if (bIntersection.status() == unreachable) {
      return true;
    }
    return std::abs(aIntersection.pathLength()) <
           std::abs(bIntersection.pathLength());
  }

  /// Comparison operator by closest distance to the reference point in the
  /// @param aIntersection First navigation target
  /// @param bIntersection Second navigation target
  /// @return true if aIntersection is closer than bIntersection
  constexpr static bool closestForwardOrder(
      const NavigationTarget& aIntersection,
      const NavigationTarget& bIntersection) noexcept {
    const double a = aIntersection.pathLength();
    const double b = bIntersection.pathLength();
    return std::signbit(a) == std::signbit(b) ? std::abs(a) < std::abs(b)
                                              : a > b;
  }
  /// @brief Define the ostream operator to print the object
  /// @param ostr: Reference to the ostream
  /// @param target: Reference to the target to print
  friend std::ostream& operator<<(std::ostream& ostr,
                                  const NavigationTarget& target) {
    target.print(ostr);
    return ostr;
  }

 private:
  /// Alias for the target variant
  using TargetVariant =
      std::variant<std::monostate, const Surface*, const Layer*,
                   const BoundarySurface*, const Portal*>;

  /// Signed path length to the target. Only the path length and the status of
  /// an intersection are ever read here, so the position that Intersection3D
  /// also carries is not stored -- that is 24 of its 40 bytes and nothing on
  /// the navigation or stepping path reads it.
  double m_pathLength = std::numeric_limits<double>::infinity();
  /// The status of the intersection
  IntersectionStatus m_status = IntersectionStatus::unreachable;
  /// The intersection index
  IntersectionIndex m_intersectionIndex = 0;
  /// The target that was intersected
  TargetVariant m_target;
  /// The surface representation of the target
  const Surface* m_surfaceRepresentation = nullptr;
  /// The boundary tolerance used for this intersection
  BoundaryTolerance m_boundaryTolerance = BoundaryTolerance::None();

  /// Default constructor creating a none target
  constexpr NavigationTarget() = default;

  /// @brief print method
  /// @param ostr: Stream to which the object is printed
  void print(std::ostream& ostr) const;
};

static_assert(std::is_trivially_copy_constructible_v<NavigationTarget>);
static_assert(std::is_trivially_move_constructible_v<NavigationTarget>);
static_assert(std::is_trivially_move_assignable_v<NavigationTarget>);

}  // namespace Acts

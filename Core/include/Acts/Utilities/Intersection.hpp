// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <type_traits>

namespace Acts {

/// Status enum
enum class IntersectionStatus : int {
  unreachable = 0,
  reachable = 1,
  onSurface = 2
};

/// Ostream-operator for the IntersectionStatus enum
/// @param os Output stream
/// @param status IntersectionStatus to output
/// @return Reference to output stream
inline std::ostream& operator<<(std::ostream& os, IntersectionStatus status) {
  constexpr static std::array<const char*, 3> names = {
      {"missed/unreachable", "reachable", "onSurface"}};

  os << names[static_cast<std::size_t>(status)];
  return os;
}

namespace detail {
/// Storage for the intersection position, present only where it is read. The
/// 3D intersection is on the navigation and stepping hot path, where only the
/// path length and the status are ever used; carrying a position there made
/// the type 40 bytes and MultiIntersection 88, returned by value out of every
/// virtual Surface::intersect. The 2D intersection is a geometry helper whose
/// consumers (Fatras digitization and masking) want the point itself.
template <unsigned int DIM, bool Store>
struct IntersectionPosition {};

template <unsigned int DIM>
struct IntersectionPosition<DIM, true> {
  /// Position of the intersection
  std::array<double, DIM> storedPosition{};
};
}  // namespace detail

/// Intersection struct containing the path length and status of an
/// intersection, plus its position where that is used.
template <unsigned int DIM, bool StorePosition = (DIM == 2)>
class Intersection : private detail::IntersectionPosition<DIM, StorePosition> {
 public:
  /// Whether this intersection carries its position
  static constexpr bool hasPosition = StorePosition;

  /// Position type
  using Position = Eigen::Map<const Vector<DIM>>;

  /// Constructor with arguments
  /// @param position is the position of the intersection, stored only if @c hasPosition
  /// @param pathLength is the path length to the intersection
  /// @param status is an enum indicating the status of the intersection
  constexpr Intersection(const Vector<DIM>& position, double pathLength,
                         IntersectionStatus status) noexcept
      : m_pathLength(pathLength), m_status(status) {
    if constexpr (hasPosition) {
      std::ranges::copy(std::span<const double, DIM>{position.data(), DIM},
                        this->storedPosition.begin());
    } else {
      static_cast<void>(position);
    }
  }

  /// Constructor from a mapped position, path length, and status
  /// @param position The intersection position
  /// @param pathLength The path length to the intersection
  /// @param status The intersection status
  constexpr Intersection(const Position& position, double pathLength,
                         IntersectionStatus status) noexcept
      : m_pathLength(pathLength), m_status(status) {
    if constexpr (hasPosition) {
      std::ranges::copy(std::span<const double, DIM>{position.data(), DIM},
                        this->storedPosition.begin());
    } else {
      static_cast<void>(position);
    }
  }

  /// Returns the position of the intersection
  /// @return Position vector of the intersection point
  Position position() const noexcept
    requires(hasPosition)
  {
    return Position{this->storedPosition.data()};
  }

  /// Constructor from path length and status
  /// @param pathLength The path length to the intersection
  /// @param status The intersection status
  constexpr Intersection(double pathLength, IntersectionStatus status) noexcept
      : m_pathLength(pathLength), m_status(status) {}

  /// Copy constructor
  constexpr Intersection(const Intersection&) noexcept = default;
  /// Move constructor
  constexpr Intersection(Intersection&&) noexcept = default;
  /// Copy assignment operator
  /// @return Reference to this intersection for chaining
  constexpr Intersection& operator=(const Intersection&) noexcept = default;
  /// Move assignment operator
  /// @return Reference to this intersection for chaining
  constexpr Intersection& operator=(Intersection&&) noexcept = default;

  /// Returns whether the intersection was successful or not
  /// @return True if intersection is reachable or on surface, false if unreachable
  constexpr bool isValid() const noexcept {
    return m_status != IntersectionStatus::unreachable;
  }

  /// Returns the path length to the intersection
  /// @return Signed path length from origin to intersection point
  constexpr double pathLength() const noexcept { return m_pathLength; }

  /// Returns the intersection status enum
  /// @return Status indicating if intersection is unreachable, reachable, or on surface
  constexpr IntersectionStatus status() const noexcept { return m_status; }

  /// Static factory to create an invalid intersection
  /// @return Invalid intersection with unreachable status
  constexpr static Intersection Invalid() noexcept { return Intersection(); }

  /// Comparison function for path length order i.e. intersection closest to
  /// -inf will be first.
  /// @param aIntersection First intersection to compare
  /// @param bIntersection Second intersection to compare
  /// @return True if first intersection has smaller path length than second
  constexpr static bool pathLengthOrder(
      const Intersection& aIntersection,
      const Intersection& bIntersection) noexcept {
    auto a = aIntersection.pathLength();
    auto b = bIntersection.pathLength();
    return a < b;
  }

  /// Comparison function for closest order i.e. intersection closest to 0 will
  /// be first.
  /// @param aIntersection First intersection to compare
  /// @param bIntersection Second intersection to compare
  /// @return True if first intersection is closer to zero path length than second
  constexpr static bool closestOrder(
      const Intersection& aIntersection,
      const Intersection& bIntersection) noexcept {
    using enum IntersectionStatus;

    if ((aIntersection.status() == unreachable) &&
        (bIntersection.status() != unreachable)) {
      return false;
    }
    if ((aIntersection.status() != unreachable) &&
        (bIntersection.status() == unreachable)) {
      return true;
    }
    // both are reachable or onSurface now
    auto a = aIntersection.pathLength();
    auto b = bIntersection.pathLength();
    return std::abs(a) < std::abs(b);
  }

  /// Comparison function for closest forward order i.e. intersection closest to
  /// 0 with positive path length will be first.
  /// @param aIntersection First intersection to compare
  /// @param bIntersection Second intersection to compare
  /// @return True if first intersection is closer to zero with preference for forward direction
  constexpr static bool closestForwardOrder(
      const Intersection& aIntersection,
      const Intersection& bIntersection) noexcept {
    auto a = aIntersection.pathLength();
    auto b = bIntersection.pathLength();
    return std::signbit(a) == std::signbit(b) ? std::abs(a) < std::abs(b)
                                              : a > b;
  }

 private:
  /// Signed path length to the intersection (if valid)
  double m_pathLength = std::numeric_limits<double>::infinity();
  /// The Status of the intersection
  IntersectionStatus m_status = IntersectionStatus::unreachable;

  constexpr Intersection() noexcept = default;
};

/// Type alias for 2D intersection
using Intersection2D = Intersection<2>;
/// Type alias for 3D intersection. Carries no position -- see
/// @c detail::IntersectionPosition.
using Intersection3D = Intersection<3>;
/// An intersection that does carry its position, for the consumers that need
/// the point rather than the path length.
template <unsigned int DIM>
using PositionedIntersection = Intersection<DIM, true>;

static_assert(std::is_trivially_copy_constructible_v<Intersection2D>);
static_assert(std::is_trivially_move_constructible_v<Intersection2D>);
static_assert(std::is_trivially_move_assignable_v<Intersection2D>);

/// Index type for intersections
using IntersectionIndex = std::uint8_t;
/// Maximum number of intersections that can be stored
static constexpr IntersectionIndex s_maximumNumberOfIntersections = 2;

/// Container for up to two intersections in a given dimension.
template <unsigned int DIM>
class MultiIntersection {
 public:
  /// Intersection type for this dimension
  using IntersectionType = Intersection<DIM>;
  /// Pair of intersection and its index
  using IndexedIntersection = std::pair<IntersectionType, IntersectionIndex>;

  /// Container type for storing intersections
  using Container =
      std::array<IntersectionType, s_maximumNumberOfIntersections>;

  /// Size type for indexing
  using size_type = IntersectionIndex;

  /// Construct from single intersection
  /// @param intersection The intersection
  constexpr explicit MultiIntersection(
      const IntersectionType& intersection) noexcept
      : m_intersections{intersection, IntersectionType::Invalid()}, m_size{1} {}
  /// Construct from two intersections
  /// @param intersection1 The first intersection
  /// @param intersection2 The second intersection
  constexpr MultiIntersection(const IntersectionType& intersection1,
                              const IntersectionType& intersection2) noexcept
      : m_intersections{intersection1, intersection2}, m_size{2} {}

  /// Copy constructor
  constexpr MultiIntersection(const MultiIntersection&) noexcept = default;
  /// Move constructor
  constexpr MultiIntersection(MultiIntersection&&) noexcept = default;
  /// Copy assignment operator
  /// @return Reference to this object
  constexpr MultiIntersection& operator=(const MultiIntersection&) noexcept =
      default;
  /// Move assignment operator
  /// @return Reference to this object
  constexpr MultiIntersection& operator=(MultiIntersection&&) noexcept =
      default;

  /// Access intersection by index
  /// @param index The index of the intersection
  /// @return Reference to the intersection
  constexpr const IntersectionType& operator[](IntersectionIndex index) const {
    return m_intersections[index];
  }

  /// Access intersection at index with bounds checking
  /// @param index The index of the intersection
  /// @return Reference to the intersection
  constexpr const IntersectionType& at(IntersectionIndex index) const {
    return m_intersections.at(index);
  }

  /// Get the number of intersections
  /// @return The number of intersections
  constexpr IntersectionIndex size() const noexcept { return m_size; }

  /// Get begin iterator
  /// @return Iterator to the beginning
  constexpr auto begin() const noexcept {
    return std::span(m_intersections.data(), m_size).begin();
  }
  /// Get end iterator
  /// @return Iterator to the end
  constexpr auto end() const noexcept {
    return std::span(m_intersections.data(), m_size).end();
  }

  /// Get closest intersection
  /// @return The closest intersection
  constexpr IntersectionType closest() const noexcept {
    return closestWithIndex().first;
  }
  /// Get closest intersection with its index
  /// @return Pair of intersection and its index
  constexpr IndexedIntersection closestWithIndex() const noexcept {
    // Only the first m_size slots are filled. Planar surfaces fill one, so
    // ranging over the whole array costs an extra closestOrder comparison on
    // every intersection the navigator does.
    const auto begin = m_intersections.begin();
    auto min =
        std::min_element(begin, begin + m_size, IntersectionType::closestOrder);
    return {*min, static_cast<IntersectionIndex>(std::distance(begin, min))};
  }

  /// Get closest forward intersection
  /// @return The closest forward intersection
  constexpr IntersectionType closestForward() const noexcept {
    return closestForwardWithIndex().first;
  }
  /// Get closest forward intersection with its index
  /// @return Pair of intersection and its index
  constexpr IndexedIntersection closestForwardWithIndex() const noexcept {
    auto min = std::ranges::min_element(m_intersections,
                                        IntersectionType::closestForwardOrder);
    return {*min, static_cast<IntersectionIndex>(
                      std::distance(m_intersections.begin(), min))};
  }

 private:
  Container m_intersections{};
  IntersectionIndex m_size{};
};

/// Container for up to two 2D intersections
using MultiIntersection2D = MultiIntersection<2>;
/// Container for up to two 3D intersections
using MultiIntersection3D = MultiIntersection<3>;

static_assert(std::is_trivially_copy_constructible_v<MultiIntersection2D>);
static_assert(std::is_trivially_move_constructible_v<MultiIntersection2D>);
static_assert(std::is_trivially_move_assignable_v<MultiIntersection2D>);

namespace detail {

/// Verbose-logging companion of checkPathLength(): prints why a path length
/// is (not) within the limits. Split out of line so the inline fast path
/// below stays free of the log-message formatting.
///
/// @param pathLength The path length of the intersection
/// @param nearLimit The minimum path length for an intersection to be considered
/// @param farLimit The maximum path length for an intersection to be considered
/// @param logger The logger to print to (at VERBOSE level)
void printCheckPathLength(double pathLength, double nearLimit, double farLimit,
                          const Logger& logger);

/// This function checks if an intersection path length is valid for the
/// specified near-limit and far-limit
///
/// This is called per candidate on the navigation hot paths, so the two
/// comparisons are inline and the (rarely enabled) verbose logging is
/// delegated to an out-of-line helper.
///
/// @param pathLength The path length of the intersection
/// @param nearLimit The minimum path length for an intersection to be considered
/// @param farLimit The maximum path length for an intersection to be considered
/// @param logger A optionally supplied logger which prints out a lot of infos
///               at VERBOSE level
inline bool checkPathLength(double pathLength, double nearLimit,
                            double farLimit) {
  // TODO why?
  const double tolerance = s_onSurfaceTolerance;
  return pathLength > nearLimit && pathLength < farLimit + tolerance;
}

/// Logging overload of the above. Kept separate so the common no-logger call
/// does not pay for materialising the dummy logger and gating on it.
///
/// @param pathLength The path length of the intersection
/// @param nearLimit The minimum path length for an intersection to be considered
/// @param farLimit The maximum path length for an intersection to be considered
/// @param logger A optionally supplied logger which prints out a lot of infos
///               at VERBOSE level
/// @return true if the path length is within the limits
inline bool checkPathLength(double pathLength, double nearLimit,
                            double farLimit, const Logger& logger) {
  if (logger.doPrint(Logging::VERBOSE)) [[unlikely]] {
    printCheckPathLength(pathLength, nearLimit, farLimit, logger);
  }
  return checkPathLength(pathLength, nearLimit, farLimit);
}

}  // namespace detail

}  // namespace Acts

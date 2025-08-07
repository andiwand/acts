// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <map>
#include <vector>

namespace Acts {

/// @brief A description of a triplet candidate.
/// @tparam external_space_point_t  The external spacepoint type.
template <typename external_space_point_t>
struct TripletCandidate {
  /// @brief Default Constructor
  TripletCandidate() = default;

  /// @brief constructor
  /// @param b The bottom space point
  /// @param m The middle space point
  /// @param t The top space point
  /// @param w The quality of the candidate
  /// @param z The z coordinate of the origin
  /// @param q Whether the candidate is high or low quality
  TripletCandidate(external_space_point_t& b, external_space_point_t& m,
                   external_space_point_t& t, float w, float z, bool q)
      : bottom(&b), middle(&m), top(&t), weight(w), zOrigin(z), isQuality(q) {}

  external_space_point_t* bottom{nullptr};
  external_space_point_t* middle{nullptr};
  external_space_point_t* top{nullptr};
  float weight{0.};
  float zOrigin{0.};
  bool isQuality{false};
};

/// @class CandidatesForMiddleSp
/// The CandidatesForMiddleSp collects the triplet candidates given a
/// fixed middle spacepoint. It internally stores the triplet candidates
/// keeping only those with the higher quality.
///
/// @tparam external_space_point_t The external spacepoint type.

template <typename external_space_point_t>
concept SatisfyCandidateConcept = requires(external_space_point_t spacePoint) {
  { spacePoint.x() } -> std::convertible_to<float>;
  { spacePoint.y() } -> std::convertible_to<float>;
  { spacePoint.z() } -> std::convertible_to<float>;
};

template <SatisfyCandidateConcept external_space_point_t>
class CandidatesForMiddleSp {
 public:
  using value_type = TripletCandidate<external_space_point_t>;

  /// @brief Setting maximum number of candidates to keep
  /// @param nLow Maximum number of candidates in the low-quality collection
  /// @param nHigh Maximum number of candidates in the high-quality collection
  void setMaxElements(std::size_t nLow, std::size_t nHigh) {
    m_maxSizeLow = nLow;
    m_maxSizeHigh = nHigh;
  }

  /// @brief Retrieve the triplet candidates, the resulting vector is already sorted,
  /// elements with higher quality first
  /// @returns Vector of triplet candidates
  std::vector<value_type> storage() {
    std::vector<value_type> result;
    result.reserve(m_storage.size());
    for (const auto& [weight, index] : m_weightIndexHigh) {
      result.push_back(m_storage[index]);
    }
    for (const auto& [weight, index] : m_weightIndexLow) {
      result.push_back(m_storage[index]);
    }
    clear();
    return result;
  }

  /// @brief Adding a new triplet candidate to the collection, should it satisfy the
  /// selection criteria
  /// @param spB Bottom space point
  /// @param spM Medium space point
  /// @param spT Top space point
  /// @param weight The quality of the triplet candidate
  /// @param zOrigin The z-coordinate of the origin
  /// @param isQuality Whether the triplet candidate is high or low quality
  /// @returns whether the triplet candidate has been added or not to the collection
  bool push(external_space_point_t& spB, external_space_point_t& spM,
            external_space_point_t& spT, float weight, float zOrigin,
            bool isQuality) {
    if (isQuality) {
      return push(spB, spM, spT, weight, zOrigin, isQuality, m_weightIndexHigh);
    }
    return push(spB, spM, spT, weight, zOrigin, isQuality, m_weightIndexLow);
  }

  /// @brief Clear the internal storage
  void clear() {
    m_storage.clear();
    m_weightIndexHigh.clear();
    m_weightIndexLow.clear();
  }

  /// @brief Retrieve the number of Low quality candidates
  /// @returns The number of Low quality candidates
  std::size_t nLowQualityCandidates() const { return m_weightIndexLow.size(); }

  /// @brief Retrieve the number of High quality candidates
  /// @returns The number of High quality candidates
  std::size_t nHighQualityCandidates() const {
    return m_weightIndexHigh.size();
  }

  const std::multimap<float, std::uint32_t, std::greater<float>>&
  weightIndexHigh() const {
    return m_weightIndexHigh;
  }

  const std::multimap<float, std::uint32_t, std::greater<float>>&
  weightIndexLow() const {
    return m_weightIndexLow;
  }

 private:
  bool push(
      external_space_point_t& spB, external_space_point_t& spM,
      external_space_point_t& spT, float weight, float zOrigin, bool isQuality,
      std::multimap<float, std::uint32_t, std::greater<float>>& weightIndex) {
    if (weightIndex.size() < m_maxSizeLow) {
      m_storage.emplace_back(spB, spM, spT, weight, zOrigin, isQuality);
      weightIndex.emplace(weight, m_storage.size() - 1);
      return true;
    }
    auto it = weightIndex.begin();
    if (weight < it->first) {
      return false;
    }
    m_storage[it->second] =
        value_type(spB, spM, spT, weight, zOrigin, isQuality);
    weightIndex.erase(it);
    weightIndex.emplace(weight, m_storage.size() - 1);
    return true;
  }

  // sizes
  // m_maxSize* is the maximum size of the indices collections. These values
  // are set by the user once
  std::size_t m_maxSizeHigh{0};
  std::size_t m_maxSizeLow{0};

  // storage contains the collection of the candidates
  std::vector<value_type> m_storage;

  std::multimap<float, std::uint32_t, std::greater<float>> m_weightIndexHigh;
  std::multimap<float, std::uint32_t, std::greater<float>> m_weightIndexLow;
};

}  // namespace Acts

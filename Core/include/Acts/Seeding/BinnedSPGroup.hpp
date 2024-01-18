// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"
#include "Acts/Utilities/GridIterator.hpp"
#include "Acts/Utilities/Holders.hpp"

#include <memory>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts {
template <typename external_spacepoint_t>
class BinnedSPGroup;

/// @c BinnedSPGroupIterator Allows to iterate over all groups of bins
/// a provided GridBinFinder can generate for each bin of a provided SPGrid

/// SpacePointGrid is a very specific structure.
/// We know it is 2D and what it contains
/// No need to be too general with this class
template <typename external_spacepoint_t>
class BinnedSPGroupIterator {
 private:
  enum INDEX : int { PHI = 0, Z = 1 };

 public:
  // Never take ownerships
  BinnedSPGroupIterator(BinnedSPGroup<external_spacepoint_t>&& group,
                        std::array<std::size_t, 2> index,
                        std::array<std::vector<std::size_t>, 2> navigation) =
      delete;
  BinnedSPGroupIterator(BinnedSPGroup<external_spacepoint_t>& group,
                        std::array<std::size_t, 2> index,
                        std::array<std::vector<std::size_t>, 2> navigation);

  BinnedSPGroupIterator(const BinnedSPGroupIterator&) = delete;
  BinnedSPGroupIterator& operator=(const BinnedSPGroupIterator&) = delete;

  BinnedSPGroupIterator(BinnedSPGroupIterator&&) noexcept = default;
  BinnedSPGroupIterator& operator=(BinnedSPGroupIterator&&) noexcept = default;

  ~BinnedSPGroupIterator() = default;

  BinnedSPGroupIterator& operator++();

  bool operator==(const BinnedSPGroupIterator& other) const;
  bool operator!=(const BinnedSPGroupIterator& other) const;

  std::tuple<boost::container::small_vector<std::size_t, 9>, std::size_t,
             boost::container::small_vector<std::size_t, 9>>
  operator*() const;

 private:
  void findNotEmptyBin();

 private:
  // The group, it contains the grid and the bin finders
  Acts::detail::RefHolder<BinnedSPGroup<external_spacepoint_t>> m_group;
  // Current grid iterator
  typename Acts::SpacePointGrid<external_spacepoint_t>::local_iterator_t
      m_gridItr;
  // End iterator
  typename Acts::SpacePointGrid<external_spacepoint_t>::local_iterator_t
      m_gridItrEnd;
};

/// @c BinnedSPGroup Provides access to begin and end BinnedSPGroupIterator
/// for given GridBinFinders and SpacePointGrid.
/// Fulfills the range_expression interface.
template <typename external_spacepoint_t>
class BinnedSPGroup {
#ifndef DOXYGEN
  friend BinnedSPGroupIterator<external_spacepoint_t>;
#endif

  enum INDEX : int { PHI = 0, Z = 1 };

 public:
  BinnedSPGroup() = delete;

  template <typename spacepoint_iterator_t, typename callable_t>
  BinnedSPGroup(spacepoint_iterator_t spBegin, spacepoint_iterator_t spEnd,
                callable_t&& toGlobal,
                std::shared_ptr<const Acts::GridBinFinder<2ul>> botBinFinder,
                std::shared_ptr<const Acts::GridBinFinder<2ul>> tBinFinder,
                std::unique_ptr<SpacePointGrid<external_spacepoint_t>> grid,
                Acts::Extent& rRangeSPExtent,
                const SeedFinderConfig<external_spacepoint_t>& _config,
                const SeedFinderOptions& _options);

  BinnedSPGroup(const BinnedSPGroup&) = delete;
  BinnedSPGroup& operator=(const BinnedSPGroup&) = delete;

  BinnedSPGroup(BinnedSPGroup&&) noexcept = default;
  BinnedSPGroup& operator=(BinnedSPGroup&&) noexcept = default;

  ~BinnedSPGroup() = default;

  std::size_t size() const;

  BinnedSPGroupIterator<external_spacepoint_t> begin();
  BinnedSPGroupIterator<external_spacepoint_t> end();

  Acts::SpacePointGrid<external_spacepoint_t>& grid() {
    return *m_grid.get();
  }

 private:
  // grid with ownership of all InternalSpacePoint
  std::unique_ptr<Acts::SpacePointGrid<external_spacepoint_t>> m_grid{nullptr};

  // GridBinFinder must return std::vector<Acts::Seeding::Bin> with content of
  // each bin sorted in r (ascending)
  std::shared_ptr<const Acts::GridBinFinder<2ul>> m_topBinFinder{nullptr};
  std::shared_ptr<const Acts::GridBinFinder<2ul>> m_bottomBinFinder{nullptr};

  // Order of z bins to loop over when searching for SPs
  std::array<std::vector<std::size_t>, 2> m_bins{};
};

}  // namespace Acts
#include "Acts/Seeding/BinnedSPGroup.ipp"

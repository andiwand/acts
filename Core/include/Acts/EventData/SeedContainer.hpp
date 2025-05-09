// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer.hpp"

#include <limits>

namespace Acts {

using SeedIndex = std::size_t;

class SeedContainer;

template <bool read_only>
class SeedProxy {
 public:
  /// Indicates whether this spacepoint proxy is read-only or if it can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  using IndexType = SeedIndex;

  using ContainerType = const_if_t<ReadOnly, SeedContainer>;

  SeedProxy(ContainerType &container, IndexType index)
      : m_container{&container}, m_index{index} {}

  SeedProxy(const SeedProxy &other) = default;

  explicit SeedProxy(const SeedProxy<false> &other)
    requires(ReadOnly)
      : m_container(&other.container()), m_index(other.index()) {}

  SeedContainer &container() { return *m_container; }

  const SeedContainer &container() const { return *m_container; }
  IndexType index() const { return m_index; }

  std::size_t size() const { return m_container->seedSize(m_index); }
  bool empty() const { return size() == 0; }

  float &quality()
    requires(!ReadOnly)
  {
    return m_container->quality(m_index);
  }
  float &vertexZ()
    requires(!ReadOnly)
  {
    return m_container->vertexZ(m_index);
  }

  class SpacePointIterator {
   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = ConstSpacePointProxy;
    using difference_type = std::ptrdiff_t;

    SpacePointIterator() = default;
    SpacePointIterator(const SpacePointContainer &spacePointContainer,
                       SpacePointIndex index)
        : m_spacePointContainer{&spacePointContainer}, m_index{index} {}

    SpacePointIterator &operator++() {
      ++m_index;
      return *this;
    }

    SpacePointIterator operator++(int) {
      SpacePointIterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool operator==(const SpacePointIterator &other) const {
      return m_index == other.m_index;
    }
    bool operator!=(const SpacePointIterator &other) const {
      return !(*this == other);
    }

    value_type operator*() const { return m_spacePointContainer->at(m_index); }

   private:
    const SpacePointContainer *m_spacePointContainer{nullptr};
    SpacePointIndex m_index{0};
  };

  class SpacePointRange {
   public:
    SpacePointRange(const SpacePointContainer &spacePointContainer,
                    SpacePointIndex begin, SpacePointIndex end)
        : m_spacePointContainer{&spacePointContainer},
          m_begin{begin},
          m_end{end} {}

    std::size_t size() const { return m_end - m_begin; }
    bool empty() const { return size() == 0; }

    ConstSpacePointProxy operator[](std::size_t index) const {
      return m_spacePointContainer->at(m_begin + index);
    }

    SpacePointIterator begin() const {
      return SpacePointIterator(*m_spacePointContainer, m_begin);
    }
    SpacePointIterator end() const {
      return SpacePointIterator(*m_spacePointContainer, m_end);
    }

   private:
    const SpacePointContainer *m_spacePointContainer{nullptr};
    SpacePointIndex m_begin;
    SpacePointIndex m_end;
  };

  SpacePointRange spacePoints(
      const SpacePointContainer &spacePointContainer) const {
    return SpacePointRange(spacePointContainer,
                           m_container->spacePointOffset(m_index),
                           m_container->spacePointOffset(m_index) +
                               m_container->seedSize(m_index));
  }

  float quality() const { return m_container->quality(m_index); }
  float vertexZ() const { return m_container->vertexZ(m_index); }

 private:
  ContainerType *m_container{nullptr};
  IndexType m_index{0};
};

using MutableSeedProxy = SeedProxy<false>;
using ConstSeedProxy = SeedProxy<true>;

class SeedContainer {
 public:
  using IndexType = SeedIndex;
  using MutableProxyType = MutableSeedProxy;
  using ConstProxyType = ConstSeedProxy;

  std::size_t size() const { return m_entries.size(); }
  bool empty() const { return size() == 0; }

  void reserve(std::size_t size) {
    m_entries.reserve(size);
    m_spacePoints.reserve(size * 3);
  }
  void clear() {
    m_entries.clear();
    m_spacePoints.clear();
  }

  IndexType addSeed(SpacePointIndex bottom, SpacePointIndex middle,
                    SpacePointIndex top) {
    m_entries.emplace_back(3, m_spacePoints.size(),
                           -std::numeric_limits<float>::infinity(), 0.f);
    m_spacePoints.push_back(bottom);
    m_spacePoints.push_back(middle);
    m_spacePoints.push_back(top);
    return m_entries.size() - 1;
  }

  MutableSeedProxy makeSeed(SpacePointIndex bottom, SpacePointIndex middle,
                            SpacePointIndex top) {
    return at(addSeed(bottom, middle, top));
  }

  MutableProxyType at(IndexType index) {
    return MutableProxyType(*this, index);
  }

  float &quality(IndexType index) { return m_entries[index].quality; }
  float &vertexZ(IndexType index) { return m_entries[index].vertexZ; }

  ConstProxyType at(IndexType index) const {
    return ConstProxyType(*this, index);
  }

  std::size_t seedSize(IndexType index) const {
    return m_entries[index].seedSize;
  }
  std::size_t spacePointOffset(IndexType index) const {
    return m_entries[index].spacePointOffset;
  }
  float quality(IndexType index) const { return m_entries[index].quality; }
  float vertexZ(IndexType index) const { return m_entries[index].vertexZ; }

  template <bool read_only>
  class SeedIterator {
   public:
    static constexpr bool ReadOnly = read_only;

    using ContainerType = const_if_t<ReadOnly, SeedContainer>;

    using iterator_category = std::forward_iterator_tag;
    using value_type = SeedProxy<ReadOnly>;
    using difference_type = std::ptrdiff_t;

    SeedIterator() = default;
    SeedIterator(ContainerType &container, IndexType index)
        : m_container(&container), m_index(index) {}

    SeedIterator &operator++() {
      ++m_index;
      return *this;
    }
    SeedIterator operator++(int) {
      SeedIterator tmp(*this);
      ++(*this);
      return tmp;
    }

    bool operator==(const SeedIterator &other) const {
      return m_index == other.m_index && m_container == other.m_container;
    }
    bool operator!=(const SeedIterator &other) const {
      return !(*this == other);
    }

    value_type operator*() const { return value_type(*m_container, m_index); }

   private:
    ContainerType *m_container{};
    IndexType m_index{};
  };
  using iterator = SeedIterator<false>;
  using const_iterator = SeedIterator<true>;

  iterator begin() { return iterator(*this, 0); }
  iterator end() { return iterator(*this, size()); }

  const_iterator begin() const { return const_iterator(*this, 0); }
  const_iterator end() const { return const_iterator(*this, size()); }

 private:
  struct Entry {
    std::size_t seedSize{};
    std::size_t spacePointOffset{};
    float quality{-std::numeric_limits<float>::infinity()};
    float vertexZ{};
  };

  std::vector<Entry> m_entries{};
  std::vector<SpacePointIndex> m_spacePoints{};
};

}  // namespace Acts

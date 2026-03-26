// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/EventData/VertexContainer.hpp"

namespace {

template <typename Tuple>
using tuple_indices =
    std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<Tuple>>>;

}

namespace ActsFatras {

static_assert(std::random_access_iterator<VertexContainer::iterator>);
static_assert(std::random_access_iterator<VertexContainer::const_iterator>);

VertexContainer::VertexContainer() noexcept = default;

VertexContainer::VertexContainer(const VertexContainer &other) noexcept
    : m_size(other.m_size), m_particleIndices(other.m_particleIndices) {
  copyColumns(other);
}

VertexContainer::VertexContainer(VertexContainer &&other) noexcept
    : m_size(other.m_size),
      m_particleIndices(std::move(other.m_particleIndices)) {
  moveColumns(other);

  other.m_size = 0;
}

VertexContainer &VertexContainer::operator=(
    const VertexContainer &other) noexcept {
  if (this == &other) {
    return *this;
  }

  copyColumns(other);
  m_particleIndices = other.m_particleIndices;
  m_size = other.m_size;

  return *this;
}

VertexContainer &VertexContainer::operator=(VertexContainer &&other) noexcept {
  if (this == &other) {
    return *this;
  }

  moveColumns(other);
  m_particleIndices = std::move(other.m_particleIndices);
  m_size = other.m_size;

  other.m_size = 0;

  return *this;
}

void VertexContainer::copyColumns(const VertexContainer &other) {
  knownColumns() = other.knownColumns();
}

void VertexContainer::moveColumns(VertexContainer &other) noexcept {
  knownColumns() = std::move(other).knownColumns();
}

void VertexContainer::reserve(std::uint32_t size,
                              float averageParticleIndices) noexcept {
  m_particleIndices.reserve(
      static_cast<std::size_t>(averageParticleIndices * size));

  const auto reserveKnownColumns = [&]<typename T>(std::vector<T> &column) {
    column.reserve(size);
  };

  [&]<std::size_t... Is>(std::index_sequence<Is...>) {
    ((reserveKnownColumns(std::get<Is>(knownColumns()))), ...);
  }(tuple_indices<decltype(knownColumns())>{});
}

void VertexContainer::clear() noexcept {
  m_size = 0;
  m_particleIndices.clear();
}

MutableVertexProxy VertexContainer::createVertex() noexcept {
  ++m_size;

  return MutableProxy(*this, size() - 1);
}

}  // namespace ActsFatras

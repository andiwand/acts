// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Utilities/detail/ContainerIterator.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <cstdint>
#include <span>

namespace Acts {
class Surface;
}

namespace ActsFatras {

using ParticleIndex = std::uint32_t;
using ParticleIndexSubset = std::span<const ParticleIndex>;

class ParticleContainer;

using VertexIndex = std::uint32_t;
using VertexIndexSubset = std::span<const VertexIndex>;

template <bool>
class VertexProxy;
using MutableVertexProxy = VertexProxy<false>;
using ConstVertexProxy = VertexProxy<true>;

class VertexContainer final {
 public:
  /// Type alias for vertex index in container
  using Index = VertexIndex;
  /// Type alias for subset of vertex indices
  using IndexSubset = VertexIndexSubset;
  /// Type alias for mutable vertex proxy
  using MutableProxy = MutableVertexProxy;
  /// Type alias for const vertex proxy
  using ConstProxy = ConstVertexProxy;

  VertexContainer() noexcept;

  /// Constructs a copy of the given vertex container.
  /// @param other The vertex container to copy.
  VertexContainer(const VertexContainer &other) noexcept;

  /// Move constructs a vertex container.
  /// @param other The vertex container to move.
  VertexContainer(VertexContainer &&other) noexcept;

  /// Detructs the vertex container.
  ~VertexContainer() noexcept = default;

  /// Assignment operator for copying a vertex container.
  /// @param other The vertex container to copy.
  /// @return A reference to this vertex container.
  VertexContainer &operator=(const VertexContainer &other) noexcept;

  /// Move assignment operator for a vertex container.
  /// @param other The vertex container to move.
  /// @return A reference to this vertex container.
  VertexContainer &operator=(VertexContainer &&other) noexcept;

  /// Assigns a particle container to this vertex container.
  /// @param particleContainer The particle container to assign.
  /// @throws std::logic_error if a particle container has already been assigned.
  void assignParticleContainer(const ParticleContainer &particleContainer);

  /// Checks if a particle container has been assigned to this vertex container.
  /// @return True if a particle container has been assigned.
  bool hasParticleContainer() const noexcept;

  /// Returns a const reference to the assigned particle container.
  /// @return A const reference to the assigned particle container.
  /// @throws std::logic_error if no particle container has been assigned.
  const ParticleContainer &particleContainer() const;

  /// Returns the number of vertices in the container.
  /// @return The number of vertices in the container.
  [[nodiscard]] std::uint32_t size() const noexcept { return m_size; }
  /// Checks if the container is empty.
  /// @return True if the container is empty, false otherwise.
  [[nodiscard]] bool empty() const noexcept { return size() == 0; }

  /// Reserves space for the given number of vertices.
  /// @param size The number of vertices to reserve space for.
  /// @param averageParticleIndices The average number of particle indices per vertex.
  void reserve(std::uint32_t size, float averageParticleIndices = 1) noexcept;

  /// Clears the container, removing all vertices and columns.
  void clear() noexcept;

  /// Creates a new vertex at the end of the container.
  /// @return A mutable proxy to the newly created vertex.
  MutableProxy createVertex() noexcept;

  /// Returns a mutable proxy to the vertex at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the vertex to access.
  /// @return A mutable proxy to the vertex at the given index.
  /// @throws std::out_of_range if the index is out of range.
  MutableProxy at(Index index);
  /// Returns a const proxy to the vertex at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the vertex to access.
  /// @return A const proxy to the vertex at the given index.
  /// @throws std::out_of_range if the index is out of range.
  ConstProxy at(Index index) const;

  /// Returns a mutable proxy to the vertex at the given index.
  /// @param index The index of the vertex to access.
  /// @return A mutable proxy to the vertex at the given index.
  MutableProxy operator[](Index index) noexcept;
  /// Returns a const proxy to the vertex at the given index.
  /// @param index The index of the vertex to access.
  /// @return A const proxy to the vertex at the given index.
  ConstProxy operator[](Index index) const noexcept;

  /// Type alias for template iterator over particles in container
  template <bool read_only>
  using Iterator = Acts::detail::ContainerIterator<
      VertexContainer, std::conditional_t<read_only, ConstProxy, MutableProxy>,
      Index, read_only>;

  /// Type alias for mutable iterator over particles
  using iterator = Iterator<false>;
  /// Type alias for const iterator over particles
  using const_iterator = Iterator<true>;

  /// @brief Returns mutable iterator to the beginning of the container
  /// @return Mutable iterator pointing to the first particle
  iterator begin() noexcept { return iterator(*this, 0); }
  /// @brief Returns mutable iterator to the end of the container
  /// @return Mutable iterator pointing past the last particle
  iterator end() noexcept { return iterator(*this, size()); }

  /// @brief Returns const iterator to the beginning of the container
  /// @return Const iterator pointing to the first particle
  const_iterator begin() const noexcept { return const_iterator(*this, 0); }
  /// @brief Returns const iterator to the end of the container
  /// @return Const iterator pointing past the last particle
  const_iterator end() const noexcept { return const_iterator(*this, size()); }

 private:
  std::uint32_t m_size{0};

  std::vector<ParticleIndex> m_particleIndices;

  std::vector<std::uint32_t> m_incomingParticleIndicesOffset;
  std::vector<std::uint8_t> m_incomingParticleIndicesCount;
  std::vector<std::uint32_t> m_outgoingParticleIndicesOffset;
  std::vector<std::uint8_t> m_outgoingParticleIndicesCount;
  std::vector<Acts::Vector4> m_positionColumn;
  std::vector<const Acts::Surface *> m_surfacePointers;
  std::vector<ProcessType> m_processType;

  const ParticleContainer *m_particleContainer{nullptr};

  template <typename Self>
  static auto knownColumns(Self &&self) noexcept {
    return std::tie(self.m_incomingParticleIndicesOffset,
                    self.m_incomingParticleIndicesCount,
                    self.m_outgoingParticleIndicesOffset,
                    self.m_outgoingParticleIndicesCount, self.m_positionColumn,
                    self.m_surfacePointers, self.m_processType);
  }
  auto knownColumns() & noexcept { return knownColumns(*this); }
  auto knownColumns() const & noexcept { return knownColumns(*this); }
  auto knownColumns() && noexcept { return knownColumns(*this); }

  void copyColumns(const VertexContainer &other);
  void moveColumns(VertexContainer &other) noexcept;
};

/// A proxy class for accessing individual vertices.
template <bool read_only>
class VertexProxy final {
 public:
  /// Indicates whether this vertex proxy is read-only or data can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  /// Type alias for vertex index type
  using Index = VertexIndex;

  /// Type alias for container type (const if read-only)
  using Container = Acts::const_if_t<ReadOnly, VertexContainer>;

  /// Constructs a vertex proxy for the given container and index.
  /// @param container The container holding the vertex.
  /// @param index The index of the vertex in the container.
  VertexProxy(Container &container, Index index) noexcept
      : m_container(&container), m_index(index) {}

  /// Copy construct a vertex proxy.
  /// @param other The vertex proxy to copy.
  VertexProxy(const VertexProxy &other) noexcept = default;

  /// Copy construct a mutable vertex proxy.
  /// @param other The mutable vertex proxy to copy.
  explicit VertexProxy(const VertexProxy<false> &other) noexcept
    requires ReadOnly
      : m_container(&other.container()), m_index(other.index()) {}

  /// Copy assign a vertex proxy.
  /// @param other The vertex proxy to copy.
  /// @return Reference to this vertex proxy after assignment.
  VertexProxy &operator=(const VertexProxy &other) noexcept = default;

  /// Copy assign a mutable vertex proxy.
  /// @param other The mutable vertex proxy to copy.
  /// @return Reference to this vertex proxy after assignment.
  VertexProxy &operator=(const VertexProxy<false> &other) noexcept
    requires ReadOnly
  {
    m_container = &other.container();
    m_index = other.index();
    return *this;
  }

  /// Move assign a vertex proxy.
  /// @param other The vertex proxy to move.
  /// @return Reference to this vertex proxy after assignment.
  VertexProxy &operator=(VertexProxy &&other) noexcept = default;

  /// Move assign a mutable vertex proxy.
  /// @param other The mutable vertex proxy to move.
  /// @return Reference to this vertex proxy after assignment.
  VertexProxy &operator=(VertexProxy<false> &&other) noexcept
    requires ReadOnly
  {
    m_container = &other.container();
    m_index = other.index();
    return *this;
  }

  /// Returns a const proxy of the vertex.
  /// @return A const proxy of the vertex.
  VertexProxy<true> asConst() const noexcept
    requires(!ReadOnly)
  {
    return {*m_container, m_index};
  }

  /// Gets the container holding the vertex.
  /// @return A reference to the container holding the vertex.
  VertexContainer &container() const noexcept
    requires(!ReadOnly)
  {
    return *m_container;
  }
  /// Gets the container holding the vertex.
  /// @return A const reference to the container holding the vertex.
  const VertexContainer &container() const noexcept { return *m_container; }
  /// Gets the index of the vertex in the container.
  /// @return The index of the vertex in the container.
  Index index() const noexcept { return m_index; }

  void assignIncomingParticleIndices(
      std::span<const ParticleIndex> particleIndices)
    requires(!ReadOnly)
  {
    assignParticleIndices(particleIndices,
                          container().m_incomingParticleIndicesOffset,
                          container().m_incomingParticleIndicesCount);
  }

  void assignOutgoingParticleIndices(
      std::span<const ParticleIndex> particleIndices)
    requires(!ReadOnly)
  {
    assignParticleIndices(particleIndices,
                          container().m_outgoingParticleIndicesOffset,
                          container().m_outgoingParticleIndicesCount);
  }

 private:
  Container *m_container{nullptr};
  Index m_index{0};

  void assignParticleIndices(std::span<const ParticleIndex> particleIndices,
                             std::vector<std::uint32_t> &offsetColumn,
                             std::vector<std::uint8_t> &countColumn)
    requires(!ReadOnly)
  {
    if (countColumn[m_index] != 0) {
      throw std::logic_error("Particle indices already assigned to the vertex");
    }

    offsetColumn[m_index] =
        static_cast<std::uint32_t>(m_container->m_parentIndices.size());
    countColumn[m_index] = static_cast<std::uint8_t>(particleIndices.size());
    m_container->m_parentIndices.insert(m_container->m_parentIndices.end(),
                                        particleIndices.begin(),
                                        particleIndices.end());
  }
};

}  // namespace ActsFatras

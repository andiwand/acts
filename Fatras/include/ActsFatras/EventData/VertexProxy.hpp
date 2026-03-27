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
#include "ActsFatras/EventData/ForwardDeclare.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <cstdint>
#include <span>

namespace Acts {
class Surface;
}

namespace ActsFatras {

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

  Acts::Vector4 &fourPosition() noexcept {
    return container().m_fourPositionColumn[m_index];
  }

  const Acts::Surface *&referenceSurface() noexcept {
    return container().m_surfacePointers[m_index];
  }

  ProcessType &generationProcess() noexcept {
    return container().m_processType[m_index];
  }

  std::span<const ParticleIndex> incomingParticleIndices() const noexcept {
    const auto &offsetColumn = container().m_incomingParticleIndicesOffset;
    const auto &countColumn = container().m_incomingParticleIndicesCount;
    const auto offset = offsetColumn[m_index];
    const auto count = countColumn[m_index];
    return std::span<const ParticleIndex>(
        container().m_parentIndices.data() + offset, count);
  }

  ConstParticleSubset incomingParticles() const noexcept;

  std::span<const ParticleIndex> outgoingParticleIndices() const noexcept {
    const auto &offsetColumn = container().m_outgoingParticleIndicesOffset;
    const auto &countColumn = container().m_outgoingParticleIndicesCount;
    const auto offset = offsetColumn[m_index];
    const auto count = countColumn[m_index];
    return std::span<const ParticleIndex>(
        container().m_parentIndices.data() + offset, count);
  }

  ConstParticleSubset outgoingParticles() const noexcept;

  const Acts::Vector4 &fourPosition() const noexcept {
    return container().m_fourPositionColumn[m_index];
  }

  const Acts::Surface *referenceSurface() const noexcept {
    return container().m_surfacePointers[m_index];
  }

  ProcessType generationProcess() const noexcept {
    return container().m_processType[m_index];
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

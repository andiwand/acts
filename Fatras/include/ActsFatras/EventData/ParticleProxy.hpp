// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/ForwardDeclare.hpp"
#include "ActsFatras/EventData/ParticleContainer.hpp"
#include "ActsFatras/EventData/ParticleOutcome.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"
#include "ActsFatras/EventData/VertexProxy.hpp"

#include <cstdint>
#include <span>

namespace Acts {
class Surface;
}

namespace ActsFatras {

/// A proxy class for accessing individual particles.
template <bool read_only>
class ParticleProxy final {
 public:
  /// Indicates whether this particle proxy is read-only or data can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  /// Type alias for particle index type
  using Index = ParticleIndex;

  /// Type alias for container type (const if read-only)
  using Container = Acts::const_if_t<ReadOnly, ParticleContainer>;

  /// Constructs a particle proxy for the given container and index.
  /// @param container The container holding the particle.
  /// @param index The index of the particle in the container.
  ParticleProxy(Container &container, Index index) noexcept
      : m_container(&container), m_index(index) {}

  /// Copy construct a particle proxy.
  /// @param other The particle proxy to copy.
  ParticleProxy(const ParticleProxy &other) noexcept = default;

  /// Copy construct a mutable particle proxy.
  /// @param other The mutable particle proxy to copy.
  explicit ParticleProxy(const ParticleProxy<false> &other) noexcept
    requires ReadOnly
      : m_container(&other.container()), m_index(other.index()) {}

  /// Copy assign a particle proxy.
  /// @param other The particle proxy to copy.
  /// @return Reference to this particle proxy after assignment.
  ParticleProxy &operator=(const ParticleProxy &other) noexcept = default;

  /// Copy assign a mutable particle proxy.
  /// @param other The mutable particle proxy to copy.
  /// @return Reference to this particle proxy after assignment.
  ParticleProxy &operator=(const ParticleProxy<false> &other) noexcept
    requires ReadOnly
  {
    m_container = &other.container();
    m_index = other.index();
    return *this;
  }

  /// Move assign a particle proxy.
  /// @param other The particle proxy to move.
  /// @return Reference to this particle proxy after assignment.
  ParticleProxy &operator=(ParticleProxy &&other) noexcept = default;

  /// Move assign a mutable particle proxy.
  /// @param other The mutable particle proxy to move.
  /// @return Reference to this particle proxy after assignment.
  ParticleProxy &operator=(ParticleProxy<false> &&other) noexcept
    requires ReadOnly
  {
    m_container = &other.container();
    m_index = other.index();
    return *this;
  }

  /// Returns a const proxy of the particle.
  /// @return A const proxy of the particle.
  ParticleProxy<true> asConst() const noexcept
    requires(!ReadOnly)
  {
    return {*m_container, m_index};
  }

  /// Gets the container holding the particle.
  /// @return A reference to the container holding the particle.
  ParticleContainer &container() const noexcept
    requires(!ReadOnly)
  {
    return *m_container;
  }

  /// Gets the container holding the particle.
  /// @return A const reference to the container holding the particle.
  const ParticleContainer &container() const noexcept { return *m_container; }

  /// Gets the index of the particle in the container.
  /// @return The index of the particle in the container.
  Index index() const noexcept { return m_index; }

  Barcode &barcode() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_barcodeColumn);
  }

  Acts::PdgParticle &pdg() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_pdgColumn);
  }

  float &charge() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_chargeColumn);
  }

  double &mass() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_massColumn);
  }

  ProcessType &generationProcess() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_generationProcessColumn);
  }

  VertexIndex &productionVertexIndex() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_productionVertexIndexColumn);
  }

  const Acts::Surface *&endReferenceSurface() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_endReferenceSurfaceColumn);
  }

  Acts::Vector4 &endFourPosition() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_endFourPositionColumn);
  }

  Acts::Vector3 &endMomentum() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_endMomentumColumn);
  }

  double &endProperTime() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_properTimeColumn);
  }

  double &endPathInX0() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_pathInX0Column);
  }

  double &endPathInL0() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_pathInL0Column);
  }

  ParticleOutcome &endOutcome() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_endOutcomeColumn);
  }

  void assignHitIndices(std::span<const HitIndex> hitIndices) noexcept
    requires(!ReadOnly)
  {
    if (accessImpl(m_container->m_hitIndicesCountColumn) != 0) {
      throw std::logic_error("Hits already assigned to the particle");
    }

    accessImpl(m_container->m_hitIndicesOffsetColumn[m_index]) =
        static_cast<std::uint32_t>(m_container->m_spacePoints.size());
    accessImpl(m_container->m_hitIndicesCountColumn[m_index]) =
        static_cast<std::uint8_t>(hitIndices.size());
    m_container->m_hitIndices.insert(m_container->m_spacePoints.end(),
                                     hitIndices.begin(), hitIndices.end());
  }

  const Barcode &barcode() const noexcept {
    return accessImpl(m_container->m_barcodeColumn);
  }

  Acts::PdgParticle pdg() const noexcept {
    return accessImpl(m_container->m_pdgColumn);
  }

  float charge() const noexcept {
    return accessImpl(m_container->m_chargeColumn);
  }

  double mass() const noexcept { return accessImpl(m_container->m_massColumn); }

  VertexIndex productionVertexIndex() const noexcept {
    return accessImpl(m_container->m_productionVertexIndexColumn);
  }

  ConstVertexProxy productionVertex() const noexcept {
    return m_container->vertexContainer().at(
        accessImpl(m_container->m_productionVertexIndexColumn));
  }

  ProcessType productionProcess() const noexcept {
    return m_container->vertexContainer()
        .at(accessImpl(m_container->m_productionVertexIndexColumn))
        .process();
  }

  const Acts::Vector4 &productionFourPosition() const noexcept {
    return accessImpl(m_container->m_initialFourPositionColumn);
  }

  const Acts::Vector3 &productionMomentum() const noexcept {
    return accessImpl(m_container->m_productionMomentumColumn);
  }

  Acts::Vector3 productionDirection() const noexcept {
    return accessImpl(m_container->m_productionMomentumColumn).normalized();
  }

  double initialAbsoluteMomentum() const noexcept {
    return accessImpl(m_container->m_productionMomentumColumn).norm();
  }

  const Acts::Surface *productionReferenceSurface() const noexcept {
    return m_container->vertexContainer()
        .at(accessImpl(m_container->m_productionVertexIndexColumn))
        .referenceSurface();
  }

  const Acts::Vector4 &endFourPosition() const noexcept {
    return accessImpl(m_container->m_endFourPositionColumn);
  }

  const Acts::Vector3 &endMomentum() const noexcept {
    return accessImpl(m_container->m_endMomentumColumn);
  }

  Acts::Vector3 endDirection() const noexcept {
    return accessImpl(m_container->m_endMomentumColumn).normalized();
  }

  double endAbsoluteMomentum() const noexcept {
    return accessImpl(m_container->m_endMomentumColumn).norm();
  }

  double endProperTime() const noexcept {
    return accessImpl(m_container->m_properTimeColumn);
  }

  double endPathInX0() const noexcept {
    return accessImpl(m_container->m_pathInX0Column);
  }

  double endPathInL0() const noexcept {
    return accessImpl(m_container->m_pathInL0Column);
  }

  const ParticleOutcome &endOutcome() const noexcept {
    return accessImpl(m_container->m_endOutcomeColumn);
  }

  std::uint32_t numberOfHits() const noexcept {
    return accessImpl(m_container->m_numberOfHitsColumn);
  }

  /// Const access to an extra column of data for the particle.
  /// @param column The extra column proxy to access.
  /// @return A const reference to the value in the extra column for the particle.
  template <typename T>
  const T &extra(const ConstParticleColumnProxy<T> &column) const noexcept {
    return column[m_index];
  }

  /// Check if this is a secondary particle.
  /// @return True if particle is a secondary (has non-zero vertex secondary, generation, or sub-particle), false otherwise
  bool isSecondary() const {
    return barcode().vertexSecondary() != 0 || barcode().generation() != 0 ||
           barcode().subParticle() != 0;
  }

  Acts::PdgParticle absolutePdg() const {
    return Acts::makeAbsolutePdgParticle(pdg());
  }

  /// Particle absolute charge.
  /// @return The absolute particle charge (positive value)
  double absoluteCharge() const { return std::abs(charge()); }

  /// Particle hypothesis.
  /// @return Particle hypothesis containing PDG, mass, and charge information
  Acts::ParticleHypothesis hypothesis() const {
    return Acts::ParticleHypothesis(absolutePdg(), static_cast<float>(mass()),
                                    static_cast<float>(absoluteCharge()));
  }

 private:
  Container *m_container{nullptr};
  Index m_index{0};

  template <typename OptColumn>
  auto &accessImpl(OptColumn &&column) const {
    assert(column.has_value() && "Column does not exist");
    assert(m_index < column->size() && "Index out of bounds");
    return column->proxy(*m_container)[m_index];
  }
};

}  // namespace ActsFatras

inline std::ostream &operator<<(std::ostream &os,
                                ActsFatras::ConstParticleProxy particle) {
  // compact format w/ only identity information but no kinematics
  os << "id=" << particle.index();
  os << "|barcode=" << "(" << particle.barcode() << ")";
  os << "|pdg=" << particle.pdg();
  os << "|q=" << particle.charge();
  os << "|m=" << particle.mass();
  os << "|p=" << particle.initialAbsoluteMomentum();
  return os;
}

inline std::ostream &operator<<(std::ostream &os,
                                ActsFatras::MutableParticleProxy particle) {
  return os << particle.asConst();
}

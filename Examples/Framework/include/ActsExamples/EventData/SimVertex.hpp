// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsFatras/EventData/GenerationProcessType.hpp"

#include <boost/container/flat_set.hpp>

namespace ActsExamples {

using SimVertexIndex = std::uint32_t;

class SimVertexBarcode {
 public:
  using PrimaryVertexId = SimBarcode::PrimaryVertexId;
  using SecondaryVertexId = SimBarcode::SecondaryVertexId;
  using ParticleId = SimBarcode::ParticleId;
  using GenerationId = SimBarcode::GenerationId;
  using SubParticleId = SimBarcode::SubParticleId;

  explicit constexpr SimVertexBarcode(SimBarcode barcode)
      : m_id(barcode.vertexId()) {}

  constexpr SimVertexBarcode() = default;

  /// Return the barcode.
  constexpr SimBarcode barcode() const { return m_id; }

  /// Return the primary vertex identifier.
  constexpr PrimaryVertexId vertexPrimary() const {
    return m_id.vertexPrimary();
  }
  /// Return the secondary vertex identifier.
  constexpr SecondaryVertexId vertexSecondary() const {
    return m_id.vertexSecondary();
  }
  /// Return the generation identifier.
  constexpr GenerationId generation() const { return m_id.generation(); }

  /// Create a new barcode with a different primary vertex identifier.
  [[nodiscard]]
  constexpr SimVertexBarcode withVertexPrimary(PrimaryVertexId id) const {
    return SimVertexBarcode(m_id.withVertexPrimary(id));
  }
  /// Create a new barcode with a different secondary vertex identifier.
  [[nodiscard]]
  constexpr SimVertexBarcode withVertexSecondary(SecondaryVertexId id) const {
    return SimVertexBarcode(m_id.withVertexSecondary(id));
  }
  /// Create a new barcode with a different generation identifier.
  [[nodiscard]]
  constexpr SimVertexBarcode withGeneration(GenerationId id) const {
    return SimVertexBarcode(m_id.withGeneration(id));
  }

  std::size_t hash() const { return m_id.hash(); }

 private:
  /// The vertex ID
  /// Note that only primary, secondary and generation should be set
  SimBarcode m_id;

  friend constexpr bool operator<(SimVertexBarcode lhs, SimVertexBarcode rhs) {
    return lhs.m_id < rhs.m_id;
  }

  friend constexpr bool operator==(SimVertexBarcode lhs, SimVertexBarcode rhs) {
    return lhs.m_id == rhs.m_id;
  }

  friend inline std::ostream& operator<<(std::ostream& os,
                                         SimVertexBarcode idx) {
    return os << idx.m_id;
  }
};

/// A simulated vertex e.g. from a physics process.
struct SimVertex {
  /// The vertex ID
  SimVertexBarcode barcode = SimVertexBarcode(SimBarcode::Invalid());
  /// The vertex four-position
  Acts::Vector4 position4 = Acts::Vector4::Zero();
  /// The vertex process type
  ActsFatras::GenerationProcessType process =
      ActsFatras::GenerationProcessType::eUndefined;
  /// The incoming particles into the vertex
  std::vector<SimParticleIndex> incoming;
  /// The outgoing particles from the vertex
  std::vector<SimParticleIndex> outgoing;

  /// Construct the vertex from a position and optional process type.
  ///
  /// @param position4_ the vertex four-position
  /// @param process_ the process type that generated this vertex
  ///
  /// Associated particles are left empty by default and must be filled by the
  /// user after construction.
  SimVertex(SimVertexBarcode barcode_, const Acts::Vector4& position4_,
            ActsFatras::GenerationProcessType process_ =
                ActsFatras::GenerationProcessType::eUndefined)
      : barcode(barcode_), position4(position4_), process(process_) {}
  // explicitly default rule-of-five.
  SimVertex() = default;
  SimVertex(const SimVertex&) = default;
  SimVertex(SimVertex&&) = default;
  SimVertex& operator=(const SimVertex&) = default;
  SimVertex& operator=(SimVertex&&) = default;

  /// The vertex three-position.
  auto position() const { return position4.head<3>(); }
  /// The vertex time.
  double time() const { return position4[3]; }
};

using SimVertexContainer = std::vector<SimVertex>;

}  // namespace ActsExamples

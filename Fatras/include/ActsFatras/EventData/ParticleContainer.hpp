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
#include "Acts/Utilities/EnumBitwiseOperators.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Utilities/detail/ContainerIterator.hpp"
#include "Acts/Utilities/detail/ContainerSubset.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/ForwardDeclare.hpp"
#include "ActsFatras/EventData/SimulationOutcome.hpp"

#include <cstdint>
#include <limits>
#include <span>

namespace Acts {
class Surface;
}

namespace ActsFatras {

enum class ParticleColumns : std::uint32_t {
  None = 0,

  Barcode = 1 << 0,
  Pdg = 1 << 1,
  Charge = 1 << 2,
  Mass = 1 << 3,
  ProductionVertex = 1 << 4,
  ProductionMomentum = 1 << 5,

  EndVertex = 1 << 6,
  EndReferenceSurface = 1 << 7,
  EndFourPosition = 1 << 8,
  EndMomentum = 1 << 9,
  EndProperTime = 1 << 10,
  EndPathInX0 = 1 << 11,
  EndPathInL0 = 1 << 12,
  EndOutcome = 1 << 13,

  Vertices = 1 << 14,
  Hits = 1 << 15,

  Generated =
      Barcode | Pdg | Charge | Mass | ProductionVertex | ProductionMomentum,
  Simulated = Generated | EndVertex | EndReferenceSurface | EndFourPosition |
              EndMomentum | EndProperTime | EndPathInX0 | EndPathInL0 |
              EndOutcome | Vertices | Hits,
  All = Barcode | Pdg | Charge | Mass | ProductionVertex | EndVertex |
        EndReferenceSurface | EndFourPosition | EndMomentum | EndProperTime |
        EndPathInX0 | EndPathInL0 | EndOutcome | Vertices | Hits,
};

/// Enable bitwise operators for ParticleColumns enum
ACTS_DEFINE_ENUM_BITWISE_OPERATORS(ParticleColumns);

namespace detail {

class ColumnHolderBase {
 public:
  virtual ~ColumnHolderBase() = default;

  virtual std::unique_ptr<ColumnHolderBase> copy() const = 0;

  virtual std::size_t size() const = 0;
  virtual void reserve(std::size_t size) = 0;
  virtual void resize(std::size_t size) = 0;
  virtual void clear() = 0;
  virtual void emplace_back() = 0;
};

template <typename T>
class ColumnHolder final : public ColumnHolderBase {
 public:
  using Value = T;
  using Container = std::vector<Value>;
  using MutableProxy = MutableParticleColumnProxy<Value>;
  using ConstProxy = ConstParticleColumnProxy<Value>;

  ColumnHolder() = default;
  explicit ColumnHolder(Value defaultValue)
      : m_default(std::move(defaultValue)) {}

  MutableProxy proxy(ParticleContainer &container) {
    return MutableProxy(container, m_data);
  }
  ConstProxy proxy(const ParticleContainer &container) const {
    return ConstProxy(container, m_data);
  }

  std::unique_ptr<ColumnHolderBase> copy() const override {
    return std::make_unique<ColumnHolder<T>>(*this);
  }

  std::size_t size() const override { return m_data.size(); }
  void reserve(std::size_t size) override { m_data.reserve(size); }
  void clear() override { m_data.clear(); }
  void resize(std::size_t size) override { m_data.resize(size, m_default); }
  void emplace_back() override { m_data.emplace_back(m_default); }

 private:
  Value m_default{};
  Container m_data;
};

}  // namespace detail

class ParticleContainer final {
 public:
  /// Type alias for particle index in container
  using Index = ParticleIndex;
  /// Type alias for subset of particle indices
  using IndexSubset = ParticleIndexSubset;
  /// Type alias for mutable particle proxy
  using MutableProxy = MutableParticleProxy;
  /// Type alias for const particle proxy
  using ConstProxy = ConstParticleProxy;
  /// Type alias for mutable particle subset
  using MutableSubset = MutableParticleSubset;
  /// Type alias for const particle subset
  using ConstSubset = ConstParticleSubset;

  explicit ParticleContainer(
      ParticleColumns columns = ParticleColumns::None) noexcept;

  /// Constructs a copy of the given particle container.
  /// @param other The particle container to copy.
  ParticleContainer(const ParticleContainer &other) noexcept;

  /// Move constructs a particle container.
  /// @param other The particle container to move.
  ParticleContainer(ParticleContainer &&other) noexcept;

  /// Detructs the particle container.
  ~ParticleContainer() noexcept = default;

  /// Assignment operator for copying a particle container.
  /// @param other The particle container to copy.
  /// @return A reference to this particle container.
  ParticleContainer &operator=(const ParticleContainer &other) noexcept;

  /// Move assignment operator for a particle container.
  /// @param other The particle container to move.
  /// @return A reference to this particle container.
  ParticleContainer &operator=(ParticleContainer &&other) noexcept;

  /// Assigns a vertex container to this vertex container.
  /// @param vertexContainer The vertex container to assign.
  /// @throws std::logic_error if a vertex container has already been assigned.
  void assignVertexContainer(const VertexContainer &vertexContainer);

  /// Assigns a hit container to this particle container.
  /// @param hitContainer The hit container to assign.
  /// @throws std::logic_error if a hit container has already been assigned.
  void assignHitContainer(const HitContainer &hitContainer);

  /// Checks if a vertex container has been assigned to this vertex container.
  /// @return True if a vertex container has been assigned.
  bool hasVertexContainer() const noexcept;

  /// Checks if a hit container has been assigned to this particle container.
  /// @return True if a hit container has been assigned.
  bool hasHitContainer() const noexcept;

  /// Returns a const reference to the assigned vertex container.
  /// @return A const reference to the assigned vertex container.
  /// @throws std::logic_error if no vertex container has been assigned.
  const VertexContainer &vertexContainer() const;

  /// Returns a const reference to the assigned hit container.
  /// @return A const reference to the assigned hit container.
  /// @throws std::logic_error if no hit container has been assigned.
  const HitContainer &hitContainer() const;

  /// Returns the number of particles in the container.
  /// @return The number of particles in the container.
  [[nodiscard]] std::uint32_t size() const noexcept { return m_size; }
  /// Checks if the container is empty.
  /// @return True if the container is empty, false otherwise.
  [[nodiscard]] bool empty() const noexcept { return size() == 0; }

  /// Reserves space for the given number of particles.
  /// @param size The number of particles to reserve space for.
  /// @param averageVerticesPerParticle The average number of vertices per particle.
  /// @param averageHitsPerParticle The average number of hits per particle.
  void reserve(std::uint32_t size, float averageVerticesPerParticle = 3,
               float averageHitsPerParticle = 10) noexcept;

  /// Clears the container, removing all particles and columns.
  void clear() noexcept;

  /// Creates a new particle at the end of the container.
  /// @return A mutable proxy to the newly created particle.
  MutableProxy createParticle() noexcept;

  /// Creates additional columns. This will create the columns if they do not
  /// already exist.
  /// @param columns The columns to create.
  void createColumns(ParticleColumns columns) noexcept;

  /// Drops the specified columns from the container.
  /// This will only drop columns if they exist.
  /// @param columns The columns to drop.
  void dropColumns(ParticleColumns columns) noexcept;

  /// Checks if the container has the given Columns.
  /// @param columns The Columns to check for.
  /// @return True if the container has all the specified Columns, false
  ///         otherwise.
  bool hasColumns(ParticleColumns columns) const noexcept {
    return (m_knownColumns & columns) == columns;
  }

  /// Creates a new column with the given name.
  /// If a column with the same name already exists, an exception is thrown.
  /// @param name The name of the column.
  /// @return A reference to the newly created column.
  /// @throws std::runtime_error if a column with the same name already exists.
  /// @throws std::runtime_error if the column name is reserved.
  template <typename T>
  MutableParticleColumnProxy<T> createColumn(const std::string &name) {
    return createColumnImpl<ColumnHolder<T>>(name);
  }

  /// Drops the column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @throws std::runtime_error if the column does not exist.
  /// @throws std::runtime_error if the column name is reserved.
  void dropColumn(const std::string &name);

  /// Checks if an Column with the given name exists.
  /// @param name The name of the column.
  /// @return True if the column exists, false otherwise.
  bool hasColumn(const std::string &name) const noexcept {
    return m_allColumns.contains(name);
  }

  /// Returns a mutable reference to the column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @return A mutable reference to the column.
  /// @throws std::runtime_error if the column does not exist.
  template <typename T>
  MutableParticleColumnProxy<T> column(const std::string &name) {
    return columnImpl<ColumnHolder<T>>(name);
  }

  /// Returns a const reference to the column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @return A const reference to the column.
  /// @throws std::runtime_error if the column does not exist.
  template <typename T>
  ConstParticleColumnProxy<T> column(const std::string &name) const {
    return columnImpl<ColumnHolder<T>>(name);
  }

  /// Returns a mutable proxy to the particle at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the particle to access.
  /// @return A mutable proxy to the particle at the given index.
  /// @throws std::out_of_range if the index is out of range.
  MutableProxy at(Index index);
  /// Returns a const proxy to the particle at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the particle to access.
  /// @return A const proxy to the particle at the given index.
  /// @throws std::out_of_range if the index is out of range.
  ConstProxy at(Index index) const;

  /// Returns a mutable proxy to the particle at the given index.
  /// @param index The index of the particle to access.
  /// @return A mutable proxy to the particle at the given index.
  MutableProxy operator[](Index index) noexcept;
  /// Returns a const proxy to the particle at the given index.
  /// @param index The index of the particle to access.
  /// @return A const proxy to the particle at the given index.
  ConstProxy operator[](Index index) const noexcept;

  /// Type alias for template iterator over particles in container
  template <bool read_only>
  using Iterator = Acts::detail::ContainerIterator<
      ParticleContainer,
      std::conditional_t<read_only, ConstProxy, MutableProxy>, Index,
      read_only>;

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

  /// Creates a mutable subset of particles from the given index subset.
  /// @param subset The index subset to create the subset from.
  /// @return A mutable subset of particles.
  MutableSubset subset(const IndexSubset &subset) noexcept;
  /// Creates a const subset of particles from the given index subset.
  /// @param subset The index subset to create the subset from.
  /// @return A const subset of particles.
  ConstSubset subset(const IndexSubset &subset) const noexcept;

 public:
  using ColumnHolderBase = detail::ColumnHolderBase;
  template <typename T>
  using ColumnHolder = detail::ColumnHolder<T>;

  template <bool>
  friend class ParticleProxy;

  std::size_t m_size{0};

  std::unordered_map<std::string, ColumnHolderBase *> m_allColumns;
  ParticleColumns m_knownColumns{ParticleColumns::None};
  std::unordered_map<std::string, std::unique_ptr<ColumnHolderBase>>
      m_dynamicColumns;

  std::optional<ColumnHolder<Barcode>> m_barcodeColumn;
  std::optional<ColumnHolder<Acts::PdgParticle>> m_pdgColumn;
  std::optional<ColumnHolder<float>> m_chargeColumn;
  std::optional<ColumnHolder<double>> m_massColumn;
  std::optional<ColumnHolder<VertexIndex>> m_productionVertexIndexColumn;
  std::optional<ColumnHolder<Acts::Vector3>> m_productionMomentumColumn;

  std::optional<ColumnHolder<VertexIndex>> m_endVertexIndexColumn;
  std::optional<ColumnHolder<const Acts::Surface *>>
      m_endReferenceSurfaceColumn;
  std::optional<ColumnHolder<Acts::Vector4>> m_endFourPositionColumn;
  std::optional<ColumnHolder<Acts::Vector3>> m_endMomentumColumn;
  std::optional<ColumnHolder<double>> m_endProperTimeColumn;
  std::optional<ColumnHolder<double>> m_endPathInX0Column;
  std::optional<ColumnHolder<double>> m_endPathInL0Column;
  std::optional<ColumnHolder<SimulationOutcome>> m_endOutcomeColumn;

  std::optional<ColumnHolder<std::uint32_t>> m_vertexIndicesOffsetColumn;
  std::optional<ColumnHolder<std::uint8_t>> m_vertexIndicesCountColumn;

  std::optional<ColumnHolder<std::uint32_t>> m_hitIndicesOffsetColumn;
  std::optional<ColumnHolder<std::uint8_t>> m_hitIndicesCountColumn;

  std::vector<VertexIndex> m_vertexIndices;
  std::vector<HitIndex> m_hitIndices;

  const VertexContainer *m_vertexContainer{nullptr};
  const HitContainer *m_hitContainer{nullptr};

  static auto knownColumnMasks() noexcept {
    using enum ParticleColumns;
    return std::tuple(Barcode, Pdg, Charge, Mass, ProductionVertex,
                      ProductionMomentum, EndVertex, EndReferenceSurface,
                      EndFourPosition, EndMomentum, EndProperTime, EndPathInX0,
                      EndPathInL0, EndOutcome, Vertices, Vertices, Hits, Hits);
  }

  static auto knownColumnNames() noexcept {
    return std::tuple(
        "barcode", "pdg", "charge", "mass", "productionVertexIndex",
        "productionMomentum", "endVertexIndex", "endReferenceSurface",
        "endFourPosition", "endMomentum", "endProperTime", "endPathInX0",
        "endPathInL0", "endOutcome", "vertexIndicesOffset",
        "vertexIndicesCount", "hitIndicesOffset", "hitIndicesCount");
  }

  static auto knownColumnDefaults() noexcept {
    return std::tuple(Barcode{}, Acts::PdgParticle::eInvalid, float{0},
                      double{0}, std::numeric_limits<VertexIndex>::max(),
                      Acts::Vector3(Acts::Vector3::Zero()),
                      std::numeric_limits<VertexIndex>::max(),
                      static_cast<const Acts::Surface *>(nullptr),
                      Acts::Vector4(Acts::Vector4::Zero()),
                      Acts::Vector3(Acts::Vector3::Zero()), double{0},
                      double{0}, double{0}, SimulationOutcome::Alive,
                      std::uint32_t{0}, std::uint8_t{0}, std::uint32_t{0},
                      std::uint8_t{0});
  }

  template <typename Self>
  static auto knownColumns(Self &&self) noexcept {
    return std::tie(
        self.m_barcodeColumn, self.m_pdgColumn, self.m_chargeColumn,
        self.m_massColumn, self.m_productionVertexIndexColumn,
        self.m_productionMomentumColumn, self.m_endVertexIndexColumn,
        self.m_endReferenceSurfaceColumn, self.m_endFourPositionColumn,
        self.m_endMomentumColumn, self.m_endProperTimeColumn,
        self.m_endPathInX0Column, self.m_endPathInL0Column,
        self.m_endOutcomeColumn, self.m_vertexIndicesOffsetColumn,
        self.m_vertexIndicesCountColumn, self.m_hitIndicesOffsetColumn,
        self.m_hitIndicesCountColumn);
  }
  auto knownColumns() & noexcept { return knownColumns(*this); }
  auto knownColumns() const & noexcept { return knownColumns(*this); }
  auto knownColumns() && noexcept { return knownColumns(*this); }

  void copyColumns(const ParticleContainer &other);
  void moveColumns(ParticleContainer &other) noexcept;

  static bool reservedColumn(const std::string &name) noexcept;

  template <typename Holder>
  MutableParticleColumnProxy<typename Holder::Value> createColumnImpl(
      const std::string &name) {
    if (reservedColumn(name)) {
      throw std::runtime_error("Column name is reserved: " + name);
    }
    if (hasColumn(name)) {
      throw std::runtime_error("Column already exists: " + name);
    }
    auto holder = std::make_unique<Holder>();
    holder->resize(size());
    auto proxy = holder->proxy(*this);
    m_allColumns.try_emplace(name, holder.get());
    m_dynamicColumns.try_emplace(name, std::move(holder));
    return proxy;
  }

  template <typename Holder, typename Self>
  static auto columnImpl(Self &&self, const std::string &name) {
    const auto it = self.m_allColumns.find(name);
    if (it == self.m_allColumns.end()) {
      throw std::runtime_error("Column not found: " + name);
    }
    auto &holder = dynamic_cast<Holder &>(*it->second);
    return holder.proxy(self);
  }

  template <typename Holder>
  MutableParticleColumnProxy<typename Holder::Value> columnImpl(
      const std::string &name) {
    return columnImpl<Holder>(*this, name);
  }

  template <typename Holder>
  ConstParticleColumnProxy<typename Holder::Value> columnImpl(
      const std::string &name) const {
    return columnImpl<Holder>(*this, name);
  }
};

/// Subset facade over arbitrary index sets.
template <bool read_only>
class ParticleSubset final
    : public Acts::detail::ContainerSubset<
          ParticleSubset<read_only>, ParticleSubset<true>, ParticleContainer,
          std::conditional_t<read_only, ConstParticleProxy,
                             MutableParticleProxy>,
          std::span<const ParticleIndex>, read_only> {
 public:
  /// Base class type
  using Base = Acts::detail::ContainerSubset<
      ParticleSubset<read_only>, ParticleSubset<true>, ParticleContainer,
      std::conditional_t<read_only, ConstParticleProxy, MutableParticleProxy>,
      std::span<const ParticleIndex>, read_only>;

  using Base::Base;
};

/// Additional column of data that can be added to the particle container.
/// The column is indexed by the particle index.
template <typename T, bool read_only>
class ParticleColumnProxy final {
 public:
  /// Flag indicating whether this particle column proxy is read-only
  constexpr static bool ReadOnly = read_only;
  /// Type alias for particle index type
  using Index = ParticleIndex;
  /// Type alias for particle index subset type
  using IndexSubset = ParticleIndexSubset;
  /// Type alias for column value type
  using Value = T;
  /// Type alias for container type (const if read-only)
  using Container = Acts::const_if_t<ReadOnly, ParticleContainer>;
  /// Type alias for column container type (const if read-only)
  using Column = Acts::const_if_t<ReadOnly, std::vector<Value>>;

  /// Constructs a particle column proxy for the given container and column.
  /// @param container The container holding the particle.
  /// @param column The column of data to access.
  ParticleColumnProxy(Container &container, Column &column)
      : m_container{&container}, m_column(&column) {}

  /// Copy construct a particle column proxy.
  /// @param other The particle column proxy to copy.
  ParticleColumnProxy(const ParticleColumnProxy &other) noexcept = default;

  /// Copy construct a mutable particle column proxy.
  /// @param other The mutable particle column proxy to copy.
  explicit ParticleColumnProxy(
      const ParticleColumnProxy<T, false> &other) noexcept
    requires ReadOnly
      : m_container(&other.container()), m_column(&other.column()) {}

  /// Returns the number of entries in the particle column.
  /// @return The size of the particle column.
  std::uint32_t size() const noexcept { return column().size(); }

  /// Returns a const proxy of the particle column.
  /// @return A const proxy of the particle column.
  ParticleColumnProxy<T, true> asConst() const noexcept
    requires(!ReadOnly)
  {
    return {*m_container, *m_column};
  }

  /// Gets the container holding the particle.
  /// @return A reference to the container holding the particle.
  ParticleContainer &container() noexcept
    requires(!ReadOnly)
  {
    return *m_container;
  }
  /// Gets the container holding the particle.
  /// @return A const reference to the container holding the particle.
  const ParticleContainer &container() const noexcept { return *m_container; }

  /// Returns a const reference to the column container.
  /// @return A const reference to the column container.
  const std::vector<Value> &column() const noexcept { return *m_column; }

  /// Returns a mutable span to the column data.
  /// @return A mutable span to the column data.
  std::span<Value> data() noexcept
    requires(!ReadOnly)
  {
    return std::span<Value>(column().data(), column().size());
  }
  /// Returns a const span to the column data.
  /// @return A const span to the column data.
  std::span<const Value> data() const noexcept {
    return std::span<const Value>(column().data(), column().size());
  }

  /// Returns a mutable reference to the column entry at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the particle to access.
  /// @return A mutable reference to the column entry at the given index.
  /// @throws std::out_of_range if the index is out of range.
  Value &at(Index index)
    requires(!ReadOnly)
  {
    if (index >= column().size()) {
      throw std::out_of_range(
          "Index out of range in ParticleContainer: " + std::to_string(index) +
          " >= " + std::to_string(size()));
    }
    return data()[index];
  }
  /// Returns a const reference to the column entry at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the particle to access.
  /// @return A const reference to the column entry at the given index.
  /// @throws std::out_of_range if the index is out of range.
  const Value &at(Index index) const {
    if (index >= column().size()) {
      throw std::out_of_range(
          "Index out of range in ParticleContainer: " + std::to_string(index) +
          " >= " + std::to_string(size()));
    }
    return data()[index];
  }

  /// Returns a mutable reference to the column entry at the given index.
  /// @param index The index of the particle to access.
  /// @return A mutable reference to the column entry at the given index.
  Value &operator[](Index index) noexcept
    requires(!ReadOnly)
  {
    assert(index < column().size() && "Index out of bounds");
    return data()[index];
  }
  /// Returns a const reference to the column entry at the given index.
  /// @param index The index of the particle to access.
  /// @return A const reference to the column entry at the given index.
  const Value &operator[](Index index) const noexcept {
    assert(index < column().size() && "Index out of bounds");
    return data()[index];
  }

  /// Subset view over selected column entries.
  class Subset final
      : public Acts::detail::ContainerSubset<Subset, Subset, Column, Value,
                                             Index, ReadOnly> {
   public:
    /// Base class type
    using Base = Acts::detail::ContainerSubset<Subset, Subset, Column, Value,
                                               Index, ReadOnly>;

    using Base::Base;
  };

  /// Creates a subset view of this particle column based on provided
  /// indices.
  ///
  /// This method creates a subset proxy that provides access to only the
  /// particles at the indices specified in the IndexSubset. The subset
  /// maintains a reference to the original column data without copying,
  /// enabling efficient access to selected particles for filtering, clustering,
  /// or other operations.
  ///
  /// @param subset The index subset specifying which particles to include
  /// @return A subset proxy providing access to the selected particles
  ///
  /// @note The returned subset shares data with the original column
  /// @note The subset remains valid only as long as the original column exists
  /// @note This operation does not copy data, providing efficient subset access
  Subset subset(const IndexSubset &subset) const noexcept {
    return Subset(*m_column, subset);
  }

 private:
  Container *m_container{};
  Column *m_column{};

  std::vector<Value> &column() noexcept
    requires(!ReadOnly)
  {
    return *m_column;
  }
};

}  // namespace ActsFatras

#include "ActsFatras/EventData/ParticleContainer.ipp"

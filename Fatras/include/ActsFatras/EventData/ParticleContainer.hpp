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
#include "ActsFatras/EventData/ParticleOutcome.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <cstdint>
#include <span>

namespace Acts {
class Surface;
}

namespace ActsFatras {

using ParticleIndex = std::uint32_t;
using ParticleIndexSubset = std::span<const ParticleIndex>;

using HitIndex = std::uint32_t;

class ParticleContainer;

template <bool>
class ParticleProxy;
using MutableParticleProxy = ParticleProxy<false>;
using ConstParticleProxy = ParticleProxy<true>;

template <typename T, bool>
class ParticleColumnProxy;
template <typename T>
using MutableParticleColumnProxy = ParticleColumnProxy<T, false>;
template <typename T>
using ConstParticleColumnProxy = ParticleColumnProxy<T, true>;

enum class ParticleColumns : std::uint32_t {
  None = 0,
  Parents = 1 << 0,
  Barcode = 1 << 1,
  Process = 1 << 2,
  Pdg = 1 << 3,
  Charge = 1 << 4,
  Mass = 1 << 5,
  Direction = 1 << 6,
  AbsoluteMomentum = 1 << 7,
  Position4 = 1 << 8,
  ProperTime = 1 << 9,
  PathInX0 = 1 << 10,
  PathInL0 = 1 << 11,
  NumberOfHits = 1 << 12,
  ReferenceSurface = 1 << 13,
  Outcome = 1 << 14,

  Generated = Barcode | Parents | Process | Pdg | Charge | Mass | Direction |
              AbsoluteMomentum | Position4,
  Simulated = Generated | ProperTime | PathInX0 | PathInL0 | NumberOfHits |
              ReferenceSurface | Outcome,
  All = Barcode | Process | Pdg | Charge | Mass | Direction | AbsoluteMomentum |
        Position4 | ProperTime | PathInX0 | PathInL0 | NumberOfHits |
        ReferenceSurface | Outcome,
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

  explicit ParticleContainer(
      ParticleColumns columns = ParticleColumns::None) noexcept;

  /// Constructs a copy of the given space point container.
  /// @param other The space point container to copy.
  ParticleContainer(const ParticleContainer &other) noexcept;

  /// Move constructs a space point container.
  /// @param other The space point container to move.
  ParticleContainer(ParticleContainer &&other) noexcept;

  /// Detructs the space point container.
  ~ParticleContainer() noexcept = default;

  /// Assignment operator for copying a space point container.
  /// @param other The space point container to copy.
  /// @return A reference to this space point container.
  ParticleContainer &operator=(const ParticleContainer &other) noexcept;

  /// Move assignment operator for a space point container.
  /// @param other The space point container to move.
  /// @return A reference to this space point container.
  ParticleContainer &operator=(ParticleContainer &&other) noexcept;

  /// Returns the number of particles in the container.
  /// @return The number of particles in the container.
  std::uint32_t size() const noexcept { return m_size; }
  /// Checks if the container is empty.
  /// @return True if the container is empty, false otherwise.
  [[nodiscard]] bool empty() const noexcept { return size() == 0; }

  /// Reserves space for the given number of particles.
  /// @param size The number of particles to reserve space for.
  /// @param averageParentIndices The average number of parent indices per particle.
  void reserve(std::uint32_t size, float averageParentIndices = 1) noexcept;

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

  /// Returns a mutable proxy to the space point at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the space point to access.
  /// @return A mutable proxy to the space point at the given index.
  /// @throws std::out_of_range if the index is out of range.
  MutableProxy at(Index index);
  /// Returns a const proxy to the space point at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the space point to access.
  /// @return A const proxy to the space point at the given index.
  /// @throws std::out_of_range if the index is out of range.
  ConstProxy at(Index index) const;

  /// Returns a mutable proxy to the space point at the given index.
  /// @param index The index of the space point to access.
  /// @return A mutable proxy to the space point at the given index.
  MutableProxy operator[](Index index) noexcept;
  /// Returns a const proxy to the space point at the given index.
  /// @param index The index of the space point to access.
  /// @return A const proxy to the space point at the given index.
  ConstProxy operator[](Index index) const noexcept;

  /// Type alias for template iterator over space points in container
  template <bool read_only>
  using Iterator = Acts::detail::ContainerIterator<
      ParticleContainer,
      std::conditional_t<read_only, ConstProxy, MutableProxy>, Index,
      read_only>;

  /// Type alias for mutable iterator over space points
  using iterator = Iterator<false>;
  /// Type alias for const iterator over space points
  using const_iterator = Iterator<true>;

  /// @brief Returns mutable iterator to the beginning of the container
  /// @return Mutable iterator pointing to the first space point
  iterator begin() noexcept { return iterator(*this, 0); }
  /// @brief Returns mutable iterator to the end of the container
  /// @return Mutable iterator pointing past the last space point
  iterator end() noexcept { return iterator(*this, size()); }

  /// @brief Returns const iterator to the beginning of the container
  /// @return Const iterator pointing to the first space point
  const_iterator begin() const noexcept { return const_iterator(*this, 0); }
  /// @brief Returns const iterator to the end of the container
  /// @return Const iterator pointing past the last space point
  const_iterator end() const noexcept { return const_iterator(*this, size()); }

  /// Subset facade over arbitrary index sets.
  template <bool read_only>
  class Subset : public Acts::detail::ContainerSubset<
                     Subset<read_only>, Subset<true>, ParticleContainer,
                     std::conditional_t<read_only, ConstProxy, MutableProxy>,
                     Index, read_only> {
   public:
    /// Base class type
    using Base = Acts::detail::ContainerSubset<
        Subset<read_only>, Subset<true>, ParticleContainer,
        std::conditional_t<read_only, ConstProxy, MutableProxy>, Index,
        read_only>;

    using Base::Base;
  };

  /// Type alias for mutable subset of space points
  using MutableSubset = Subset<false>;
  /// Type alias for const subset of space points
  using ConstSubset = Subset<true>;

  /// Creates a mutable subset of space points from the given index subset.
  /// @param subset The index subset to create the subset from.
  /// @return A mutable subset of space points.
  MutableSubset subset(const IndexSubset &subset) noexcept {
    return MutableSubset(*this, subset);
  }
  /// Creates a const subset of space points from the given index subset.
  /// @param subset The index subset to create the subset from.
  /// @return A const subset of space points.
  ConstSubset subset(const IndexSubset &subset) const noexcept {
    return ConstSubset(*this, subset);
  }

 public:
  using ColumnHolderBase = detail::ColumnHolderBase;
  template <typename T>
  using ColumnHolder = detail::ColumnHolder<T>;

  std::size_t m_size{0};

  std::unordered_map<std::string, ColumnHolderBase *> m_allColumns;
  ParticleColumns m_knownColumns{ParticleColumns::None};
  std::unordered_map<std::string, std::unique_ptr<ColumnHolderBase>>
      m_dynamicColumns;

  std::vector<ParticleIndex> m_parentIndices;

  std::optional<ColumnHolder<std::uint32_t>> m_parentIndicesOffsetColumn;
  std::optional<ColumnHolder<std::uint8_t>> m_parentIndicesCountColumn;
  std::optional<ColumnHolder<Barcode>> m_barcodeColumn;
  std::optional<ColumnHolder<ProcessType>> m_processColumn;
  std::optional<ColumnHolder<Acts::PdgParticle>> m_pdgColumn;
  std::optional<ColumnHolder<double>> m_chargeColumn;
  std::optional<ColumnHolder<double>> m_massColumn;
  std::optional<ColumnHolder<Acts::Vector3>> m_directionColumn;
  std::optional<ColumnHolder<double>> m_absMomentumColumn;
  std::optional<ColumnHolder<Acts::Vector4>> m_position4Column;
  std::optional<ColumnHolder<double>> m_properTimeColumn;
  std::optional<ColumnHolder<double>> m_pathInX0Column;
  std::optional<ColumnHolder<double>> m_pathInL0Column;
  std::optional<ColumnHolder<std::uint32_t>> m_numberOfHitsColumn;
  std::optional<ColumnHolder<const Acts::Surface *>> m_referenceSurfaceColumn;
  std::optional<ColumnHolder<ParticleOutcome>> m_outcomeColumn;

  static auto knownColumnMasks() noexcept {
    using enum ParticleColumns;
    return std::tuple(Parents, Parents, Barcode, Process, Pdg, Charge, Mass,
                      Direction, AbsoluteMomentum, Position4, ProperTime,
                      PathInX0, PathInL0, NumberOfHits, ReferenceSurface,
                      Outcome);
  }

  static auto knownColumnNames() noexcept {
    return std::tuple("parentsOffset", "parentsCount", "barcode", "process",
                      "pdg", "charge", "mass", "direction", "absoluteMomentum",
                      "position4", "properTime", "pathInX0", "pathInL0",
                      "numberOfHits", "referenceSurface", "outcome");
  }

  static auto knownColumnDefaults() noexcept {
    return std::tuple(
        std::uint32_t{0}, std::uint8_t{0}, Barcode{}, ProcessType::eUndefined,
        Acts::PdgParticle::eInvalid, double{0}, double{0},
        Acts::Vector3(Acts::Vector3::Zero()), double{0},
        Acts::Vector4(Acts::Vector4::Zero()), double{0}, double{0}, double{0},
        std::uint32_t{0}, static_cast<const Acts::Surface *>(nullptr),
        ParticleOutcome::Alive);
  }

  template <typename Self>
  static auto knownColumns(Self &&self) noexcept {
    return std::tie(self.m_parentIndicesOffsetColumn,
                    self.m_parentIndicesCountColumn, self.m_barcodeColumn,
                    self.m_processColumn, self.m_pdgColumn, self.m_chargeColumn,
                    self.m_massColumn, self.m_directionColumn,
                    self.m_absMomentumColumn, self.m_position4Column,
                    self.m_properTimeColumn, self.m_pathInX0Column,
                    self.m_pathInL0Column, self.m_numberOfHitsColumn,
                    self.m_referenceSurfaceColumn, self.m_outcomeColumn);
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

/// A proxy class for accessing individual particles.
template <bool read_only>
class ParticleProxy {
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

  void assignParentIndices(std::span<const ParticleIndex> parentIndices)
    requires(!ReadOnly)
  {
    if (!m_container->m_sourceLinkOffsetColumn.has_value() ||
        !m_container->m_sourceLinkCountColumn.has_value()) {
      throw std::logic_error("No source links column available");
    }
    if (m_container->m_sourceLinkCountColumn->proxy(*this)[m_index] != 0) {
      throw std::logic_error(
          "Source links already assigned to the space point");
    }

    m_container->m_sourceLinkOffsetColumn->proxy(*this)[m_index] =
        static_cast<std::uint32_t>(m_container->m_sourceLinks.size());
    m_container->m_sourceLinkCountColumn->proxy(*this)[m_index] =
        static_cast<std::uint8_t>(parentIndices.size());
    m_container->m_sourceLinks.insert(m_container->m_sourceLinks.end(),
                                      parentIndices.begin(),
                                      parentIndices.end());
  }

  std::span<ParticleIndex> parentIndices() noexcept
    requires(!ReadOnly)
  {
    return {m_container->m_parentIndices.data() +
                m_container->m_parentIndicesOffsetColumn->at(m_index),
            m_container->m_parentIndicesCountColumn->at(m_index)};
  }

  ProcessType &process() noexcept
    requires(!ReadOnly)
  {
    return m_container->m_processColumn[m_index];
  }

  Acts::PdgParticle &pdg() noexcept
    requires(!ReadOnly)
  {
    return m_container->m_pdgColumn[m_index];
  }

  double &charge() noexcept
    requires(!ReadOnly)
  {
    return m_container->m_chargeColumn[m_index];
  }

  double &mass() noexcept
    requires(!ReadOnly)
  {
    return m_container->m_massColumn[m_index];
  }

  Acts::Vector3 &direction() noexcept
    requires(!ReadOnly)
  {
    return m_container->m_directionColumn[m_index];
  }

  double &absoluteMomentum() noexcept
    requires(!ReadOnly)
  {
    return m_container->m_absMomentumColumn[m_index];
  }

  Acts::Vector4 &position4() noexcept
    requires(!ReadOnly)
  {
    return m_container->m_position4Column[m_index];
  }

  double &properTime() noexcept
    requires(!ReadOnly)
  {
    return m_container->m_properTimeColumn[m_index];
  }

  double &pathInX0() noexcept
    requires(!ReadOnly)
  {
    return m_container->m_pathInX0Column[m_index];
  }

  double &pathInL0() noexcept
    requires(!ReadOnly)
  {
    return m_container->m_pathInL0Column[m_index];
  }

  std::uint32_t &numberOfHits() noexcept
    requires(!ReadOnly)
  {
    return m_container->m_numberOfHitsColumn[m_index];
  }

  const Acts::Surface *&referenceSurface() noexcept
    requires(!ReadOnly)
  {
    return m_container->m_referenceSurfaceColumn[m_index];
  }

  ParticleOutcome &outcome() noexcept
    requires(!ReadOnly)
  {
    return m_container->m_outcomeColumn[m_index];
  }

  std::span<const ParticleIndex> parentIndices() const noexcept {
    return {m_container->m_parentIndices.data() +
                m_container->m_parentIndicesOffsetColumn[m_index],
            m_container->m_parentIndicesCountColumn[m_index]};
  }

  const ProcessType &process() const noexcept {
    return m_container->m_processColumn[m_index];
  }

  const Acts::PdgParticle &pdg() const noexcept {
    return m_container->m_pdgColumn[m_index];
  }

  const double &charge() const noexcept {
    return m_container->m_chargeColumn[m_index];
  }

  const double &mass() const noexcept {
    return m_container->m_massColumn[m_index];
  }

  const Acts::Vector3 &direction() const noexcept {
    return m_container->m_directionColumn[m_index];
  }

  const double &absoluteMomentum() const noexcept {
    return m_container->m_absMomentumColumn[m_index];
  }

  const Acts::Vector4 &position4() const noexcept {
    return m_container->m_position4Column[m_index];
  }

  const double &properTime() const noexcept {
    return m_container->m_properTimeColumn[m_index];
  }

  const double &pathInX0() const noexcept {
    return m_container->m_pathInX0Column[m_index];
  }

  const double &pathInL0() const noexcept {
    return m_container->m_pathInL0Column[m_index];
  }

  const std::uint32_t &numberOfHits() const noexcept {
    return m_container->m_numberOfHitsColumn[m_index];
  }

  const Acts::Surface *const &referenceSurface() const noexcept {
    return m_container->m_referenceSurfaceColumn[m_index];
  }

  const ParticleOutcome &outcome() const noexcept {
    return m_container->m_outcomeColumn[m_index];
  }

  /// Const access to an extra column of data for the space point.
  /// @param column The extra column proxy to access.
  /// @return A const reference to the value in the extra column for the space point.
  template <typename T>
  const T &extra(const ConstParticleColumnProxy<T> &column) const noexcept {
    return column[m_index];
  }

 private:
  Container *m_container{nullptr};
  Index m_index{0};
};

/// Additional column of data that can be added to the space point container.
/// The column is indexed by the space point index.
template <typename T, bool read_only>
class ParticleColumnProxy {
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
  /// @param index The index of the space point to access.
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
  class Subset : public Acts::detail::ContainerSubset<Subset, Subset, Column,
                                                      Value, Index, ReadOnly> {
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

  friend class ParticleContainer;
};

}  // namespace ActsFatras

#include "ActsFatras/EventData/ParticleContainer.ipp"

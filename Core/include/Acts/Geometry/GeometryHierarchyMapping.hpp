#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <iostream>
#include <string_view>

#include <boost/container/flat_map.hpp>
#include <boost/container/small_vector.hpp>

namespace Acts {

using GeometryLevelIdentifier = std::uint32_t;
inline static constexpr std::size_t s_geometryHierarchyDepth = 10;
template <typename T>
using GeometryHierarchyVector =
    boost::container::small_vector<T, s_geometryHierarchyDepth>;

using GeometryHierarchyMapping = GeometryHierarchyVector<
    std::pair<GeometryLevelIdentifier, std::string_view>>;

using GeometryHierarchyMapper =
    Delegate<GeometryHierarchyMapping(GeometryIdentifier)>;

class SimpleGeometryHierarchyMapper {
 public:
  struct Level {
    std::uint8_t bits;
    std::string name;
  };

  explicit SimpleGeometryHierarchyMapper(GeometryHierarchyVector<Level> levels)
      : m_levels(std::move(levels)) {}

  GeometryHierarchyMapping operator()(GeometryIdentifier id) const {
    GeometryHierarchyMapping result;
    std::uint64_t value = id.value();
    for (const Level& level : m_levels) {
      GeometryLevelIdentifier levelId =
          (value >> (64 - level.bits)) & ((1 << level.bits) - 1);
      if (levelId == 0) {
        break;
      }
      result.push_back({levelId, level.name});
      value <<= level.bits;
    }
    return result;
  }

 private:
  GeometryHierarchyVector<Level> m_levels;
};

class NestedGeometryHierarchyMapper {
 public:
  struct Level {
    std::uint8_t bits;
    std::string name;
    boost::container::flat_map<GeometryLevelIdentifier, GeometryHierarchyMapper>
        mappers;
  };

  explicit NestedGeometryHierarchyMapper(Level level)
      : m_level(std::move(level)) {}

  GeometryHierarchyMapping operator()(GeometryIdentifier id) const {
    std::uint64_t value = id.value();
    GeometryLevelIdentifier levelId =
        (value >> (64 - m_level.bits)) & ((1 << m_level.bits) - 1);
    GeometryHierarchyMapping result;
    result.push_back({levelId, m_level.name});
    if (auto it = m_level.mappers.find(levelId); it != m_level.mappers.end()) {
      GeometryHierarchyMapping subResult = it->second(id);
      result.insert(result.end(), subResult.begin(), subResult.end());
    }
    return result;
  }

 private:
  Level m_level;
};

std::ostream& operator<<(std::ostream& os,
                         const GeometryHierarchyMapping& mapping) {
  for (std::size_t i = 0; i < mapping.size(); ++i) {
    const auto& [id, name] = mapping[i];
    os << (i == 0 ? "" : ",") << name << "=" << id;
  }
  return os;
}

}  // namespace Acts

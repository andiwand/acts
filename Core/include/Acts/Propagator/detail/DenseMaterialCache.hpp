// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"

namespace Acts::detail {

struct DenseMaterialCache {
  struct Entry {
    Material material;
    double pathLength{};
    double momentumInitial{};
    double momentumFinal{};
  };

  std::vector<Entry> entries;
  MaterialSlab accumulatedMaterial;

  void addEntry(const Material& material, double pathLength,
                double momentumInitial, double momentumFinal) {
    const auto& entry = entries.emplace_back(material, pathLength,
                                             momentumInitial, momentumFinal);
    accumulatedMaterial = MaterialSlab::combineLayers(
        accumulatedMaterial, MaterialSlab(entry.material, entry.pathLength));
  }

  bool empty() const { return entries.empty(); }

  double pathLength() const { return accumulatedMaterial.thickness(); }

  void clear() {
    entries.clear();
    accumulatedMaterial = MaterialSlab();
  }

  struct SampleIterator {
    const DenseMaterialCache* cache{};
    double maxXOverX0Step{};

    std::optional<std::size_t> entryIndex;
    std::size_t stepIndex{};
    std::size_t stepCount{};

    explicit SampleIterator(const DenseMaterialCache& cache_)
        : entryIndex(cache_.entries.size()) {}

    SampleIterator(const DenseMaterialCache& cache_, double maxXOverX0Step_)
        : cache(&cache_), maxXOverX0Step(maxXOverX0Step_) {}

    SampleIterator& operator++() {
      if (entryIndex >= cache->entries.size()) {
        return *this;
      }

      if (!entryIndex.has_value()) {
        entryIndex = 0;
      } else if (stepIndex + 1 < stepCount) {
        ++stepIndex;
        return *this;
      } else {
        ++(*entryIndex);
      }

      const auto& entry = cache->entries[*entryIndex];
      MaterialSlab slab(entry.material, entry.pathLength);

      stepIndex = 0;
      stepCount = entry.material.isValid()
                      ? static_cast<std::size_t>(
                            std::ceil(slab.thicknessInX0() / maxXOverX0Step))
                      : 1;

      return *this;
    }

    bool operator!=(const SampleIterator& other) const {
      return entryIndex != other.entryIndex || stepIndex != other.stepIndex;
    }

    Entry operator*() const {
      const Entry& entry = cache->entries[entryIndex.value()];

      if (stepCount == 1) {
        return entry;
      }

      Entry result;
      result.material = entry.material;
      result.pathLength = entry.pathLength * (1.0 / stepCount);
      result.momentumInitial = entry.momentumInitial +
                               (entry.momentumFinal - entry.momentumInitial) *
                                   (static_cast<double>(stepIndex) / stepCount);
      result.momentumFinal =
          entry.momentumFinal -
          (entry.momentumFinal - entry.momentumInitial) *
              (1.0 - static_cast<double>(stepIndex + 1) / stepCount);

      return result;
    }
  };

  struct SampleRange {
    const DenseMaterialCache* cache{};
    double maxXOverX0Step{};

    SampleRange(const DenseMaterialCache& cache_, double maxXOverX0Step_)
        : cache(&cache_), maxXOverX0Step(maxXOverX0Step_) {}

    SampleIterator begin() const {
      SampleIterator preBegin(*cache, maxXOverX0Step);
      ++preBegin;
      return preBegin;
    }
    SampleIterator end() const { return SampleIterator(*cache); }
  };

  SampleRange samples(double maxXOverX0Step) const {
    return SampleRange(*this, maxXOverX0Step);
  }
};

}  // namespace Acts::detail

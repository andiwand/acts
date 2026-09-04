// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/Rz/RzMeasurementGrid.hpp"

#include <algorithm>
#include <cmath>

namespace Acts::Experimental {

RzMeasurementGrid::RzMeasurementGrid(const RzLayout& layout)
    : m_layout(&layout) {
  m_binStart.assign(layout.totalBins + 1, 0);
  m_layerHasStrips.assign(layout.layers.size(), false);
}

void RzMeasurementGrid::clear() {
  m_entries.clear();
  m_binOf.clear();
  m_order.clear();
  std::ranges::fill(m_binStart, 0u);
  m_layerHasStrips.assign(m_layerHasStrips.size(), false);
}

void RzMeasurementGrid::add(std::uint32_t module, std::uint8_t dim,
                            std::span<const std::uint8_t> localIndices,
                            std::span<const double> localParams,
                            std::span<const double> localCov,
                            std::uint32_t source) {
  const RzModule& m = m_layout->modules[module];
  RzMeasurement e;
  e.u = m.u;
  e.v = m.v;
  e.normal = m.normal;
  e.halfV = m.halfV;
  e.maxDistance = m_layout->layers[m.layer].moduleDistance;
  e.module = module;
  e.source = source;
  e.dim = dim;

  double lu = 0.;
  double lv = 0.;
  if (dim == 2) {
    // whichever order the caller measures in, the frame is (u, v)
    const bool swapped = localIndices[0] == 1;
    lu = localParams[swapped ? 1 : 0];
    lv = localParams[swapped ? 0 : 1];
    e.cov00 = localCov[swapped ? 3 : 0];
    e.cov11 = localCov[swapped ? 0 : 3];
    e.cov01 = localCov[1];
  } else {
    if (localIndices[0] == 1) {
      // a strip measuring v: swap the frame so that u is what it measures
      std::swap(e.u, e.v);
      e.normal = -e.normal;
    }
    lu = localParams[0];
    e.cov00 = localCov[0];
  }
  e.position = m.center + lu * e.u + lv * e.v;

  const RzLayer& layer = m_layout->layers[m.layer];
  if (dim == 1) {
    m_layerHasStrips[m.layer] = true;
  }
  const RzSurface& surface = m_layout->surfaces[layer.surface];
  const double phi = std::atan2(e.position.y(), e.position.x());
  const double along = surface.shape == RzShape::Cylinder
                           ? e.position.z()
                           : std::hypot(e.position.x(), e.position.y());
  const std::uint32_t bin = m_layout->bin(m.layer, phi, along);

  m_entries.push_back(e);
  m_binOf.push_back(bin);
}

void RzMeasurementGrid::add(std::uint32_t module, std::uint8_t dim,
                            const Vector3& position, const Vector3& u,
                            const Vector3& v, double cov00, double cov01,
                            double cov11, std::uint32_t source) {
  const RzModule& m = m_layout->modules[module];
  RzMeasurement e;
  e.position = position;
  e.u = u;
  e.v = v;
  e.normal = m.normal;
  // the room along a strip: the module's box in either direction, the
  // strip's own direction not being a module axis in general
  e.halfV = std::max(m.halfU, m.halfV);
  e.maxDistance = m_layout->layers[m.layer].moduleDistance;
  e.module = module;
  e.source = source;
  e.dim = dim;
  e.cov00 = cov00;
  e.cov01 = cov01;
  e.cov11 = cov11;

  const RzLayer& layer = m_layout->layers[m.layer];
  if (dim == 1) {
    m_layerHasStrips[m.layer] = true;
  }
  const RzSurface& surface = m_layout->surfaces[layer.surface];
  const double phi = std::atan2(e.position.y(), e.position.x());
  const double along = surface.shape == RzShape::Cylinder
                           ? e.position.z()
                           : std::hypot(e.position.x(), e.position.y());
  const std::uint32_t bin = m_layout->bin(m.layer, phi, along);
  m_entries.push_back(e);
  m_binOf.push_back(bin);
}

void RzMeasurementGrid::finalize() {
  // counting sort into the bins
  std::ranges::fill(m_binStart, 0u);
  for (const std::uint32_t b : m_binOf) {
    ++m_binStart[b + 1];
  }
  for (std::size_t i = 1; i < m_binStart.size(); ++i) {
    m_binStart[i] += m_binStart[i - 1];
  }
  // the entries themselves go into bin order, so that a window walks memory
  // linearly; the indices handed out from here on are the sorted ones
  std::vector<RzMeasurement> sorted(m_entries.size());
  std::vector<std::uint32_t> fill(m_binStart.begin(), m_binStart.end() - 1);
  for (std::uint32_t i = 0; i < m_binOf.size(); ++i) {
    sorted[fill[m_binOf[i]]++] = m_entries[i];
  }
  m_entries.swap(sorted);
  m_order.resize(m_entries.size());
  for (std::uint32_t i = 0; i < m_order.size(); ++i) {
    m_order[i] = i;
  }
}

}  // namespace Acts::Experimental

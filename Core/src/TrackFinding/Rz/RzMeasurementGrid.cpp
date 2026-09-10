// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/Rz/RzMeasurementGrid.hpp"

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <algorithm>
#include <cmath>

namespace Acts::Experimental {

RzMeasurementGrid::RzMeasurementGrid(const RzLayout& layout)
    : m_layout(&layout) {
  m_moduleStart.assign(layout.modules.size() + 1, 0);
}

void RzMeasurementGrid::clear() {
  m_entries.clear();
  m_frames.clear();
  m_moduleOf.clear();
  std::ranges::fill(m_moduleStart, 0u);
  m_grouped = true;
  m_lastModule = 0;
}

void RzMeasurementGrid::reserve(std::size_t n) {
  m_entries.reserve(n);
  m_moduleOf.reserve(n);
}

/// Grow the frame table to match the entries, so that a polar module's frames
/// stay parallel to them however few polar modules the event has
void RzMeasurementGrid::padFrames() {
  if (!m_frames.empty()) {
    m_frames.resize(m_entries.size());
  }
}

std::uint32_t RzMeasurementGrid::add(std::uint32_t module,
                                     const RzMeasurement& measurement,
                                     const RzMeasurementFrame& frame) {
  if (module < m_lastModule) {
    m_grouped = false;
  }
  m_lastModule = module;
  const std::uint32_t index = m_moduleStart[module + 1]++;
  m_entries.push_back(measurement);
  m_moduleOf.push_back(module);
  if (m_layout->modules[module].polar) {
    m_frames.resize(m_entries.size());
    m_frames.back() = frame;
  }
  return index;
}

void RzMeasurementGrid::addRange(std::uint32_t module,
                                 std::span<const RzMeasurement> measurements,
                                 std::span<const RzMeasurementFrame> frames) {
  if (module < m_lastModule) {
    m_grouped = false;
  }
  m_lastModule = module;
  m_moduleStart[module + 1] += static_cast<std::uint32_t>(measurements.size());
  m_entries.insert(m_entries.end(), measurements.begin(), measurements.end());
  m_moduleOf.insert(m_moduleOf.end(), measurements.size(), module);
  if (!frames.empty()) {
    m_frames.resize(m_entries.size() - measurements.size());
    m_frames.insert(m_frames.end(), frames.begin(), frames.end());
  }
}

RzMeasurement RzMeasurementGrid::fromBound(
    const RzLayout& layout, std::uint32_t module, const Surface& surface,
    const GeometryContext& gctx, std::uint8_t dim,
    std::span<const std::uint8_t> boundIndices,
    std::span<const double> boundParams, std::span<const double> boundCov,
    std::uint32_t source, RzMeasurementFrame& frame) {
  const RzModule& m = layout.modules[module];

  RzMeasurement e;
  e.source = source;
  if (dim == 2) {
    e.projector = RzProjector::Both;
  } else {
    e.projector = boundIndices[0] == 1 ? RzProjector::Loc1 : RzProjector::Loc0;
  }

  // what was measured where it was measured, the rest at the module centre
  Vector2 local = m.boundCenter;
  for (std::uint8_t i = 0; i < dim; ++i) {
    local[boundIndices[i]] = boundParams[i];
  }
  // whichever order the caller measures in, the variance of a coordinate
  // sits in that coordinate's own slot
  const bool measuresLoc1 = dim == 1 && boundIndices[0] == 1;
  const bool swapped = dim == 2 && boundIndices[0] == 1;
  double var0 = 0.;
  double var1 = 0.;
  double cov01 = 0.;
  if (dim == 2) {
    var0 = boundCov[swapped ? 3 : 0];
    var1 = boundCov[swapped ? 0 : 3];
    cov01 = boundCov[1];
  } else if (measuresLoc1) {
    var1 = boundCov[0];
  } else {
    var0 = boundCov[0];
  }

  if (!m.polar) {
    // A cartesian frame is the module's own, which the layout read off the
    // surface once: the offset from the centre is the difference of the bound
    // coordinates, the map to the module axes is the identity, and the
    // variance is already a length squared. Nothing here touches the surface,
    // which is what makes this path a copy — and it is every pixel and every
    // barrel strip.
    const Vector2 offset = local - m.boundCenter;
    e.loc0 = offset.x();
    e.loc1 = offset.y();
    e.cov00 = var0;
    e.cov01 = cov01;
    e.cov11 = var1;
    return e;
  }

  // d(global) / d(bound local): the surface's own rotation, composed with the
  // bounds' map to the cartesian frame. Its columns are the directions the two
  // bound coordinates move the point in, and their lengths are what turns a
  // variance in the bound coordinates into one in length units. The frame is
  // the measurement's own, because it turns with the strip.
  const Vector3 position = surface.localToGlobal(gctx, local, m.normal);
  Eigen::Matrix<double, 3, 2> jac =
      surface.localToGlobalTransform(gctx).rotation().leftCols<2>();
  jac *= surface.bounds().boundToCartesianJacobian(local);
  const double scale0 = jac.col(0).norm();
  const double scale1 = jac.col(1).norm();
  frame.u = jac.col(0) / scale0;
  frame.v = jac.col(1) / scale1;
  frame.normal = frame.u.cross(frame.v);
  const Vector3 d = position - m.center;
  e.loc0 = frame.u.dot(d);
  e.loc1 = frame.v.dot(d);
  e.cov00 = var0 * scale0 * scale0;
  e.cov01 = cov01 * scale0 * scale1;
  e.cov11 = var1 * scale1 * scale1;
  // The lever arm the azimuth was converted with: the distance from the polar
  // frame's origin, which is where the entry sits. Only a strip measuring the
  // azimuth has one — `loc0` is a radius, already a length, and a pixel that
  // measures both is placed by them together.
  // CONTROL: reproduce the old expression exactly
  const double oldScaleU = measuresLoc1 ? scale1 : scale0;
  e.invLever = oldScaleU > 0. ? 1. / oldScaleU : 0.;
  return e;
}

std::uint32_t RzMeasurementGrid::addBound(
    std::uint32_t module, const Surface& surface, const GeometryContext& gctx,
    std::uint8_t dim, std::span<const std::uint8_t> boundIndices,
    std::span<const double> boundParams, std::span<const double> boundCov,
    std::uint32_t source) {
  RzMeasurementFrame frame;
  const RzMeasurement e =
      fromBound(*m_layout, module, surface, gctx, dim, boundIndices,
                boundParams, boundCov, source, frame);
  return add(module, e, frame);
}

void RzMeasurementGrid::finalize() {
  // `m_moduleStart[i + 1]` counts what module `i` holds; the prefix sum turns
  // it into where module `i` starts
  for (std::size_t i = 1; i < m_moduleStart.size(); ++i) {
    m_moduleStart[i] += m_moduleStart[i - 1];
  }
  if (m_grouped) {
    // the caller added in module order, so the entries already sit where the
    // offsets say and there is nothing to move
    padFrames();
    return;
  }
  // counting sort by module, permuting the entries themselves so that a
  // module's measurements are contiguous: the search reads one module at a
  // time and nothing else, and an entry is small enough that moving it costs
  // less than the indirection would
  padFrames();
  std::vector<RzMeasurement> sorted(m_entries.size());
  std::vector<RzMeasurementFrame> sortedFrames(m_frames.size());
  std::vector<std::uint32_t> fill(m_moduleStart.begin(),
                                  m_moduleStart.end() - 1);
  const bool withFrames = !m_frames.empty();
  for (std::uint32_t i = 0; i < m_moduleOf.size(); ++i) {
    const std::uint32_t to = fill[m_moduleOf[i]]++;
    sorted[to] = m_entries[i];
    if (withFrames) {
      sortedFrames[to] = m_frames[i];
    }
  }
  m_entries.swap(sorted);
  m_frames.swap(sortedFrames);
}

}  // namespace Acts::Experimental

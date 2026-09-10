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
  m_blocks.assign(layout.modules.size(), ModuleBlock{});
}

void RzMeasurementGrid::clear() {
  m_entries.clear();
  m_frames.clear();
  m_moduleOf.clear();
  std::ranges::fill(m_blocks, ModuleBlock{});
  m_contiguous = true;
}

void RzMeasurementGrid::reserve(std::size_t n) {
  m_entries.reserve(n);
  m_moduleOf.reserve(n);
}

/// Grow the frame table to match the entries, so that a polar module's frames
/// stay parallel to them however few polar modules the event has
std::uint32_t RzMeasurementGrid::add(std::uint32_t module,
                                     const RzMeasurement& measurement,
                                     const RzMeasurementFrame& frame) {
  ModuleBlock& block = m_blocks[module];
  const std::uint32_t at = static_cast<std::uint32_t>(m_entries.size());
  const bool polar = m_layout->modules[module].polar;
  if (block.size == 0) {
    block.begin = at;
    if (polar) {
      block.frame = static_cast<std::uint32_t>(m_frames.size());
    }
  } else if (block.begin + block.size != at) {
    // this module's measurements were interrupted by another module's, so
    // one run no longer holds them and `finalize` has to group them
    m_contiguous = false;
  }
  const std::uint32_t index = block.size++;
  m_entries.push_back(measurement);
  m_moduleOf.push_back(module);
  if (polar) {
    m_frames.push_back(frame);
  }
  return index;
}

void RzMeasurementGrid::addRange(std::uint32_t module,
                                 std::span<const RzMeasurement> measurements,
                                 std::span<const RzMeasurementFrame> frames) {
  ModuleBlock& block = m_blocks[module];
  const std::uint32_t at = static_cast<std::uint32_t>(m_entries.size());
  const bool polar = m_layout->modules[module].polar;
  if (block.size == 0) {
    block.begin = at;
    if (polar) {
      block.frame = static_cast<std::uint32_t>(m_frames.size());
    }
  } else if (block.begin + block.size != at) {
    m_contiguous = false;
  }
  block.size += static_cast<std::uint32_t>(measurements.size());
  m_entries.insert(m_entries.end(), measurements.begin(), measurements.end());
  m_moduleOf.insert(m_moduleOf.end(), measurements.size(), module);
  if (polar) {
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
  if (m_contiguous) {
    // every module's measurements arrived in one run, so the blocks already
    // point at them and there is nothing to move — which is the case for any
    // caller whose container is grouped by module, in whatever order it
    // visits the modules
    return;
  }
  // Counting sort by module, permuting the entries themselves so that a
  // module's measurements are contiguous: the search reads one module at a
  // time and nothing else, and an entry is small enough that moving it costs
  // less than the indirection would. The blocks are laid out in module order
  // here, which the search does not need but nothing forbids.
  std::vector<std::uint32_t> fill(m_blocks.size());
  std::vector<std::uint32_t> fillFrame(m_blocks.size());
  std::uint32_t at = 0;
  std::uint32_t atFrame = 0;
  for (std::size_t m = 0; m < m_blocks.size(); ++m) {
    ModuleBlock& block = m_blocks[m];
    block.begin = at;
    at += block.size;
    if (m_layout->modules[m].polar && block.size != 0) {
      block.frame = atFrame;
      atFrame += block.size;
    }
    fill[m] = block.begin;
    fillFrame[m] = block.frame;
  }
  std::vector<RzMeasurement> sorted(m_entries.size());
  std::vector<RzMeasurementFrame> sortedFrames(m_frames.size());
  // the frames were added in the order the polar entries were, so walking the
  // entries in that order walks the frames too
  std::uint32_t from = 0;
  for (std::uint32_t i = 0; i < m_moduleOf.size(); ++i) {
    const std::uint32_t m = m_moduleOf[i];
    sorted[fill[m]++] = m_entries[i];
    if (m_layout->modules[m].polar) {
      sortedFrames[fillFrame[m]++] = m_frames[from++];
    }
  }
  m_entries.swap(sorted);
  m_frames.swap(sortedFrames);
  m_contiguous = true;
}

}  // namespace Acts::Experimental

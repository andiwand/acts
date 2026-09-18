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

RzMeasurementGrid::ModuleBlock& RzMeasurementGrid::prepareModule(
    std::uint32_t module) {
  ModuleBlock& block = m_blocks[module];
  const std::uint32_t at = static_cast<std::uint32_t>(m_entries.size());
  if (block.size == 0) {
    block.begin = at;
    if (m_layout->modules[module].polar) {
      block.frame = static_cast<std::uint32_t>(m_frames.size());
    }
  } else if (block.begin + block.size != at) {
    // Interleaved modules need regrouping in finalize().
    m_contiguous = false;
  }
  return block;
}

std::uint32_t RzMeasurementGrid::add(std::uint32_t module,
                                     const RzMeasurement& measurement,
                                     const RzMeasurementFrame& frame) {
  ModuleBlock& block = prepareModule(module);
  const std::uint32_t index = block.size++;
  m_entries.push_back(measurement);
  m_moduleOf.push_back(module);
  if (m_layout->modules[module].polar) {
    m_frames.push_back(frame);
  }
  return index;
}

void RzMeasurementGrid::addRange(std::uint32_t module,
                                 std::span<const RzMeasurement> measurements,
                                 std::span<const RzMeasurementFrame> frames) {
  ModuleBlock& block = prepareModule(module);
  block.size += static_cast<std::uint32_t>(measurements.size());
  m_entries.insert(m_entries.end(), measurements.begin(), measurements.end());
  m_moduleOf.insert(m_moduleOf.end(), measurements.size(), module);
  if (m_layout->modules[module].polar) {
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
  // Map covariance entries to their coordinate slots, regardless of order.
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
    // Cartesian measurements already use the module axes and length units.
    const Vector2 offset = local - m.boundCenter;
    e.loc0 = offset.x();
    e.loc1 = offset.y();
    e.cov00 = var0;
    e.cov01 = cov01;
    e.cov11 = var1;
    return e;
  }

  // Polar axes depend on the hit; convert covariance to length units.
  double scale0 = 0.;
  double scale1 = 0.;
  if (m.polarIsPlain) {
    // The layout verified this plain polar map against the surface.
    const double r = local[0];
    const double c = std::cos(local[1]);
    const double sn = std::sin(local[1]);
    frame.u = c * m.u + sn * m.v;
    frame.v = -sn * m.u + c * m.v;
    frame.uU = c;
    frame.uV = sn;
    frame.vU = -sn;
    frame.vV = c;
    scale0 = 1.;
    scale1 = r;
    // the offset from the module centre, taken on the frame's own axes
    e.loc0 = r - (c * m.localCenter.x() + sn * m.localCenter.y());
    e.loc1 = sn * m.localCenter.x() - c * m.localCenter.y();
  } else {
    // Jacobian columns give the measurement axes and their length scales.
    const Vector3 position = surface.localToGlobal(gctx, local, m.normal);
    Eigen::Matrix<double, 3, 2> jac =
        surface.localToGlobalTransform(gctx).rotation().leftCols<2>();
    jac *= surface.bounds().boundToCartesianJacobian(local);
    scale0 = jac.col(0).norm();
    scale1 = jac.col(1).norm();
    frame.u = jac.col(0) / scale0;
    frame.v = jac.col(1) / scale1;
    frame.uU = frame.u.dot(m.u);
    frame.uV = frame.u.dot(m.v);
    frame.vU = frame.v.dot(m.u);
    frame.vV = frame.v.dot(m.v);
    const Vector3 d = position - m.center;
    e.loc0 = frame.u.dot(d);
    e.loc1 = frame.v.dot(d);
  }
  e.cov00 = var0 * scale0 * scale0;
  e.cov01 = cov01 * scale0 * scale1;
  e.cov11 = var1 * scale1 * scale1;
  // Inverse angular scale used to convert variance to length units.
  e.invLever =
      e.projector != RzProjector::Loc0 && scale1 > 0. ? 1. / scale1 : 0.;
  return e;
}

std::uint32_t RzMeasurementGrid::addBound(
    std::uint32_t module, const Surface& surface, const GeometryContext& gctx,
    std::uint8_t dim, std::span<const std::uint8_t> boundIndices,
    std::span<const double> boundParams, std::span<const double> boundCov,
    std::uint32_t source, double time, double timeVariance) {
  RzMeasurementFrame frame;
  RzMeasurement e =
      fromBound(*m_layout, module, surface, gctx, dim, boundIndices,
                boundParams, boundCov, source, frame);
  e.time = time;
  e.timeVariance = timeVariance;
  return add(module, e, frame);
}

void RzMeasurementGrid::finalize() {
  if (m_contiguous) {
    // Each module already occupies one contiguous run.
    return;
  }
  // Counting-sort entries and frames into contiguous module blocks.
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
  // Frames follow the insertion order of polar entries.
  std::uint32_t from = 0;
  for (std::uint32_t i = 0; i < m_moduleOf.size(); ++i) {
    const std::uint32_t m = m_moduleOf[i];
    sorted[fill[m]++] = m_entries[i];
    if (m_layout->modules[m].polar) {
      sortedFrames[fillFrame[m]++] = m_frames[from++];
    }
  }
  // Keep the module tags in the same order for a subsequent append/finalize.
  for (std::uint32_t m = 0; m < m_blocks.size(); ++m) {
    const ModuleBlock& block = m_blocks[m];
    std::fill_n(m_moduleOf.begin() + block.begin, block.size, m);
  }
  m_entries.swap(sorted);
  m_frames.swap(sortedFrames);
  m_contiguous = true;
}

}  // namespace Acts::Experimental

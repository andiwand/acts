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
  m_moduleOf.clear();
  m_order.clear();
  std::ranges::fill(m_moduleStart, 0u);
}

void RzMeasurementGrid::reserve(std::size_t n) {
  m_entries.reserve(n);
  m_moduleOf.reserve(n);
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

  m_entries.push_back(e);
  m_moduleOf.push_back(module);
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

  m_entries.push_back(e);
  m_moduleOf.push_back(module);
}

void RzMeasurementGrid::addBound(std::uint32_t module, const Surface& surface,
                                 const GeometryContext& gctx, std::uint8_t dim,
                                 std::span<const std::uint8_t> boundIndices,
                                 std::span<const double> boundParams,
                                 std::span<const double> boundCov,
                                 std::uint32_t source) {
  const RzModule& m = m_layout->modules[module];
  // what was measured where it was measured, the rest at the module centre
  Vector2 local = m.boundCenter;
  for (std::uint8_t i = 0; i < dim; ++i) {
    local[boundIndices[i]] = boundParams[i];
  }

  const bool measuresLoc1 = dim == 1 && boundIndices[0] == 1;
  Vector3 position;
  Vector3 u;
  Vector3 v;
  double scaleU = 1.;
  double scaleV = 1.;
  if (!m.polar) {
    // A cartesian frame is the module's own, which was read off the surface
    // once when the layout was built: the point is the centre plus the
    // measured offsets along the module axes, the map to the global frame is
    // those axes, and it is unit. Worth the branch — this is every pixel and
    // every barrel strip, and the general form below costs four virtual calls
    // into the surface for each of them.
    const Vector2 offset = local - m.boundCenter;
    position = m.center + offset.x() * m.u + offset.y() * m.v;
    u = measuresLoc1 ? m.v : m.u;
    v = measuresLoc1 ? m.u : m.v;
  } else {
    position = surface.localToGlobal(gctx, local, m.normal);
    // d(global) / d(bound local): the surface's own rotation, composed with
    // the bounds' map to the cartesian frame. Its columns are the directions
    // the two bound coordinates move the point in, and their lengths are what
    // turns a covariance in the bound coordinates into one in length units.
    Eigen::Matrix<double, 3, 2> jac =
        surface.localToGlobalTransform(gctx).rotation().leftCols<2>();
    jac *= surface.bounds().boundToCartesianJacobian(local);
    const double scale0 = jac.col(0).norm();
    const double scale1 = jac.col(1).norm();
    const std::uint8_t iu = measuresLoc1 ? 1 : 0;
    scaleU = measuresLoc1 ? scale1 : scale0;
    scaleV = measuresLoc1 ? scale0 : scale1;
    u = jac.col(iu) / scaleU;
    v = jac.col(1 - iu) / scaleV;
  }

  double cov00 = 0.;
  double cov01 = 0.;
  double cov11 = 0.;
  if (dim == 2) {
    // whichever order the caller measures in, the frame is (u, v)
    const bool swapped = boundIndices[0] == 1;
    cov00 = boundCov[swapped ? 3 : 0] * scaleU * scaleU;
    cov11 = boundCov[swapped ? 0 : 3] * scaleV * scaleV;
    cov01 = boundCov[1] * scaleU * scaleV;
  } else {
    cov00 = boundCov[0] * scaleU * scaleU;
  }

  RzMeasurement e;
  e.position = position;
  e.u = u;
  e.v = v;
  e.normal = u.cross(v);
  e.invLever = m.polar && scaleU > 0. ? 1. / scaleU : 0.;
  // the room a search opens along a strip: the module's extent along the
  // coordinate it does not measure, and for a polar frame, where neither
  // bound coordinate is a module axis, the box in either direction
  e.halfV = !m.polar ? (measuresLoc1 ? m.halfU : m.halfV)
                     : std::max(m.halfU, m.halfV);
  e.maxDistance = m_layout->layers[m.layer].moduleDistance;
  e.module = module;
  e.source = source;
  e.dim = dim;
  e.cov00 = cov00;
  e.cov01 = cov01;
  e.cov11 = cov11;

  m_entries.push_back(e);
  m_moduleOf.push_back(module);
}

void RzMeasurementGrid::finalize() {
  // counting sort by module. Only the index is sorted: the entries stay where
  // they were added, which on a detector with a million measurements an event
  // is worth more than the locality a full copy would buy — a search touches
  // one module at a time, and a module's measurements arrive together anyway
  // because the caller's container is grouped by module.
  std::ranges::fill(m_moduleStart, 0u);
  for (const std::uint32_t module : m_moduleOf) {
    ++m_moduleStart[module + 1];
  }
  for (std::size_t i = 1; i < m_moduleStart.size(); ++i) {
    m_moduleStart[i] += m_moduleStart[i - 1];
  }
  m_order.resize(m_entries.size());
  std::vector<std::uint32_t> fill(m_moduleStart.begin(),
                                  m_moduleStart.end() - 1);
  for (std::uint32_t i = 0; i < m_moduleOf.size(); ++i) {
    m_order[fill[m_moduleOf[i]]++] = i;
  }
}

}  // namespace Acts::Experimental

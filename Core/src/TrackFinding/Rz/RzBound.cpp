// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/Rz/RzBound.hpp"

#include <cmath>

namespace Acts::Experimental {

BoundMatrix rzBoundCovariance(const RzFreeToBoundMatrix& j, const RzMatrix& c) {
  // Coefficient-based products, asked for by name: Eigen's size heuristics
  // can send a 6x8 by 8x8 through its blocked GEMM, with the packing that
  // costs more than the product, and a scalar loop of length seven is a
  // latency chain. The time row of J is zero, so its row and column of the
  // result are zero and the time variance is set by hand.
  const Eigen::Matrix<double, eBoundSize, eRzSize> jc = j.lazyProduct(c);
  BoundMatrix out = jc.lazyProduct(j.transpose());
  out(eBoundTime, eBoundTime) = 1.;
  return out;
}

void rzFillDirectionRows(const Vector3& direction, RzFreeToBoundMatrix& j) {
  // the same derivatives the surfaces use: phi = atan2(dy, dx) and
  // theta = atan2(perp, dz), written against the unit direction
  const double dx = direction.x();
  const double dy = direction.y();
  const double dz = direction.z();
  const double sinTheta = std::sqrt(std::max(0., 1. - dz * dz));
  const double invSinTheta = sinTheta > 0. ? 1. / sinTheta : 0.;
  const double cosPhi = dx * invSinTheta;
  const double sinPhi = dy * invSinTheta;
  j(eBoundPhi, eRzDir0) = -sinPhi * invSinTheta;
  j(eBoundPhi, eRzDir1) = cosPhi * invSinTheta;
  j(eBoundPhi, eRzDir2) = 0.;
  j(eBoundTheta, eRzDir0) = cosPhi * dz;
  j(eBoundTheta, eRzDir1) = sinPhi * dz;
  j(eBoundTheta, eRzDir2) = -sinTheta;
  j(eBoundQOverP, eRzQOverP) = 1.;
}

std::optional<RzBoundState> rzBoundOnModule(const RzModule& module,
                                            const RzVector& v,
                                            const RzMatrix& c) {
  if (module.polar && !module.polarIsPlain) {
    return std::nullopt;
  }
  const Vector3 position = v.segment<3>(eRzPos0);
  const Vector3 direction = v.segment<3>(eRzDir0);
  RzBoundState out;
  RzFreeToBoundMatrix j = RzFreeToBoundMatrix::Zero();
  // the two position rows: the module axes for a cartesian module; for a
  // polar one the radial and the azimuthal direction about the surface
  // origin, the latter over the radius, as the disc surface writes them
  Vector3 row0;
  Vector3 row1;
  if (!module.polar) {
    const Vector3 d = position - module.center;
    out.parameters[eBoundLoc0] = module.u.dot(d);
    out.parameters[eBoundLoc1] = module.v.dot(d);
    row0 = module.u;
    row1 = module.v;
  } else {
    const Vector3 origin = module.center - module.localCenter.x() * module.u -
                           module.localCenter.y() * module.v;
    const Vector3 d = position - origin;
    const double x = module.u.dot(d);
    const double y = module.v.dot(d);
    const double r = std::sqrt(x * x + y * y);
    out.parameters[eBoundLoc0] = r;
    out.parameters[eBoundLoc1] = std::atan2(y, x);
    const double invR = r > 0. ? 1. / r : 0.;
    const double cosPhi = x * invR;
    const double sinPhi = y * invR;
    row0 = cosPhi * module.u + sinPhi * module.v;
    row1 = (cosPhi * module.v - sinPhi * module.u) * invR;
  }
  for (unsigned int k = 0; k < 3; ++k) {
    j(eBoundLoc0, eRzPos0 + k) = row0[k];
    j(eBoundLoc1, eRzPos0 + k) = row1[k];
  }
  rzFillDirectionRows(direction, j);
  out.parameters[eBoundTime] = 0.;
  out.parameters[eBoundPhi] = std::atan2(direction.y(), direction.x());
  out.parameters[eBoundTheta] = std::atan2(
      std::sqrt(direction.x() * direction.x() + direction.y() * direction.y()),
      direction.z());
  out.parameters[eBoundQOverP] = v[eRzQOverP];
  out.covariance = rzBoundCovariance(j, c);
  return out;
}

}  // namespace Acts::Experimental

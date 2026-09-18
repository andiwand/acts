// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/Rz/RzBound.hpp"

#include "Acts/Utilities/MathHelpers.hpp"

#include <cmath>

namespace Acts::Experimental {

BoundMatrix rzBoundCovariance(const RzFreeToBoundMatrix& j, const RzMatrix& c) {
  // Avoid blocked GEMM packing for these small matrices.
  // The spatial Jacobian has a zero time row; supply its variance below.
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
  const double sinTheta = fastHypot(dx, dy);
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

std::optional<RzBoundState> rzBoundOnModule(
    const RzModule& module, const RzVector& v, const RzMatrix& c,
    const RzHelix::StepJacobian* transport) {
  if (module.polar && !module.polarIsPlain) {
    return std::nullopt;
  }
  const Vector3 position = v.segment<3>(eRzPos0);
  const Vector3 direction = v.segment<3>(eRzDir0);
  RzBoundState out;
  RzFreeToBoundMatrix j = RzFreeToBoundMatrix::Zero();
  // Position rows use Cartesian axes or the local polar Jacobian.
  Vector3 row0;
  Vector3 row1;
  if (!module.polar) {
    const Vector3 d = position - module.center;
    out.parameters[eBoundLoc0] = module.localCenter.x() + module.u.dot(d);
    out.parameters[eBoundLoc1] = module.localCenter.y() + module.v.dot(d);
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
  for (std::uint32_t k = 0; k < 3; ++k) {
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
  if (transport != nullptr) {
    // Compose the sparse helix Jacobian with the bound rows before applying
    // the covariance. There is no intermediate 7x7 transported covariance.
    const auto& t = *transport;
    RzFreeToBoundMatrix composed;
    for (std::uint32_t r = 0; r < eBoundSize; ++r) {
      const double w = j.row(r).dot(t.d);
      composed(r, eRzPos0) = j(r, eRzPos0);
      composed(r, eRzPos1) = j(r, eRzPos1);
      composed(r, eRzPos2) = j(r, eRzPos2);
      composed(r, eRzDir0) = t.f1 * j(r, eRzPos0) - t.f2 * j(r, eRzPos1) +
                             t.cs * j(r, eRzDir0) - t.sn * j(r, eRzDir1);
      composed(r, eRzDir1) = t.f2 * j(r, eRzPos0) + t.f1 * j(r, eRzPos1) +
                             t.sn * j(r, eRzDir0) + t.cs * j(r, eRzDir1);
      composed(r, eRzDir2) = t.s * j(r, eRzPos2) + j(r, eRzDir2);
      composed(r, eRzQOverP) = t.a1 * j(r, eRzPos0) + t.a2 * j(r, eRzPos1) +
                               t.a3 * j(r, eRzPos2) + t.b1 * j(r, eRzDir0) +
                               t.b2 * j(r, eRzDir1) + t.b3 * j(r, eRzDir2) +
                               j(r, eRzQOverP);
      composed.row(r) += w * t.dsdv.transpose();
    }
    j = composed;
  }
  out.covariance = rzBoundCovariance(j, c);
  return out;
}

}  // namespace Acts::Experimental

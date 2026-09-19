// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFinding/Rz/RzTransport.hpp"

#include <cmath>

#include "RzMaterial.hpp"

namespace Acts::Experimental::detail::rz {

struct State {
  RzVector v;
  RzMatrix c;
  double time{};
  double timeVariance{};
  double timeCorrection{};
  double massOverCharge{};
  ProcessNoise pending;
  double turned{};
  /// Covariance anchor and path since its last transport.
  RzVector anchor;
  double pathSince{};
  /// `Bz` the state moves in from here, and the one at the anchor
  double bz{};
  double anchorBz{};
  /// The radial field where the state stands
  double br{};
  /// Accumulated radial-field correction to the Jacobian q/p column.
  Vector3 brQopPosition{Vector3::Zero()};
  Vector3 brQopDirection{Vector3::Zero()};

  void materialise(const Vector3& normal) {
    pending.apply(c, v.segment<3>(eRzDir0), normal);
  }

  /// Walk on without the covariance
  void travel(double s) {
    pathSince += s;
    const double mOverP = massOverCharge * v[eRzQOverP];
    time += s * std::sqrt(1. + mOverP * mOverP);
  }
  /// Bring the covariance to the state, on a surface with the given normal
  void moveCovariance(const RzHelix& helix, const Vector3& normal) {
    if (pathSince == 0.) {
      return;
    }
    c = helix
            .stepJacobianOnto(anchor, pathSince, v, normal,
                              detail::stepTrig(helix.kappa(anchor) * pathSince),
                              brQopPosition, brQopDirection)
            .transport(c);
    anchor = v;
    anchorBz = bz;
    pathSince = 0.;
    brQopPosition.setZero();
    brQopDirection.setZero();
  }
};

}  // namespace Acts::Experimental::detail::rz

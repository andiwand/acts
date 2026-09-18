// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Convert an RZ free state to bound parameters on its module.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/TrackFinding/Rz/RzLayout.hpp"
#include "Acts/TrackFinding/Rz/RzTransport.hpp"

#include <optional>

namespace Acts::Experimental {

/// Bound parameters and covariance of an RZ state on a surface
struct RzBoundState {
  BoundVector parameters{BoundVector::Zero()};
  BoundMatrix covariance{BoundMatrix::Zero()};
};

/// The Jacobian of the bound parameters with respect to the RZ free state:
/// six rows, the time row zero, seven columns
using RzFreeToBoundMatrix = Eigen::Matrix<double, eBoundSize, eRzSize>;

/// `J C J^T` for a free-to-bound Jacobian and an RZ covariance, with the
/// time variance set to one so that the result stays invertible.
/// @param j the Jacobian
/// @param c the RZ covariance
/// @return the bound covariance
BoundMatrix rzBoundCovariance(const RzFreeToBoundMatrix& j, const RzMatrix& c);

/// Fill the phi, theta and q/p rows of the free-to-bound Jacobian.
/// @param direction the unit direction
/// @param j the Jacobian, position rows left as they are
void rzFillDirectionRows(const Vector3& direction, RzFreeToBoundMatrix& j);

/// Convert a state already on the module plane to bound parameters.
/// @param module the module
/// @param v the RZ state on the module plane
/// @param c covariance before transport, or on the module if transport is null
/// @param transport optional transport from the covariance to the module
/// @return the bound state, or nothing if the polar frame is unsupported
std::optional<RzBoundState> rzBoundOnModule(
    const RzModule& module, const RzVector& v, const RzMatrix& c,
    const RzHelix::StepJacobian* transport = nullptr);

}  // namespace Acts::Experimental

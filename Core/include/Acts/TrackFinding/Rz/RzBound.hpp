// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// An RZ free state expressed as bound parameters on the module it sits on,
/// in closed form from the module's own frame. The surface would give the
/// same answer, through two virtual calls, a general transform inverse and a
/// product on an 8x8 matrix whose time row and column are zero; a track
/// state is written once per measurement, and this is most of what writing
/// it costs.

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

/// The rows of the free-to-bound Jacobian that every surface shares: the
/// direction to azimuth and polar angle, and q/p to itself. The position rows
/// are the caller's.
/// @param direction the unit direction
/// @param j the Jacobian, position rows left as they are
void rzFillDirectionRows(const Vector3& direction, RzFreeToBoundMatrix& j);

/// The state as bound parameters on a module, from the module's frame. A
/// cartesian module measures along its axes from its centre; a polar module
/// measures the radius and azimuth about its surface origin, which the layout
/// has checked to be plain polar. The state is taken to sit on the module
/// plane already.
/// @param module the module
/// @param v the RZ state on the module plane
/// @param c its covariance
/// @return the bound state, or nothing for a polar module the layout could
///         not confirm as plain polar, where the surface has to be asked
std::optional<RzBoundState> rzBoundOnModule(const RzModule& module,
                                            const RzVector& v,
                                            const RzMatrix& c);

}  // namespace Acts::Experimental

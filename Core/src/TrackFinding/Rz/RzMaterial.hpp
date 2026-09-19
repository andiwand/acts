// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFinding/Rz/RzLayout.hpp"
#include "Acts/TrackFinding/Rz/RzTypes.hpp"

namespace Acts::Experimental::detail::rz {

struct State;

/// Process noise accumulated between covariance materialisations.
struct ProcessNoise {
  double varAngle{};
  double varPosition{};
  double covAnglePosition{};
  double varQOverP{};

  void apply(RzMatrix& covariance, const Vector3& direction,
             const Vector3& normal);

  bool empty() const { return varAngle == 0. && varQOverP == 0.; }
  // Signed path: position-direction correlations reverse when walking inward.
  void advance(double s) {
    varPosition += 2. * covAnglePosition * s + varAngle * s * s;
    covAnglePosition += varAngle * s;
  }
};

bool applyMaterial(State& state, const ParticleHypothesis& hypothesis,
                   const MaterialSlab& slab, const Vector3& normal,
                   double direction = 1.);
bool applyMaterial(State& state, const ParticleHypothesis& hypothesis,
                   const RzSurface& surface, std::int32_t band,
                   const Vector3& normal, double direction = 1.);
void regainEnergy(State& state, const ParticleHypothesis& hypothesis,
                  const RzSurface& surface, std::int32_t band,
                  const Vector3& normal);

}  // namespace Acts::Experimental::detail::rz

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/detail/GaussianMixture.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/detail/GeneralMixture.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/detail/Highland.hpp"

#include <array>
#include <cmath>
#include <numbers>
#include <random>
#include <vector>

namespace ActsFatras {

/// Prototype hook: one large scatter of one particle.
struct ScatteringRecord {
  Barcode particleId;
  Acts::Vector3 position;
  /// largest projected angle in units of the Highland theta0
  double ratio = 0;
  double thicknessInX0 = 0;
  double absMomentum = 0;
};

/// Prototype hook: the simulation sets this per event to collect the scatters
/// above `recordThreshold`.
inline thread_local std::vector<ScatteringRecord> *scatteringRecorder = nullptr;
/// Prototype hook: record threshold in units of the Highland theta0.
inline constexpr double recordThreshold = 3.;

/// Simulate (multiple) scattering using a configurable scattering model.
///
/// @tparam scattering_model_t Model implementation to draw a scattering angle.
template <typename scattering_model_t>
struct GenericScattering {
  /// The scattering formula
  scattering_model_t angle;

  /// Simulate scattering and update the particle parameters.
  ///
  /// @param[in]     generator is the random number generator
  /// @param[in]     slab      defines the passed material
  /// @param[in,out] particle  is the particle being updated
  /// @return Empty secondaries containers.
  ///
  /// @tparam generator_t is a RandomNumberEngine
  template <typename generator_t>
  std::array<Particle, 0> operator()(generator_t &generator,
                                     const Acts::MaterialSlab &slab,
                                     Particle &particle) const {
    // the scattered direction can be computed by rotating the initial
    // direction around a vector orthogonal to the initial direction, i.e. the
    // scattering deflector, by the scattering angle. there are an infinite
    // number of vectors orthogonal to the initial direction. the deflector is
    // rotated by some angle relative to some fixpoint.
    //
    // thus two random angles are required: the random deflector orientation
    // angle drawn uniformly from the [-pi,pi) range and the scattering angle
    // drawn from the specific scattering model distribution.

    // draw the random orientation angle
    const double psi = std::uniform_real_distribution<double>(
        -std::numbers::pi, std::numbers::pi)(generator);
    // draw the scattering angle
    const double theta = angle(generator, slab, particle);

    // construct the combined rotation to the scattered direction
    const Acts::RotationMatrix3 rotation(
        // rotation of the scattering deflector axis relative to the reference
        Acts::AngleAxis3(psi, particle.direction()) *
        // rotation by the scattering angle around the deflector axis
        Acts::AngleAxis3(theta,
                         Acts::createCurvilinearUnitU(particle.direction())));
    if (scatteringRecorder != nullptr) {
      const double theta0 = Acts::computeMultipleScatteringTheta0(
          slab, particle.absolutePdg(), particle.mass(), particle.qOverP(),
          particle.absoluteCharge());
      const double ratio =
          std::abs(theta) *
          std::max(std::abs(std::cos(psi)), std::abs(std::sin(psi))) / theta0;
      if (ratio > recordThreshold) {
        scatteringRecorder->push_back(
            {particle.particleId(), particle.position(), ratio,
             slab.thicknessInX0(), particle.absoluteMomentum()});
      }
    }

    particle.setDirection(rotation * particle.direction());

    // scattering is non-destructive and produces no secondaries
    return {};
  }
};

/// Scattering with Gaussian mixture model
using GaussianMixtureScattering = GenericScattering<detail::GaussianMixture>;
/// Scattering with general mixture model
using GeneralMixtureScattering = GenericScattering<detail::GeneralMixture>;
/// Scattering with Highland model
using HighlandScattering = GenericScattering<detail::Highland>;

}  // namespace ActsFatras

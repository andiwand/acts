// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialSlab.hpp"
#include "ActsFatras/EventData/ParticleContainer.hpp"
#include "ActsFatras/EventData/ParticleSimulationQueue.hpp"
#include "ActsFatras/Utilities/detail/FpeSafeGammaDistribution.hpp"

#include <array>
#include <cmath>
#include <numbers>
#include <random>

namespace ActsFatras {

/// Simulate electron energy loss using the Bethe-Heitler description.
///
/// Bethe-Heitler for electron bremsstrahlung description as described here:
/// "A Gaussian-mixture approximation of the Bethe–Heitler model of electron
/// energy loss by bremsstrahlung" R. Frühwirth
struct BetheHeitler {
  /// A scaling factor to
  double scaleFactor = 1.;

  /// Flag for simplified uniform Hertz dipole angle evaluation
  bool uniformHertzDipoleAngle = false;

  /// Simulate the photon emission
  ///
  /// @param [in] particle The unmodified electron
  /// @param [in] gammaE Energy of the photon
  /// @param [in] rndPsi Random number for the azimuthal angle
  /// @param [in] rndTheta1 Random number for the polar angle
  /// @param [in] rndTheta2 Random number for the polar angle
  /// @param [in] rndTheta3 Random number for the polar angle
  /// @param [in] queue The particle simulation queue for creating the photon
  /// @return Particle representing the emitted photon
  MutableParticleProxy bremPhoton(ConstParticleProxy particle, double gammaE,
                                  double rndPsi, double rndTheta1,
                                  double rndTheta2, double rndTheta3,
                                  const ParticleSimulationQueue &queue) const;

  /// Simulate energy loss and update the particle parameters.
  ///
  /// @param[in]     generator is the random number generator
  /// @param[in]     slab      defines the passed material
  /// @param[in,out] particle  is the particle being updated
  /// @return Produced photon.
  ///
  /// @tparam generator_t is a RandomNumberEngine
  template <typename generator_t>
  std::array<MutableParticleProxy, 1> operator()(
      generator_t &generator, const Acts::MaterialSlab &slab,
      MutableParticleProxy particle,
      const ParticleSimulationQueue &queue) const {
    // Take a random gamma-distributed value - depending on t/X0
    detail::FpeSafeGammaDistribution gDist(
        slab.thicknessInX0() / std::numbers::ln2, 1.);

    const double u = gDist(generator);
    const double z = std::exp(-u);
    const double sampledEnergyLoss =
        std::abs(scaleFactor * particle.simulationState().energy() * (z - 1.));

    std::uniform_real_distribution<double> uDist(0., 1.);
    // Build the produced photon
    const double rndPsi = uDist(generator);
    const double rndTheta1 = uDist(generator);
    const double rndTheta2 = uDist(generator);
    const double rndTheta3 = uDist(generator);
    const MutableParticleProxy photon =
        bremPhoton(particle.asConst(), sampledEnergyLoss, rndPsi, rndTheta1,
                   rndTheta2, rndTheta3, queue);

    // Recoil input momentum and apply the energy loss to the electron
    const Acts::Vector3 newMomentum =
        particle.simulationState().direction() *
            particle.simulationState().absoluteMomentum() -
        photon.simulationState().energy() *
            photon.simulationState().direction();
    particle.simulationState().direction() = newMomentum.normalized();
    particle.simulationState().absoluteMomentum() = newMomentum.norm();

    return {photon};
  }
};

}  // namespace ActsFatras

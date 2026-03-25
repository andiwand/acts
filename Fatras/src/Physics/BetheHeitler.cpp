// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Physics/ElectroMagnetic/BetheHeitler.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/ParticleContainer.hpp"
#include "ActsFatras/EventData/ParticleSimulationQueue.hpp"

#include <cmath>
#include <numbers>

namespace ActsFatras {

MutableParticleProxy ActsFatras::BetheHeitler::bremPhoton(
    ConstParticleProxy particle, double gammaE, double rndPsi, double rndTheta1,
    double rndTheta2, double rndTheta3,
    const ParticleSimulationQueue &queue) const {
  // ------------------------------------------------------
  // simple approach
  // (a) simulate theta uniform within the opening angle of the relativistic
  // Hertz dipole
  //      theta_max = 1/gamma
  // (b)Following the Geant4 approximation from L. Urban -> encapsulate that
  // later
  //      the azimutal angle

  const double psi = 2. * std::numbers::pi * rndPsi;

  // the start of the equation
  double theta = 0.;
  if (uniformHertzDipoleAngle) {
    // the simplest simulation
    theta = particle.mass() / particle.simulationState().energy() * rndTheta1;
  } else {
    // ----->
    theta = particle.mass() / particle.simulationState().energy();
    // follow
    constexpr double a = 0.625;  // 5/8
    const double u = -std::log(rndTheta2 * rndTheta3) / a;
    theta *= (rndTheta1 < 0.25) ? u : u / 3.;  // 9./(9.+27) = 0.25
  }

  const Acts::Vector3 particleDirection =
      particle.simulationState().direction();

  // construct the combined rotation to the scattered direction
  const Acts::RotationMatrix3 rotation(
      // rotation of the scattering deflector axis relative to the reference
      Acts::AngleAxis3(psi, particleDirection) *
      // rotation by the scattering angle around the deflector axis
      Acts::AngleAxis3(theta, Acts::createCurvilinearUnitU(particleDirection)));
  const Acts::Vector3 photonDirection = rotation * particleDirection;

  MutableParticleProxy photon = queue.createParticle();
  photon.barcode() = particle.barcode().makeDescendant(0);
  photon.pdg() = Acts::PdgParticle::eGamma;
  photon.generationProcess() = GenerationProcessType::eBremsstrahlung;
  photon.generationState().fourPosition() =
      particle.simulationState().fourPosition();
  photon.generationState().direction() = photonDirection;
  photon.generationState().absoluteMomentum() = gammaE;
  photon.generationState().referenceSurface() =
      particle.simulationState().referenceSurface();
  return photon;
}

}  // namespace ActsFatras

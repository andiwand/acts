// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Material/Interactions.hpp"
#include "ActsFatras/EventData/ParticleContainer.hpp"

#include <numbers>
#include <random>

namespace ActsFatras::detail {

/// Generate scattering angles using a general mixture model.
///
/// Emulates core and tail scattering as described in
///
///     General mixture model Fruehwirth, M. Liendl.
///     Comp. Phys. Comm. 141 (2001) 230-246
///
struct GeneralMixture {
  /// Steering parameter
  bool log_include = true;
  /// Scale the mixture level
  double genMixtureScalor = 1.;

  /// Generate a single 3D scattering angle.
  ///
  /// @param[in]     generator is the random number generator
  /// @param[in]     slab      defines the passed material
  /// @param[in,out] particle  is the particle being scattered
  /// @return a 3d scattering angle
  ///
  /// @tparam generator_t is a RandomNumberEngine
  template <typename generator_t>
  double operator()(generator_t &generator, const Acts::MaterialSlab &slab,
                    MutableParticleProxy particle) const {
    double theta = 0.0;

    if (particle.absolutePdg() != Acts::PdgParticle::eElectron) {
      //----------------------------------------------------------------------------
      // see Mixture models of multiple scattering: computation and simulation.
      // -
      // R.Fruehwirth, M. Liendl. -
      // Computer Physics Communications 141 (2001) 230â246
      //----------------------------------------------------------------------------
      std::array<double, 4> scattering_params{};
      // Decide which mixture is best
      //   beta² = (p/E)² = p²/(p² + m²) = 1/(1 + (m/p)²)
      // 1/beta² = 1 + (m/p)²
      //    beta = 1/sqrt(1 + (m/p)²)
      const double mOverP =
          particle.mass() / particle.simulationState().absoluteMomentum();
      const double beta2Inv = 1 + mOverP * mOverP;
      const double beta = 1 / std::sqrt(beta2Inv);
      const double tInX0 = slab.thicknessInX0();
      const double tob2 = tInX0 * beta2Inv;
      if (tob2 > 0.6 / std::pow(slab.material().Z(), 0.6)) {
        // Gaussian mixture or pure Gaussian
        if (tob2 > 10) {
          scattering_params =
              getGaussian(beta, particle.simulationState().absoluteMomentum(),
                          tInX0, genMixtureScalor);
        } else {
          scattering_params =
              getGaussmix(beta, particle.simulationState().absoluteMomentum(),
                          tInX0, slab.material().Z(), genMixtureScalor);
        }
        // Simulate
        theta = gaussmix(generator, scattering_params);
      } else {
        // Semigaussian mixture - get parameters
        auto scattering_params_sg =
            getSemigauss(beta, particle.simulationState().absoluteMomentum(),
                         tInX0, slab.material().Z(), genMixtureScalor);
        // Simulate
        theta = semigauss(generator, scattering_params_sg);
      }
    } else {
      // for electrons we fall back to the Highland (extension)
      // return projection factor times sigma times gauss random
      const auto theta0 = Acts::computeMultipleScatteringTheta0(
          slab, particle.absolutePdg(), particle.mass(),
          particle.simulationState().qOverP(), particle.absoluteCharge());
      theta = std::normal_distribution<double>(0.0, theta0)(generator);
    }
    // scale from planar to 3d angle
    return std::numbers::sqrt2 * theta;
  }

  // helper methods for getting parameters and simulating

  std::array<double, 4> getGaussian(double beta, double p, double tInX0,
                                    double scale) const {
    std::array<double, 4> scattering_params{};
    // Total standard deviation of mixture
    scattering_params[0] = 15. / beta / p * std::sqrt(tInX0) * scale;
    scattering_params[1] = 1.0;  // Variance of core
    scattering_params[2] = 1.0;  // Variance of tails
    scattering_params[3] = 0.5;  // Mixture weight of tail component
    return scattering_params;
  }

  std::array<double, 4> getGaussmix(double beta, double p, double tInX0,
                                    double Z, double scale) const {
    std::array<double, 4> scattering_params{};
    scattering_params[0] = 15. / beta / p * std::sqrt(tInX0) *
                           scale;  // Total standard deviation of mixture
    const double d1 = std::log(tInX0 / (beta * beta));
    const double d2 = std::log(std::pow(Z, 2.0 / 3.0) * tInX0 / (beta * beta));
    double epsi = 0;
    const double var1 =
        (-1.843e-3 * d1 + 3.347e-2) * d1 + 8.471e-1;  // Variance
                                                      // of core
    if (d2 < 0.5) {
      epsi = (6.096e-4 * d2 + 6.348e-3) * d2 + 4.841e-2;
    } else {
      epsi = (-5.729e-3 * d2 + 1.106e-1) * d2 - 1.908e-2;
    }
    scattering_params[1] = var1;                            // Variance of core
    scattering_params[2] = (1 - (1 - epsi) * var1) / epsi;  // Variance of tails
    scattering_params[3] = epsi;  // Mixture weight of tail component
    return scattering_params;
  }

  std::array<double, 6> getSemigauss(double beta, double p, double tInX0,
                                     double Z, double scale) const {
    std::array<double, 6> scattering_params{};
    const double N = tInX0 * 1.587E7 * std::pow(Z, 1.0 / 3.0) / (beta * beta) /
                     (Z + 1) / std::log(287 / std::sqrt(Z));
    scattering_params[4] = 15. / beta / p * std::sqrt(tInX0) *
                           scale;  // Total standard deviation of mixture
    const double rho = 41000 / std::pow(Z, 2.0 / 3.0);
    const double b = rho / std::sqrt(N * (std::log(rho) - 0.5));
    const double n = std::pow(Z, 0.1) * std::log(N);
    const double var1 = (5.783E-4 * n + 3.803E-2) * n + 1.827E-1;
    const double a =
        (((-4.590E-5 * n + 1.330E-3) * n - 1.355E-2) * n + 9.828E-2) * n +
        2.822E-1;
    const double epsi = (1 - var1) / (a * a * (std::log(b / a) - 0.5) - var1);
    scattering_params[3] =
        (epsi > 0) ? epsi : 0.0;  // Mixture weight of tail component
    scattering_params[0] = a;     // Parameter 1 of tails
    scattering_params[1] = b;     // Parameter 2 of tails
    scattering_params[2] = var1;  // Variance of core
    scattering_params[5] = N;     // Average number of scattering processes
    return scattering_params;
  }

  /// @brief Retrieve the gaussian mixture
  ///
  /// @tparam generator_t Type of the generator
  ///
  /// @param udist The uniform distribution handed over by the call operator
  /// @param scattering_params the tuned parameters for the generation
  ///
  /// @return a double value that represents the gaussian mixture
  template <typename generator_t>
  double gaussmix(generator_t &generator,
                  const std::array<double, 4> &scattering_params) const {
    std::uniform_real_distribution<double> udist(0.0, 1.0);
    const double sigma_tot = scattering_params[0];
    const double var1 = scattering_params[1];
    const double var2 = scattering_params[2];
    const double epsi = scattering_params[3];
    const bool ind = udist(generator) > epsi;
    const double u = udist(generator);
    if (ind) {
      return std::sqrt(var1) * std::sqrt(-2 * std::log(u)) * sigma_tot;
    } else {
      return std::sqrt(var2) * std::sqrt(-2 * std::log(u)) * sigma_tot;
    }
  }

  /// @brief Retrieve the semi-gaussian mixture
  ///
  /// @tparam generator_t Type of the generator
  ///
  /// @param udist The uniform distribution handed over by the call operator
  /// @param scattering_params the tuned parameters for the generation
  ///
  /// @return a double value that represents the gaussian mixture
  template <typename generator_t>
  double semigauss(generator_t &generator,
                   const std::array<double, 6> &scattering_params) const {
    std::uniform_real_distribution<double> udist(0.0, 1.0);
    const double a = scattering_params[0];
    const double b = scattering_params[1];
    const double var1 = scattering_params[2];
    const double epsi = scattering_params[3];
    const double sigma_tot = scattering_params[4];
    const bool ind = udist(generator) > epsi;
    const double u = udist(generator);
    if (ind) {
      return std::sqrt(var1) * std::sqrt(-2 * std::log(u)) * sigma_tot;
    } else {
      return a * b * std::sqrt((1 - u) / (u * b * b + a * a)) * sigma_tot;
    }
  }
};

}  // namespace ActsFatras::detail

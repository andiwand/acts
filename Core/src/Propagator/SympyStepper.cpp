// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/SympyStepper.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/detail/SympyCovarianceEngine.hpp"
#include "Acts/Propagator/detail/SympyJacobianEngine.hpp"

#include <cmath>

#include "codegen/sympy_stepper_math.hpp"

namespace Acts {

namespace {

std::optional<FreeMatrix> computeAdditionalFreeCovariance(
    const SympyStepper& stepper, SympyStepper::State& state) {
  if (!state.denseCache.accumulatedMaterial.isValid()) {
    return std::nullopt;
  }

  FreeMatrix additionalFreeCovariance = FreeMatrix::Zero();

  auto particleHypothesis = stepper.particleHypothesis(state);
  float absQ = particleHypothesis.absoluteCharge();
  float mass = particleHypothesis.mass();

  // handle multiple scattering
  {
    // for derivation see
    // https://github.com/andiwand/cern-scripts/blob/5f0ebf1bef35db65322f28c2e840c1db1aaaf9a7/notebooks/2023-12-07_qp-dense-nav.ipynb
    //
    Vector3 direction = stepper.direction(state);
    SquareMatrix3 directionProjection =
        (ActsSquareMatrix<3>::Identity() - direction * direction.transpose());

    additionalFreeCovariance.block<3, 3>(eFreeDir0, eFreeDir0) =
        state.denseCache.varAngle * directionProjection;
    additionalFreeCovariance.block<3, 3>(eFreePos0, eFreePos0) =
        state.denseCache.varPosition * directionProjection;
  }

  // handle energy loss covariance
  {
    double qOverP = absQ / state.denseCache.initialMomentum;

    float qOverPSigma = computeEnergyLossLandauSigmaQOverP(
        state.denseCache.accumulatedMaterial, mass, qOverP, absQ);

    additionalFreeCovariance(eFreeQOverP, eFreeQOverP) =
        qOverPSigma * qOverPSigma;
  }

  state.denseCache.initialMomentum = 0;
  state.denseCache.accumulatedMaterial = MaterialSlab();
  state.denseCache.varPosition = 0;
  state.denseCache.varAngle = 0;
  state.denseCache.covAnglePosition = 0;

  return additionalFreeCovariance;
}

}  // namespace

SympyStepper::SympyStepper(std::shared_ptr<const MagneticFieldProvider> bField)
    : m_bField(std::move(bField)) {}

SympyStepper::SympyStepper(const Config& config) : m_bField(config.bField) {}

SympyStepper::State SympyStepper::makeState(
    const Options& options, const BoundTrackParameters& par) const {
  State state{options, m_bField->makeCache(options.magFieldContext)};

  state.particleHypothesis = par.particleHypothesis();

  Vector3 position = par.position(options.geoContext);
  Vector3 direction = par.direction();
  state.pars.template segment<3>(eFreePos0) = position;
  state.pars.template segment<3>(eFreeDir0) = direction;
  state.pars[eFreeTime] = par.time();
  state.pars[eFreeQOverP] = par.parameters()[eBoundQOverP];

  // Init the jacobian matrix if needed
  if (par.covariance()) {
    // Get the reference surface for navigation
    const auto& surface = par.referenceSurface();
    // set the covariance transport flag to true and copy
    state.covTransport = true;
    state.cov = BoundSquareMatrix(*par.covariance());
    state.jacToGlobal =
        surface.boundToFreeJacobian(options.geoContext, position, direction);
  }

  state.stepSize = ConstrainedStep(options.maxStepSize);

  return state;
}

void SympyStepper::resetState(State& state, const BoundVector& boundParams,
                              const BoundSquareMatrix& cov,
                              const Surface& surface,
                              const double stepSize) const {
  FreeVector freeParams = transformBoundToFreeParameters(
      surface, state.options.geoContext, boundParams);

  // Update the stepping state
  state.pars = freeParams;
  state.cov = cov;
  state.stepSize = ConstrainedStep(stepSize);
  state.pathAccumulated = 0.;

  // Reinitialize the stepping jacobian
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.options.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
  state.jacobian = BoundMatrix::Identity();
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

Result<std::tuple<BoundTrackParameters, BoundMatrix, double>>
SympyStepper::boundState(
    State& state, const Surface& surface, bool transportCov,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  std::optional<FreeMatrix> additionalFreeCovariance =
      computeAdditionalFreeCovariance(*this, state);
  return detail::sympy::boundState(
      state.options.geoContext, surface, state.cov, state.jacobian,
      state.jacTransport, state.derivative, state.jacToGlobal,
      additionalFreeCovariance, state.pars, state.particleHypothesis,
      state.covTransport && transportCov, state.pathAccumulated,
      freeToBoundCorrection);
}

std::tuple<CurvilinearTrackParameters, BoundMatrix, double>
SympyStepper::curvilinearState(State& state, bool transportCov) const {
  std::optional<FreeMatrix> additionalFreeCovariance =
      computeAdditionalFreeCovariance(*this, state);
  return detail::sympy::curvilinearState(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, additionalFreeCovariance, state.pars,
      state.particleHypothesis, state.covTransport && transportCov,
      state.pathAccumulated);
}

void SympyStepper::update(State& state, const FreeVector& freeParams,
                          const BoundVector& /*boundParams*/,
                          const Covariance& covariance,
                          const Surface& surface) const {
  state.pars = freeParams;
  state.cov = covariance;
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.options.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
}

void SympyStepper::update(State& state, const Vector3& uposition,
                          const Vector3& udirection, double qOverP,
                          double time) const {
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = qOverP;
}

void SympyStepper::transportCovarianceToCurvilinear(State& state) const {
  detail::sympy::transportCovarianceToCurvilinear(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, std::nullopt,
      state.pars.template segment<3>(eFreeDir0));
}

void SympyStepper::transportCovarianceToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  detail::sympy::transportCovarianceToBound(
      state.options.geoContext, surface, state.cov, state.jacobian,
      state.jacTransport, state.derivative, state.jacToGlobal, std::nullopt,
      state.pars, freeToBoundCorrection);
}

Result<double> SympyStepper::stepImpl(State& state, Direction stepDirection,
                                      double stepTolerance,
                                      double stepSizeCutOff,
                                      std::size_t maxRungeKuttaStepTrials,
                                      const IVolumeMaterial* material) const {
  double m = particleHypothesis(state).mass();
  PdgParticle absPdg = particleHypothesis(state).absolutePdg();
  double q = charge(state);
  Vector3 pos = position(state);
  Vector3 dir = direction(state);
  double t = time(state);
  double qop = qOverP(state);
  double pabs = absoluteMomentum(state);

  const auto getB = [&](const double* p) -> Result<Vector3> {
    return getField(state, {p[0], p[1], p[2]});
  };

  const auto getG = [&](const double* p, double l) -> double {
    return computeEnergyLossMean(
        MaterialSlab(material->material({p[0], p[1], p[2]}),
                     1.0f * UnitConstants::mm),
        absPdg, m, l, q);
  };

  const auto calcStepSizeScaling = [&](const double errorEstimate_) -> double {
    // For details about these values see ATL-SOFT-PUB-2009-001
    constexpr double lower = 0.25;
    constexpr double upper = 4.0;
    // This is given by the order of the Runge-Kutta method
    constexpr double exponent = 0.25;

    double x = stepTolerance / errorEstimate_;

    if constexpr (exponent == 0.25) {
      // This is 3x faster than std::pow
      x = std::sqrt(std::sqrt(x));
    } else {
      x = std::pow(x, exponent);
    }

    return std::clamp(x, lower, upper);
  };

  double h = state.stepSize.value() * stepDirection;
  double initialH = h;
  std::size_t nStepTrials = 0;
  double errorEstimate = 0.;

  while (true) {
    ++nStepTrials;
    ++state.statistics.nAttemptedSteps;

    // For details about the factor 4 see ATL-SOFT-PUB-2009-001
    Result<bool> res = Result<bool>::success(false);
    if (material == nullptr) {
      res = rk4_vacuum(
          pos.data(), dir.data(), t, h, qop, m, pabs, getB, &errorEstimate,
          4 * stepTolerance, state.pars.template segment<3>(eFreePos0).data(),
          state.pars.template segment<1>(eFreeTime).data(),
          state.pars.template segment<3>(eFreeDir0).data(),
          state.derivative.data(),
          state.covTransport ? state.jacTransport.data() : nullptr);
    } else {
      res = rk4_dense(pos.data(), dir.data(), t, h, qop, m, q, pabs, getB, getG,
                      &errorEstimate, 4 * stepTolerance,
                      state.pars.template segment<3>(eFreePos0).data(),
                      state.pars.template segment<1>(eFreeTime).data(),
                      state.pars.template segment<3>(eFreeDir0).data(),
                      state.pars.template segment<1>(eFreeQOverP).data(),
                      state.derivative.data(),
                      state.covTransport ? state.jacTransport.data() : nullptr);
    }
    if (!res.ok()) {
      return res.error();
    }
    // Protect against division by zero
    errorEstimate = std::max(1e-20, errorEstimate);

    if (*res) {
      break;
    }

    ++state.statistics.nRejectedSteps;

    const double stepSizeScaling = calcStepSizeScaling(errorEstimate);
    h *= stepSizeScaling;

    // If step size becomes too small the particle remains at the initial
    // place
    if (std::abs(h) < std::abs(stepSizeCutOff)) {
      // Not moving due to too low momentum needs an aborter
      return EigenStepperError::StepSizeStalled;
    }

    // If the parameter is off track too much or given stepSize is not
    // appropriate
    if (nStepTrials > maxRungeKuttaStepTrials) {
      // Too many trials, have to abort
      return EigenStepperError::StepSizeAdjustmentFailed;
    }
  }

  state.pathAccumulated += h;
  ++state.nSteps;
  state.nStepTrials += nStepTrials;

  ++state.statistics.nSuccessfulSteps;
  if (stepDirection != Direction::fromScalarZeroAsPositive(initialH)) {
    ++state.statistics.nReverseSteps;
  }
  state.statistics.pathLength += h;
  state.statistics.absolutePathLength += std::abs(h);

  const double stepSizeScaling = calcStepSizeScaling(errorEstimate);
  const double nextAccuracy = std::abs(h * stepSizeScaling);
  const double previousAccuracy = std::abs(state.stepSize.accuracy());
  const double initialStepLength = std::abs(initialH);
  if (nextAccuracy < initialStepLength || nextAccuracy > previousAccuracy) {
    state.stepSize.setAccuracy(nextAccuracy);
  }

  if (material != nullptr || state.denseCache.accumulatedMaterial.isValid()) {
    if (!state.denseCache.accumulatedMaterial.isValid()) {
      state.denseCache.initialMomentum = pabs;
    }

    Material mat = material != nullptr ? material->material(pos) : Material();
    MaterialSlab slab(mat, h);

    std::size_t substepCount =
        mat.isValid()
            ? static_cast<std::size_t>(std::ceil(slab.thicknessInX0() /
                                                 state.options.maxXOverX0Step))
            : 1;
    double substep = h * 1.0 / substepCount;
    MaterialSlab subslab(mat, substep);

    MaterialSlab accumulatedMaterial = state.denseCache.accumulatedMaterial;
    for (std::size_t i = 0; i < substepCount; ++i) {
      float theta0in = computeMultipleScatteringTheta0(accumulatedMaterial,
                                                       absPdg, m, qop, q);

      accumulatedMaterial =
          MaterialSlab::combineLayers(accumulatedMaterial, subslab);

      float theta0out = computeMultipleScatteringTheta0(accumulatedMaterial,
                                                        absPdg, m, qop, q);

      double deltaVarTheta = square(theta0out) - square(theta0in);
      double deltaVarPos = state.denseCache.varAngle * square(substep) +
                           2 * state.denseCache.covAnglePosition * substep +
                           deltaVarTheta * (square(substep) / 3);
      state.denseCache.varPosition += deltaVarPos;
      state.denseCache.covAnglePosition += state.denseCache.varAngle * substep;
      state.denseCache.varAngle += deltaVarTheta;
    }
    state.denseCache.accumulatedMaterial = accumulatedMaterial;
  }

  return h;
}

void SympyStepper::setIdentityJacobian(State& state) const {
  state.jacobian = BoundMatrix::Identity();
}

}  // namespace Acts

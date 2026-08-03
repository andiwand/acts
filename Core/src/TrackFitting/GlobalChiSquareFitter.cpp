// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GlobalChiSquareFitter.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <cmath>

void Acts::Experimental::updateGx2fParams(
    BoundTrackParameters& params, const Eigen::VectorXd& deltaParamsExtended,
    const Gx2fParameterLayout& layout,
    std::unordered_map<GeometryIdentifier, Gx2fMaterialProperties>& materialMap,
    const std::vector<GeometryIdentifier>& geoIdVector) {
  assert(geoIdVector.size() == layout.nMaterialSurfaces() &&
         "Number of visited material surfaces does not match the layout.");

  // update params
  params.parameters() +=
      deltaParamsExtended.topLeftCorner<eBoundSize, 1>().eval();

  // update the material parameters
  for (std::size_t matSurface = 0; matSurface < layout.nMaterialSurfaces();
       matSurface++) {
    const GeometryIdentifier geoId = geoIdVector[matSurface];
    const auto materialMapId = materialMap.find(geoId);
    assert(materialMapId != materialMap.end() &&
           "No material properties found for material surface.");

    if (layout.fitScattering()) {
      materialMapId->second.scatteringAngles().segment<2>(eBoundPhi) +=
          deltaParamsExtended.segment<2>(layout.scatteringOffset(matSurface))
              .eval();
    }

    if (layout.fitEnergyLoss()) {
      materialMapId->second.qOverPOffset() +=
          deltaParamsExtended(layout.energyLossOffset(matSurface));
    }
  }

  return;
}

void Acts::Experimental::updateGx2fCovarianceParams(
    BoundMatrix& fullCovariancePredicted, Gx2fSystem& extendedSystem) {
  // Regularize unmeasured parameters without changing the fitted system.
  Eigen::MatrixXd aMatrix = extendedSystem.aMatrix();
  for (Eigen::Index i = 0; i < aMatrix.rows(); ++i) {
    if (aMatrix(i, i) == 0.) {
      aMatrix(i, i) = 1.;
    }
  }

  const Eigen::VectorXd scale = aMatrix.diagonal().cwiseSqrt().cwiseInverse();
  const Eigen::MatrixXd aScaled =
      scale.asDiagonal() * aMatrix * scale.asDiagonal();
  const Eigen::MatrixXd covariance =
      scale.asDiagonal() * aScaled.inverse() * scale.asDiagonal();

  for (Eigen::Index i = 0; i < eBoundSize; ++i) {
    for (Eigen::Index j = 0; j < eBoundSize; ++j) {
      if (extendedSystem.aMatrix()(i, i) != 0. &&
          extendedSystem.aMatrix()(j, j) != 0.) {
        fullCovariancePredicted(i, j) = covariance(i, j);
      }
    }
  }
}

void Acts::Experimental::addMeasurementToGx2fSumsBackend(
    Gx2fSystem& extendedSystem,
    const std::vector<BoundMatrix>& jacobianFromStart,
    const Eigen::MatrixXd& covarianceMeasurement, const BoundVector& predicted,
    const Eigen::VectorXd& measurement, const Eigen::MatrixXd& projector,
    const Logger& logger) {
  // Skip measurements with a non-invertible covariance.
  const auto safeInvCovMeasurement = safeInverse(covarianceMeasurement);
  if (!safeInvCovMeasurement) {
    ACTS_WARNING("addMeasurementToGx2fSums: safeInvCovMeasurement failed.");
    ACTS_VERBOSE("    covarianceMeasurement:\n" << covarianceMeasurement);
    return;
  }

  // Map all fitted parameters to this state.
  Eigen::MatrixXd extendedJacobian =
      Eigen::MatrixXd::Zero(eBoundSize, extendedSystem.nDims());

  // Transport from the reference surface.
  extendedJacobian.topLeftCorner<eBoundSize, eBoundSize>() =
      jacobianFromStart[0];

  // Include only material surfaces reached before or at this measurement.
  const Gx2fParameterLayout& layout = extendedSystem.layout();
  assert(jacobianFromStart.size() <= layout.nMaterialSurfaces() + 1 &&
         "More Jacobians than fitted material surfaces.");

  for (std::size_t matSurface = 1; matSurface < jacobianFromStart.size();
       matSurface++) {
    const BoundMatrix& jac = jacobianFromStart[matSurface];

    // The index of the material surface this Jacobian starts from
    const std::size_t k = matSurface - 1;

    if (layout.fitScattering()) {
      extendedJacobian.template block<eBoundSize, 2>(
          0, layout.scatteringOffset(k)) =
          jac * Gx2fConstants::phiThetaProjector;
    }

    if (layout.fitEnergyLoss()) {
      // Projecting onto q/p is just picking that column of the Jacobian
      extendedJacobian.template block<eBoundSize, 1>(
          0, layout.energyLossOffset(k)) = jac.col(eBoundQOverP);
    }
  }

  const Eigen::MatrixXd projJacobian = projector * extendedJacobian;

  const Eigen::VectorXd projPredicted = projector * predicted;

  const Eigen::VectorXd residual = measurement - projPredicted;

  // Finally contribute to chi2sum, aMatrix, and bVector
  extendedSystem.chi2() +=
      (residual.transpose() * (*safeInvCovMeasurement) * residual)(0, 0);

  extendedSystem.aMatrix() +=
      (projJacobian.transpose() * (*safeInvCovMeasurement) * projJacobian)
          .eval();

  extendedSystem.bVector() +=
      (residual.transpose() * (*safeInvCovMeasurement) * projJacobian)
          .eval()
          .transpose();

  ACTS_VERBOSE(
      "Contributions in addMeasurementToGx2fSums:\n"
      << "    predicted:   " << predicted.transpose() << "\n"
      << "    measurement: " << measurement.transpose() << "\n"
      << "    covarianceMeasurement:\n"
      << covarianceMeasurement << "\n"
      << "    projector:\n"
      << projector.eval() << "\n"
      << "    projJacobian:\n"
      << projJacobian.eval() << "\n"
      << "    projPredicted: " << (projPredicted.transpose()).eval() << "\n"
      << "    residual: " << (residual.transpose()).eval() << "\n"
      << "    extendedJacobian:\n"
      << extendedJacobian << "\n"
      << "    aMatrix contribution:\n"
      << (projJacobian.transpose() * (*safeInvCovMeasurement) * projJacobian)
             .eval()
      << "\n"
      << "    bVector contribution: "
      << (residual.transpose() * (*safeInvCovMeasurement) * projJacobian).eval()
      << "\n"
      << "    chi2sum contribution: "
      << (residual.transpose() * (*safeInvCovMeasurement) * residual)(0, 0)
      << "\n"
      << "    safeInvCovMeasurement:\n"
      << (*safeInvCovMeasurement));
}

double Acts::Experimental::computeGx2fQOverPOffset(
    const MaterialSlab& slab, const ParticleHypothesis& particleHypothesis,
    const double qOverP, const Direction direction,
    const Gx2fEnergyLossMode mode) {
  const PdgParticle absPdg = particleHypothesis.absolutePdg();
  const double mass = particleHypothesis.mass();
  const double absQ = particleHypothesis.absoluteCharge();

  if (slab.isVacuum() || particleHypothesis.hasMomentumHypothesis() ||
      absQ <= 0. || mass <= 0. || qOverP == 0.) {
    return 0.;
  }

  // Keep the momentum update in double precision.
  const double eLoss =
      (mode == Gx2fEnergyLossMode::Mean)
          ? static_cast<double>(computeEnergyLossMean(
                slab, absPdg, static_cast<float>(mass),
                static_cast<float>(qOverP), static_cast<float>(absQ)))
          : static_cast<double>(computeEnergyLossMode(
                slab, absPdg, static_cast<float>(mass),
                static_cast<float>(qOverP), static_cast<float>(absQ)));

  const double momentum = particleHypothesis.extractMomentum(qOverP);

  // in forward(backward) propagation, energy decreases(increases)
  const double nextE = fastHypot(mass, momentum) - eLoss * direction;
  // put the particle at rest if the energy loss is too large
  double nextP = (mass < nextE) ? fastCathetus(nextE, mass) : 0.;

  // minimum momentum below which we will not push particles via material update
  static constexpr double minP = 10 * UnitConstants::MeV;
  nextP = std::max(minP, nextP);

  const double nextQOverP =
      particleHypothesis.qOverP(nextP, std::copysign(absQ, qOverP));

  return nextQOverP - qOverP;
}

Eigen::VectorXd Acts::Experimental::computeGx2fDeltaParams(
    const Acts::Experimental::Gx2fSystem& extendedSystem) {
  // Equilibrate the mixed parameter units, leaving zero diagonals unscaled.
  Eigen::VectorXd scale = extendedSystem.aMatrix().diagonal().cwiseAbs();
  scale = (scale.array() > 0.).select(scale, 1.);
  scale = scale.cwiseSqrt().cwiseInverse();

  const Eigen::MatrixXd aScaled =
      scale.asDiagonal() * extendedSystem.aMatrix() * scale.asDiagonal();
  const Eigen::VectorXd bScaled = scale.asDiagonal() * extendedSystem.bVector();

  return scale.asDiagonal() * aScaled.colPivHouseholderQr().solve(bScaled);
}

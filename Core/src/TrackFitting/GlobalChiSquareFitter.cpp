// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GlobalChiSquareFitter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"

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
  }

  return;
}

void Acts::Experimental::updateGx2fCovarianceParams(
    BoundMatrix& fullCovariancePredicted, Gx2fSystem& extendedSystem) {
  // make invertible
  for (std::size_t i = 0; i < extendedSystem.nDims(); ++i) {
    if (extendedSystem.aMatrix()(i, i) == 0.) {
      extendedSystem.aMatrix()(i, i) = 1.;
    }
  }

  visit_measurement(extendedSystem.findRequiredNdf(), [&](auto N) {
    fullCovariancePredicted.topLeftCorner<N, N>() =
        extendedSystem.aMatrix().inverse().topLeftCorner<N, N>();
  });

  return;
}

void Acts::Experimental::addMeasurementToGx2fSumsBackend(
    Gx2fSystem& extendedSystem,
    const std::vector<BoundMatrix>& jacobianFromStart,
    const Eigen::MatrixXd& covarianceMeasurement, const BoundVector& predicted,
    const Eigen::VectorXd& measurement, const Eigen::MatrixXd& projector,
    const Logger& logger) {
  // First, w try to invert the covariance matrix. If the inversion fails, we
  // can already abort.
  const auto safeInvCovMeasurement = safeInverse(covarianceMeasurement);
  if (!safeInvCovMeasurement) {
    ACTS_WARNING("addMeasurementToGx2fSums: safeInvCovMeasurement failed.");
    ACTS_VERBOSE("    covarianceMeasurement:\n" << covarianceMeasurement);
    return;
  }

  // Create an extended Jacobian. This one contains only eBoundSize rows,
  // because the rest is irrelevant. We fill it in the next steps.
  // TODO make dimsExtendedParams template with unrolling
  Eigen::MatrixXd extendedJacobian =
      Eigen::MatrixXd::Zero(eBoundSize, extendedSystem.nDims());
  assert(extendedJacobian.cols() ==
             static_cast<Eigen::Index>(extendedSystem.nDims()) &&
         "Extended Jacobian does not match the system dimensions.");

  // This part of the Jacobian comes from the material-less propagation
  extendedJacobian.topLeftCorner<eBoundSize, eBoundSize>() =
      jacobianFromStart[0];

  // If we have material, loop here over all Jacobians. We add extra columns for
  // the parameters attached to each material surface. These parts account for
  // the propagation of the scattering angles.
  // We hold one Jacobian per material surface passed so far, plus the one from
  // the start of the track. Material surfaces downstream of this measurement
  // have not been reached yet, hence the inequality.
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

Eigen::VectorXd Acts::Experimental::computeGx2fDeltaParams(
    const Acts::Experimental::Gx2fSystem& extendedSystem) {
  return extendedSystem.aMatrix().colPivHouseholderQr().solve(
      extendedSystem.bVector());
}

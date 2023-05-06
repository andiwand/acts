// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PerigeeSurface.hpp"

#include <iostream>

template <typename propagator_t, typename propagator_options_t>
Acts::Result<Acts::LinearizedTrack>
Acts::HelicalTrackLinearizerWithTime<propagator_t, propagator_options_t>::
    linearizeTrack(const BoundTrackParameters& params, const Vector4& linPoint,
                   const Acts::GeometryContext& gctx,
                   const Acts::MagneticFieldContext& mctx, State& state) const {
  // estimated vertex position
  Vector3 linPointPos = VectorHelpers::position(linPoint);
  // coordinate system of perigee has axes parallel to gloab axes
  // origin at (linPointPos, linPoint.time())
  const std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(linPointPos);

  auto intersection = perigeeSurface->intersect(gctx, params.position(gctx),
                                                params.unitDirection(), false);

  // Create propagator options
  propagator_options_t pOptions(gctx, mctx);
  pOptions.direction = intersection.intersection.pathLength >= 0
                           ? Direction::Forward
                           : Direction::Backward;

  // Do the propagation to the PCA of linPointPos
  auto result = m_cfg.propagator->propagate(params, *perigeeSurface, pOptions);
  if (not result.ok()) {
    return result.error();
  }

  const auto& endParams = *result->endParameters;

  // parameters of the original track at the linearization point (i.e., the PCA
  // to the vertex, i.e., the point V)
  BoundVector paramsAtPCA = endParams.parameters();
  // this needs to be included in the linearization point in the future
  double ref_time_point =
      0.;  // I don't think we need this here since R = linPointPos

  // time needs to be relative to the new perigee (just like its spatial
  // parameters)
  paramsAtPCA[eBoundTime] = endParams.time() + ref_time_point - linPoint[eTime];
  /*
  if (m_cfg.verbose) {
    std::cout<<"PARAMS AT PCA (LIN POINT)"<<std::endl;
    std::cout<<"\n"<<paramsAtPCA<<std::endl;
    std::cout<<"ADDING TIME"<<std::endl;
    std::cout<<"\n"<<paramsAtPCA<<std::endl;
  }
  */
  // global coordinates of the PCA (i.e., the point V)

  Vector4 positionAtPCA = Vector4::Zero();
  {
    auto pos = endParams.position(gctx);
    positionAtPCA[ePos0] = pos[ePos0];
    positionAtPCA[ePos1] = pos[ePos1];
    positionAtPCA[ePos2] = pos[ePos2];

    // The time in global frame
    positionAtPCA[eTime] = endParams.time() + ref_time_point;
    // positionAtPCA[eTime] = endParams.time() - linPoint[eTime];
  }
  BoundSymMatrix parCovarianceAtPCA = endParams.covariance().value();

  if (parCovarianceAtPCA.determinant() <= 0) {
    // Use the original parameters
    paramsAtPCA = params.parameters();
    auto pos = endParams.position(gctx);
    // the following 3 lines are the same in the if confition above - shouldn't
    // they be outside the if then?
    positionAtPCA[ePos0] = pos[ePos0];
    positionAtPCA[ePos1] = pos[ePos1];
    positionAtPCA[ePos2] = pos[ePos2];
    // Add timing
    // positionAtPCA[eTime] = params.time() - linPoint[eTime];

    // The time in global frame
    positionAtPCA[eTime] = params.time() + ref_time_point;
    parCovarianceAtPCA = params.covariance().value();
  }
  /*
    if (m_cfg.verbose) {
      std::cout<<"PF::PositionAtPCA - Linearizer"<<std::endl;
      std::cout<<positionAtPCA<<std::endl;
      std::cout<<"Parameter covariance at PCA"<<std::endl;
      std::cout<<parCovarianceAtPCA<<std::endl;
    }
    */

  // phiV and functions
  double phiV = paramsAtPCA(BoundIndices::eBoundPhi);
  double sinPhiV = std::sin(phiV);
  double cosPhiV = std::cos(phiV);

  // theta and functions
  double th = paramsAtPCA(BoundIndices::eBoundTheta);
  const double sinTh = std::sin(th);
  const double tanTh = std::tan(th);

  // q over p
  double qOvP = paramsAtPCA(BoundIndices::eBoundQOverP);
  // if the charge is negative the particle moves counterclockwise (when looking
  // in negative z-direction) therefore the sign of the curvature is negative
  double sgnH = (qOvP < 0.) ? -1 : 1;

  Vector3 momentumAtPCA(phiV, th, qOvP);

  // get B-field z-component at current position
  auto field = m_cfg.bField->getField(VectorHelpers::position(positionAtPCA),
                                      state.fieldCache);
  if (!field.ok()) {
    return field.error();
  }
  double Bz = (*field)[eZ];
  double rho = 0;
  // Curvature is infinite w/o b field
  if (Bz == 0. || std::abs(qOvP) < m_cfg.minQoP) {
    rho = sgnH * m_cfg.maxRho;
  } else {
    rho = sinTh * (1. / qOvP) / Bz;
  }

  // Eq. 5.34 in Ref(1) (see .hpp)
  // the reference point R corresponds to linearization point
  double X = positionAtPCA(0) - linPointPos.x() + rho * sinPhiV;
  double Y = positionAtPCA(1) - linPointPos.y() - rho * cosPhiV;
  const double S2 = (X * X + Y * Y);
  const double S = std::sqrt(S2);

  /// F(V, p_i) at PCA in Billoir paper
  /// (see FullBilloirVertexFitter.hpp for paper reference,
  /// Page 140, Eq. (2) )
  BoundVector predParamsAtPCA;

  int sgnX = (X < 0.) ? -1 : 1;
  int sgnY = (Y < 0.) ? -1 : 1;

  double phiAtPCA = 0;  // I would change the variable name to predPhiAtPCA
  // why do we evaluate the arctan here so strangely?
  if (std::abs(X) > std::abs(Y)) {
    phiAtPCA = sgnH * sgnX * std::acos(-sgnH * Y / S);
  } else {
    phiAtPCA = std::asin(sgnH * X / S);
    if ((sgnH * sgnY) > 0) {
      phiAtPCA = sgnH * sgnX * M_PI - phiAtPCA;
    }
  }

  // Eq. 5.33 in Ref(1) (see .hpp)
  // sgnH corresponds to the sign of rho in he reference above
  // Do we even need the equations below?
  predParamsAtPCA[0] = rho - sgnH * S;
  predParamsAtPCA[1] =
      positionAtPCA[eZ] - linPointPos.z() +
      rho * (phiV - phiAtPCA) /
          tanTh;  // sign mistake in Ref(1) - here the equation is correct
  predParamsAtPCA[2] = phiAtPCA;
  predParamsAtPCA[3] = th;
  predParamsAtPCA[4] = qOvP;
  // predParamsAtPCA[5] = 0.;

  // Adding timing to the predicted parameters

  // The sign of the DeltaT is positive if the track has to move in the
  // forward direction in order to reach the new point.

  // int sgnt = intersection.intersection.pathLength >= 0 ? 1 : -1;
  ActsScalar pInGeV = std::abs(1.0 / qOvP);
  ActsScalar massInGeV =
      0.1057;  // in GeV - pion, mass hypothesis! TODO: find out why
  ActsScalar beta = pInGeV / std::hypot(pInGeV, massInGeV);

  double extrap_term = 0.;

  // If we neglect correlation between vtx position and vtx time
  // we should get the weighted avg.
  // DeltaPhi = phi_p - phi_V
  // if (m_cfg.TimeAndPosFit)
  extrap_term = -rho * (phiAtPCA - phiV) / (beta * sinTh);
  // I think we should write paramsAtPCA[eBoundTime] instead of
  // positionAtPCA[eTime] - linPoint[eTime]
  predParamsAtPCA[5] = positionAtPCA[eTime] - linPoint[eTime] +
                       extrap_term;  // should be relative to linPoint

  std::cout << "Parameters of V:\n"
            << paramsAtPCA << "\n"
            << "Parameters of P:\n"
            << predParamsAtPCA << "\n";

  // Fill position jacobian (D_k matrix), Eq. 5.36 in Ref(1)
  ActsMatrix<eBoundSize, 4> positionJacobian;
  positionJacobian.setZero();
  // First row
  positionJacobian(0, 0) = -sgnH * X / S;
  positionJacobian(0, 1) = -sgnH * Y / S;

  const double S2tanTh = S2 * tanTh;
  const double S2sinTh = S2 * sinTh;

  // Second row
  positionJacobian(1, 0) = rho * Y / S2tanTh;
  positionJacobian(1, 1) = -rho * X / S2tanTh;
  positionJacobian(1, 2) = 1.;

  // Third row
  positionJacobian(2, 0) = -Y / S2;
  positionJacobian(2, 1) = X / S2;

  // x Last row
  /*
    if (m_cfg.TimeAndPosFit) {

      if (m_cfg.verbose) {
        std::cout<<"dz0/dXV="<<rho * Y / S2tanTh<<std::endl;
        std::cout<<"beta="<<beta<<std::endl;
        std::cout<<"S2tanTh="<<S2tanTh<<std::endl;
        std::cout<<"S2sinTh="<<S2sinTh<<std::endl;
      }

      positionJacobian(5, 0) = (rho / beta) * Y / S2sinTh;
      positionJacobian(5, 1) = (-rho / beta) * X / S2sinTh;
    }
    */
  positionJacobian(5, 0) = (rho / beta) * Y / S2sinTh;
  positionJacobian(5, 1) = (-rho / beta) * X / S2sinTh;
  positionJacobian(5, 3) = 1.;

  // Fill momentum jacobian (E_k matrix), Eq. 5.37 in Ref(1)
  ActsMatrix<eBoundSize, 3> momentumJacobian;
  momentumJacobian.setZero();

  double R = X * cosPhiV + Y * sinPhiV;
  double Q = X * sinPhiV - Y * cosPhiV;
  double dPhi = phiAtPCA - phiV;  // I think we can neglect all dPhi terms

  // First row
  momentumJacobian(0, 0) = -sgnH * rho * R / S;

  double qOvSred = 1 - sgnH * Q / S;

  momentumJacobian(0, 1) = qOvSred * rho / tanTh;
  momentumJacobian(0, 2) = -qOvSred * rho / qOvP;

  const double rhoOverS2 = rho / S2;

  // Second row
  momentumJacobian(1, 0) = (1 - rhoOverS2 * Q) * rho / tanTh;
  momentumJacobian(1, 1) = (dPhi + rho * R / (S2tanTh * tanTh)) * rho;
  momentumJacobian(1, 2) = (dPhi - rhoOverS2 * R) * rho / (qOvP * tanTh);

  // Third row
  momentumJacobian(2, 0) = rhoOverS2 * Q;
  momentumJacobian(2, 1) = -rho * R / S2tanTh;
  momentumJacobian(2, 2) = rhoOverS2 * R / qOvP;

  // Fourth and fifth row:
  momentumJacobian(3, 1) = 1.;
  momentumJacobian(4, 2) = 1.;

  // TODO: fill last row
  double oneOverOmega = rho / (sinTh * beta);
  momentumJacobian(5, 0) = oneOverOmega * (1 - rhoOverS2 * Q);
  momentumJacobian(5, 1) = oneOverOmega * rho * R / S2 * tanTh;
  momentumJacobian(5, 2) =
      oneOverOmega / qOvP * (beta * beta * dPhi - R * rhoOverS2);
  // const term F(V_0, p_0) in Talyor expansion
  // F(V_0,p_0) = F(V,p) - D * dV - E * dP

  BoundVector constTerm =
      predParamsAtPCA - positionJacobian * positionAtPCA -
      momentumJacobian * momentumAtPCA;  // should it not be paramsAtPCA instead
                                         // of predParamsAtPCA?

  // The parameter weight
  // ActsSymMatrix<5> parWeight = (parCovarianceAtPCA.block<5, 5>(0,
  // 0)).inverse(); BoundSymMatrix weightAtPCA{BoundSymMatrix::Identity()};
  // weightAtPCA.block<5, 5>(0, 0) = parWeight;

  BoundSymMatrix weightAtPCA = parCovarianceAtPCA.inverse();

  return LinearizedTrack(paramsAtPCA, parCovarianceAtPCA, weightAtPCA, linPoint,
                         positionJacobian, momentumJacobian, positionAtPCA,
                         momentumAtPCA, constTerm);
}

template <typename propagator_t, typename propagator_options_t>
void Acts::HelicalTrackLinearizerWithTime<propagator_t, propagator_options_t>::
    calculateNumericalJacobians(const BoundTrackParameters& params,
                                const Vector4& linPoint,
                                const Acts::GeometryContext& gctx,
                                const Acts::MagneticFieldContext& mctx,
                                State& state) const {
  // estimated vertex position
  Vector3 linPointPos = VectorHelpers::position(linPoint);
  // coordinate system of perigee has axes parallel to gloab axes
  // origin at (linPointPos, linPoint.time())
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(linPointPos);

  auto intersection = perigeeSurface->intersect(gctx, params.position(gctx),
                                                params.unitDirection(), false);

  // Create propagator options
  propagator_options_t pOptions(gctx, mctx);
  pOptions.targetTolerance = 1e-10;
  pOptions.direction = intersection.intersection.pathLength >= 0
                           ? Direction::Forward
                           : Direction::Backward;

  // Do the propagation to the PCA of linPointPos
  auto result = m_cfg.propagator->propagate(params, *perigeeSurface, pOptions);

  // parameters of the original track at the linearization point (i.e., the PCA
  // to the vertex, i.e., the point V)
  auto endParams = *result->endParameters;
  BoundVector paramsAtPCA = endParams.parameters();

  // calculate the params for a slightly shifted vertex
  bool shiftRef = false;
  double delta = 0.0001;

  Vector3 globalCoords = endParams.position(gctx);

  for (int i = 0; i < 3; i++) {
    std::cout << "param number " << i << std::endl;

    Vector3 globalCoordsShifted = globalCoords;
    SingleBoundTrackParameters endParamsShifted = endParams;

    if (shiftRef) {
      linPointPos = VectorHelpers::position(linPoint);
      linPointPos(i) -= delta;
      perigeeSurface = Surface::makeShared<PerigeeSurface>(linPointPos);
    } else {
      globalCoordsShifted(i) += delta;
      auto shiftedD0Z0 = *endParams.referenceSurface().globalToLocal(
          gctx, globalCoordsShifted, endParams.momentum());
      endParamsShifted.parameters().template head<2>() = shiftedD0Z0;
    }

    std::cout << "endParams " << paramsAtPCA << std::endl;
    std::cout << "endParams after wiggle: " << endParamsShifted << std::endl;

    intersection = perigeeSurface->intersect(
        gctx, globalCoordsShifted, endParamsShifted.unitDirection(), false);

    std::cout << "distance to perigee " << intersection.intersection.pathLength
              << std::endl;

    // Create propagator options
    pOptions.direction = intersection.intersection.pathLength >= 0
                             ? Direction::Forward
                             : Direction::Backward;

    // Do the propagation to the PCA of linPointPos
    auto result = m_cfg.propagator->propagate(endParamsShifted, *perigeeSurface,
                                              pOptions);

    // parameters of the original track at the linearization point (i.e., the
    // PCA to the vertex, i.e., the point V)
    BoundVector paramsAtPCAShifted = (*result->endParameters).parameters();

    std::cout << "endParams " << paramsAtPCA << std::endl;
    std::cout << "endParams after prop " << paramsAtPCAShifted << std::endl;

    for (int j = 0; j < 6; j++) {
      std::cout << i << " " << j << " numerical diff "
                << (paramsAtPCAShifted[j] - paramsAtPCA[j]) / delta
                << std::endl;
    }
  }
}

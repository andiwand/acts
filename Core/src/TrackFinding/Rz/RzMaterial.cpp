// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RzMaterial.hpp"

#include "Acts/Material/Interactions.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <cmath>

#include "RzState.hpp"

namespace Acts::Experimental::detail::rz {
namespace {

auto materialInterpolator(double momentum) {
  const double x = (std::log(momentum) - RzMaterialTable::logMinP()) /
                   RzMaterialTable::logStep();
  const double xc =
      std::clamp(x, 0., static_cast<double>(RzMaterialTable::kBins - 1));
  const std::uint32_t i =
      std::min(static_cast<std::uint32_t>(xc), RzMaterialTable::kBins - 2);
  const double w = xc - i;
  return [i, w](const std::array<float, RzMaterialTable::kBins>& a) {
    return (1. - w) * a[i] + w * a[i + 1];
  };
}

}  // namespace

bool applyMaterial(State& state, const ParticleHypothesis& hyp,
                   const MaterialSlab& slab, const Vector3& normal,
                   double direction) {
  const Vector3 dir = state.v.segment<3>(eRzDir0);
  const double cosIncidence = std::max(std::abs(normal.dot(dir)), 1e-3);
  const MaterialSlab crossed(
      slab.material(), static_cast<float>(slab.thickness() / cosIncidence));

  const float mass = hyp.mass();
  const float absQ = hyp.absoluteCharge();
  const PdgParticle pdg = hyp.absolutePdg();
  const double qOverP = state.v[eRzQOverP];

  const double theta0 = computeMultipleScatteringTheta0(
      crossed, pdg, mass, static_cast<float>(qOverP), absQ);
  state.pending.varAngle += theta0 * theta0;
  const double sigmaQOverP = computeEnergyLossLandauSigmaQOverP(
      crossed, mass, static_cast<float>(qOverP), absQ);
  state.pending.varQOverP += sigmaQOverP * sigmaQOverP;

  const double dE = computeEnergyLossMean(crossed, pdg, mass,
                                          static_cast<float>(qOverP), absQ);
  const double p = hyp.extractMomentum(qOverP);
  const double e = fastHypot(mass, p) - direction * dE;
  if (e <= mass) {
    return false;
  }
  const double pNew = std::sqrt(e * e - mass * mass);
  state.v[eRzQOverP] = hyp.qOverP(pNew, hyp.extractCharge(qOverP));
  return true;
}

bool applyMaterial(State& state, const ParticleHypothesis& hyp,
                   const RzSurface& surface, std::int32_t band,
                   const Vector3& normal, double direction) {
  if (surface.materialTables.empty()) {
    return applyMaterial(state, hyp, surface.materialBands[band], normal,
                         direction);
  }
  const RzMaterialTable& t = surface.materialTables[band];
  const double qOverP = state.v[eRzQOverP];
  const double p = hyp.extractMomentum(qOverP);
  const auto lerp = materialInterpolator(p);

  const Vector3 dir = state.v.segment<3>(eRzDir0);
  const double factor = 1. / std::max(std::abs(normal.dot(dir)), 1e-3);
  // Highland: theta0^2 ~ t (1 + 0.038 ln(t/X0))^2, so the path factor enters
  // the logarithm as well as the thickness
  const double lnT = t.logThicknessInX0;
  const double highland =
      (1. + 0.038 * (lnT + std::log(factor))) / (1. + 0.038 * lnT);
  state.pending.varAngle += lerp(t.theta0Sq) * factor * highland * highland;
  state.pending.varQOverP += lerp(t.sigmaQOverPSq) * factor;

  const double dE = lerp(t.energyLoss) * factor;
  const double mass = hyp.mass();
  const double e = fastHypot(mass, p) - direction * dE;
  if (e <= mass) {
    return false;
  }
  const double pNew = std::sqrt(e * e - mass * mass);
  state.v[eRzQOverP] = hyp.qOverP(pNew, hyp.extractCharge(qOverP));
  return true;
}

void regainEnergy(State& state, const ParticleHypothesis& hyp,
                  const RzSurface& surface, std::int32_t band,
                  const Vector3& normal) {
  const double qOverP = state.v[eRzQOverP];
  const double p = hyp.extractMomentum(qOverP);
  const Vector3 dir = state.v.segment<3>(eRzDir0);
  const double factor = 1. / std::max(std::abs(normal.dot(dir)), 1e-3);
  double dE = 0.;
  if (!surface.materialTables.empty()) {
    const RzMaterialTable& t = surface.materialTables[band];
    dE = materialInterpolator(p)(t.energyLoss) * factor;
  } else {
    const MaterialSlab& slab = surface.materialBands[band];
    const MaterialSlab crossed(slab.material(),
                               static_cast<float>(slab.thickness() * factor));
    dE =
        computeEnergyLossMean(crossed, hyp.absolutePdg(), hyp.mass(),
                              static_cast<float>(qOverP), hyp.absoluteCharge());
  }
  const double mass = hyp.mass();
  const double e = fastHypot(mass, p) + dE;
  const double pNew = std::sqrt(e * e - mass * mass);
  state.v[eRzQOverP] = hyp.qOverP(pNew, hyp.extractCharge(qOverP));
}

void ProcessNoise::apply(RzMatrix& covariance, const Vector3& d,
                         const Vector3& normal) {
  ProcessNoise& p = *this;
  if (p.empty()) {
    return;
  }
  const SquareMatrix3 transverse =
      SquareMatrix3::Identity() - d * d.transpose();
  covariance.block<3, 3>(eRzDir0, eRzDir0) += p.varAngle * transverse;
  covariance.block<3, 3>(eRzPos0, eRzPos0) += p.varPosition * transverse;
  covariance.block<3, 3>(eRzPos0, eRzDir0) += p.covAnglePosition * transverse;
  covariance.block<3, 3>(eRzDir0, eRzPos0) += p.covAnglePosition * transverse;
  covariance(eRzQOverP, eRzQOverP) += p.varQOverP;
  // Project position noise onto the surface along the track:
  // apply I - d n^T / (n.d) to the position rows and columns.
  const Vector3 dOver = d / normal.dot(d);
  for (std::uint32_t c = 0; c < eRzSize; ++c) {
    const double nc = normal.x() * covariance(eRzPos0, c) +
                      normal.y() * covariance(eRzPos1, c) +
                      normal.z() * covariance(eRzPos2, c);
    covariance(eRzPos0, c) -= dOver.x() * nc;
    covariance(eRzPos1, c) -= dOver.y() * nc;
    covariance(eRzPos2, c) -= dOver.z() * nc;
  }
  for (std::uint32_t r = 0; r < eRzSize; ++r) {
    const double nr = normal.x() * covariance(r, eRzPos0) +
                      normal.y() * covariance(r, eRzPos1) +
                      normal.z() * covariance(r, eRzPos2);
    covariance(r, eRzPos0) -= dOver.x() * nr;
    covariance(r, eRzPos1) -= dOver.y() * nr;
    covariance(r, eRzPos2) -= dOver.z() * nr;
  }
  p = ProcessNoise{};
}

}  // namespace Acts::Experimental::detail::rz

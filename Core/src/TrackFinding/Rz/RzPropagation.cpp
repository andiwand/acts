// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RzPropagation.hpp"

#include "Acts/Material/Interactions.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <cmath>

namespace Acts::Experimental::detail::rz {

/// The normal of an RZ surface at a position on it
Vector3 surfaceNormal(const RzSurface& surface, const RzVector& v) {
  if (surface.shape == RzShape::Disc) {
    return Vector3::UnitZ();
  }
  const double r = fastHypot(v[eRzPos0], v[eRzPos1]);
  return Vector3(v[eRzPos0] / r, v[eRzPos1] / r, 0.);
}

/// The coordinate an RZ surface extends in, at a position
double alongCoordinate(const RzSurface& surface, const RzVector& v) {
  return surface.shape == RzShape::Cylinder ? v[eRzPos2]
                                            : fastHypot(v[eRzPos0], v[eRzPos1]);
}
/// A surface's field at a crossing, or the constant one
double bzAt(const RzSurface& surface, double along, double fallback) {
  return surface.bzAt(along).value_or(fallback);
}

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

/// Apply the radial kick and project disc crossings back onto their plane.
void radialKick(RzVector& v, const RzVector& from, double s, double br0,
                double br1, const RzSurface& surface, Vector3& qopPosition,
                Vector3& qopDirection) {
  rzRadialKick(v, from, s, br0, br1, qopPosition, qopDirection);
  if (surface.shape == RzShape::Disc && v[eRzDir2] != 0.) {
    const double back = (surface.refCoord - v[eRzPos2]) / v[eRzDir2];
    v.segment<3>(eRzPos0) += back * v.segment<3>(eRzDir0);
  }
}
}  // namespace

bool Stepper::applyMaterial(State& state, const MaterialSlab& slab,
                            const Vector3& normal, double direction) const {
  const Vector3 dir = state.v.segment<3>(eRzDir0);
  const double cosIncidence = std::max(std::abs(normal.dot(dir)), 1e-3);
  const MaterialSlab crossed(
      slab.material(), static_cast<float>(slab.thickness() / cosIncidence));

  const ParticleHypothesis& hyp = m_cfg.particleHypothesis;
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

bool Stepper::applyMaterial(State& state, const RzSurface& surface,
                            std::int32_t band, const Vector3& normal,
                            double direction) const {
  if (surface.materialTables.empty()) {
    return applyMaterial(state, surface.materialBands[band], normal, direction);
  }
  const RzMaterialTable& t = surface.materialTables[band];
  const ParticleHypothesis& hyp = m_cfg.particleHypothesis;
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

void Stepper::regainEnergy(State& state, const RzSurface& surface,
                           std::int32_t band, const Vector3& normal) const {
  const ParticleHypothesis& hyp = m_cfg.particleHypothesis;
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

void Stepper::materialise(State& state, const Vector3& normal) const {
  Pending& p = state.pending;
  if (p.empty()) {
    return;
  }
  const Vector3 d = state.v.segment<3>(eRzDir0);
  const SquareMatrix3 transverse =
      SquareMatrix3::Identity() - d * d.transpose();
  state.c.block<3, 3>(eRzDir0, eRzDir0) += p.varAngle * transverse;
  state.c.block<3, 3>(eRzPos0, eRzPos0) += p.varPosition * transverse;
  state.c.block<3, 3>(eRzPos0, eRzDir0) += p.covAnglePosition * transverse;
  state.c.block<3, 3>(eRzDir0, eRzPos0) += p.covAnglePosition * transverse;
  state.c(eRzQOverP, eRzQOverP) += p.varQOverP;
  // Project position noise onto the surface along the track:
  // apply I - d n^T / (n.d) to the position rows and columns.
  const Vector3 dOver = d / normal.dot(d);
  for (std::uint32_t c = 0; c < eRzSize; ++c) {
    const double nc = normal.x() * state.c(eRzPos0, c) +
                      normal.y() * state.c(eRzPos1, c) +
                      normal.z() * state.c(eRzPos2, c);
    state.c(eRzPos0, c) -= dOver.x() * nc;
    state.c(eRzPos1, c) -= dOver.y() * nc;
    state.c(eRzPos2, c) -= dOver.z() * nc;
  }
  for (std::uint32_t r = 0; r < eRzSize; ++r) {
    const double nr = normal.x() * state.c(r, eRzPos0) +
                      normal.y() * state.c(r, eRzPos1) +
                      normal.z() * state.c(r, eRzPos2);
    state.c(r, eRzPos0) -= dOver.x() * nr;
    state.c(r, eRzPos1) -= dOver.y() * nr;
    state.c(r, eRzPos2) -= dOver.z() * nr;
  }
  p = Pending{};
}

Vector3 Stepper::land(State& state, RzVector& landed, double path,
                      const RzSurface& surface, double along) const {
  const double brLanded = surface.brAt(along);
  if (m_cfg.radialField) {
    radialKick(landed, state.v, path, state.br, brLanded, surface,
               state.brQopPosition, state.brQopDirection);
  }
  state.v = landed;
  const Vector3 normal = surfaceNormal(surface, state.v);
  state.travel(path);
  state.pending.advance(path);
  state.bz = bzAt(surface, along, m_bz);
  state.br = brLanded;
  return normal;
}

void Navigator::initialize(NavigationState& navigation,
                           const RzVector& start) const {
  navigation = NavigationState{};
  const RzLayout& layout = m_layout;
  const double r0 = fastHypot(start[eRzPos0], start[eRzPos1]);
  while (navigation.cyl < layout.cylinders.size() &&
         layout.surfaces[layout.cylinders[navigation.cyl]].refCoord <= r0) {
    ++navigation.cyl;
  }
  const bool forward = start[eRzDir2] >= 0.;
  navigation.discStep = forward ? 1 : -1;
  // The first disc ahead, using the sorted disc positions.
  navigation.disc =
      forward ? std::ranges::upper_bound(layout.discCoord, start[eRzPos2]) -
                    layout.discCoord.begin()
              : std::ranges::lower_bound(layout.discCoord, start[eRzPos2]) -
                    layout.discCoord.begin() - 1;
}

std::optional<Target> Navigator::next(NavigationState& navigation,
                                      const State& state,
                                      double maxPath) const {
  const RzLayout& layout = m_layout;
  const RzHelix helix{state.bz};
  const double kappa = std::abs(helix.kappa(state.v));
  const auto discValid = [&]() {
    return navigation.disc >= 0 &&
           navigation.disc < static_cast<std::int32_t>(layout.discs.size());
  };
  // Cache state-dependent quantities across navigation probes at this stop.
  const double dTransverse = fastHypot(state.v[eRzDir0], state.v[eRzDir1]);
  const double pz = state.v[eRzPos2];
  const double dz = state.v[eRzDir2];
  const double invDz = dz != 0. ? 1. / dz : 0.;
  const double px = state.v[eRzPos0];
  const double py = state.v[eRzPos1];
  const double dxDir = state.v[eRzDir0];
  const double dyDir = state.v[eRzDir1];
  const double halfKappaT = 0.5 * kappa * dTransverse;
  std::optional<double> sDisc;
  while (navigation.discsLeft && discValid()) {
    const std::size_t di = static_cast<std::size_t>(navigation.disc);
    const double sTry = (layout.discCoord[di] - pz) * invDz;
    if (sTry <= 0.) {
      // behind us
      sDisc.reset();
      navigation.disc += navigation.discStep;
      continue;
    }
    if (sTry > maxPath) {
      // z grows monotonically, so every disc beyond is out of reach too
      sDisc.reset();
      navigation.discsLeft = false;
      break;
    }
    // Reject discs outside the straight-line radius plus sagitta.
    // Compare squared radii; check exact bounds after landing.
    const double xs = px + dxDir * sTry;
    const double ys = py + dyDir * sTry;
    const double r2 = xs * xs + ys * ys;
    const double sagitta = halfKappaT * sTry * sTry;
    const double lo = layout.discMin[di] - sagitta;
    const double hi = layout.discMax[di] + sagitta;
    if ((lo > 0. && r2 < lo * lo) || r2 > hi * hi) {
      sDisc.reset();
      navigation.disc += navigation.discStep;
      continue;
    }
    sDisc = sTry;
    break;
  }

  std::optional<double> sCyl;
  if (navigation.cylindersLeft && navigation.cyl < layout.cylinders.size()) {
    const double rCyl = layout.cylCoord[navigation.cyl];
    bool tryCylinder = true;
    if (navigation.inEndcap && sDisc.has_value()) {
      // Skip the cylinder solve if the disc lies inside it, allowing sagitta.
      const double rAtDisc =
          fastHypot(px + dxDir * *sDisc, py + dyDir * *sDisc);
      const double sagitta = halfKappaT * *sDisc * *sDisc;
      tryCylinder = rAtDisc + sagitta + 1. >= rCyl;
    }
    if (tryCylinder && navigation.cylCached == navigation.cyl) {
      // Rejecting a disc leaves the cached cylinder intersection valid.
      sCyl = navigation.cylCachedPath;
    } else if (tryCylinder) {
      sCyl = helix.pathToCylinder(state.v, rCyl);
      if (!sCyl.has_value() || *sCyl > maxPath) {
        // the helix never reaches this radius, so none beyond it either
        sCyl.reset();
        navigation.cylindersLeft = false;
      }
      navigation.cylCached = static_cast<std::uint32_t>(navigation.cyl);
      navigation.cylCachedPath = sCyl;
    }
  }
  const bool takeCyl =
      sCyl.has_value() && (!sDisc.has_value() || *sCyl <= *sDisc);
  if (!sCyl.has_value() && !sDisc.has_value()) {
    return std::nullopt;
  }
  const double s = takeCyl ? *sCyl : *sDisc;
  const std::uint32_t surfaceIndex = takeCyl ? layout.cylinders[navigation.cyl]
                                             : layout.discs[navigation.disc];
  if (takeCyl) {
    ++navigation.cyl;
  } else {
    navigation.disc += navigation.discStep;
  }

  return Target{surfaceIndex, s, takeCyl};
}

void Propagator::initialize(PropagationState& propagation,
                            const RzVector& start,
                            std::uint32_t startSurface) const {
  propagation.startSurface = startSurface;
  propagation.status = PropagationStatus::Active;
  m_navigator.initialize(propagation.navigation, start);
}

std::optional<Crossing> Propagator::advance(State& state,
                                            PropagationState& propagation,
                                            RzTrackCandidate& candidate) const {
  if (propagation.status != PropagationStatus::Active) {
    return std::nullopt;
  }
  NavigationState& navigation = propagation.navigation;
  const RzLayout& layout = m_layout;
  while (true) {
    const RzHelix helix{state.bz};
    const double kappa = std::abs(helix.kappa(state.v));
    const double maxPath =
        kappa > 0. ? (m_cfg.maxTurningAngle - state.turned) / kappa
                   : 2. * (layout.escapeRadius + layout.escapeHalfZ);
    if (maxPath <= 0.) {
      propagation.status = PropagationStatus::TurningLimit;
      return std::nullopt;
    }
    const auto target = m_navigator.next(navigation, state, maxPath);
    if (!target) {
      propagation.status = PropagationStatus::NoTarget;
      return std::nullopt;
    }
    const double s = target->path;
    const std::uint32_t surfaceIndex = target->surface;
    const RzSurface& surface = layout.surfaces[surfaceIndex];
    // Check surface bounds before committing transport or covariance work.
    RzVector landed = state.v;
    helix.step(landed, s);
    const double along = alongCoordinate(surface, landed);
    if (!surface.contains(along)) {
      // the state itself stays put; the next candidate is measured from here
      continue;
    }
    ++candidate.stops;
    const std::uint32_t stop =
        static_cast<std::uint32_t>(candidate.stopSurfaces.size());
    candidate.stopSurfaces.push_back(surfaceIndex);
    candidate.stopPaths.push_back(s);
    candidate.stopAlong.push_back(along);

    const Vector3 normal = m_stepper.land(state, landed, s, surface, along);
    Navigator::moved(navigation, target->cylinder);
    state.turned += std::abs(helix.kappa(state.v)) * s;

    const double r = fastHypot(state.v[eRzPos0], state.v[eRzPos1]);
    if (r > layout.escapeRadius ||
        std::abs(state.v[eRzPos2]) > layout.escapeHalfZ ||
        state.turned > m_cfg.maxTurningAngle) {
      propagation.status = state.turned > m_cfg.maxTurningAngle
                               ? PropagationStatus::TurningLimit
                               : PropagationStatus::Escape;
      return std::nullopt;
    }

    if (surfaceIndex == propagation.startSurface) {
      continue;
    }

    if (m_cfg.applyMaterial) {
      if (const std::int32_t band = surface.materialBandAt(along);
          band >= 0 && !m_stepper.applyMaterial(state, surface, band, normal)) {
        propagation.status = PropagationStatus::MaterialFailure;
        return std::nullopt;
      }
    }

    if (surface.layer == kRzNone) {
      continue;
    }
    state.moveCovariance(RzHelix{state.anchorBz}, normal);
    m_stepper.materialise(state, normal);
    return Crossing{surface.layer, stop};
  }
}

}  // namespace Acts::Experimental::detail::rz

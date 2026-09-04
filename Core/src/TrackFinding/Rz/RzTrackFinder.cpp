// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/Rz/RzTrackFinder.hpp"

#include "Acts/Material/Interactions.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Acts::Experimental {

namespace {

/// `sqrt(x^2 + y^2)` without hypot's range handling, which the finder does
/// not need and pays for at every stop
double norm2(double x, double y) {
  return std::sqrt(x * x + y * y);
}

/// The normal of an RZ surface at a position on it
Vector3 surfaceNormal(const RzSurface& surface, const RzVector& v) {
  if (surface.shape == RzShape::Disc) {
    return Vector3::UnitZ();
  }
  const double r = norm2(v[eRzPos0], v[eRzPos1]);
  return Vector3(v[eRzPos0] / r, v[eRzPos1] / r, 0.);
}

/// The coordinate an RZ surface extends in, at a position
double alongCoordinate(const RzSurface& surface, const RzVector& v) {
  return surface.shape == RzShape::Cylinder ? v[eRzPos2]
                                            : norm2(v[eRzPos0], v[eRzPos1]);
}

}  // namespace

RzTrackFinder::RzTrackFinder(const RzTrackFinderConfig& config,
                             const RzLayout& layout, double bz)
    : m_cfg(config), m_layout(&layout), m_helix{bz} {}

bool RzTrackFinder::applyMaterial(State& state, const MaterialSlab& slab,
                                  const Vector3& normal,
                                  double direction) const {
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
  const double e = norm2(mass, p) - direction * dE;
  if (e <= mass) {
    return false;
  }
  const double pNew = std::sqrt(e * e - mass * mass);
  state.v[eRzQOverP] = hyp.qOverP(pNew, hyp.extractCharge(qOverP));
  return true;
}

bool RzTrackFinder::applyMaterial(State& state, const RzSurface& surface,
                                  int band, const Vector3& normal,
                                  double direction) const {
  if (surface.materialTables.empty()) {
    return applyMaterial(state, surface.materialBands[band], normal, direction);
  }
  const RzMaterialTable& t = surface.materialTables[band];
  const ParticleHypothesis& hyp = m_cfg.particleHypothesis;
  const double qOverP = state.v[eRzQOverP];
  const double p = hyp.extractMomentum(qOverP);
  const double x =
      (std::log(p) - RzMaterialTable::logMinP()) / RzMaterialTable::logStep();
  const double xc =
      std::clamp(x, 0., static_cast<double>(RzMaterialTable::kBins - 1));
  const std::uint32_t i =
      std::min(static_cast<std::uint32_t>(xc), RzMaterialTable::kBins - 2);
  const double w = xc - i;
  auto lerp = [&](const std::array<float, RzMaterialTable::kBins>& a) {
    return (1. - w) * a[i] + w * a[i + 1];
  };

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
  const double e = norm2(mass, p) - direction * dE;
  if (e <= mass) {
    return false;
  }
  const double pNew = std::sqrt(e * e - mass * mass);
  state.v[eRzQOverP] = hyp.qOverP(pNew, hyp.extractCharge(qOverP));
  return true;
}

void RzTrackFinder::materialise(State& state, const Vector3& normal) const {
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
  // a displacement normal to the surface is, moved along the track, one in
  // its plane; the covariance has to say so before the next transport
  const SquareMatrix3 project =
      SquareMatrix3::Identity() - d * normal.transpose() / normal.dot(d);
  state.c.block<3, eRzSize>(eRzPos0, 0) =
      project * state.c.block<3, eRzSize>(eRzPos0, 0);
  state.c.block<eRzSize, 3>(0, eRzPos0) =
      state.c.block<eRzSize, 3>(0, eRzPos0) * project.transpose();
  p = Pending{};
}

std::optional<RzTrackFinder::Evaluation> RzTrackFinder::evaluate(
    const State& state, const RzMeasurement& m, bool gate) const {
  // the straight-line crossing of the module plane first: it is the Newton
  // start anyway, and enough for the gate
  const Vector3 p0 = state.v.segment<3>(eRzPos0);
  const Vector3 d0 = state.v.segment<3>(eRzDir0);
  const double along = m.normal.dot(d0);
  if (std::abs(along) < 1e-9) {
    return std::nullopt;
  }
  const double s0 = m.normal.dot(m.position - p0) / along;
  if (std::abs(s0) > m_cfg.maxModuleDistance) {
    return std::nullopt;
  }
  // Most candidates in the window are nowhere near: a straight-line residual
  // against the covariance at the stop, widened by what the direction
  // uncertainty does over the module distance, is enough to drop them before
  // the helix is solved and the transport Jacobian is built
  if (gate) {
    const Vector3 d = m.position - (p0 + s0 * d0);
    const double ru = m.u.dot(d);
    const double rv = m.v.dot(d);
    const auto cPos = state.c.block<3, 3>(eRzPos0, eRzPos0);
    const double spread =
        state.c.block<3, 3>(eRzDir0, eRzDir0).trace() * s0 * s0 +
        state.pending.varPosition;
    const double su = m.u.dot(cPos * m.u) + m.cov00 + spread;
    double chi2 = ru * ru / su;
    if (m.dim == 2) {
      const double sv = m.v.dot(cPos * m.v) + m.cov11 + spread;
      chi2 += rv * rv / sv;
    }
    if (chi2 > m_cfg.gateFactor * m_cfg.chi2Cut) {
      return std::nullopt;
    }
  }
  const std::optional<double> s =
      m_helix.pathToPlane(state.v, m.position, m.normal);
  if (!s.has_value() || std::abs(*s) > m_cfg.maxModuleDistance) {
    return std::nullopt;
  }
  RzVector w = state.v;
  m_helix.step(w, *s);
  const Vector3 d = m.position - w.segment<3>(eRzPos0);
  const double ru = m.u.dot(d);
  const double rv = m.v.dot(d);
  if (m.dim == 1 && std::abs(rv) > m.halfV + m_cfg.stripMargin) {
    return std::nullopt;
  }
  // S = H C H^T + R with H the two frame axes on the position block of the
  // covariance moved to the module: the rows of H J, and only they, are
  // formed, and C (H J)^T is what the update needs too
  Evaluation e;
  const Eigen::Matrix<double, 3, eRzSize> jPos =
      m_helix.stepJacobianOnto(state.v, *s, w, m.normal).positionRows();
  const Eigen::Matrix<double, 1, eRzSize> hu = m.u.transpose() * jPos;
  e.ch.col(0) = state.c * hu.transpose();
  const double s00 = hu.dot(e.ch.col(0)) + m.cov00;
  e.sInv.setZero();
  if (m.dim == 1) {
    e.ch.col(1).setZero();
    e.sInv(0, 0) = 1. / s00;
    e.chi2 = ru * ru * e.sInv(0, 0);
  } else {
    const Eigen::Matrix<double, 1, eRzSize> hv = m.v.transpose() * jPos;
    e.ch.col(1) = state.c * hv.transpose();
    const double s01 = hv.dot(e.ch.col(0)) + m.cov01;
    const double s11 = hv.dot(e.ch.col(1)) + m.cov11;
    const double det = s00 * s11 - s01 * s01;
    if (det <= 0.) {
      return std::nullopt;
    }
    e.sInv(0, 0) = s11 / det;
    e.sInv(0, 1) = -s01 / det;
    e.sInv(1, 0) = e.sInv(0, 1);
    e.sInv(1, 1) = s00 / det;
    e.chi2 = ru * ru * e.sInv(0, 0) + 2. * ru * rv * e.sInv(0, 1) +
             rv * rv * e.sInv(1, 1);
  }
  e.residual << ru, rv;
  return e;
}

void RzTrackFinder::update(State& state, const Evaluation& e) const {
  const Eigen::Matrix<double, eRzSize, 2> k = e.ch * e.sInv;
  state.v += k * e.residual;
  state.v.segment<3>(eRzDir0).normalize();
  state.anchor = state.v;
  state.c -= k * e.ch.transpose();
  state.c = 0.5 * (state.c + state.c.transpose()).eval();
}

std::optional<double> RzTrackFinder::pathBackward(const RzVector& v,
                                                  const RzSurface& surface,
                                                  double guess) const {
  if (surface.shape == RzShape::Disc) {
    const double dz = v[eRzDir2];
    return dz != 0. ? std::optional((surface.refCoord - v[eRzPos2]) / dz)
                    : std::nullopt;
  }
  // Newton from the forward path: the state has moved by an update or two
  // since, so the root is next to the guess
  double s = guess;
  for (int i = 0; i < 6; ++i) {
    RzVector w = v;
    m_helix.step(w, s);
    const double f = w[eRzPos0] * w[eRzPos0] + w[eRzPos1] * w[eRzPos1] -
                     surface.refCoord * surface.refCoord;
    const double df = 2. * (w[eRzPos0] * w[eRzDir0] + w[eRzPos1] * w[eRzDir1]);
    if (df == 0.) {
      break;
    }
    const double ds = f / df;
    s -= ds;
    if (std::abs(ds) < 1e-9) {
      // a root far from the guess is the other crossing of the circle
      if (std::abs(s - guess) < std::max(20., 0.2 * std::abs(guess))) {
        return s;
      }
      break;
    }
  }
  // the closed form the other way round, as a fallback
  const RzVector r = RzHelix::reversed(v);
  const std::optional<double> back =
      m_helix.pathToCylinder(r, surface.refCoord);
  if (!back.has_value()) {
    return std::nullopt;
  }
  return -*back;
}

bool RzTrackFinder::onModule(std::uint32_t layerIndex,
                             const State& state) const {
  const RzLayout& layout = *m_layout;
  const RzLayer& layer = layout.layers[layerIndex];
  const RzSurface& surface = layout.surfaces[layer.surface];
  const RzVector& v = state.v;
  const double r = norm2(v[eRzPos0], v[eRzPos1]);
  const double phi = std::atan2(v[eRzPos1], v[eRzPos0]);
  const double along = alongCoordinate(surface, v);
  const double room = layer.maxHalfExtent + m_cfg.moduleEdgeTolerance;
  // the straight-line crossing of each module plane, with the sagitta over
  // the module distance as the allowance on top of the edge tolerance
  const Vector3 p = v.segment<3>(eRzPos0);
  const Vector3 dir = v.segment<3>(eRzDir0);
  const double kappa = std::abs(m_helix.kappa(v));
  bool hit = false;
  RzMeasurementGrid::visitBins(
      layout, layerIndex, phi, along, room / r, room, [&](std::uint32_t b) {
        if (hit) {
          return;
        }
        for (std::uint32_t i = layout.moduleBinStart[b];
             i < layout.moduleBinStart[b + 1]; ++i) {
          const RzModule& m = layout.modules[layout.moduleOrder[i]];
          const double alongNormal = m.normal.dot(dir);
          if (std::abs(alongNormal) < 1e-9) {
            continue;
          }
          const double s = m.normal.dot(m.center - p) / alongNormal;
          if (std::abs(s) > m_cfg.maxModuleDistance) {
            continue;
          }
          const Vector3 d = p + s * dir - m.center;
          const double tolerance =
              m_cfg.moduleEdgeTolerance + 0.5 * kappa * s * s;
          if (std::abs(m.u.dot(d)) <= m.halfU + tolerance &&
              std::abs(m.v.dot(d)) <= m.halfV + tolerance) {
            hit = true;
            return;
          }
        }
      });
  return hit;
}

std::uint32_t RzTrackFinder::searchLayer(const RzMeasurementGrid& grid,
                                         std::uint32_t layerIndex,
                                         std::uint32_t stop, State& state,
                                         RzTrackCandidate& candidate) const {
  const RzLayer& layer = m_layout->layers[layerIndex];
  const RzSurface& surface = m_layout->surfaces[layer.surface];
  const RzVector& v = state.v;
  const double x = v[eRzPos0];
  const double y = v[eRzPos1];
  const double r = norm2(x, y);
  const double phi = std::atan2(y, x);

  // the window from the position uncertainty in the surface's coordinates
  const auto cPos = state.c.block<3, 3>(eRzPos0, eRzPos0);
  const Vector3 tangent(-y / r, x / r, 0.);
  const Vector3 radial(x / r, y / r, 0.);
  const double varPending = state.pending.varPosition;
  const double sigmaPhi =
      std::sqrt(std::max(0., tangent.dot(cPos * tangent) + varPending));
  const Vector3 alongDir =
      surface.shape == RzShape::Cylinder ? Vector3::UnitZ() : radial;
  const double sigmaAlong =
      std::sqrt(std::max(0., alongDir.dot(cPos * alongDir) + varPending));
  // modules off the RZ surface are met at a different (phi, along) than the
  // stop: by the layer's thickness times the track's slope in each
  const Vector3 d = v.segment<3>(eRzDir0);
  const double dRadial = std::max(std::abs(radial.dot(d)), 1e-6);
  const double dTangent = std::abs(tangent.dot(d));
  double thicknessPhi = 0.;
  double thicknessAlong = 0.;
  if (surface.shape == RzShape::Cylinder) {
    thicknessAlong = layer.halfThickness * std::abs(d.z()) / dRadial;
    thicknessPhi = layer.halfThickness * dTangent / dRadial;
  } else {
    const double dz = std::max(std::abs(d.z()), 1e-6);
    thicknessAlong = layer.halfThickness * dRadial / dz;
    thicknessPhi = layer.halfThickness * dTangent / dz;
  }
  const double halfPhi =
      (m_cfg.windowSigmas * sigmaPhi + m_cfg.windowMin + thicknessPhi) / r;
  const double halfAlong = m_cfg.windowSigmas * sigmaAlong + m_cfg.windowMin +
                           grid.stripHalfV(layerIndex) + thicknessAlong;
  const double along = alongCoordinate(surface, v);

  std::uint32_t accepted = 0;
  std::vector<std::uint32_t> usedModules;
  for (std::uint32_t round = 0; round < m_cfg.maxMeasurementsPerLayer;
       ++round) {
    std::uint32_t bestIndex = kRzNone;
    Evaluation best;
    best.chi2 = std::numeric_limits<double>::max();
    grid.visit(
        layerIndex, phi, along, halfPhi, halfAlong, [&](std::uint32_t i) {
          const RzMeasurement& m = grid.entry(i);
          if (std::ranges::find(usedModules, m.module) != usedModules.end()) {
            return;
          }
          ++candidate.candidatesTested;
          const std::optional<Evaluation> e = evaluate(state, m);
          if (e.has_value() && e->chi2 < best.chi2) {
            bestIndex = i;
            best = *e;
          }
        });
    if (bestIndex == kRzNone || best.chi2 > m_cfg.chi2Cut) {
      break;
    }
    const RzMeasurement& m = grid.entry(bestIndex);
    update(state, best);
    const std::uint32_t forwardState =
        static_cast<std::uint32_t>(candidate.forwardStates.size());
    candidate.forwardStates.emplace_back(state.v, state.c);
    candidate.hits.push_back(
        {layerIndex, bestIndex, stop, forwardState, best.chi2});
    candidate.chi2 += best.chi2;
    usedModules.push_back(m.module);
    ++accepted;
  }
  return accepted;
}

void RzTrackFinder::backwardPass(const RzMeasurementGrid& grid,
                                 const State& forward,
                                 RzTrackCandidate& candidate) const {
  // where to start: the last measurement, or the outermost of the innermost
  // `backwardLayers` of them
  auto hit = candidate.hits.rbegin();
  while (hit != candidate.hits.rend() && hit->isHole()) {
    ++hit;
  }
  if (hit == candidate.hits.rend()) {
    candidate.backwardFailure = 4;
    return;
  }
  const bool partial =
      m_cfg.backwardLayers > 0 && candidate.measurements > m_cfg.backwardLayers;
  if (partial) {
    std::uint32_t seen = 0;
    for (auto it = candidate.hits.begin(); it != candidate.hits.end(); ++it) {
      if (it->isHole()) {
        continue;
      }
      if (++seen == m_cfg.backwardLayers) {
        hit = std::make_reverse_iterator(it + 1);
        break;
      }
    }
  }
  const RzVector& startV =
      partial ? candidate.forwardStates[hit->forwardState].first : forward.v;
  const RzMatrix& startC =
      partial ? candidate.forwardStates[hit->forwardState].second : forward.c;

  // Forget what the forward filter knew, so that the inner end carries the
  // hits' own precision: uncorrelated, or a transport turns the correlations
  // into q/p shifts from the first residual on, and with the null space a
  // free state has, no variance along the direction and none normal to the
  // plane the state sits on. Loose rather than infinite: with a radian of
  // direction uncertainty two strips millimetres apart would set the
  // direction.
  State state;
  state.v = startV;
  state.anchor = startV;
  {
    const Vector3 d = startV.segment<3>(eRzDir0);
    const Vector3 n = grid.entry(hit->measurement).normal;
    const double varPos =
        startC.block<3, 3>(eRzPos0, eRzPos0).trace() * m_cfg.backwardInflation;
    const double varDir =
        startC.block<3, 3>(eRzDir0, eRzDir0).trace() * m_cfg.backwardInflation;
    state.c.setZero();
    state.c.block<3, 3>(eRzPos0, eRzPos0) =
        varPos * (SquareMatrix3::Identity() - n * n.transpose());
    state.c.block<3, 3>(eRzDir0, eRzDir0) =
        varDir * (SquareMatrix3::Identity() - d * d.transpose());
    if (partial) {
      // the momentum is the whole track's: the forward filter's final word,
      // brought in to this hit by regaining what the stops in between took
      state.v[eRzQOverP] = forward.v[eRzQOverP];
      state.c(eRzQOverP, eRzQOverP) =
          forward.c(eRzQOverP, eRzQOverP) * m_cfg.backwardQOverPScale;
      if (m_cfg.applyMaterial) {
        auto last = candidate.hits.rbegin();
        while (last->isHole()) {
          ++last;
        }
        const std::uint32_t outerStop = last->stop;
        const std::uint32_t innerStop = hit->stop;
        for (std::uint32_t j = outerStop; j != innerStop && j != kRzNone; --j) {
          const RzSurface& surface =
              m_layout->surfaces[candidate.stopSurfaces[j]];
          const int band = surface.materialBandAt(candidate.stopAlong[j]);
          if (band < 0) {
            continue;
          }
          const Vector3 normal = surfaceNormal(surface, state.v);
          State only;
          only.v = state.v;
          applyMaterial(only, surface, band, normal, -1.);
          state.v[eRzQOverP] = only.v[eRzQOverP];
        }
      }
    } else {
      state.c(eRzQOverP, eRzQOverP) =
          startC(eRzQOverP, eRzQOverP) * m_cfg.backwardInflation;
    }
  }

  // the forward pass at a stop went: transport, material, materialise,
  // update; replayed from the last hit inwards that is: update the stop's
  // hits, its material, transport to the stop before, materialise
  auto updateHitsAt = [&](std::uint32_t stop) {
    // hits are stored outward, so this stop's hits are the next ones inward
    while (hit != candidate.hits.rend() && hit->stop == stop) {
      if (!hit->isHole()) {
        const RzMeasurement& m = grid.entry(hit->measurement);
        const std::optional<Evaluation> e = evaluate(state, m, false);
        if (!e.has_value()) {
          return false;
        }
        update(state, *e);
      }
      ++hit;
    }
    return true;
  };

  std::uint32_t stop = hit->stop;
  if (stop != kRzNone) {
    for (std::ptrdiff_t j = static_cast<std::ptrdiff_t>(stop); j >= 0; --j) {
      const RzSurface& surface = m_layout->surfaces[candidate.stopSurfaces[j]];
      if (static_cast<std::uint32_t>(j) != stop) {
        const std::optional<double> s =
            pathBackward(state.v, surface, -candidate.stopPaths[j + 1]);
        if (!s.has_value()) {
          candidate.backwardFailure = 1;
          return;
        }
        m_helix.step(state.v, *s);
        const Vector3 normal = surfaceNormal(surface, state.v);
        state.travel(*s);
        state.pending.advance(-*s);
        // The covariance is needed where there is something to update, and
        // where scattering is waiting to be put in: materialising it at the
        // stop it belongs to, and letting the Jacobians carry it from there,
        // keeps the lever arms exact for the parameters this pass is for.
        if (!state.pending.empty() ||
            (hit != candidate.hits.rend() &&
             hit->stop == static_cast<std::uint32_t>(j))) {
          state.moveCovariance(m_helix, normal);
          materialise(state, normal);
        }
      }
      if (!updateHitsAt(static_cast<std::uint32_t>(j))) {
        candidate.backwardFailure = 3;
        return;
      }
      if (m_cfg.applyMaterial) {
        const Vector3 normal = surfaceNormal(surface, state.v);
        if (const int band =
                surface.materialBandAt(alongCoordinate(surface, state.v));
            band >= 0 && !applyMaterial(state, surface, band, normal, -1.)) {
          candidate.backwardFailure = 2;
          return;
        }
      }
    }
  }
  // what the track started with, found before any stop
  while (hit != candidate.hits.rend()) {
    if (!hit->isHole()) {
      const RzMeasurement& m = grid.entry(hit->measurement);
      const std::optional<double> s =
          m_helix.pathToPlane(state.v, m.position, m.normal);
      if (!s.has_value()) {
        candidate.backwardFailure = 1;
        return;
      }
      m_helix.step(state.v, *s);
      state.travel(*s);
      state.moveCovariance(m_helix, m.normal);
      state.pending.advance(std::abs(*s));
      materialise(state, m.normal);
      const std::optional<Evaluation> e = evaluate(state, m, false);
      if (!e.has_value()) {
        candidate.backwardFailure = 3;
        return;
      }
      update(state, *e);
    }
    ++hit;
  }
  if (partial && m_cfg.backwardQOverPScale == 0.) {
    state.c(eRzQOverP, eRzQOverP) = forward.c(eRzQOverP, eRzQOverP);
  }
  candidate.innerParameters = state.v;
  candidate.innerCovariance = state.c;
  candidate.hasInner = true;
}

bool RzTrackFinder::findTrack(const RzMeasurementGrid& grid,
                              const RzVector& start,
                              const RzMatrix& startCovariance,
                              std::uint32_t startModule,
                              RzTrackCandidate& candidate) const {
  candidate.clear();
  State state;
  state.v = start;
  state.anchor = start;
  state.c = startCovariance;

  const RzLayout& layout = *m_layout;
  std::uint32_t startSurface = kRzNone;
  if (startModule != kRzNone) {
    const std::uint32_t layer = layout.modules[startModule].layer;
    startSurface = layout.layers[layer].surface;
    if (searchLayer(grid, layer, kRzNone, state, candidate) == 0 &&
        onModule(layer, state)) {
      candidate.hits.push_back({layer, kRzNone, kRzNone, kRzNone, 0.});
    }
  }

  // navigation cursors: the next cylinder outward and the next disc along z
  const double r0 = norm2(state.v[eRzPos0], state.v[eRzPos1]);
  std::size_t cyl = 0;
  while (cyl < layout.cylinders.size() &&
         layout.surfaces[layout.cylinders[cyl]].refCoord <= r0) {
    ++cyl;
  }
  const bool forward = state.v[eRzDir2] >= 0.;
  const int discStep = forward ? 1 : -1;
  std::ptrdiff_t disc =
      forward ? 0 : static_cast<std::ptrdiff_t>(layout.discs.size()) - 1;
  auto discValid = [&]() {
    return disc >= 0 && disc < static_cast<std::ptrdiff_t>(layout.discs.size());
  };
  while (discValid()) {
    const double z = layout.surfaces[layout.discs[disc]].refCoord;
    if (forward ? z > state.v[eRzPos2] : z < state.v[eRzPos2]) {
      break;
    }
    disc += discStep;
  }

  std::uint32_t holes = static_cast<std::uint32_t>(candidate.hits.size());
  std::uint32_t consecutiveHoles = holes;
  bool cylindersLeft = true;
  bool discsLeft = true;
  // the state at the last accepted measurement is what the track keeps; the
  // last stop may be the escape
  State lastHit = state;
  // what the last stop was: a track in the barrel stays there until a disc
  // comes first, one in the endcap until a cylinder does, so the other
  // kind's stop is looked at only once it can be nearer
  bool inEndcap = false;

  while (true) {
    // how far the track may still go before it has turned the whole budget;
    // a stop beyond that is not reached, whatever the geometry says
    const double kappa = std::abs(m_helix.kappa(state.v));
    const double maxPath =
        kappa > 0. ? (m_cfg.maxTurningAngle - state.turned) / kappa
                   : 2. * (layout.escapeRadius + layout.escapeHalfZ);
    if (maxPath <= 0.) {
      break;
    }

    std::optional<double> sDisc;
    while (discsLeft && discValid()) {
      sDisc = m_helix.pathToDisc(state.v,
                                 layout.surfaces[layout.discs[disc]].refCoord);
      if (sDisc.has_value() && *sDisc > maxPath) {
        // z grows monotonically, so every disc beyond is out of reach too
        sDisc.reset();
        discsLeft = false;
        break;
      }
      if (sDisc.has_value()) {
        break;
      }
      disc += discStep;
    }

    std::optional<double> sCyl;
    if (cylindersLeft && cyl < layout.cylinders.size()) {
      const double rCyl = layout.surfaces[layout.cylinders[cyl]].refCoord;
      bool tryCylinder = true;
      if (inEndcap && sDisc.has_value()) {
        // the radius the track has reached at the disc, with the sagitta
        // over that path as the margin: short of the cylinder means the
        // disc comes first and the cylinder need not be solved
        const double dT = norm2(state.v[eRzDir0], state.v[eRzDir1]);
        const double rAtDisc =
            norm2(state.v[eRzPos0] + state.v[eRzDir0] * *sDisc,
                  state.v[eRzPos1] + state.v[eRzDir1] * *sDisc);
        const double sagitta = 0.5 * kappa * dT * *sDisc * *sDisc;
        tryCylinder = rAtDisc + sagitta + 1. >= rCyl;
      }
      if (tryCylinder) {
        sCyl = m_helix.pathToCylinder(state.v, rCyl);
        if (!sCyl.has_value() || *sCyl > maxPath) {
          // the helix never reaches this radius, so none beyond it either
          sCyl.reset();
          cylindersLeft = false;
        }
      }
    }
    const bool takeCyl =
        sCyl.has_value() && (!sDisc.has_value() || *sCyl <= *sDisc);
    if (!sCyl.has_value() && !sDisc.has_value()) {
      break;
    }
    const double s = takeCyl ? *sCyl : *sDisc;
    const std::uint32_t surfaceIndex =
        takeCyl ? layout.cylinders[cyl] : layout.discs[disc];
    const RzSurface& surface = layout.surfaces[surfaceIndex];
    if (takeCyl) {
      ++cyl;
    } else {
      disc += discStep;
    }

    // land there: a stop off the surface's extent costs no covariance work
    RzVector landed = state.v;
    m_helix.step(landed, s);
    const double along = alongCoordinate(surface, landed);
    if (!surface.contains(along)) {
      // the state itself stays put; the next candidate is measured from here
      continue;
    }
    inEndcap = !takeCyl;
    ++candidate.stops;
    const std::uint32_t stop =
        static_cast<std::uint32_t>(candidate.stopSurfaces.size());
    candidate.stopSurfaces.push_back(surfaceIndex);
    candidate.stopPaths.push_back(s);
    candidate.stopAlong.push_back(along);

    state.v = landed;
    const Vector3 normal = surfaceNormal(surface, state.v);
    state.travel(s);
    state.pending.advance(s);
    state.turned += std::abs(m_helix.kappa(state.v)) * s;
    candidate.pathLength += s;

    const double r = norm2(state.v[eRzPos0], state.v[eRzPos1]);
    if (r > layout.escapeRadius ||
        std::abs(state.v[eRzPos2]) > layout.escapeHalfZ ||
        state.turned > m_cfg.maxTurningAngle) {
      break;
    }

    if (surfaceIndex == startSurface) {
      continue;
    }

    if (m_cfg.applyMaterial) {
      if (const int band = surface.materialBandAt(along);
          band >= 0 && !applyMaterial(state, surface, band, normal)) {
        break;
      }
    }

    if (surface.layer == kRzNone) {
      continue;
    }
    state.moveCovariance(m_helix, normal);
    materialise(state, normal);
    if (searchLayer(grid, surface.layer, stop, state, candidate) > 0) {
      consecutiveHoles = 0;
      lastHit = state;
    } else if (!onModule(surface.layer, state)) {
      // between modules: nothing was expected here
      continue;
    } else {
      candidate.hits.push_back({surface.layer, kRzNone, stop, kRzNone, 0.});
      ++holes;
      ++consecutiveHoles;
      if (holes > m_cfg.maxHoles ||
          consecutiveHoles > m_cfg.maxConsecutiveHoles) {
        break;
      }
    }
  }

  while (!candidate.hits.empty() && candidate.hits.back().isHole()) {
    candidate.hits.pop_back();
  }
  for (const RzTrackHit& hit : candidate.hits) {
    (hit.isHole() ? candidate.holes : candidate.measurements) += 1;
  }
  candidate.parameters = lastHit.v;
  candidate.covariance = lastHit.c;
  if (candidate.measurements < m_cfg.minMeasurements) {
    return false;
  }
  if (m_cfg.backwardPass) {
    backwardPass(grid, lastHit, candidate);
  }
  return true;
}

}  // namespace Acts::Experimental

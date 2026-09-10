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
    : m_cfg(config), m_layout(&layout), m_bz(bz) {}

namespace {
/// A surface's field at a crossing, or the constant one
double bzAt(const RzSurface& surface, double along, double fallback) {
  return surface.bzAt(along).value_or(fallback);
}
}  // namespace

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
  const double maxDistance = std::max(m_cfg.maxModuleDistance, m.maxDistance);
  if (std::abs(s0) > maxDistance) {
    return std::nullopt;
  }
  const RzHelix helix = helixAt(state.bz);
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
      helix.pathToPlane(state.v, m.position, m.normal);
  if (!s.has_value() || std::abs(*s) > maxDistance) {
    return std::nullopt;
  }
  RzVector w = state.v;
  helix.step(w, *s);
  const Vector3 d = m.position - w.segment<3>(eRzPos0);
  const double ru = m.u.dot(d);
  const double rv = m.v.dot(d);
  if (m.dim == 1 && std::abs(rv) > m.halfV + m_cfg.stripMargin) {
    return std::nullopt;
  }
  // In a polar frame the measurement is an angle, so the length it stands for
  // grows with the distance from the frame's origin. The entry carries the
  // variance at the module's own radius; the crossing is `rv` further along
  // the radial direction, and that is where the variance belongs.
  const double lever = 1. - rv * m.invLever;
  const double cov00 = m.cov00 * lever * lever;
  // S = H C H^T + R with H the two frame axes on the position block of the
  // covariance moved to the module: the rows of H J, and only they, are
  // formed, and C (H J)^T is what the update needs too
  Evaluation e;
  const Eigen::Matrix<double, 3, eRzSize> jPos =
      helix.stepJacobianOnto(state.v, *s, w, m.normal).positionRows();
  const Eigen::Matrix<double, 1, eRzSize> hu = m.u.transpose() * jPos;
  e.ch.col(0) = state.c * hu.transpose();
  const double s00 = hu.dot(e.ch.col(0)) + cov00;
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

std::optional<double> RzTrackFinder::pathBackward(const RzHelix& helix,
                                                  const RzVector& v,
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
    helix.step(w, s);
    const double f = w[eRzPos0] * w[eRzPos0] + w[eRzPos1] * w[eRzPos1] -
                     surface.refCoord * surface.refCoord;
    const double df = 2. * (w[eRzPos0] * w[eRzDir0] + w[eRzPos1] * w[eRzDir1]);
    if (df == 0.) {
      break;
    }
    const double ds = f / df;
    s -= ds;
    // a nanometre: every physical scale here is millimetres, and each further
    // iteration is a full helix step
    if (std::abs(ds) < 1e-6) {
      // a root far from the guess is the other crossing of the circle
      if (std::abs(s - guess) < std::max(20., 0.2 * std::abs(guess))) {
        return s;
      }
      break;
    }
  }
  // the closed form the other way round, as a fallback
  const RzVector r = RzHelix::reversed(v);
  const std::optional<double> back = helix.pathToCylinder(r, surface.refCoord);
  if (!back.has_value()) {
    return std::nullopt;
  }
  return -*back;
}

void RzTrackFinder::modulesAt(std::uint32_t layerIndex, const State& state,
                              ModuleList& modules, bool& onModule) const {
  modules.clear();
  onModule = false;
  const RzLayout& layout = *m_layout;
  const RzLayer& layer = layout.layers[layerIndex];
  const RzSurface& surface = layout.surfaces[layer.surface];
  const RzVector& v = state.v;
  const double r = norm2(v[eRzPos0], v[eRzPos1]);
  const double phi = std::atan2(v[eRzPos1], v[eRzPos0]);
  const double along = alongCoordinate(surface, v);

  // Where the state could be, not just where it is: the module test has to
  // open by the same amount the layer search used to, or a module just past
  // the crossing point is never looked at and its measurement is lost.
  const auto cPos = state.c.block<3, 3>(eRzPos0, eRzPos0);
  const double varPending = state.pending.varPosition;
  // the largest variance along any direction is at most the trace
  const double sigmaMax = std::sqrt(std::max(0., cPos.trace() + varPending));
  const double window = m_cfg.windowSigmas * sigmaMax + m_cfg.windowMin;

  const double room = layer.maxHalfExtent + m_cfg.moduleEdgeTolerance + window;
  const Vector3 p = v.segment<3>(eRzPos0);
  const Vector3 dir = v.segment<3>(eRzDir0);
  const double kappa = std::abs(helixAt(state.bz).kappa(v));
  const double maxDistance =
      std::max(m_cfg.maxModuleDistance, layer.moduleDistance);
  const double sigmas2 = m_cfg.windowSigmas * m_cfg.windowSigmas;
  RzMeasurementGrid::visitBins(
      layout, layerIndex, phi, along, room / r, room, [&](std::uint32_t b) {
        for (std::uint32_t i = layout.moduleBinStart[b];
             i < layout.moduleBinStart[b + 1]; ++i) {
          if (modules.size() == modules.capacity()) {
            return;
          }
          const std::uint32_t index = layout.moduleOrder[i];
          const RzModule& m = layout.modules[index];
          const double alongNormal = m.normal.dot(dir);
          if (std::abs(alongNormal) < 1e-9) {
            continue;
          }
          const double s = m.normal.dot(m.center - p) / alongNormal;
          if (std::abs(s) > maxDistance) {
            continue;
          }
          const Vector3 d = p + s * dir - m.center;
          // the sagitta over the distance to the module, on top of the edge
          // tolerance and the state's own spread along each module axis.
          // Most modules the bins hand over are rejected, so the spread is
          // asked for only where the fixed part of the tolerance is already
          // exceeded, and compared squared rather than rooted.
          const double sagitta = 0.5 * kappa * s * s;
          const double fixed =
              m_cfg.moduleEdgeTolerance + sagitta + m_cfg.windowMin;
          const double du = std::abs(m.u.dot(d)) - m.halfU - fixed;
          if (du > 0.) {
            const double varU = std::max(0., m.u.dot(cPos * m.u) + varPending);
            if (du * du > sigmas2 * varU) {
              continue;
            }
          }
          const double dv = std::abs(m.v.dot(d)) - m.halfV - fixed;
          if (dv > 0.) {
            const double varV = std::max(0., m.v.dot(cPos * m.v) + varPending);
            if (dv * dv > sigmas2 * varV) {
              continue;
            }
          }
          if (std::ranges::find(modules, index) == modules.end()) {
            modules.push_back(index);
          }
          // the hole decision, on the crossing itself rather than on where
          // the state might be
          const double tight = m_cfg.moduleEdgeTolerance + sagitta;
          if (std::abs(m.u.dot(d)) <= m.halfU + tight &&
              std::abs(m.v.dot(d)) <= m.halfV + tight) {
            onModule = true;
          }
        }
      });
}

std::uint32_t RzTrackFinder::searchLayer(
    const RzMeasurementGrid& grid, std::uint32_t layerIndex,
    std::uint32_t stop, const ModuleList& modules, State& state,
    RzTrackCandidate& candidate, std::uint32_t skipRounds,
    std::uint32_t usedModule) const {
  // The modules are what the crossing landed on; which of their measurements
  // is worth the full transport is `evaluate`'s gate, which takes the
  // straight-line crossing of the module plane and a diagonal chi2 against
  // the covariance widened over the module distance. A window in the stop's
  // own (phi, along) cannot do that job: the module is offset from the RZ
  // surface, so the track meets it somewhere else entirely.
  std::uint32_t accepted = 0;
  ModuleList usedModules;
  if (usedModule != kRzNone) {
    // a measurement the caller named is already on the track; the layer may
    // still hold the overlap hit, on one of its other modules
    usedModules.push_back(usedModule);
  }
  for (std::uint32_t round = skipRounds;
       round < m_cfg.maxMeasurementsPerLayer; ++round) {
    std::uint32_t bestIndex = kRzNone;
    Evaluation best;
    best.chi2 = std::numeric_limits<double>::max();
    // the best of the cartesian candidates by the straight-line chi2, which
    // is transported exactly once the layer has been walked
    std::uint32_t bestGateIndex = kRzNone;
    double bestGateChi2 = std::numeric_limits<double>::max();
    // Everything the gate needs except the measurement's own position and
    // variance belongs to the module and to the state, not to the
    // measurement: the plane crossing, the residual's origin along the
    // module axes, and the covariance projected on them and widened over
    // the distance to the module. Formed once per module here, the gate
    // costs a subtraction and two multiplies per measurement instead of a
    // division and two quadratic forms. Polar modules carry their own axes
    // per measurement, so they take the general path.
    const Vector3 p0 = state.v.segment<3>(eRzPos0);
    const Vector3 d0 = state.v.segment<3>(eRzDir0);
    const SquareMatrix3 cPos = state.c.block<3, 3>(eRzPos0, eRzPos0);
    const double dirTrace = state.c.block<3, 3>(eRzDir0, eRzDir0).trace();
    const double varPending = state.pending.varPosition;
    const double gate2 = m_cfg.gateFactor * m_cfg.chi2Cut;
    for (const std::uint32_t module : modules) {
      if (std::ranges::find(usedModules, module) != usedModules.end()) {
        continue;
      }
      const RzModule& mod = m_layout->modules[module];
      bool hoisted = false;
      double cu = 0.;
      double cv = 0.;
      double su0 = 0.;
      double sv0 = 0.;
      if (!mod.polar) {
        const double alongNormal = mod.normal.dot(d0);
        if (std::abs(alongNormal) > 1e-9) {
          const double s0 = mod.normal.dot(mod.center - p0) / alongNormal;
          const Vector3 crossing = p0 + s0 * d0;
          cu = mod.u.dot(crossing);
          cv = mod.v.dot(crossing);
          const double spread = dirTrace * s0 * s0 + varPending;
          su0 = mod.u.dot(cPos * mod.u) + spread;
          sv0 = mod.v.dot(cPos * mod.v) + spread;
          hoisted = su0 > 0. && sv0 > 0.;
        }
      }
      for (const std::uint32_t i : grid.moduleRange(module)) {
        const RzMeasurement& m = grid.entry(i);
        ++candidate.candidatesTested;
        if (hoisted) {
          // chi2 = ru^2/su (+ rv^2/sv), tested without the divisions. Only
          // the best of these is worth a full transport: the straight-line
          // chi2 is the same quantity with the covariance widened over the
          // distance to the module, so it orders the candidates, and the one
          // it picks is then evaluated exactly and has to pass the cut.
          const double ru = mod.u.dot(m.position) - cu;
          const double su = su0 + m.cov00;
          if (ru * ru > gate2 * su) {
            continue;
          }
          double chi2Gate = ru * ru / su;
          if (m.dim == 2) {
            const double rv = mod.v.dot(m.position) - cv;
            const double sv = sv0 + m.cov11;
            if (ru * ru * sv + rv * rv * su > gate2 * su * sv) {
              continue;
            }
            chi2Gate += rv * rv / sv;
          }
          if (chi2Gate < bestGateChi2) {
            bestGateChi2 = chi2Gate;
            bestGateIndex = i;
          }
          continue;
        }
        const std::optional<Evaluation> e =
            evaluate(state, m, true);
        if (!e.has_value()) {
          continue;
        }
        if (e->chi2 < best.chi2) {
          bestIndex = i;
          best = *e;
        }
      }
    }
    if (bestGateIndex != kRzNone) {
      if (const std::optional<Evaluation> e =
              evaluate(state, grid.entry(bestGateIndex), false);
          e.has_value() && e->chi2 < best.chi2) {
        bestIndex = bestGateIndex;
        best = *e;
      }
    }
    if (bestIndex == kRzNone || best.chi2 > m_cfg.chi2Cut) {
      break;
    }
    const RzMeasurement& m = grid.entry(bestIndex);
    update(state, best);
    const std::uint32_t forwardState =
        static_cast<std::uint32_t>(candidate.forwardStates.size());
    candidate.forwardStates.emplace_back(state.v, state.c);
    candidate.hits.push_back(
        {layerIndex, bestIndex, stop, forwardState, m.module, best.chi2});
    candidate.chi2 += best.chi2;
    if (usedModules.size() < usedModules.capacity()) {
      usedModules.push_back(m.module);
    }
    ++accepted;
  }
  return accepted;
}

bool RzTrackFinder::takeKnownHit(const RzMeasurementGrid& grid,
                                 std::uint32_t entry, std::uint32_t layerIndex,
                                 std::uint32_t stop, State& state,
                                 RzTrackCandidate& candidate) const {
  const RzMeasurement& m = grid.entry(entry);
  const std::optional<Evaluation> e = evaluate(state, m, false);
  // The caller names the measurement, so no search and no selection - but a
  // measurement the prediction cannot reach at all would be pulled onto the
  // track whatever it does to the fit, so the gate's own threshold still
  // has to hold. Where it does not, the layer is searched as any other.
  if (!e.has_value() || !std::isfinite(e->chi2) ||
      e->chi2 > m_cfg.gateFactor * m_cfg.chi2Cut) {
    return false;
  }
  update(state, *e);
  const std::uint32_t forwardState =
      static_cast<std::uint32_t>(candidate.forwardStates.size());
  candidate.forwardStates.emplace_back(state.v, state.c);
  candidate.hits.push_back(
      {layerIndex, entry, stop, forwardState, m.module, e->chi2});
  candidate.chi2 += e->chi2;
  return true;
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
  state.bz = forward.bz;
  state.anchorBz = forward.bz;
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
        const RzHelix helix = helixAt(state.bz);
        const std::optional<double> s =
            pathBackward(helix, state.v, surface, -candidate.stopPaths[j + 1]);
        if (!s.has_value()) {
          candidate.backwardFailure = 1;
          return;
        }
        helix.step(state.v, *s);
        const Vector3 normal = surfaceNormal(surface, state.v);
        state.travel(*s);
        state.bz = bzAt(surface, candidate.stopAlong[j], m_bz);
        state.pending.advance(-*s);
        // The covariance is needed where there is something to update, and
        // where scattering is waiting to be put in: materialising it at the
        // stop it belongs to, and letting the Jacobians carry it from there,
        // keeps the lever arms exact for the parameters this pass is for.
        if (!state.pending.empty() ||
            (hit != candidate.hits.rend() &&
             hit->stop == static_cast<std::uint32_t>(j))) {
          state.moveCovariance(helixAt(state.anchorBz), normal);
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
      const RzHelix helix = helixAt(state.bz);
      const std::optional<double> s =
          helix.pathToPlane(state.v, m.position, m.normal);
      if (!s.has_value()) {
        candidate.backwardFailure = 1;
        return;
      }
      helix.step(state.v, *s);
      state.travel(*s);
      state.moveCovariance(helixAt(state.anchorBz), m.normal);
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

  if (m_cfg.inwardSearch) {
    // Everything found from here on is inside what the forward pass saw, and
    // is found outward to inward, so it has to be turned around and put in
    // front to keep the hits ordered from the beam line out.
    const std::size_t before = candidate.hits.size();
    const bool reached = inwardSearch(grid, state, candidate);
    if (candidate.hits.size() > before) {
      const auto first =
          candidate.hits.begin() + static_cast<std::ptrdiff_t>(before);
      std::reverse(first, candidate.hits.end());
      std::rotate(candidate.hits.begin(), first, candidate.hits.end());
      // the counts the forward pass took no longer describe the track
      candidate.measurements = 0;
      candidate.holes = 0;
      for (const RzTrackHit& found : candidate.hits) {
        (found.isHole() ? candidate.holes : candidate.measurements) += 1;
      }
    }
    if (!reached) {
      candidate.backwardFailure = 5;
    }
    candidate.innerAtPerigee = reached;
  }

  candidate.innerParameters = state.v;
  candidate.innerCovariance = state.c;
  candidate.hasInner = true;
}

bool RzTrackFinder::inwardSearch(const RzMeasurementGrid& grid, State& state,
                                 RzTrackCandidate& candidate) const {
  const RzLayout& layout = *m_layout;

  // Path lengths inward are negative, the convention the backward replay
  // already uses. A disc is linear in the path; a cylinder comes from the
  // reversed state, which runs the same helix the other way.
  auto pathInwardToDisc = [](const RzVector& v, double z) {
    const double dz = v[eRzDir2];
    return dz != 0. ? std::optional((z - v[eRzPos2]) / dz) : std::nullopt;
  };
  auto pathInwardToCylinder = [](const RzHelix& helix, const RzVector& v,
                                 double radius) -> std::optional<double> {
    const RzVector reversed = RzHelix::reversed(v);
    const std::optional<double> forward =
        helix.pathToCylinder(reversed, radius);
    return forward.has_value() ? std::optional(-*forward) : std::nullopt;
  };

  // cursors, mirrored: cylinders inward from the current radius, discs back
  // toward z = 0 against the direction of travel
  const double r0 = norm2(state.v[eRzPos0], state.v[eRzPos1]);
  std::ptrdiff_t cyl = static_cast<std::ptrdiff_t>(layout.cylinders.size()) - 1;
  while (cyl >= 0 && layout.surfaces[layout.cylinders[cyl]].refCoord >= r0) {
    --cyl;
  }
  const bool travellingForward = state.v[eRzDir2] >= 0.;
  const int discStep = travellingForward ? -1 : 1;
  std::ptrdiff_t disc =
      travellingForward ? static_cast<std::ptrdiff_t>(layout.discs.size()) - 1
                        : 0;
  auto discValid = [&]() {
    return disc >= 0 && disc < static_cast<std::ptrdiff_t>(layout.discs.size());
  };
  while (discValid()) {
    const double z = layout.surfaces[layout.discs[disc]].refCoord;
    if (travellingForward ? z < state.v[eRzPos2] : z > state.v[eRzPos2]) {
      break;
    }
    disc += discStep;
  }

  ModuleList crossedModules;
  // The same three things the outward walk had to stop wasting: the closest
  // approach and the state's own constants are only worth recomputing once
  // the state has moved, a cylinder solve holds while it stands still, and a
  // disc whose radial extent the track cannot reach costs nothing to skip.
  bool stateMoved = true;
  RzHelix helix = helixAt(state.bz);
  double sPerigee = 0.;
  double pzIn = 0.;
  double invDzIn = 0.;
  double pxIn = 0.;
  double pyIn = 0.;
  double dxIn = 0.;
  double dyIn = 0.;
  double halfKappaTIn = 0.;
  std::ptrdiff_t cylCached = -1;
  std::optional<double> cylCachedPath;
  while (true) {
    if (stateMoved) {
      helix = helixAt(state.bz);
      // where the track is closest to the beam axis; nothing inside that is
      // still on the way in
      sPerigee = helix.pathToPerigee(state.v);
      pzIn = state.v[eRzPos2];
      pxIn = state.v[eRzPos0];
      pyIn = state.v[eRzPos1];
      dxIn = state.v[eRzDir0];
      dyIn = state.v[eRzDir1];
      const double dzIn = state.v[eRzDir2];
      invDzIn = dzIn != 0. ? 1. / dzIn : 0.;
      halfKappaTIn =
          0.5 * std::abs(helix.kappa(state.v)) * norm2(dxIn, dyIn);
      cylCached = -1;
      stateMoved = false;
    }
    if (sPerigee >= 0.) {
      break;
    }

    std::optional<double> sCyl;
    std::optional<double> sDisc;
    if (cyl >= 0) {
      if (cylCached == cyl) {
        sCyl = cylCachedPath;
      } else {
        sCyl = pathInwardToCylinder(helix, state.v, layout.cylCoord[cyl]);
        cylCached = cyl;
        cylCachedPath = sCyl;
      }
    }
    if (discValid()) {
      const std::size_t di = static_cast<std::size_t>(disc);
      const double sTry = (layout.discCoord[di] - pzIn) * invDzIn;
      if (sTry < 0. && sTry > sPerigee) {
        const double xs = pxIn + dxIn * sTry;
        const double ys = pyIn + dyIn * sTry;
        const double r2 = xs * xs + ys * ys;
        const double sagitta = halfKappaTIn * sTry * sTry;
        const double lo = layout.discMin[di] - sagitta;
        const double hi = layout.discMax[di] + sagitta;
        if ((lo > 0. && r2 < lo * lo) || r2 > hi * hi) {
          disc += discStep;
          continue;
        }
      }
      sDisc = pathInwardToDisc(state.v, layout.discCoord[di]);
    }
    // inward is negative, so the nearer stop is the larger of the two
    const bool takeCyl =
        sCyl.has_value() && (!sDisc.has_value() || *sCyl >= *sDisc);
    if (!sCyl.has_value() && !sDisc.has_value()) {
      break;
    }
    const double step = takeCyl ? *sCyl : *sDisc;
    if (step > 0. || step <= sPerigee) {
      // behind us, or beyond the closest approach
      break;
    }
    const std::uint32_t surfaceIndex =
        takeCyl ? layout.cylinders[cyl] : layout.discs[disc];
    const RzSurface& surface = layout.surfaces[surfaceIndex];
    if (takeCyl) {
      --cyl;
    } else {
      disc += discStep;
    }

    RzVector landed = state.v;
    helix.step(landed, step);
    const double along = alongCoordinate(surface, landed);
    if (!surface.contains(along)) {
      continue;
    }
    const std::uint32_t stop =
        static_cast<std::uint32_t>(candidate.stopSurfaces.size());
    candidate.stopSurfaces.push_back(surfaceIndex);
    candidate.stopPaths.push_back(step);
    candidate.stopAlong.push_back(along);
    ++candidate.stops;

    state.v = landed;
    stateMoved = true;
    const Vector3 normal = surfaceNormal(surface, state.v);
    state.travel(step);
    state.pending.advance(std::abs(step));
    state.bz = bzAt(surface, along, m_bz);

    // going inward the particle gains back what it lost on the way out
    if (m_cfg.applyMaterial) {
      if (const int band = surface.materialBandAt(along);
          band >= 0 && !applyMaterial(state, surface, band, normal, -1.)) {
        return false;
      }
    }
    if (surface.layer == kRzNone) {
      continue;
    }
    state.moveCovariance(helixAt(state.anchorBz), normal);
    materialise(state, normal);
    // holes are the forward pass's business: this pass is here to pick up
    // what the seed's own layers hold, not to judge what is missing
    bool onModule = false;
    modulesAt(surface.layer, state, crossedModules, onModule);
    if (crossedModules.empty()) {
      continue;
    }
    searchLayer(grid, surface.layer, stop, crossedModules, state, candidate);
  }

  // finish on the beam line: the parameters a caller wants are here, with the
  // material of everything crossed already in the covariance
  const RzHelix endHelix = helixAt(state.bz);
  const double sEnd = endHelix.pathToPerigee(state.v);
  RzVector end = state.v;
  endHelix.step(end, sEnd);
  const double dt = std::hypot(end[eRzDir0], end[eRzDir1]);
  if (dt <= 0.) {
    return false;
  }
  const Vector3 normal(end[eRzDir0] / dt, end[eRzDir1] / dt, 0.);
  state.v = end;
  state.travel(sEnd);
  state.pending.advance(std::abs(sEnd));
  state.moveCovariance(helixAt(state.anchorBz), normal);
  materialise(state, normal);
  return true;
}

bool RzTrackFinder::findTrack(const RzMeasurementGrid& grid,
                              const RzVector& start,
                              const RzMatrix& startCovariance,
                              std::uint32_t startModule,
                              RzTrackCandidate& candidate,
                              std::span<const std::uint32_t> seedEntries)
    const {
  candidate.clear();
  State state;
  state.v = start;
  state.anchor = start;
  state.c = startCovariance;
  state.bz = m_bz;
  state.anchorBz = m_bz;

  const RzLayout& layout = *m_layout;
  // the layer each seed measurement sits on, so a crossing can ask in a few
  // comparisons whether it is one the caller already knows the answer to
  boost::container::static_vector<std::pair<std::uint32_t, std::uint32_t>, 8>
      knownHits;
  for (const std::uint32_t entry : seedEntries) {
    if (entry == kRzNone || knownHits.size() == knownHits.capacity()) {
      continue;
    }
    const std::uint32_t layer = layout.modules[grid.entry(entry).module].layer;
    if (layer != kRzNone) {
      knownHits.emplace_back(layer, entry);
    }
  }
  const auto knownAt = [&](std::uint32_t layer) {
    for (const auto& [l, entry] : knownHits) {
      if (l == layer) {
        return entry;
      }
    }
    return kRzNone;
  };
  std::uint32_t startSurface = kRzNone;
  if (startModule != kRzNone) {
    const std::uint32_t layer = layout.modules[startModule].layer;
    startSurface = layout.layers[layer].surface;
    state.bz =
        bzAt(layout.surfaces[startSurface],
             alongCoordinate(layout.surfaces[startSurface], state.v), m_bz);
    state.anchorBz = state.bz;
    const std::uint32_t known = knownAt(layer);
    const bool took =
        known != kRzNone &&
        takeKnownHit(grid, known, layer, kRzNone, state, candidate);
    ModuleList startModules;
    bool startOnModule = false;
    modulesAt(layer, state, startModules, startOnModule);
    if (!startModules.empty() &&
        searchLayer(grid, layer, kRzNone, startModules, state, candidate,
                    took ? 1u : 0u,
                    took ? grid.entry(known).module : kRzNone) == 0 &&
        !took && startOnModule) {
      candidate.hits.push_back(
          {layer, kRzNone, kRzNone, kRzNone, startModules.front(), 0.});
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


  // the start layer has already run, and what it left is either measurements
  // or one hole; counting its hits as holes would spend the budget on them
  std::uint32_t holes = 0;
  for (const RzTrackHit& hit : candidate.hits) {
    holes += hit.isHole() ? 1 : 0;
  }
  std::uint32_t consecutiveHoles = holes;
  std::uint32_t layersCrossed = 0;
  std::uint32_t measurementsFound =
      static_cast<std::uint32_t>(candidate.hits.size()) - holes;
  ModuleList crossedModules;
  bool cylindersLeft = true;
  bool discsLeft = true;
  // the state at the last accepted measurement is what the track keeps; the
  // last stop may be the escape
  State lastHit = state;
  // what the last stop was: a track in the barrel stays there until a disc
  // comes first, one in the endcap until a cylinder does, so the other
  // kind's stop is looked at only once it can be nearer
  bool inEndcap = false;
  // the cylinder solve, kept while the state has not moved
  std::uint32_t cylCached = kRzNone;
  std::optional<double> cylCachedPath;

  while (true) {
    // how far the track may still go before it has turned the whole budget;
    // a stop beyond that is not reached, whatever the geometry says
    const RzHelix helix = helixAt(state.bz);
    const double kappa = std::abs(helix.kappa(state.v));
    const double maxPath =
        kappa > 0. ? (m_cfg.maxTurningAngle - state.turned) / kappa
                   : 2. * (layout.escapeRadius + layout.escapeHalfZ);
    if (maxPath <= 0.) {
      break;
    }

    // The probe below runs several times per stop and keeps almost none of
    // what it looks at, so everything that does not change while the state
    // stands still is formed here: the reciprocal of dz, the position and
    // direction, and half the transverse curvature for the sagitta.
    const double dTransverse = norm2(state.v[eRzDir0], state.v[eRzDir1]);
    const double pz = state.v[eRzPos2];
    const double dz = state.v[eRzDir2];
    const double invDz = dz != 0. ? 1. / dz : 0.;
    const double px = state.v[eRzPos0];
    const double py = state.v[eRzPos1];
    const double dxDir = state.v[eRzDir0];
    const double dyDir = state.v[eRzDir1];
    const double halfKappaT = 0.5 * kappa * dTransverse;
    std::optional<double> sDisc;
    while (discsLeft && discValid()) {
      const std::size_t di = static_cast<std::size_t>(disc);
      const double sTry = (layout.discCoord[di] - pz) * invDz;
      if (sTry <= 0.) {
        // behind us
        sDisc.reset();
        disc += discStep;
        continue;
      }
      if (sTry > maxPath) {
        // z grows monotonically, so every disc beyond is out of reach too
        sDisc.reset();
        discsLeft = false;
        break;
      }
      // the radius the straight line reaches in the disc's plane, with the
      // sagitta over that path as the margin: a disc whose extent that cannot
      // touch is not crossed, and skipping it here saves the step and the
      // trigonometry the landing would cost. Squared, so the probe needs no
      // root - the landing's own `contains` is what decides either way.
      const double xs = px + dxDir * sTry;
      const double ys = py + dyDir * sTry;
      const double r2 = xs * xs + ys * ys;
      const double sagitta = halfKappaT * sTry * sTry;
      const double lo = layout.discMin[di] - sagitta;
      const double hi = layout.discMax[di] + sagitta;
      if ((lo > 0. && r2 < lo * lo) || r2 > hi * hi) {
        sDisc.reset();
        disc += discStep;
        continue;
      }
      sDisc = sTry;
      break;
    }

    std::optional<double> sCyl;
    if (cylindersLeft && cyl < layout.cylinders.size()) {
      const double rCyl = layout.cylCoord[cyl];
      bool tryCylinder = true;
      if (inEndcap && sDisc.has_value()) {
        // the radius the track has reached at the disc, with the sagitta
        // over that path as the margin: short of the cylinder means the
        // disc comes first and the cylinder need not be solved
        const double rAtDisc = norm2(px + dxDir * *sDisc, py + dyDir * *sDisc);
        const double sagitta = halfKappaT * *sDisc * *sDisc;
        tryCylinder = rAtDisc + sagitta + 1. >= rCyl;
      }
      if (tryCylinder && cylCached == cyl) {
        // an iteration that only rejected a disc left the state where it
        // was, so the Newton solve for this cylinder still holds
        sCyl = cylCachedPath;
      } else if (tryCylinder) {
        sCyl = helix.pathToCylinder(state.v, rCyl);
        if (!sCyl.has_value() || *sCyl > maxPath) {
          // the helix never reaches this radius, so none beyond it either
          sCyl.reset();
          cylindersLeft = false;
        }
        cylCached = cyl;
        cylCachedPath = sCyl;
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
    helix.step(landed, s);
    const double along = alongCoordinate(surface, landed);
    if (!surface.contains(along)) {
      // the state itself stays put; the next candidate is measured from here
      continue;
    }
    inEndcap = !takeCyl;
    ++candidate.stops;
    {
      const bool sens = surface.layer != kRzNone;
      if (takeCyl && sens) {
      } else if (takeCyl) {
      } else if (sens) {
      } else {
      }
    }
    const std::uint32_t stop =
        static_cast<std::uint32_t>(candidate.stopSurfaces.size());
    candidate.stopSurfaces.push_back(surfaceIndex);
    candidate.stopPaths.push_back(s);
    candidate.stopAlong.push_back(along);

    state.v = landed;
    cylCached = kRzNone;
    const Vector3 normal = surfaceNormal(surface, state.v);
    state.travel(s);
    state.pending.advance(s);
    state.turned += std::abs(helix.kappa(state.v)) * s;
    state.bz = bzAt(surface, along, m_bz);
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
    state.moveCovariance(helixAt(state.anchorBz), normal);
    materialise(state, normal);
    ++layersCrossed;
    // a layer the seed names does not have to be searched for that hit
    const std::uint32_t known = knownAt(surface.layer);
    const bool took =
        known != kRzNone &&
        takeKnownHit(grid, known, surface.layer, stop, state, candidate);
    // one geometry pass: the modules the crossing landed on are what the
    // search looks at, and whether there were any is the hole decision
    bool onModule = false;
    modulesAt(surface.layer, state, crossedModules, onModule);
    if (crossedModules.empty()) {
      if (took) {
        // the caller's own measurement is on the track even where the window
        // found no module to search for a second one
        ++measurementsFound;
        consecutiveHoles = 0;
        lastHit = state;
      } else {
        // nothing to look at here
      }
      continue;
    }
    const std::uint32_t accepted =
        searchLayer(grid, surface.layer, stop, crossedModules, state,
                    candidate, took ? 1u : 0u,
                    took ? grid.entry(known).module : kRzNone) +
        (took ? 1u : 0u);
    if (accepted > 0) {
      measurementsFound += accepted;
      consecutiveHoles = 0;
      lastHit = state;
      // Branch stopper. A candidate that has picked up a measurement and is
      // now soft was following noise: the transverse momentum comes out of
      // the filter, so it is only meaningful once a measurement has moved it.
      if (m_cfg.ptMin > 0.) {
        const double qOverP = std::abs(state.v[eRzQOverP]);
        const double sinTheta = norm2(state.v[eRzDir0], state.v[eRzDir1]);
        if (qOverP > 0. && sinTheta / qOverP < m_cfg.ptMin) {
          break;
        }
      }
    } else if (!onModule) {
      // passed between the modules: nothing was expected here
      for (const std::uint32_t module : crossedModules) {
        if (!grid.moduleRange(module).empty()) {
          break;
        }
      }
      continue;
    } else {
      std::size_t onCrossed = 0;
      for (const std::uint32_t module : crossedModules) {
        onCrossed += grid.moduleRange(module).size();
      }
      if (onCrossed == 0) {
      } else {
      }
      candidate.hits.push_back(
          {surface.layer, kRzNone, stop, kRzNone, crossedModules.front(), 0.});
      ++holes;
      ++consecutiveHoles;
      if (holes > m_cfg.maxHoles ||
          consecutiveHoles > m_cfg.maxConsecutiveHoles) {
        break;
      }
    }
    // A candidate that has crossed this many layers and has too little to
    // show for it will not reach `minMeasurements` either
    if (m_cfg.layersForMinMeasurements > 0 &&
        layersCrossed >= m_cfg.layersForMinMeasurements &&
        measurementsFound < m_cfg.minMeasurementsAtLayer) {
      break;
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

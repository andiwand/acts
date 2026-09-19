// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "RzTrackFinderImpl.hpp"

namespace Acts::Experimental::detail::rz {

Placed Finder::placeHit(const RzMeasurementAccessor& measurements,
                        const RzTrackHit& hit) const {
  const RzModuleMeasurements group = measurements(hit.module);
  const RzMeasurementFrame* frame =
      group.frames.empty() ? nullptr : &group.frames[hit.measurement];
  return place(m_layout->modules[hit.module], group.entries[hit.measurement],
               frame);
}

Placed Finder::place(const RzModule& mod, const RzMeasurement& m,
                     const RzMeasurementFrame* frame) const {
  // Polar measurements may provide their own local frame.
  const Vector3& au = frame != nullptr ? frame->u : mod.u;
  const Vector3& av = frame != nullptr ? frame->v : mod.v;
  Placed p;
  p.position = mod.center + m.loc0 * au + m.loc1 * av;
  const bool swapped = m.projector == RzProjector::Loc1 ||
                       (m.projector == RzProjector::Both && m.invLever != 0.);
  // Put a strip's measured coordinate first. For a polar pixel, put the
  // angular coordinate first too, so only its variance gets the lever factor.
  p.u = swapped ? av : au;
  p.v = swapped ? au : av;
  p.normal = p.u.cross(p.v);
  p.cov00 = swapped ? m.cov11 : m.cov00;
  p.cov11 = swapped ? m.cov00 : m.cov11;
  p.cov01 = m.cov01;
  p.invLever = m.invLever;
  p.time = m.time;
  p.timeVariance = m.timeVariance;
  // Bound the unmeasured coordinate; polar frames may span either module axis.
  p.halfV = frame != nullptr ? std::max(mod.halfU, mod.halfV)
                             : (swapped ? mod.halfU : mod.halfV);
  p.maxDistance = m_layout->layers[mod.layer].moduleDistance;
  p.pixel = m.projector == RzProjector::Both;
  return p;
}

template <bool Cache>
std::optional<Evaluation> Finder::evaluate(
    const State& state, const Placed& m, bool gate, bool useTime,
    std::optional<Prediction>* prediction) const {
  // Use the straight-line plane crossing for the inexpensive pre-gate.
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
  const RzHelix helix = RzHelix{state.bz};
  // Reject distant hits using a straight-line residual and covariance widened
  // by direction uncertainty before computing the exact transport.
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
    if (m.pixel) {
      const double sv = m.v.dot(cPos * m.v) + m.cov11 + spread;
      chi2 += rv * rv / sv;
    }
    if (chi2 > m_cfg.gateFactor * m_cfg.chi2Cut) {
      return std::nullopt;
    }
  }
  // Reuse the converged crossing and trigonometry for the Jacobian.
  std::optional<RzHelix::PlaneStep> crossing;
  if constexpr (Cache) {
    // Rotating measurement frames can change the plane within a module.
    if (*prediction && (*prediction)->normal == m.normal &&
        std::abs(m.normal.dot(m.position - (*prediction)->planePosition)) <
            1e-12) {
      crossing = (*prediction)->crossing;
    }
  }
  if (!crossing) {
    crossing = helix.stepToPlane(state.v, m.position, m.normal);
    if constexpr (Cache) {
      if (crossing) {
        prediction->emplace(
            Prediction{m.position, m.normal, *crossing, {}, false});
      }
    }
  }
  if (!crossing.has_value() || std::abs(crossing->s) > maxDistance) {
    return std::nullopt;
  }
  const detail::StepTrig& trig = crossing->trig;
  const RzVector& w = crossing->state;
  double timeChi2 = 0.;
  double timeResidual = 0.;
  double timeGain = 0.;
  if (useTime && m.timeVariance > 0.) {
    const double mOverP = state.massOverCharge * state.v[eRzQOverP];
    const double predictedTime =
        state.time + crossing->s * std::sqrt(1. + mOverP * mOverP);
    timeResidual = m.time - predictedTime;
    const double variance = state.timeVariance + m.timeVariance;
    timeChi2 = timeResidual * timeResidual / variance;
    timeGain = state.timeVariance / variance;
  }
  const Vector3 d = m.position - w.segment<3>(eRzPos0);
  const double ru = m.u.dot(d);
  const double rv = m.v.dot(d);
  if (!m.pixel && std::abs(rv) > m.halfV + m_cfg.stripMargin) {
    return std::nullopt;
  }
  // Rescale angular variance from the measured radius to the crossing radius.
  const double lever = 1. - rv * m.invLever;
  const double cov00 = m.cov00 * lever * lever;
  // Form only the measured rows of H J; reuse C (H J)^T in the update.
  Evaluation e;
  Eigen::Matrix<double, 3, eRzSize> jPos;
  if constexpr (Cache) {
    if ((*prediction)->hasJacobian) {
      jPos = (*prediction)->jPos;
    } else {
      jPos = helix.stepJacobianOnto(state.v, crossing->s, w, m.normal, trig)
                 .positionRows();
      (*prediction)->jPos = jPos;
      (*prediction)->hasJacobian = true;
    }
  } else {
    jPos = helix.stepJacobianOnto(state.v, crossing->s, w, m.normal, trig)
               .positionRows();
  }
  // the two products by hand: Eigen takes a 1x3 by 3x7 and a 7x7 by 7x1
  // through its general kernels, out of line, for 21 and 49 multiplies
  const auto projectRows = [&](const Vector3& axis,
                               Eigen::Matrix<double, 1, eRzSize>& h) {
    for (std::uint32_t c = 0; c < eRzSize; ++c) {
      h[c] =
          axis.x() * jPos(0, c) + axis.y() * jPos(1, c) + axis.z() * jPos(2, c);
    }
  };
  const auto covarianceTimes = [&](const Eigen::Matrix<double, 1, eRzSize>& h,
                                   std::uint32_t col) {
    for (std::uint32_t r = 0; r < eRzSize; ++r) {
      double acc = 0.;
      for (std::uint32_t c = 0; c < eRzSize; ++c) {
        acc += state.c(r, c) * h[c];
      }
      e.ch(r, col) = acc;
    }
  };
  Eigen::Matrix<double, 1, eRzSize> hu;
  projectRows(m.u, hu);
  covarianceTimes(hu, 0);
  const double s00 = hu.dot(e.ch.col(0)) + cov00;
  if (!(s00 > 0.)) {
    return std::nullopt;
  }
  e.sInv.setZero();
  if (!m.pixel) {
    e.ch.col(1).setZero();
    e.sInv(0, 0) = 1. / s00;
    e.chi2 = ru * ru * e.sInv(0, 0);
  } else {
    Eigen::Matrix<double, 1, eRzSize> hv;
    projectRows(m.v, hv);
    covarianceTimes(hv, 1);
    const double s01 = hv.dot(e.ch.col(0)) + m.cov01 * lever;
    const double s11 = hv.dot(e.ch.col(1)) + m.cov11;
    const double det = s00 * s11 - s01 * s01;
    if (!(det > 0.)) {
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
  e.chi2 += timeChi2;
  e.pixel = m.pixel;
  e.timeResidual = timeResidual;
  e.timeGain = timeGain;
  e.hasTime = useTime && m.timeVariance > 0.;
  return e;
}

void Finder::update(State& state, const Evaluation& e) const {
  if (e.hasTime) {
    const double correction = e.timeGain * e.timeResidual;
    state.time += correction;
    state.timeCorrection += correction;
    state.timeVariance *= 1. - e.timeGain;
  }
  if (!e.pixel) {
    const RzVector k = e.ch.col(0) * e.sInv(0, 0);
    state.v += k * e.residual[0];
    state.v.segment<3>(eRzDir0).normalize();
    state.anchor = state.v;
    for (std::uint32_t r = 0; r < eRzSize; ++r) {
      for (std::uint32_t c = 0; c <= r; ++c) {
        state.c(r, c) -= k[r] * e.ch(c, 0);
        state.c(c, r) = state.c(r, c);
      }
    }
    return;
  }
  const Eigen::Matrix<double, eRzSize, 2> k = e.ch * e.sInv;
  state.v += k * e.residual;
  state.v.segment<3>(eRzDir0).normalize();
  state.anchor = state.v;
  // Compute the symmetric covariance update once per lower-triangle entry.
  for (std::uint32_t r = 0; r < eRzSize; ++r) {
    for (std::uint32_t c = 0; c <= r; ++c) {
      const double d = k(r, 0) * e.ch(c, 0) + k(r, 1) * e.ch(c, 1);
      state.c(r, c) -= d;
      state.c(c, r) = state.c(r, c);
    }
  }
}

void Finder::modulesAt(std::uint32_t layerIndex, const State& state,
                       ModuleList& modules, bool& onModule,
                       RzTrackCandidate& candidate) const {
  modules.clear();
  onModule = false;
  const RzLayout& layout = *m_layout;
  const RzLayer& layer = layout.layers[layerIndex];
  const RzSurface& surface = layout.surfaces[layer.surface];
  const RzVector& v = state.v;
  const double r = fastHypot(v[eRzPos0], v[eRzPos1]);
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
  const double kappa = std::abs(RzHelix{state.bz}.kappa(v));
  const double maxDistance =
      std::max(m_cfg.maxModuleDistance, layer.moduleDistance);
  const double sigmas2 = m_cfg.windowSigmas * m_cfg.windowSigmas;
  layout.visitBins(
      layerIndex, phi, along, room / r, room, [&](std::uint32_t b) {
        ++candidate.binsVisited;
        for (std::uint32_t i = layout.moduleBinStart[b];
             i < layout.moduleBinStart[b + 1]; ++i) {
          ++candidate.modulesTested;
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
          // Add sagitta and edge tolerance first; compute projected variance
          // only for modules outside that margin, comparing squared distances.
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

std::uint32_t Finder::searchLayer(const RzMeasurementAccessor& measurements,
                                  std::uint32_t layerIndex, std::uint32_t stop,
                                  const ModuleList& modules, State& state,
                                  RzTrackCandidate& candidate,
                                  std::uint32_t skipRounds,
                                  std::uint32_t usedModule) const {
  // Gate at each crossed module before the exact transport: an RZ-stop window
  // would miss hits on modules offset from that stop.
  std::uint32_t accepted = 0;
  ModuleList usedModules;
  if (usedModule != kRzNone) {
    // a measurement the caller named is already on the track; the layer may
    // still hold the overlap hit, on one of its other modules
    usedModules.push_back(usedModule);
  }
  for (std::uint32_t round = skipRounds; round < m_cfg.maxMeasurementsPerLayer;
       ++round) {
    std::uint32_t bestIndex = kRzNone;
    std::uint32_t bestModule = kRzNone;
    Evaluation best;
    best.chi2 = std::numeric_limits<double>::max();
    // Reuse the straight-line crossing and covariance projection across a
    // module's hits; rotate them for each polar measurement frame.
    const Vector3 p0 = state.v.segment<3>(eRzPos0);
    const Vector3 d0 = state.v.segment<3>(eRzDir0);
    const SquareMatrix3 cPos = state.c.block<3, 3>(eRzPos0, eRzPos0);
    const double dirTrace = state.c.block<3, 3>(eRzDir0, eRzDir0).trace();
    const double varPending = state.pending.varPosition;
    const double gate2 = m_cfg.gateFactor * m_cfg.chi2Cut;
    const double mOverP = state.massOverCharge * state.v[eRzQOverP];
    const double timePerPath = std::sqrt(1. + mOverP * mOverP);
    // Evaluate survivors exactly in their original module and hit order.
    struct Survivor {
      std::uint32_t module;
      std::uint32_t index;
      bool gated;
    };
    boost::container::small_vector<Survivor, 16> survivors;
    for (const std::uint32_t module : modules) {
      if (std::ranges::find(usedModules, module) != usedModules.end()) {
        continue;
      }
      const RzModule& mod = m_layout->modules[module];
      const RzModuleMeasurements group = measurements(module);
      if (group.entries.empty()) {
        continue;
      }
      // Project once per module; rotate these quantities for polar hits.
      bool crossed = false;
      double c0 = 0.;
      double c1 = 0.;
      double su0 = 0.;
      double sv0 = 0.;
      double cUU = 0.;
      double cUV = 0.;
      double cVV = 0.;
      double spread = 0.;
      double moduleStep = 0.;
      {
        const double alongNormal = mod.normal.dot(d0);
        if (std::abs(alongNormal) > 1e-9) {
          const double s0 = mod.normal.dot(mod.center - p0) / alongNormal;
          moduleStep = s0;
          // the crossing on the module's own axes, which is the frame the
          // measurements are already in, so a residual is a subtraction
          const Vector3 delta = p0 + s0 * d0 - mod.center;
          c0 = mod.u.dot(delta);
          c1 = mod.v.dot(delta);
          spread = dirTrace * s0 * s0 + varPending;
          cUU = mod.u.dot(cPos * mod.u);
          cVV = mod.v.dot(cPos * mod.v);
          su0 = cUU + spread;
          sv0 = cVV + spread;
          if (mod.polar) {
            cUV = mod.u.dot(cPos * mod.v);
          }
          crossed = true;
        }
      }
      const bool hoisted = crossed && !mod.polar && su0 > 0. && sv0 > 0.;
      const bool polarGate = crossed && mod.polar && !group.frames.empty();
      for (std::uint32_t i = 0; i < group.entries.size(); ++i) {
        const RzMeasurement& m = group.entries[i];
        ++candidate.candidatesTested;
        if (crossed && m.timeVariance > 0.) {
          const double predictedTime = state.time + moduleStep * timePerPath;
          const double residual = m.time - predictedTime;
          if (residual * residual >
              gate2 * (state.timeVariance + m.timeVariance)) {
            continue;
          }
        }
        bool gated = hoisted;
        if (polarGate) {
          // the module's crossing and spread turned onto this measurement's
          // axes: the same gate `evaluate` would take, before the placement
          // and the exact transport it would then go on to
          const RzMeasurementFrame& f = group.frames[i];
          const double r0 = f.uU * c0 + f.uV * c1;
          const double r1 = f.vU * c0 + f.vV * c1;
          const double s0v = f.uU * f.uU * cUU + 2. * f.uU * f.uV * cUV +
                             f.uV * f.uV * cVV + spread;
          const double s1v = f.vU * f.vU * cUU + 2. * f.vU * f.vV * cUV +
                             f.vV * f.vV * cVV + spread;
          if (s0v > 0. && s1v > 0.) {
            double chi2Gate = 0.;
            if (m.projector != RzProjector::Loc1) {
              const double d = m.loc0 - r0;
              chi2Gate += d * d / (s0v + m.cov00);
            }
            if (m.projector != RzProjector::Loc0) {
              const double d = m.loc1 - r1;
              chi2Gate += d * d / (s1v + m.cov11);
            }
            if (chi2Gate > gate2) {
              continue;
            }
            gated = true;
            survivors.push_back({module, i, true});
            continue;
          }
        }
        if (hoisted) {
          // The diagonal straight-line chi2 rejects distant measurements;
          // it does not preserve the exact chi2 ordering or strip acceptance.
          double chi2Gate = 0.;
          if (m.projector != RzProjector::Loc1) {
            const double r0 = m.loc0 - c0;
            const double s0v = su0 + m.cov00;
            if (r0 * r0 > gate2 * s0v) {
              continue;
            }
            chi2Gate = r0 * r0 / s0v;
          }
          if (m.projector != RzProjector::Loc0) {
            const double r1 = m.loc1 - c1;
            const double s1v = sv0 + m.cov11;
            if (r1 * r1 > gate2 * s1v) {
              continue;
            }
            chi2Gate += r1 * r1 / s1v;
            if (chi2Gate > gate2) {
              continue;
            }
          }
          survivors.push_back({module, i, true});
          continue;
        }
        // no gate of its own here: `evaluate` gates it, and it is always
        // looked at
        survivors.push_back({module, i, gated});
      }
    }
    std::optional<Prediction> prediction;
    std::uint32_t predictionModule = kRzNone;
    for (std::size_t i = 0; i < survivors.size(); ++i) {
      const Survivor& c = survivors[i];
      const RzModule& mod = m_layout->modules[c.module];
      const RzModuleMeasurements group = measurements(c.module);
      const RzMeasurementFrame* frame =
          group.frames.empty() ? nullptr : &group.frames[c.index];
      if (predictionModule != c.module) {
        prediction.reset();
        predictionModule = c.module;
      }
      const Placed placed = place(mod, group.entries[c.index], frame);
      const bool sharedPlane =
          (i + 1 < survivors.size() && survivors[i + 1].module == c.module) ||
          prediction.has_value();
      const std::optional<Evaluation> e =
          sharedPlane
              ? evaluate<true>(state, placed, !c.gated, true, &prediction)
              : evaluate(state, placed, !c.gated);
      if (!e.has_value()) {
        continue;
      }
      ++candidate.exactEvaluated;
      if (e->chi2 < best.chi2) {
        bestIndex = c.index;
        bestModule = c.module;
        best = *e;
      }
    }
    if (bestIndex == kRzNone || best.chi2 > m_cfg.chi2Cut) {
      break;
    }
    update(state, best);
    const std::uint32_t forwardState = saveForwardState(state, candidate);
    candidate.hits.push_back(
        {layerIndex, bestIndex, stop, forwardState, bestModule, best.chi2});
    candidate.chi2 += best.chi2;
    usedModules.push_back(bestModule);
    ++accepted;
  }
  return accepted;
}

std::uint32_t Finder::saveForwardState(const State& state,
                                       RzTrackCandidate& candidate) const {
  ++candidate.measurements;
  if (!m_cfg.storeForwardStates &&
      (!m_cfg.backwardPass || candidate.measurements != m_cfg.backwardLayers)) {
    return kRzNone;
  }
  const auto index = static_cast<std::uint32_t>(candidate.forwardStates.size());
  candidate.forwardStates.push_back(
      {state.v, state.c, state.time, state.timeVariance, state.timeCorrection});
  return index;
}

bool Finder::takeKnownHit(const RzMeasurementAccessor& measurements,
                          const RzSeedMeasurement& seed,
                          std::uint32_t layerIndex, std::uint32_t stop,
                          State& state, RzTrackCandidate& candidate) const {
  const RzModuleMeasurements group = measurements(seed.module);
  if (seed.index >= group.entries.size()) {
    return false;
  }
  const RzModule& mod = m_layout->modules[seed.module];
  const RzMeasurementFrame* frame =
      group.frames.empty() ? nullptr : &group.frames[seed.index];
  const std::optional<Evaluation> e =
      evaluate(state, place(mod, group.entries[seed.index], frame), false);
  // Known hits must still pass the gate. Otherwise search the layer normally.
  if (!e.has_value() || !std::isfinite(e->chi2) ||
      e->chi2 > m_cfg.gateFactor * m_cfg.chi2Cut) {
    return false;
  }
  ++candidate.exactEvaluated;
  update(state, *e);
  const std::uint32_t forwardState = saveForwardState(state, candidate);
  candidate.hits.push_back(
      {layerIndex, seed.index, stop, forwardState, seed.module, e->chi2});
  candidate.chi2 += e->chi2;
  return true;
}

// The backward refit evaluates known hits without the module-plane cache.
template std::optional<Evaluation> Finder::evaluate<false>(
    const State&, const Placed&, bool, bool, std::optional<Prediction>*) const;

}  // namespace Acts::Experimental::detail::rz

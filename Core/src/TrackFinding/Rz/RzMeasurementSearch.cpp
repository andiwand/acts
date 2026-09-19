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
  const RzModule& module = m_layout->modules[hit.module];
  return placeMeasurement(module, group.entries[hit.measurement], frame,
                          m_layout->layers[module.layer].moduleDistance);
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
      const Placed placed =
          placeMeasurement(mod, group.entries[c.index], frame,
                           m_layout->layers[mod.layer].moduleDistance);
      const bool sharedPlane =
          (i + 1 < survivors.size() && survivors[i + 1].module == c.module) ||
          prediction.has_value();
      const std::optional<Evaluation> e =
          sharedPlane ? m_evaluator.evaluate<true>(state, placed, !c.gated,
                                                   true, &prediction)
                      : m_evaluator.evaluate(state, placed, !c.gated);
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
    kalmanUpdate(state, best);
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
  const std::optional<Evaluation> e = m_evaluator.evaluate(
      state,
      placeMeasurement(mod, group.entries[seed.index], frame,
                       m_layout->layers[mod.layer].moduleDistance),
      false);
  // Known hits must still pass the gate. Otherwise search the layer normally.
  if (!e.has_value() || !std::isfinite(e->chi2) ||
      e->chi2 > m_cfg.gateFactor * m_cfg.chi2Cut) {
    return false;
  }
  ++candidate.exactEvaluated;
  kalmanUpdate(state, *e);
  const std::uint32_t forwardState = saveForwardState(state, candidate);
  candidate.hits.push_back(
      {layerIndex, seed.index, stop, forwardState, seed.module, e->chi2});
  candidate.chi2 += e->chi2;
  return true;
}

}  // namespace Acts::Experimental::detail::rz

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

#include "RzTrackFinderImpl.hpp"

namespace Acts::Experimental::detail::rz {

namespace {

std::optional<double> pathBackward(const RzHelix& helix, const RzVector& v,
                                   const RzSurface& surface, double guess) {
  if (surface.shape == RzShape::Disc) {
    const double dz = v[eRzDir2];
    return dz != 0. ? std::optional((surface.refCoord - v[eRzPos2]) / dz)
                    : std::nullopt;
  }
  // Newton from the forward path: the state has moved by an update or two
  // since, so the root is next to the guess
  double s = guess;
  for (std::int32_t i = 0; i < 6; ++i) {
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

}  // namespace

void Finder::backwardPass(const RzMeasurementAccessor& measurements,
                          const State& forward,
                          RzTrackCandidate& candidate) const {
  // Start at the last measurement or the configured inner checkpoint.
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
      partial ? candidate.forwardStates[hit->forwardState].parameters
              : forward.v;
  const RzMatrix& startC =
      partial ? candidate.forwardStates[hit->forwardState].covariance
              : forward.c;

  // Refit from an inflated, uncorrelated prior in the surface and direction
  // tangent planes. Preserve their null directions.
  State state;
  state.v = startV;
  state.anchor = startV;
  state.time = forward.time;
  state.timeVariance = forward.timeVariance;
  state.massOverCharge = forward.massOverCharge;
  if (partial) {
    const auto& checkpoint = candidate.forwardStates[hit->forwardState];
    state.time =
        checkpoint.time + forward.timeCorrection - checkpoint.timeCorrection;
  }
  const RzSurface& startSurface =
      m_layout->surfaces[m_layout->layers[hit->layer].surface];
  const double startAlong = hit->stop == kRzNone
                                ? alongCoordinate(startSurface, startV)
                                : candidate.stopAlong[hit->stop];
  state.bz = bzAt(startSurface, startAlong, m_bz);
  state.br = startSurface.brAt(startAlong);
  state.anchorBz = state.bz;
  {
    const Vector3 d = startV.segment<3>(eRzDir0);
    const Vector3 n = placeHit(measurements, *hit).normal;
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
      // Undo energy loss from the final forward state back to the checkpoint.
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
          const std::int32_t band =
              surface.materialBandAt(candidate.stopAlong[j]);
          if (band < 0) {
            continue;
          }
          m_stepper.regainEnergy(state, surface, band,
                                 surfaceNormal(surface, state.v));
        }
      }
    } else {
      state.c(eRzQOverP, eRzQOverP) =
          startC(eRzQOverP, eRzQOverP) * m_cfg.backwardInflation;
    }
  }

  // Replay hits, material and transport inward in reverse stop order.
  // Keep forward timing information without applying the same hits twice.
  auto updateHitsAt = [&](std::uint32_t stop) {
    // hits are stored outward, so this stop's hits are the next ones inward
    while (hit != candidate.hits.rend() && hit->stop == stop) {
      if (!hit->isHole()) {
        const std::optional<Evaluation> e =
            evaluate(state, placeHit(measurements, *hit), false, false);
        if (!e.has_value()) {
          return false;
        }
        ++candidate.exactEvaluated;
        update(state, *e);
      }
      ++hit;
    }
    return true;
  };

  std::uint32_t stop = hit->stop;
  if (stop != kRzNone) {
    for (std::int32_t j = static_cast<std::int32_t>(stop); j >= 0; --j) {
      const RzSurface& surface = m_layout->surfaces[candidate.stopSurfaces[j]];
      if (static_cast<std::uint32_t>(j) != stop) {
        const RzHelix helix = RzHelix{state.bz};
        const std::optional<double> s =
            pathBackward(helix, state.v, surface, -candidate.stopPaths[j + 1]);
        if (!s.has_value()) {
          candidate.backwardFailure = 1;
          return;
        }
        RzVector landed = state.v;
        helix.step(landed, *s);
        const Vector3 normal =
            m_stepper.land(state, landed, *s, surface, candidate.stopAlong[j]);
        // Materialise pending noise at its stop to preserve scattering lever
        // arms.
        if (!state.pending.empty() ||
            (hit != candidate.hits.rend() &&
             hit->stop == static_cast<std::uint32_t>(j))) {
          state.moveCovariance(RzHelix{state.anchorBz}, normal);
          m_stepper.materialise(state, normal);
        }
      }
      if (!updateHitsAt(static_cast<std::uint32_t>(j))) {
        candidate.backwardFailure = 3;
        return;
      }
      if (m_cfg.applyMaterial) {
        const Vector3 normal = surfaceNormal(surface, state.v);
        if (const std::int32_t band =
                surface.materialBandAt(alongCoordinate(surface, state.v));
            band >= 0 &&
            !m_stepper.applyMaterial(state, surface, band, normal, -1.)) {
          candidate.backwardFailure = 2;
          return;
        }
      }
    }
  }
  // what the track started with, found before any stop
  while (hit != candidate.hits.rend()) {
    if (!hit->isHole()) {
      const Placed m = placeHit(measurements, *hit);
      const RzHelix helix = RzHelix{state.bz};
      const std::optional<double> s =
          helix.pathToPlane(state.v, m.position, m.normal);
      if (!s.has_value()) {
        candidate.backwardFailure = 1;
        return;
      }
      helix.step(state.v, *s);
      state.travel(*s);
      state.moveCovariance(RzHelix{state.anchorBz}, m.normal);
      state.pending.advance(*s);
      m_stepper.materialise(state, m.normal);
      const std::optional<Evaluation> e = evaluate(state, m, false, false);
      if (!e.has_value()) {
        candidate.backwardFailure = 3;
        return;
      }
      ++candidate.exactEvaluated;
      update(state, *e);
    }
    ++hit;
  }
  if (partial && m_cfg.backwardQOverPScale == 0.) {
    state.c(eRzQOverP, eRzQOverP) = forward.c(eRzQOverP, eRzQOverP);
  }

  if (m_cfg.inwardSearch) {
    // Reverse new inward hits and prepend them to preserve outward hit order.
    const std::size_t before = candidate.hits.size();
    const bool reached = inwardSearch(measurements, state, candidate);
    if (candidate.hits.size() > before) {
      const auto first =
          candidate.hits.begin() + static_cast<std::int32_t>(before);
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
  }

  candidate.innerTime = state.time;
  candidate.innerTimeVariance = state.timeVariance;
  candidate.innerParameters = state.v;
  candidate.innerCovariance = state.c;
  candidate.hasInner = true;
}

bool Finder::inwardSearch(const RzMeasurementAccessor& measurements,
                          State& state, RzTrackCandidate& candidate) const {
  const RzLayout& layout = *m_layout;

  // Inward paths are negative; solve cylinders with the reversed helix.
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
  const double r0 = fastHypot(state.v[eRzPos0], state.v[eRzPos1]);
  std::int32_t cyl = static_cast<std::int32_t>(layout.cylinders.size()) - 1;
  while (cyl >= 0 && layout.surfaces[layout.cylinders[cyl]].refCoord >= r0) {
    --cyl;
  }
  const bool travellingForward = state.v[eRzDir2] >= 0.;
  const std::int32_t discStep = travellingForward ? -1 : 1;
  // the first disc behind, by binary search on the sorted disc positions
  std::int32_t disc =
      travellingForward
          ? std::ranges::lower_bound(layout.discCoord, state.v[eRzPos2]) -
                layout.discCoord.begin() - 1
          : std::ranges::upper_bound(layout.discCoord, state.v[eRzPos2]) -
                layout.discCoord.begin();
  auto discValid = [&]() {
    return disc >= 0 && disc < static_cast<std::int32_t>(layout.discs.size());
  };

  ModuleList crossedModules;
  // Cache navigation quantities until the state changes.
  bool stateMoved = true;
  RzHelix helix = RzHelix{state.bz};
  double sPerigee = 0.;
  double pzIn = 0.;
  double invDzIn = 0.;
  double pxIn = 0.;
  double pyIn = 0.;
  double dxIn = 0.;
  double dyIn = 0.;
  double halfKappaTIn = 0.;
  std::int32_t cylCached = -1;
  std::optional<double> cylCachedPath;
  while (true) {
    if (stateMoved) {
      helix = RzHelix{state.bz};
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
          0.5 * std::abs(helix.kappa(state.v)) * fastHypot(dxIn, dyIn);
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

    const Vector3 normal = m_stepper.land(state, landed, step, surface, along);
    stateMoved = true;

    // going inward the particle gains back what it lost on the way out
    if (m_cfg.applyMaterial) {
      if (const std::int32_t band = surface.materialBandAt(along);
          band >= 0 &&
          !m_stepper.applyMaterial(state, surface, band, normal, -1.)) {
        return false;
      }
    }
    if (surface.layer == kRzNone) {
      continue;
    }
    // Skip previously searched layers. Refit rounding can place the state just
    // past its starting surface and otherwise select the same hit twice.
    if (std::ranges::any_of(candidate.hits, [&](const RzTrackHit& hit) {
          return hit.layer == surface.layer;
        })) {
      continue;
    }
    state.moveCovariance(RzHelix{state.anchorBz}, normal);
    m_stepper.materialise(state, normal);
    // Inward extension adds measurements without counting new holes.
    bool onModule = false;
    modulesAt(surface.layer, state, crossedModules, onModule, candidate);
    if (crossedModules.empty()) {
      continue;
    }
    searchLayer(measurements, surface.layer, stop, crossedModules, state,
                candidate);
  }

  // Finish at closest approach with the accumulated material covariance.
  const RzHelix endHelix = RzHelix{state.bz};
  const double sEnd = endHelix.pathToPerigee(state.v);
  RzVector end = state.v;
  endHelix.step(end, sEnd);
  const double dt = fastHypot(end[eRzDir0], end[eRzDir1]);
  if (dt <= 0.) {
    return false;
  }
  const Vector3 normal(end[eRzDir0] / dt, end[eRzDir1] / dt, 0.);
  state.v = end;
  state.travel(sEnd);
  state.pending.advance(sEnd);
  state.moveCovariance(RzHelix{state.anchorBz}, normal);
  m_stepper.materialise(state, normal);
  return true;
}

}  // namespace Acts::Experimental::detail::rz

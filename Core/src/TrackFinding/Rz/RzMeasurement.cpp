// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RzMeasurement.hpp"

#include <algorithm>
#include <cmath>

namespace Acts::Experimental::detail::rz {

Placed placeMeasurement(const RzModule& mod, const RzMeasurement& m,
                        const RzMeasurementFrame* frame, double maxDistance) {
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
  p.maxDistance = maxDistance;
  p.pixel = m.projector == RzProjector::Both;
  return p;
}

template <bool Cache>
std::optional<Evaluation> MeasurementEvaluator::evaluate(
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
  const double maxDistance = std::max(m_maxModuleDistance, m.maxDistance);
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
    if (chi2 > m_gateChi2) {
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
  if (!m.pixel && std::abs(rv) > m.halfV + m_stripMargin) {
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
  // Explicit products avoid Eigen's out-of-line kernels for these small sizes.
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

void kalmanUpdate(State& state, const Evaluation& e) {
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

template std::optional<Evaluation> MeasurementEvaluator::evaluate<false>(
    const State&, const Placed&, bool, bool, std::optional<Prediction>*) const;
template std::optional<Evaluation> MeasurementEvaluator::evaluate<true>(
    const State&, const Placed&, bool, bool, std::optional<Prediction>*) const;

}  // namespace Acts::Experimental::detail::rz

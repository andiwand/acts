// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <utility>

namespace Acts {

/// Prototype: how the CKF uses the non-Gaussian tail of multiple scattering.
enum class ScatteringTailMode {
  /// Gaussian core only
  None,
  /// Retry the track state creation with the tail covariance on a miss
  Retry,
  /// As `Retry`, but also keep the miss (hole or outlier) as a branch
  RetryBranch,
  /// As `Retry`, but update with the core covariance
  SelectOnly,
  /// Select with the mixture chi2 and update with the posterior covariance
  Mixture,
};

namespace detail {

/// Two-Gaussian approximation of the Fatras general mixture (semi-Gaussian
/// regime, Z = 14, beta = 1) for one projected scattering angle.
///
/// Returns the tail weight and the tail width in units of the Highland width.
/// Fitted by EM on 2M samples for each thickness, interpolated in ln(x/X0).
inline std::pair<double, double> scatteringTailModel(double thicknessInX0) {
  static constexpr std::array<double, 6> lnT = {-5.8091, -5.2983, -4.6052,
                                                -3.9120, -3.2189, -2.5257};
  static constexpr std::array<double, 6> eps = {0.0152, 0.0151, 0.0156,
                                                0.0171, 0.0199, 0.0276};
  static constexpr std::array<double, 6> width = {7.263, 6.945, 6.207,
                                                  5.373, 4.554, 3.628};
  const double x = std::clamp(std::log(thicknessInX0), lnT.front(), lnT.back());
  std::size_t i = 0;
  while (i + 2 < lnT.size() && x > lnT[i + 1]) {
    ++i;
  }
  const double f = (x - lnT[i]) / (lnT[i + 1] - lnT[i]);
  return {eps[i] + f * (eps[i + 1] - eps[i]),
          width[i] + f * (width[i + 1] - width[i])};
}

/// Prototype: number of track states accepted through the tail, all threads.
inline std::atomic<std::size_t>& scatteringTailCounter() {
  static std::atomic<std::size_t> counter{0};
  return counter;
}

/// The tail hypothesis since the last accepted measurement.
///
/// Each scatterer adds the excess `(k^2 - 1) Q` of its tail over the core.
/// Only one scatterer is assumed to be in the tail, so the tail covariance is
/// the tail-weighted mean of the excesses.
struct ScatteringTailState {
  /// sum of eps_i * D_i, in the frame of the last transport
  BoundMatrix weightedExcess = BoundMatrix::Zero();
  /// sum of eps_i
  double weightSum = 0;
  /// sum of ln(1 - eps_i)
  double logNoTail = 0;

  void reset() {
    weightedExcess.setZero();
    weightSum = 0;
    logNoTail = 0;
  }

  bool empty() const { return weightSum == 0; }

  /// Transport with a bound-to-bound jacobian of one segment.
  void transport(const BoundMatrix& jacobian) {
    if (!empty()) {
      weightedExcess = jacobian * weightedExcess * jacobian.transpose();
    }
  }

  /// Add one scatterer with its core variances in phi and theta.
  void add(double thicknessInX0, double variancePhi, double varianceTheta) {
    if (thicknessInX0 <= 0 || (variancePhi <= 0 && varianceTheta <= 0)) {
      return;
    }
    const auto [eps, width] = scatteringTailModel(thicknessInX0);
    const double excess = width * width - 1;
    weightedExcess(eBoundPhi, eBoundPhi) += eps * excess * variancePhi;
    weightedExcess(eBoundTheta, eBoundTheta) += eps * excess * varianceTheta;
    weightSum += eps;
    logNoTail += std::log1p(-eps);
  }

  /// The covariance excess of the tail component.
  BoundMatrix covariance() const { return weightedExcess / weightSum; }

  /// The prior weight of the tail component.
  double weight() const { return -std::expm1(logNoTail); }
};

}  // namespace detail
}  // namespace Acts

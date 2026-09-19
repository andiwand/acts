// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFinding/Rz/RzMeasurementGrid.hpp"

#include <optional>

#include "RzState.hpp"

namespace Acts::Experimental::detail::rz {

/// Residual and projected covariance at the current stop.
struct Evaluation {
  double chi2{};
  double timeResidual{};
  double timeGain{};
  bool hasTime{};
  bool pixel{};
  /// `C (H J)^T`, one column per measured coordinate
  Eigen::Matrix<double, eRzSize, 2> ch;
  Eigen::Matrix<double, 2, 1> residual;
  Eigen::Matrix<double, 2, 2> sInv;
};

/// Exact transport shared by surviving hits on the same module plane.
struct Prediction {
  Vector3 planePosition;
  Vector3 normal;
  RzHelix::PlaneStep crossing;
  Eigen::Matrix<double, 3, eRzSize> jPos;
  bool hasJacobian = false;
};

/// Global measurement frame, constructed only for exact evaluation.
struct Placed {
  Vector3 position{Vector3::Zero()};
  /// The direction the measured coordinate is taken along
  Vector3 u{Vector3::Zero()};
  /// The other one, which a strip does not measure
  Vector3 v{Vector3::Zero()};
  Vector3 normal{Vector3::Zero()};
  /// Variance along `u`
  double cov00{};
  double cov01{};
  /// Variance along `v`, unused by a strip
  double cov11{};
  double invLever{};
  double time{};
  double timeVariance{};
  /// Room along `v`, the coordinate a strip does not measure
  double halfV{};
  /// How far from the RZ stop the module may be met
  double maxDistance{};
  bool pixel{};
};

Placed placeMeasurement(const RzModule& module,
                        const RzMeasurement& measurement,
                        const RzMeasurementFrame* frame, double maxDistance);

class MeasurementEvaluator {
 public:
  MeasurementEvaluator(double maxModuleDistance, double gateChi2,
                       double stripMargin)
      : m_maxModuleDistance(maxModuleDistance),
        m_gateChi2(gateChi2),
        m_stripMargin(stripMargin) {}

  /// With Cache=true, supply a cache and reset it whenever the state changes.
  template <bool Cache = false>
  std::optional<Evaluation> evaluate(
      const State& state, const Placed& measurement, bool gate = true,
      bool useTime = true,
      std::optional<Prediction>* prediction = nullptr) const;

 private:
  double m_maxModuleDistance;
  double m_gateChi2;
  double m_stripMargin;
};

void kalmanUpdate(State& state, const Evaluation& evaluation);

}  // namespace Acts::Experimental::detail::rz

// This file is part of the Acts project.
//
// Copyright (C) 2021-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/MeasurementSelector.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"

#include <algorithm>

double Acts::calculateChi2(unsigned int calibratedSize,
                           const double* calibratedRaw,
                           const double* calibratedCovarianceRaw,
                           const BoundVector& predicted,
                           const BoundMatrix& predictedCovariance,
                           const BoundMatrix& projector) {
  return visit_measurement(calibratedSize, [&](auto N) -> double {
    constexpr std::size_t kMeasurementSize = decltype(N)::value;

    typename TrackStateTraits<kMeasurementSize, true>::Measurement calibrated{
        calibratedRaw};

    typename TrackStateTraits<kMeasurementSize, true>::MeasurementCovariance
        calibratedCovariance{calibratedCovarianceRaw};

    using ParametersVector = ActsVector<kMeasurementSize>;

    // Take the projector (measurement mapping function)
    const auto H =
        projector.template topLeftCorner<kMeasurementSize, eBoundSize>().eval();

    // Get the residuals
    ParametersVector res = calibrated - H * predicted;

    // Get the chi2
    return (res.transpose() *
            (calibratedCovariance + H * predictedCovariance * H.transpose())
                .inverse() *
            res)
        .eval()(0, 0);
  });
}

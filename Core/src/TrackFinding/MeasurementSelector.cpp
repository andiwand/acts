// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/MeasurementSelector.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>

namespace Acts {

MeasurementSelector::MeasurementSelector()
    : MeasurementSelector({{GeometryIdentifier(), MeasurementSelectorCuts{}}}) {
}

MeasurementSelector::MeasurementSelector(const MeasurementSelectorCuts& cuts)
    : MeasurementSelector({{GeometryIdentifier(), cuts}}) {}

MeasurementSelector::MeasurementSelector(const Config& config) {
  if (config.empty()) {
    throw std::invalid_argument(
        "MeasurementSelector: Configuration must not be empty");
  }

  std::vector<InternalConfig::InputElement> tmp;
  tmp.reserve(config.size());
  for (std::size_t i = 0; i < config.size(); ++i) {
    GeometryIdentifier geoID = config.idAt(i);
    MeasurementSelectorCuts cuts = config.valueAt(i);

    InternalCutBins internalCuts = convertCutBins(cuts);
    tmp.emplace_back(geoID, std::move(internalCuts));
  }
  m_config = InternalConfig(std::move(tmp));
}

MeasurementSelector::InternalCutBins MeasurementSelector::convertCutBins(
    const MeasurementSelectorCuts& config) {
  InternalCutBins cutBins;

  auto getEtaOrInf = [](const auto& vec, std::size_t bin) {
    if (bin >= vec.size()) {
      return std::numeric_limits<double>::infinity();
    }
    assert(vec[bin] >= 0 && "Eta bins must be positive");
    return vec[bin];
  };

  auto getBinOrBackOrMax = [](const auto& vec, std::size_t bin) {
    using Value = std::remove_reference_t<decltype(vec[0])>;
    static constexpr Value max = std::numeric_limits<Value>::max();
    return vec.empty() ? max : (bin < vec.size() ? vec[bin] : vec.back());
  };

  for (std::size_t bin = 0; bin < config.etaBins.size() + 1; ++bin) {
    InternalCutBin cuts;
    cuts.maxTheta = getEtaOrInf(config.etaBins, bin);
    cuts.maxNumMeasurements =
        getBinOrBackOrMax(config.numMeasurementsCutOff, bin);
    cuts.maxChi2Measurement = getBinOrBackOrMax(config.chi2CutOff, bin);
    cuts.maxChi2Outlier = getBinOrBackOrMax(config.chi2CutOffOutlier, bin);
    cutBins.push_back(cuts);
  }

  return cutBins;
}

double MeasurementSelector::calculateChi2(
    const double* fullCalibrated, const double* fullCalibratedCovariance,
    TrackStateTraits<kMeasurementSizeMax, false>::Parameters predicted,
    TrackStateTraits<kMeasurementSizeMax, false>::Covariance
        predictedCovariance,
    BoundSubspaceIndices projector, unsigned int calibratedSize) const {
  return visit_measurement(
      calibratedSize,
      [&fullCalibrated, &fullCalibratedCovariance, &predicted,
       &predictedCovariance, &projector](auto N) -> double {
        constexpr std::size_t kMeasurementSize = decltype(N)::value;

        typename TrackStateTraits<kMeasurementSize, true>::Calibrated
            calibrated{fullCalibrated};

        typename TrackStateTraits<kMeasurementSize, true>::CalibratedCovariance
            calibratedCovariance{fullCalibratedCovariance};

        using ParametersVector = Vector<kMeasurementSize>;

        std::span<const std::uint8_t, kMeasurementSize> validSubspaceIndices(
            projector.begin(), projector.begin() + kMeasurementSize);
        FixedBoundSubspaceHelper<kMeasurementSize> subspaceHelper(
            validSubspaceIndices);

        // Get the residuals
        ParametersVector res =
            calibrated - subspaceHelper.projectVector(predicted);

        // Get the chi2
        return (res.transpose() *
                (calibratedCovariance +
                 subspaceHelper.projectMatrix(predictedCovariance))
                    .inverse() *
                res)
            .eval()(0, 0);
      });
}

std::pair<double, double> MeasurementSelector::calculateMixtureChi2(
    const double* fullCalibrated, const double* fullCalibratedCovariance,
    TrackStateTraits<kMeasurementSizeMax, false>::Parameters predicted,
    TrackStateTraits<kMeasurementSizeMax, false>::Covariance
        predictedCovariance,
    const BoundMatrix& tailCov, double logPriorRatio, double residualGate,
    double rejectAbove, BoundSubspaceIndices projector,
    unsigned int calibratedSize) const {
  return visit_measurement(
      calibratedSize, [&](auto N) -> std::pair<double, double> {
        constexpr std::size_t kMeasurementSize = decltype(N)::value;

        typename TrackStateTraits<kMeasurementSize, true>::Calibrated
            calibrated{fullCalibrated};
        typename TrackStateTraits<kMeasurementSize, true>::CalibratedCovariance
            calibratedCovariance{fullCalibratedCovariance};

        std::span<const std::uint8_t, kMeasurementSize> validSubspaceIndices(
            projector.begin(), projector.begin() + kMeasurementSize);
        FixedBoundSubspaceHelper<kMeasurementSize> subspaceHelper(
            validSubspaceIndices);

        const Vector<kMeasurementSize> res =
            calibrated - subspaceHelper.projectVector(predicted);
        const SquareMatrix<kMeasurementSize> predictedProj =
            subspaceHelper.projectMatrix(predictedCovariance);
        const SquareMatrix<kMeasurementSize> tailProj =
            subspaceHelper.projectMatrix(tailCov);

        const SquareMatrix<kMeasurementSize> covCore =
            calibratedCovariance + predictedProj;
        const SquareMatrix<kMeasurementSize> invCore = covCore.inverse();
        const double chi2Core = (res.transpose() * invCore * res).eval()(0, 0);
        if (residualGate > 0 && (invCore * tailProj).trace() < residualGate) {
          return {chi2Core, -std::numeric_limits<double>::infinity()};
        }

        // same sums as `predictedCovariance + tailCov` projected
        const SquareMatrix<kMeasurementSize> covTail =
            calibratedCovariance +
            SquareMatrix<kMeasurementSize>(predictedProj + tailProj);
        const double chi2Tail =
            (res.transpose() * covTail.inverse() * res).eval()(0, 0);

        // ln det S_core - ln det S_tail <= 0, and softplus(a) <= max(a, 0) +
        // ln 2, so this bounds the mixture chi2 from below
        const double aMax = logPriorRatio + 0.5 * (chi2Core - chi2Tail);
        const double chi2Min =
            chi2Core - 2 * (std::max(aMax, 0.) + std::numbers::ln2) - 1e-6;
        if (chi2Min >= rejectAbove) {
          return {chi2Min, -std::numeric_limits<double>::infinity()};
        }
        const double logDetCore = std::log(covCore.determinant());
        const double logDetTail = std::log(covTail.determinant());

        const double a = logPriorRatio + 0.5 * (logDetCore - logDetTail) +
                         0.5 * (chi2Core - chi2Tail);
        const double softplus =
            a > 0 ? a + std::log1p(std::exp(-a)) : std::log1p(std::exp(a));
        return {chi2Core - 2 * softplus, a};
      });
}

std::pair<double, double> MeasurementSelector::calculateChi2AndLogDet(
    const double* fullCalibrated, const double* fullCalibratedCovariance,
    TrackStateTraits<kMeasurementSizeMax, false>::Parameters predicted,
    TrackStateTraits<kMeasurementSizeMax, false>::Covariance
        predictedCovariance,
    BoundSubspaceIndices projector, unsigned int calibratedSize) const {
  return visit_measurement(
      calibratedSize,
      [&fullCalibrated, &fullCalibratedCovariance, &predicted,
       &predictedCovariance, &projector](auto N) -> std::pair<double, double> {
        constexpr std::size_t kMeasurementSize = decltype(N)::value;

        typename TrackStateTraits<kMeasurementSize, true>::Calibrated
            calibrated{fullCalibrated};
        typename TrackStateTraits<kMeasurementSize, true>::CalibratedCovariance
            calibratedCovariance{fullCalibratedCovariance};

        std::span<const std::uint8_t, kMeasurementSize> validSubspaceIndices(
            projector.begin(), projector.begin() + kMeasurementSize);
        FixedBoundSubspaceHelper<kMeasurementSize> subspaceHelper(
            validSubspaceIndices);

        const Vector<kMeasurementSize> res =
            calibrated - subspaceHelper.projectVector(predicted);
        const SquareMatrix<kMeasurementSize> cov =
            calibratedCovariance +
            subspaceHelper.projectMatrix(predictedCovariance);
        return {(res.transpose() * cov.inverse() * res).eval()(0, 0),
                std::log(cov.determinant())};
      });
}

MeasurementSelector::Cuts MeasurementSelector::getCutsByTheta(
    const InternalCutBins& config, double theta) {
  // since theta is in [0, pi] and we have a symmetric cut in eta, we can just
  // look at the positive half of the Z axis
  const double constrainedTheta = std::min(theta, std::numbers::pi - theta);

  auto it = std::ranges::find_if(
      config, [constrainedTheta](const InternalCutBin& cuts) {
        return constrainedTheta < cuts.maxTheta;
      });
  assert(it != config.end());
  return {it->maxNumMeasurements, it->maxChi2Measurement, it->maxChi2Outlier};
}

Result<MeasurementSelector::Cuts> MeasurementSelector::getCuts(
    const GeometryIdentifier& geoID, double theta) const {
  // Find the appropriate cuts
  auto cuts = m_config.find(geoID);
  if (cuts == m_config.end()) {
    // for now we consider missing cuts an unrecoverable error
    // TODO consider other options e.g. do not add measurements at all (not
    // even as outliers)
    return CombinatorialKalmanFilterError::MeasurementSelectionFailed;
  }
  return getCutsByTheta(*cuts, theta);
}

}  // namespace Acts

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterError.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cassert>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace Acts {

/// Selection cuts for associating measurements with predicted track
/// parameters on a surface.
///
/// The default configuration only takes the best matching measurement without a
/// cut on the local chi2.
struct MeasurementSelectorCuts {
  /// bins in |eta| to specify variable selections
  std::vector<double> etaBins{};
  /// Maximum local chi2 contribution to classify as measurement.
  std::vector<double> chi2CutOff{15};
  /// Maximum number of associated measurements on a single surface.
  std::vector<std::size_t> numMeasurementsCutOff{1};
  /// Maximum local chi2 contribution to classify as outlier.
  std::vector<double> chi2CutOffOutlier{};
};

/// @brief Measurement selection struct selecting those measurements compatible
/// with the given track parameter against provided criteria on one surface
///
/// The selection criteria could be allowed maximum chi2
/// and allowed maximum number of measurements on one surface
///
/// If there is no compatible measurement, the measurement with the minimum
/// chi2 will be selected and the status will be tagged as an outlier
///
class MeasurementSelector {
 public:
  /// Geometry-dependent cut configuration.
  ///
  /// Different components on the geometry can require different cut settings.
  /// The configuration must either contain explicit settings for all geometry
  /// components that are used or contain a global default.
  using Config = Acts::GeometryHierarchyMap<MeasurementSelectorCuts>;

  /// @brief Default constructor
  ///
  /// This will use the default configuration for the cuts.
  MeasurementSelector();

  /// @brief Constructor with cuts
  ///
  /// @param cuts The cuts to use
  explicit MeasurementSelector(const MeasurementSelectorCuts& cuts);

  /// @brief Constructor with config
  ///
  /// @param config a config instance
  explicit MeasurementSelector(const Config& config);

  /// @brief Function that select the measurements compatible with
  /// the given track parameter on a surface
  ///
  /// @param candidates The track state candidates which already contain predicted parameters
  /// @param isOutlier The indicator for outlier or not
  /// @param logger The logger wrapper
  ///
  /// @return Pair of iterators into @a candidates marking the range of selected candidates
  ///

  template <typename traj_t>
  Result<std::pair<
      typename std::vector<typename traj_t::TrackStateProxy>::iterator,
      typename std::vector<typename traj_t::TrackStateProxy>::iterator>>
  select(std::vector<typename traj_t::TrackStateProxy>& candidates,
         bool& isOutlier, const Logger& logger) const;

  /// Prototype: select with the chi2 of a two-component prediction.
  ///
  /// The tail component has the predicted covariance plus @a tailCov and the
  /// prior weight @a tailWeight. The chi2 is -2 ln of the mixture likelihood,
  /// shifted to be the core chi2 for a small residual, so the cuts keep their
  /// meaning in the core.
  ///
  /// @param tailPosterior is filled with the posterior tail weight of each
  ///        selected candidate, in the order of the returned range
  template <typename traj_t>
  Result<std::pair<
      typename std::vector<typename traj_t::TrackStateProxy>::iterator,
      typename std::vector<typename traj_t::TrackStateProxy>::iterator>>
  selectWithTail(std::vector<typename traj_t::TrackStateProxy>& candidates,
                 const BoundMatrix& tailCov, double tailWeight,
                 std::vector<double>& tailPosterior, bool& isOutlier,
                 const Logger& logger) const;

  /// Prototype: the mixture chi2 and the posterior tail weight of a candidate.
  template <typename track_state_t>
  std::pair<double, double> mixtureChi2(track_state_t& trackState,
                                        const BoundMatrix& tailCov,
                                        double tailWeight) const;

 private:
  template <typename traj_t, typename chi2_fn_t>
  Result<std::pair<
      typename std::vector<typename traj_t::TrackStateProxy>::iterator,
      typename std::vector<typename traj_t::TrackStateProxy>::iterator>>
  selectImpl(std::vector<typename traj_t::TrackStateProxy>& candidates,
             bool& isOutlier, const Logger& logger,
             const chi2_fn_t& chi2Fn) const;

  /// chi2 and ln det of the residual covariance
  std::pair<double, double> calculateChi2AndLogDet(
      const double* fullCalibrated, const double* fullCalibratedCovariance,
      TrackStateTraits<kMeasurementSizeMax, false>::Parameters predicted,
      TrackStateTraits<kMeasurementSizeMax, false>::Covariance
          predictedCovariance,
      BoundSubspaceIndices projector, unsigned int calibratedSize) const;

  struct InternalCutBin {
    double maxTheta{};
    std::size_t maxNumMeasurements{};
    double maxChi2Measurement{};
    double maxChi2Outlier{};
  };
  using InternalCutBins = std::vector<InternalCutBin>;
  using InternalConfig = Acts::GeometryHierarchyMap<InternalCutBins>;

  struct Cuts {
    std::size_t numMeasurements{};
    double chi2Measurement{};
    double chi2Outlier{};
  };

  static InternalCutBins convertCutBins(const MeasurementSelectorCuts& config);

  static Cuts getCutsByTheta(const InternalCutBins& config, double theta);
  Result<Cuts> getCuts(const GeometryIdentifier& geoID, double theta) const;

  double calculateChi2(
      const double* fullCalibrated, const double* fullCalibratedCovariance,
      TrackStateTraits<kMeasurementSizeMax, false>::Parameters predicted,
      TrackStateTraits<kMeasurementSizeMax, false>::Covariance
          predictedCovariance,
      BoundSubspaceIndices projector, unsigned int calibratedSize) const;

  InternalConfig m_config;
};

}  // namespace Acts

#include "MeasurementSelector.ipp"

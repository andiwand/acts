// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackStateProxy.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

template <typename source_link_iterator_t>
using SourceLinkAccessorDelegate =
    Delegate<std::pair<source_link_iterator_t, source_link_iterator_t>(
        const Surface&)>;

template <typename propagator_t, typename traj_t,
          typename source_link_iterator_t>
class MyTrackFinding {
 public:
  using SourceLinkAccessor = SourceLinkAccessorDelegate<source_link_iterator_t>;

  struct Config {
    propagator_t propagator;
  };

  struct Options {
    Options(const GeometryContext& gctx, const MagneticFieldContext& mctx,
            std::reference_wrapper<const CalibrationContext> cctx,
            SourceLinkAccessor sourceLinkAccessor_)
        : geoContext(gctx),
          magFieldContext(mctx),
          calibrationContext(cctx),
          sourceLinkAccessor(std::move(sourceLinkAccessor_)) {}

    /// Context object for the geometry
    std::reference_wrapper<const GeometryContext> geoContext;
    /// Context object for the magnetic field
    std::reference_wrapper<const MagneticFieldContext> magFieldContext;
    /// context object for the calibration
    std::reference_wrapper<const CalibrationContext> calibrationContext;

    SourceLinkAccessor sourceLinkAccessor;
  };

  MyTrackFinding(const Config& cfg, std::unique_ptr<const Logger> logger =
                                        getDefaultLogger("CKF", Logging::INFO))
      : m_cfg(cfg), m_logger(std::move(logger)) {}

  template <typename start_parameters_t, typename track_container_t,
            template <typename> class holder_t,
            typename parameters_t = BoundTrackParameters>
  auto findTracks(
      const start_parameters_t& initialParameters, const Options& options,
      TrackContainer<track_container_t, traj_t, holder_t>& trackContainer) const
      -> Result<std::vector<
          typename std::decay_t<decltype(trackContainer)>::TrackProxy>> {
    using TrackContainer = typename std::decay_t<decltype(trackContainer)>;
    using TrackProxy = typename TrackContainer::TrackProxy;
    using IndexType = typename TrackContainer::IndexType;
    using PM = TrackStatePropMask;

    auto& trackStateContainer = trackContainer.trackStateContainer();

    TrackProxy track = trackContainer.makeTrack();

    PropagatorOptions<ActionList<>, AbortList<AnySurfaceReached>> pOptions(
        options.geoContext, options.magFieldContext);

    auto pState = m_cfg.propagator.makeState(initialParameters, pOptions);

    m_cfg.propagator.initialize(pState);

    IndexType tipIndex = TrackContainer::kInvalid;

    for (size_t i = 0; i < 100; ++i) {
      pState.navigation.startSurface = pState.navigation.currentSurface;
      m_cfg.propagator.propagate(pState);

      const Surface* surface = pState.navigation.currentSurface;
      const TrackingVolume* volume = pState.navigation.currentVolume;

      if (surface == nullptr || volume == nullptr) {
        break;
      }
      ACTS_VERBOSE("surface reached " << surface->geometryId());

      const DetectorElementBase* detectorElement =
          surface->associatedDetectorElement();
      const ISurfaceMaterial* material = surface->surfaceMaterial();

      if (material != nullptr) {
        ACTS_VERBOSE("surface has material");
      }

      if (detectorElement != nullptr) {
        ACTS_VERBOSE("surface has detector element");

        auto [slBegin, slEnd] = options.sourceLinkAccessor(*surface);

        ACTS_VERBOSE("found " << std::distance(slBegin, slEnd)
                              << " source links");

        if (slBegin != slEnd) {
          auto boundRes =
              m_cfg.propagator.stepper().boundState(pState.stepping, *surface);

          const auto& [boundParams, jacobian, pathLength] = *boundRes;

          PM mask = PM::Predicted | PM::Jacobian;
          auto trackState = trackStateContainer.makeTrackState(mask, tipIndex);
          trackState.predicted() = boundParams.parameters();
          trackState.predictedCovariance() = *boundParams.covariance();
          trackState.jacobian() = jacobian;
          trackState.setReferenceSurface(
              boundParams.referenceSurface().getSharedPtr());

          tipIndex = trackState.index();
        }
      }
    }
    ACTS_VERBOSE("track finding done");

    track.tipIndex() = tipIndex;

    std::vector<typename TrackContainer::TrackProxy> tracks;
    tracks.push_back(track);
    return tracks;
  }

 private:
  Config m_cfg;
  std::unique_ptr<const Logger> m_logger;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

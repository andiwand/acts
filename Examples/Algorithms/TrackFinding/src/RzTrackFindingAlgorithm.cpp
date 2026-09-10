// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/RzTrackFindingAlgorithm.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackStateType.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/TrackFinding/Rz/RzBound.hpp"
#include "Acts/TrackFinding/Rz/RzMeasurementGrid.hpp"
#include "Acts/TrackFinding/Rz/RzTransport.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <stdexcept>

namespace ActsExamples {

namespace {

using namespace Acts::Experimental;

/// The free components the RZ state carries, in the order of `RzIndices`
constexpr std::array<unsigned int, eRzSize> kFreeOf = {
    Acts::eFreePos0, Acts::eFreePos1, Acts::eFreePos2,  Acts::eFreeDir0,
    Acts::eFreeDir1, Acts::eFreeDir2, Acts::eFreeQOverP};

}  // namespace

RzTrackFindingAlgorithm::RzTrackFindingAlgorithm(
    const Config& config, std::unique_ptr<const Acts::Logger> log)
    : IAlgorithm("RzTrackFindingAlgorithm", std::move(log)), m_cfg(config) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements input collection");
  }
  if (m_cfg.inputInitialTrackParameters.empty()) {
    throw std::invalid_argument("Missing initial track parameters input");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing tracks output collection");
  }
  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (m_cfg.magneticField == nullptr) {
    throw std::invalid_argument("Missing magnetic field");
  }
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputInitialTrackParameters.initialize(m_cfg.inputInitialTrackParameters);
  m_outputTracks.initialize(m_cfg.outputTracks);

  RzLayoutOptions options;
  options.phiBins = m_cfg.phiBins;
  options.alongBinWidth = m_cfg.alongBinWidth;
  options.materialTables = m_cfg.materialTables;
  if (!m_cfg.excludeVolumes.empty()) {
    options.surfaceSelector =
        [exclude = m_cfg.excludeVolumes](const Acts::Surface& surface) {
          return std::ranges::find(exclude, surface.geometryId().volume()) ==
                 exclude.end();
        };
  }
  // the field along every RZ surface, so that a stop moves in its own Bz
  const Acts::MagneticFieldContext fieldContext;
  auto fieldCache = m_cfg.magneticField->makeCache(fieldContext);
  options.fieldSampler = [&](const Acts::Vector3& position) {
    const auto field = m_cfg.magneticField->getField(position, fieldCache);
    return field.ok() ? *field : Acts::Vector3::Zero();
  };
  m_layout = makeRzLayout(*m_cfg.trackingGeometry,
                          Acts::GeometryContext::dangerouslyDefaultConstruct(),
                          options, logger());
  m_perigee =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3::Zero());
}

ProcessCode RzTrackFindingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const MeasurementContainer& measurements = m_inputMeasurements(ctx);
  const TrackParametersContainer& initialParameters =
      m_inputInitialTrackParameters(ctx);

  // the field at the origin is the field everywhere, for now
  auto fieldCache = m_cfg.magneticField->makeCache(ctx.magFieldContext);
  const auto field =
      m_cfg.magneticField->getField(Acts::Vector3::Zero(), fieldCache);
  if (!field.ok()) {
    ACTS_ERROR("Failed to look up the magnetic field");
    return ProcessCode::ABORT;
  }
  const double bz = field->z();

  using Clock = std::chrono::steady_clock;
  const auto tFill0 = Clock::now();
  RzMeasurementGrid grid(m_layout);
  grid.reserve(measurements.size());
  for (std::uint32_t i = 0; i < measurements.size(); ++i) {
    const auto measurement = measurements.getMeasurement(i);
    const auto module = m_layout.moduleIndex.find(measurement.geometryId());
    if (module == m_layout.moduleIndex.end()) {
      continue;
    }
    // keep the two position components, whatever else was measured
    const auto subspace = measurement.subspaceIndexVector();
    const auto values = measurement.parameters();
    const auto covariance = measurement.covariance();
    std::array<std::uint8_t, 2> indices{};
    std::array<double, 2> params{};
    std::array<double, 4> cov{};
    std::array<std::uint8_t, 2> rows{};
    std::uint8_t dim = 0;
    for (std::size_t a = 0; a < measurement.size() && dim < 2; ++a) {
      if (subspace[a] != Acts::eBoundLoc0 && subspace[a] != Acts::eBoundLoc1) {
        continue;
      }
      indices[dim] = static_cast<std::uint8_t>(subspace[a]);
      params[dim] = values[a];
      rows[dim] = static_cast<std::uint8_t>(a);
      ++dim;
    }
    if (dim == 0) {
      continue;
    }
    for (std::uint8_t a = 0; a < dim; ++a) {
      for (std::uint8_t b = 0; b < dim; ++b) {
        cov[a * 2 + b] = covariance(rows[a], rows[b]);
      }
    }
    grid.addBound(module->second, *m_layout.modules[module->second].surface,
                  ctx.recoGeoContext, dim, std::span(indices.data(), dim),
                  std::span(params.data(), dim), std::span(cov.data(), 4), i);
  }
  const auto tFill1 = Clock::now();
  grid.finalize();
  const auto tFill2 = Clock::now();
  m_nsFill += static_cast<std::size_t>(
      std::chrono::duration_cast<std::chrono::nanoseconds>(tFill1 - tFill0)
          .count());
  m_nsFinalize += static_cast<std::size_t>(
      std::chrono::duration_cast<std::chrono::nanoseconds>(tFill2 - tFill1)
          .count());
  m_nMeasurementsBinned += grid.size();
  ACTS_DEBUG("Binned " << grid.size() << " of " << measurements.size()
                       << " measurements");

  RzTrackFinderConfig finderConfig;
  finderConfig.chi2Cut = m_cfg.chi2Cut;
  finderConfig.windowSigmas = m_cfg.windowSigmas;
  finderConfig.windowMin = m_cfg.windowMin;
  finderConfig.maxHoles = m_cfg.maxHoles;
  finderConfig.maxConsecutiveHoles = m_cfg.maxConsecutiveHoles;
  finderConfig.minMeasurements = m_cfg.minMeasurements;
  finderConfig.maxMeasurementsPerLayer = m_cfg.maxMeasurementsPerLayer;
  finderConfig.ptMin = m_cfg.ptMin;
  finderConfig.minMeasurementsAtLayer = m_cfg.minMeasurementsAtLayer;
  finderConfig.layersForMinMeasurements = m_cfg.layersForMinMeasurements;
  finderConfig.applyMaterial = m_cfg.applyMaterial;
  finderConfig.inwardSearch = m_cfg.inwardSearch;
  finderConfig.backwardPass = m_cfg.backwardPass;
  finderConfig.backwardInflation = m_cfg.backwardInflation;
  finderConfig.backwardLayers = m_cfg.backwardLayers;
  finderConfig.backwardQOverPScale = m_cfg.backwardQOverPScale;
  const RzTrackFinder finder(finderConfig, m_layout, bz);
  // the inner hits sit next to the beam line, where the field is the central
  // one
  const RzHelix helix = finder.helixAt(finder.bz());

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);
  PassThroughCalibrator calibrator;

  const RzMeasurementAccessor accessor = grid.accessor();
  const auto tFind0 = Clock::now();
  // room for the tracks and their states, so that writing does not grow
  // and copy: on average a track a seed here, sixteen states a track
  trackContainer->reserve(initialParameters.size());
  trackStateContainer->reserve(16 * initialParameters.size());
  RzTrackCandidate candidate;
  std::size_t nTracks = 0;
  std::size_t nStops = 0;
  std::size_t nCandidates = 0;
  std::size_t nMeasurements = 0;
  std::size_t nHoles = 0;
  std::size_t nBackwardFailures = 0;
  for (const TrackParameters& start : initialParameters) {
    // bound to the RZ free state, time dropped
    const Acts::Vector3 position = start.position(ctx.recoGeoContext);
    const Acts::Vector3 direction = start.direction();
    RzVector v;
    v.segment<3>(eRzPos0) = position;
    v.segment<3>(eRzDir0) = direction;
    v[eRzQOverP] = start.qOverP();
    RzMatrix c = RzMatrix::Zero();
    if (start.covariance().has_value()) {
      const Acts::BoundToFreeMatrix j =
          start.referenceSurface().boundToFreeJacobian(ctx.recoGeoContext,
                                                       position, direction);
      const Acts::FreeMatrix free = j * (*start.covariance()) * j.transpose();
      for (unsigned int a = 0; a < eRzSize; ++a) {
        for (unsigned int b = 0; b < eRzSize; ++b) {
          c(a, b) = free(kFreeOf[a], kFreeOf[b]);
        }
      }
    }
    std::uint32_t startModule = kRzNone;
    if (const auto it =
            m_layout.moduleIndex.find(start.referenceSurface().geometryId());
        it != m_layout.moduleIndex.end()) {
      startModule = it->second;
    }

    const bool found = finder.findTrack(accessor, v, c, startModule, candidate);
    nStops += candidate.stops;
    nCandidates += candidate.candidatesTested;
    if (!found) {
      continue;
    }
    ++nTracks;
    nMeasurements += candidate.measurements;
    nHoles += candidate.holes;
    if (candidate.backwardFailure != 0) {
      ++nBackwardFailures;
      ACTS_DEBUG("Backward pass failed with "
                 << candidate.backwardFailure << " on a track with "
                 << candidate.measurements << " measurements");
    }
    const auto tMake = Clock::now();
    // the inner state to the perigee, the closest approach to the beam axis
    RzVector w =
        candidate.hasInner ? candidate.innerParameters : candidate.parameters;
    const RzMatrix& cInner =
        candidate.hasInner ? candidate.innerCovariance : candidate.covariance;
    const double s = helix.pathToPerigee(w);
    RzMatrix j = helix.stepJacobian(w, s);
    helix.step(w, s);
    const double dt = std::hypot(w[eRzDir0], w[eRzDir1]);
    const Acts::Vector3 normal(w[eRzDir0] / dt, w[eRzDir1] / dt, 0.);
    RzHelix::constrainToSurface(j, helix.derivative(w), normal);
    const RzMatrix cPerigee = j * cInner * j.transpose();
    const Acts::Vector3 pos = w.segment<3>(eRzPos0);
    const Acts::Vector3 dir = w.segment<3>(eRzDir0);
    const auto bound = Acts::transformFreeToBoundParameters(
        pos, 0., dir, w[eRzQOverP], *m_perigee, ctx.recoGeoContext);
    if (!bound.ok()) {
      ACTS_WARNING("Perigee conversion failed: " << bound.error().message());
      continue;
    }
    // the perigee's own Jacobian, on the seven components the RZ state has;
    // the product is formed on them rather than on the 8x8 with a zero time
    const Acts::FreeToBoundMatrix jf2b =
        m_perigee->freeToBoundJacobian(ctx.recoGeoContext, pos, dir);
    RzFreeToBoundMatrix jPerigee;
    for (unsigned int r = 0; r < Acts::eBoundSize; ++r) {
      for (unsigned int a = 0; a < eRzSize; ++a) {
        jPerigee(r, a) = jf2b(r, kFreeOf[a]);
      }
    }
    const Acts::BoundMatrix boundCov = rzBoundCovariance(jPerigee, cPerigee);

    auto track = tracks.makeTrack();
    track.setReferenceSurface(m_perigee);
    track.parameters() = *bound;
    track.covariance() = boundCov;
    const auto tStates = Clock::now();
    for (const RzTrackHit& hit : candidate.hits) {
      if (hit.isHole()) {
        auto state = track.appendTrackState(Acts::TrackStatePropMask::None);
        state.typeFlags().setUnchecked(Acts::TrackStateFlag::IsHole);
        state.setReferenceSurface(
            m_layout.surfaces[m_layout.layers[hit.layer].surface].surface);
        continue;
      }
      auto state = track.appendTrackState(Acts::TrackStatePropMask::Filtered |
                                          Acts::TrackStatePropMask::Calibrated);
      const RzModuleMeasurements on = grid.moduleRange(hit.module);
      const RzMeasurement& m = on.entries[hit.measurement];
      const RzModule& module = m_layout.modules[hit.module];
      state.setReferenceSurface(module.surface);
      const IndexSourceLink sourceLink(module.geometryId, m.source);
      calibrator.calibrate(measurements, nullptr, ctx.recoGeoContext,
                           ctx.calibContext, Acts::SourceLink{sourceLink},
                           state);
      state.typeFlags().setUnchecked(Acts::TrackStateFlag::HasMeasurement);
      state.chi2() = static_cast<float>(hit.chi2);
      if (hit.forwardState == kRzNone) {
        state.filtered().setZero();
        state.filteredCovariance().setIdentity();
        continue;
      }
      // The finder updates at the RZ stop the measurement was found from,
      // not on the module, so the state is walked the last bit onto the
      // module plane, which is the surface the track state lives on
      const auto& [v0, c0] = candidate.forwardStates[hit.forwardState];
      const Acts::Vector3& au =
          on.frames.empty() ? module.u : on.frames[hit.measurement].u;
      const Acts::Vector3& av =
          on.frames.empty() ? module.v : on.frames[hit.measurement].v;
      const Acts::Vector3 measured = module.center + m.loc0 * au + m.loc1 * av;
      RzVector v = v0;
      RzMatrix c;
      if (const std::optional<double> step =
              helix.pathToPlane(v0, measured, module.normal);
          step.has_value()) {
        helix.step(v, *step);
        c = helix.stepJacobianOnto(v0, *step, v, module.normal).transport(c0);
      } else {
        c = c0;
      }
      std::optional<RzBoundState> onModule = rzBoundOnModule(module, v, c);
      if (!onModule.has_value()) {
        // through the surface, for a polar module the layout could not
        // confirm as plain polar
        const Acts::Vector3 position = v.segment<3>(eRzPos0);
        const Acts::Vector3 direction = v.segment<3>(eRzDir0);
        const auto stateBound = Acts::transformFreeToBoundParameters(
            position, 0., direction, v[eRzQOverP], *module.surface,
            ctx.recoGeoContext);
        if (!stateBound.ok()) {
          state.filtered().setZero();
          state.filteredCovariance().setIdentity();
          continue;
        }
        Acts::FreeMatrix freeCov = Acts::FreeMatrix::Zero();
        for (unsigned int a = 0; a < eRzSize; ++a) {
          for (unsigned int b = 0; b < eRzSize; ++b) {
            freeCov(kFreeOf[a], kFreeOf[b]) = c(a, b);
          }
        }
        const Acts::FreeToBoundMatrix jm = module.surface->freeToBoundJacobian(
            ctx.recoGeoContext, position, direction);
        Acts::BoundMatrix stateCov = jm * freeCov * jm.transpose();
        stateCov(Acts::eBoundTime, Acts::eBoundTime) = 1.;
        onModule = RzBoundState{*stateBound, stateCov};
      }
      state.filtered() = onModule->parameters;
      state.filteredCovariance() = onModule->covariance;
    }
    Acts::calculateTrackQuantities(track);
    const auto tEnd = Clock::now();
    m_nsMakeStates += static_cast<std::size_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tStates)
            .count());
    m_nsMake += static_cast<std::size_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tMake)
            .count());
  }

  const auto tFind1 = Clock::now();
  m_nsFind += static_cast<std::size_t>(
      std::chrono::duration_cast<std::chrono::nanoseconds>(tFind1 - tFind0)
          .count());
  m_nSeeds += initialParameters.size();
  m_nTracks += nTracks;
  m_nStops += nStops;
  m_nCandidates += nCandidates;
  m_nMeasurementsOnTracks += nMeasurements;
  m_nHolesOnTracks += nHoles;
  m_nBackwardFailures += nBackwardFailures;
  ACTS_DEBUG("Found " << nTracks << " tracks from " << initialParameters.size()
                      << " seeds");

  auto constTrackStateContainer =
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer));
  auto constTrackContainer = std::make_shared<Acts::ConstVectorTrackContainer>(
      std::move(*trackContainer));
  m_outputTracks(
      ctx, ConstTrackContainer{constTrackContainer, constTrackStateContainer});
  return ProcessCode::SUCCESS;
}

ProcessCode RzTrackFindingAlgorithm::finalize() {
  const double seeds = std::max<double>(1., m_nSeeds.load());
  const double tracks = std::max<double>(1., m_nTracks.load());
  ACTS_INFO("RzTrackFinding: "
            << m_nSeeds << " seeds, " << m_nTracks << " tracks, "
            << m_nStops / seeds << " stops and " << m_nCandidates / seeds
            << " candidates per seed, " << m_nMeasurementsOnTracks / tracks
            << " measurements and " << m_nHolesOnTracks / tracks
            << " holes per track, " << m_nBackwardFailures
            << " backward pass failures");
  const double binned = std::max<double>(1., m_nMeasurementsBinned.load());
  const double ms = 1e6;
  ACTS_INFO(
      "RzTrackFinding timing: fill "
      << m_nsFill / ms << " ms (" << m_nsFill / binned
      << " ns per measurement), finalize " << m_nsFinalize / ms << " ms, find "
      << m_nsFind / ms << " ms (" << m_nsFind / seeds / 1e3
      << " us per seed, of which " << m_nsMake / seeds / 1e3
      << " us writing: " << m_nsMake / tracks / 1e3 << " us per track, "
      << m_nsMakeStates / tracks / 1e3 << " us for its "
      << m_nMeasurementsOnTracks / tracks << " states), "
      << 100. * (m_nsFill + m_nsFinalize) /
             std::max<double>(1., m_nsFill + m_nsFinalize + m_nsFind)
      << "% spent preparing " << m_nMeasurementsBinned << " measurements");
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

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
#include "Acts/TrackFinding/Rz/RzMeasurementGrid.hpp"
#include "Acts/TrackFinding/Rz/RzTransport.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"

#include <array>
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

  RzMeasurementGrid grid(m_layout);
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
    grid.add(module->second, dim, std::span(indices.data(), dim),
             std::span(params.data(), dim), std::span(cov.data(), 4), i);
  }
  grid.finalize();
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
  finderConfig.applyMaterial = m_cfg.applyMaterial;
  finderConfig.backwardInflation = m_cfg.backwardInflation;
  finderConfig.backwardLayers = m_cfg.backwardLayers;
  finderConfig.backwardQOverPScale = m_cfg.backwardQOverPScale;
  const RzTrackFinder finder(finderConfig, m_layout, bz);
  const RzHelix& helix = finder.helix();

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);
  PassThroughCalibrator calibrator;

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

    const bool found = finder.findTrack(grid, v, c, startModule, candidate);
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
    Acts::FreeMatrix freeCov = Acts::FreeMatrix::Zero();
    for (unsigned int a = 0; a < eRzSize; ++a) {
      for (unsigned int b = 0; b < eRzSize; ++b) {
        freeCov(kFreeOf[a], kFreeOf[b]) = cPerigee(a, b);
      }
    }
    const Acts::FreeToBoundMatrix jf2b =
        m_perigee->freeToBoundJacobian(ctx.recoGeoContext, pos, dir);
    Acts::BoundMatrix boundCov = jf2b * freeCov * jf2b.transpose();
    // no time on the RZ state; a finite variance keeps the matrix invertible
    boundCov(Acts::eBoundTime, Acts::eBoundTime) = 1.;

    auto track = tracks.makeTrack();
    track.setReferenceSurface(m_perigee);
    track.parameters() = *bound;
    track.covariance() = boundCov;
    for (const RzTrackHit& hit : candidate.hits) {
      auto state = track.appendTrackState(
          hit.isHole() ? Acts::TrackStatePropMask::None
                       : Acts::TrackStatePropMask::Calibrated);
      if (hit.isHole()) {
        state.typeFlags().setUnchecked(Acts::TrackStateFlag::IsHole);
        state.setReferenceSurface(
            m_layout.surfaces[m_layout.layers[hit.layer].surface].surface);
        continue;
      }
      const RzMeasurement& m = grid.entry(hit.measurement);
      state.setReferenceSurface(m_layout.modules[m.module].surface);
      const IndexSourceLink sourceLink(m_layout.modules[m.module].geometryId,
                                       m.source);
      calibrator.calibrate(measurements, nullptr, ctx.recoGeoContext,
                           ctx.calibContext, Acts::SourceLink{sourceLink},
                           state);
      state.typeFlags().setUnchecked(Acts::TrackStateFlag::HasMeasurement);
      state.chi2() = static_cast<float>(hit.chi2);
    }
    Acts::calculateTrackQuantities(track);
  }

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
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

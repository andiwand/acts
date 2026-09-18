// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/TrackFinding/Rz/RzTrackFinder.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

using namespace Acts;
using namespace Acts::Experimental;

namespace {
RzLayout makeLayout(std::uint32_t modules = 1) {
  RzLayout layout;
  RzSurface surface;
  surface.shape = RzShape::Disc;
  surface.minBound = 1.;
  surface.maxBound = 200.;
  surface.layer = 0;
  layout.surfaces.push_back(surface);
  RzLayer layer;
  layer.surface = 0;
  layer.phiBins = 1;
  layer.alongBins = 1;
  layer.alongMin = 1.;
  layer.alongMax = 200.;
  layer.maxHalfExtent = 30.;
  layout.layers.push_back(layer);
  for (std::uint32_t i = 0; i < modules; ++i) {
    RzModule module;
    module.center = Vector3(100., 0., 0.);
    module.u = Vector3::UnitX();
    module.v = Vector3::UnitY();
    module.normal = Vector3::UnitZ();
    module.halfU = 20.;
    module.halfV = 20.;
    module.layer = 0;
    layout.modules.push_back(module);
    layout.moduleOrder.push_back(i);
  }
  layout.moduleBinStart = {0, modules};
  layout.escapeRadius = 200.;
  layout.escapeHalfZ = 200.;
  return layout;
}

RzTrackFinderConfig config() {
  RzTrackFinderConfig cfg;
  cfg.minMeasurements = 1;
  cfg.maxMeasurementsPerLayer = 1;
  cfg.applyMaterial = false;
  cfg.backwardPass = false;
  return cfg;
}

RzVector start() {
  RzVector v = RzVector::Zero();
  v[eRzPos0] = 100.;
  v[eRzDir2] = 1.;
  v[eRzQOverP] = 1.;
  return v;
}

RzMatrix covariance() {
  RzMatrix c = RzMatrix::Zero();
  c(eRzPos0, eRzPos0) = 1.;
  c(eRzPos1, eRzPos1) = 1.;
  return c;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(RzTrackFinderSuite)

BOOST_AUTO_TEST_CASE(RejectedGateWinnerDoesNotHideValidStrip) {
  const RzLayout layout = makeLayout();
  RzMeasurementGrid grid(layout);
  RzMeasurement strip;
  strip.projector = RzProjector::Loc0;
  strip.cov00 = 1.;
  strip.loc1 = 100.;  // Best measured residual, but outside the strip extent.
  grid.add(0, strip);
  strip.loc0 = 1.;
  strip.loc1 = 0.;
  grid.add(0, strip);
  grid.finalize();
  RzTrackCandidate candidate;
  BOOST_REQUIRE(
      RzTrackFinder(config(), layout, 0.)
          .findTrack(grid.accessor(), start(), covariance(), 0, candidate));
  BOOST_REQUIRE_EQUAL(candidate.hits.size(), 1u);
  BOOST_CHECK_EQUAL(candidate.hits.front().measurement, 1u);
}

BOOST_AUTO_TEST_CASE(ExactChi2OrdersCorrelatedPixels) {
  const RzLayout layout = makeLayout();
  RzMeasurementGrid grid(layout);
  RzMeasurement pixel;
  pixel.cov00 = 1.;
  pixel.cov11 = 1.;
  pixel.cov01 = 0.9;
  pixel.loc0 = 1.;
  pixel.loc1 = -1.;
  grid.add(0, pixel);
  pixel.loc0 = 1.1;
  pixel.loc1 = 1.1;
  grid.add(0, pixel);
  grid.finalize();
  RzTrackCandidate candidate;
  BOOST_REQUIRE(
      RzTrackFinder(config(), layout, 0.)
          .findTrack(grid.accessor(), start(), covariance(), 0, candidate));
  BOOST_CHECK_EQUAL(candidate.hits.front().measurement, 1u);
}

BOOST_AUTO_TEST_CASE(TimeRejectsCloserPileupHit) {
  const RzLayout layout = makeLayout();
  RzMeasurementGrid grid(layout);
  RzMeasurement pixel;
  pixel.cov00 = 1.;
  pixel.cov11 = 1.;
  pixel.time = 100.;
  pixel.timeVariance = 1.;
  grid.add(0, pixel);
  pixel.loc0 = 1.;
  pixel.time = 0.;
  grid.add(0, pixel);
  grid.finalize();
  RzTrackStart startState{start(), covariance()};
  startState.module = 0;
  startState.timeVariance = 1.;
  const std::array starts{startState};
  bool found = false;
  RzTrackFinder(config(), layout, 0.)
      .findTracks(grid.accessor(), starts, 1,
                  [&](std::size_t, bool ok, RzTrackCandidate& candidate) {
                    found = ok;
                    BOOST_REQUIRE_EQUAL(candidate.hits.size(), 1u);
                    BOOST_CHECK_EQUAL(candidate.hits.front().measurement, 1u);
                  });
  BOOST_CHECK(found);
}

BOOST_AUTO_TEST_CASE(SearchIncludesModulesBeyondInlineCapacity) {
  const RzLayout layout = makeLayout(12);
  RzMeasurementGrid grid(layout);
  RzMeasurement pixel;
  pixel.cov00 = 1.;
  pixel.cov11 = 1.;
  grid.add(11, pixel);
  grid.finalize();
  RzTrackCandidate candidate;
  BOOST_REQUIRE(
      RzTrackFinder(config(), layout, 0.)
          .findTrack(grid.accessor(), start(), covariance(), 0, candidate));
  BOOST_CHECK_EQUAL(candidate.hits.front().module, 11u);
}

BOOST_AUTO_TEST_CASE(PolarCovarianceUsesAngularLeverOnly) {
  for (bool plain : {false, true}) {
    for (std::uint8_t dim : {1, 2}) {
      RzLayout layout = makeLayout();
      RzModule& module = layout.modules.front();
      module.polar = true;
      module.polarIsPlain = plain;
      module.localCenter = Vector2(100., 0.);
      module.boundCenter = Vector2(100., 0.);
      const auto surface = Surface::makeShared<DiscSurface>(
          Transform3::Identity(), std::make_shared<RadialBounds>(1., 200.));
      const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
      RzMeasurementGrid grid(layout);
      const std::array<std::uint8_t, 2> indices{0, 1};
      const std::array<double, 2> params{100., 0.};
      const std::array<double, 4> cov{4., 0.001, 0.001, 0.0001};
      grid.addBound(0, *surface, gctx, dim, indices, params, cov, 0);
      grid.finalize();
      RzVector v = start();
      v[eRzPos0] = 110.;
      v[eRzPos1] = 1.;
      auto cfg = config();
      cfg.chi2Cut = 100.;
      RzTrackCandidate candidate;
      const std::array<RzSeedMeasurement, 1> seed{{{0, 0}}};
      BOOST_REQUIRE(
          RzTrackFinder(cfg, layout, 0.)
              .findTrack(grid.accessor(), v, covariance(), 0, candidate, seed));
      // Radial variance stays 4. Angular variance and cross covariance
      // scale from radius 100 to the radial projection 110.
      const double expected = dim == 1 ? 100. / 5.
                                       : (100. * 2.21 - 20. * 0.11 + 5.) /
                                             (5. * 2.21 - 0.11 * 0.11);
      BOOST_CHECK_CLOSE(candidate.chi2, expected, 1e-8);
    }
  }
}

BOOST_AUTO_TEST_CASE(BackwardScatteringHasSignedPositionDirectionCovariance) {
  RzLayout layout = makeLayout();
  RzSurface outer = layout.surfaces.front();
  outer.refCoord = 10.;
  outer.layer = 1;
  outer.materialEdges = {1., 200.};
  outer.materialBands.emplace_back(
      Material::fromMolarDensity(93.7f, 465.2f, 28.0855f, 14.f, 0.083f), 1.f);
  layout.surfaces.push_back(outer);
  RzLayer layer = layout.layers.front();
  layer.surface = 1;
  layer.binOffset = 1;
  layout.layers.push_back(layer);
  RzModule module = layout.modules.front();
  module.center.z() = 10.;
  module.layer = 1;
  layout.modules.push_back(module);
  layout.moduleOrder = {0, 1};
  layout.moduleBinStart = {0, 1, 2};
  layout.discs = {0, 1};
  layout.discCoord = {0., 10.};
  layout.discMin = {1., 1.};
  layout.discMax = {200., 200.};
  RzMeasurementGrid grid(layout);
  RzMeasurement pixel;
  pixel.cov00 = 1.;
  pixel.cov11 = 1.;
  grid.add(0, pixel);
  grid.add(1, pixel);
  grid.finalize();
  auto cfg = config();
  cfg.applyMaterial = true;
  cfg.backwardPass = true;
  cfg.backwardLayers = 0;
  cfg.inwardSearch = false;
  // Isolate the backward process noise from the initial covariance.
  cfg.backwardInflation = 0.;
  RzTrackCandidate candidate;
  BOOST_REQUIRE(
      RzTrackFinder(cfg, layout, 0.)
          .findTrack(grid.accessor(), start(), covariance(), 0, candidate));
  BOOST_REQUIRE_EQUAL(candidate.measurements, 2u);
  BOOST_REQUIRE(candidate.hasInner);
  BOOST_REQUIRE_EQUAL(candidate.backwardFailure, 0u);
  BOOST_CHECK_LT(candidate.innerCovariance(eRzPos0, eRzDir0), 0.);
  BOOST_CHECK_LT(candidate.innerCovariance(eRzPos1, eRzDir1), 0.);
}

BOOST_AUTO_TEST_CASE(BatchedSearchMatchesIndividualTracks) {
  RzLayout layout = makeLayout(12);
  layout.discs = {0};
  layout.discCoord = {0.};
  layout.discMin = {1.};
  layout.discMax = {200.};
  RzMeasurementGrid grid(layout);
  RzMeasurement pixel;
  pixel.cov00 = 1.;
  pixel.cov11 = 1.;
  pixel.cov01 = 0.9;
  pixel.loc0 = 1.;
  pixel.loc1 = -1.;
  grid.add(11, pixel);
  pixel.loc0 = 1.1;
  pixel.loc1 = 1.1;
  grid.add(11, pixel);
  grid.finalize();
  const RzTrackFinder finder(config(), layout, 0.);
  std::vector<RzTrackStart> starts(5);
  std::vector<RzTrackCandidate> expected(starts.size());
  std::vector<bool> found(starts.size());
  for (std::size_t i = 0; i < starts.size(); ++i) {
    auto& seed = starts[i];
    seed.parameters = start();
    seed.covariance = covariance();
    // Approach the disc from both z directions, with an empty search among
    // successful ones to check callback order and state reuse across batches.
    seed.parameters[eRzDir2] = i % 2 == 0 ? 1. : -1.;
    seed.parameters[eRzPos2] = -10. * seed.parameters[eRzDir2];
    if (i == 2) {
      seed.parameters[eRzPos0] = 500.;
    }
    found[i] = finder.findTrack(grid.accessor(), seed.parameters,
                                seed.covariance, seed.module, expected[i]);
    BOOST_REQUIRE_EQUAL(found[i], i != 2);
    if (found[i]) {
      BOOST_REQUIRE_EQUAL(expected[i].hits.front().module, 11u);
      BOOST_REQUIRE_EQUAL(expected[i].hits.front().measurement, 1u);
    }
  }
  for (std::size_t width : {0u, 1u, 2u, 3u, 8u}) {
    std::size_t called = 0;
    finder.findTracks(
        grid.accessor(), starts, width,
        [&](std::size_t index, bool ok, RzTrackCandidate& candidate) {
          BOOST_REQUIRE_EQUAL(index, called++);
          const auto& reference = expected[index];
          BOOST_CHECK_EQUAL(ok, found[index]);
          BOOST_CHECK_EQUAL(candidate.measurements, reference.measurements);
          BOOST_CHECK_EQUAL(candidate.chi2, reference.chi2);
          BOOST_CHECK_EQUAL(candidate.modulesTested, reference.modulesTested);
          BOOST_CHECK_EQUAL(candidate.binsVisited, reference.binsVisited);
          BOOST_CHECK_EQUAL(candidate.exactEvaluated, reference.exactEvaluated);
          BOOST_CHECK(candidate.parameters == reference.parameters);
          BOOST_CHECK(candidate.covariance == reference.covariance);
        });
    BOOST_CHECK_EQUAL(called, starts.size());
  }
}

BOOST_AUTO_TEST_CASE(CheckpointHistoryMatchesFullHistory) {
  RzLayout layout = makeLayout();
  for (unsigned int i = 1; i < 5; ++i) {
    auto surface = layout.surfaces.front();
    surface.refCoord = 10. * i;
    surface.layer = i;
    layout.surfaces.push_back(surface);
    auto layer = layout.layers.front();
    layer.surface = i;
    layer.binOffset = i;
    layout.layers.push_back(layer);
    auto module = layout.modules.front();
    module.center.z() = surface.refCoord;
    module.layer = i;
    layout.modules.push_back(module);
    layout.moduleOrder.push_back(i);
    layout.moduleBinStart.push_back(i + 1);
  }
  RzMeasurementGrid grid(layout);
  for (unsigned int i = 0; i < 5; ++i) {
    layout.discs.push_back(i);
    layout.discCoord.push_back(10. * i);
    layout.discMin.push_back(1.);
    layout.discMax.push_back(200.);
    RzMeasurement pixel;
    pixel.cov00 = 1.;
    pixel.cov11 = 1.;
    pixel.loc0 = 0.01 * i;
    // Include a hole before the checkpoint: it must count measurements.
    if (i != 1) {
      grid.add(i, pixel);
    }
  }
  grid.finalize();
  for (bool backward : {false, true}) {
    for (unsigned int layers : {0u, 2u, 6u}) {
      for (bool inward : {false, true}) {
        auto cfg = config();
        cfg.backwardPass = backward;
        cfg.backwardLayers = layers;
        cfg.inwardSearch = inward;
        RzTrackCandidate full, checkpoint;
        RzMatrix c = covariance();
        c.block<3, 3>(eRzDir0, eRzDir0).diagonal().setConstant(0.001);
        c(eRzQOverP, eRzQOverP) = 0.01;
        BOOST_REQUIRE(RzTrackFinder(cfg, layout, 0.)
                          .findTrack(grid.accessor(), start(), c, 0, full));
        cfg.storeForwardStates = false;
        BOOST_REQUIRE(
            RzTrackFinder(cfg, layout, 0.)
                .findTrack(grid.accessor(), start(), c, 0, checkpoint));
        BOOST_REQUIRE_EQUAL(full.measurements, 4u);
        BOOST_CHECK_EQUAL(checkpoint.measurements, full.measurements);
        BOOST_CHECK_EQUAL(checkpoint.holes, full.holes);
        BOOST_CHECK_EQUAL(checkpoint.backwardFailure, full.backwardFailure);
        BOOST_CHECK_EQUAL(checkpoint.chi2, full.chi2);
        BOOST_CHECK(checkpoint.parameters == full.parameters);
        BOOST_CHECK(checkpoint.covariance == full.covariance);
        BOOST_CHECK_EQUAL(checkpoint.hasInner, full.hasInner);
        BOOST_CHECK(checkpoint.innerParameters == full.innerParameters);
        BOOST_CHECK(checkpoint.innerCovariance == full.innerCovariance);
        BOOST_CHECK_LE(checkpoint.forwardStates.size(), 1u);
        BOOST_REQUIRE_EQUAL(checkpoint.hits.size(), full.hits.size());
        for (std::size_t i = 0; i < full.hits.size(); ++i) {
          BOOST_CHECK_EQUAL(checkpoint.hits[i].module, full.hits[i].module);
          BOOST_CHECK_EQUAL(checkpoint.hits[i].measurement,
                            full.hits[i].measurement);
        }
        if (backward && layers == 2) {
          BOOST_REQUIRE_EQUAL(checkpoint.forwardStates.size(), 1u);
          BOOST_CHECK(checkpoint.forwardStates.front() ==
                      full.forwardStates[1]);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(PartialBackwardUsesCheckpointField) {
  RzLayout layout = makeLayout();
  for (unsigned int i = 1; i < 3; ++i) {
    auto surface = layout.surfaces.front();
    surface.refCoord = 10. * i;
    surface.layer = i;
    layout.surfaces.push_back(surface);
    auto layer = layout.layers.front();
    layer.surface = i;
    layer.binOffset = i;
    layout.layers.push_back(layer);
    auto module = layout.modules.front();
    module.center.z() = surface.refCoord;
    module.layer = i;
    layout.modules.push_back(module);
    layout.moduleOrder.push_back(i);
    layout.moduleBinStart.push_back(i + 1);
    layout.discs.push_back(i);
    layout.discCoord.push_back(surface.refCoord);
    layout.discMin.push_back(1.);
    layout.discMax.push_back(200.);
  }
  layout.surfaces[0].bzTable = {0.002};
  layout.surfaces[1].bzTable = {0.003};
  for (auto& surface : layout.surfaces) {
    surface.fieldBinWidth = 200.;
  }
  RzMeasurementGrid grid(layout);
  RzMeasurement pixel;
  pixel.cov00 = 1.;
  pixel.cov11 = 1.;
  for (unsigned int i = 0; i < 3; ++i) {
    grid.add(i, pixel);
  }
  grid.finalize();
  auto cfg = config();
  cfg.backwardPass = true;
  cfg.backwardLayers = 2;
  cfg.inwardSearch = false;
  RzMatrix c = covariance();
  c.block<3, 3>(eRzDir0, eRzDir0).diagonal().setConstant(0.001);
  RzVector v = start();
  v[eRzDir0] = 0.1;
  v[eRzDir2] = std::sqrt(0.99);
  RzTrackCandidate reference, changed;
  layout.surfaces[2].bzTable = {0.004};
  BOOST_REQUIRE(RzTrackFinder(cfg, layout, 0.)
                    .findTrack(grid.accessor(), v, c, 0, reference));
  layout.surfaces[2].bzTable = {-0.004};
  BOOST_REQUIRE(RzTrackFinder(cfg, layout, 0.)
                    .findTrack(grid.accessor(), v, c, 0, changed));
  BOOST_REQUIRE_EQUAL(reference.backwardFailure, 0u);
  BOOST_REQUIRE_EQUAL(changed.backwardFailure, 0u);
  BOOST_CHECK(
      reference.innerParameters.isApprox(changed.innerParameters, 1e-12));
  BOOST_CHECK(
      reference.innerCovariance.isApprox(changed.innerCovariance, 1e-12));
}

BOOST_AUTO_TEST_SUITE_END()

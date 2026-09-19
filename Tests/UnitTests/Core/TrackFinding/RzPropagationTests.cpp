// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <ranges>

#include "RzPropagation.hpp"

using namespace Acts;
using namespace Acts::Experimental;
namespace rz = Acts::Experimental::detail::rz;

static_assert(std::ranges::input_range<rz::SensitiveCrossings>);
static_assert(!std::ranges::forward_range<rz::SensitiveCrossings>);

namespace {
RzLayout discLayout() {
  RzLayout layout;
  layout.escapeRadius = 200.;
  layout.escapeHalfZ = 200.;
  for (std::uint32_t i = 0; i < 4; ++i) {
    RzSurface surface;
    surface.shape = RzShape::Disc;
    surface.refCoord = 5. * (i + 1);
    surface.minBound = 1.;
    surface.maxBound = 200.;
    surface.layer = i % 2 == 0 ? kRzNone : i / 2;
    layout.surfaces.push_back(surface);
    layout.discs.push_back(i);
    layout.discCoord.push_back(surface.refCoord);
    layout.discMin.push_back(surface.minBound);
    layout.discMax.push_back(surface.maxBound);
  }
  return layout;
}

rz::State start() {
  rz::State state;
  state.v = RzVector::Zero();
  state.v[eRzPos0] = 10.;
  state.v[eRzDir2] = 1.;
  state.v[eRzQOverP] = 1.;
  state.anchor = state.v;
  state.c = RzMatrix::Zero();
  return state;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(RzPropagationSuite)

BOOST_AUTO_TEST_CASE(NextCrossingUsesFilteredDirection) {
  RzLayout layout = discLayout();
  RzSurface cylinder;
  cylinder.shape = RzShape::Cylinder;
  cylinder.refCoord = 25.;
  cylinder.minBound = -200.;
  cylinder.maxBound = 200.;
  cylinder.layer = 2;
  layout.surfaces.push_back(cylinder);
  layout.cylinders = {4};
  layout.cylCoord = {25.};

  RzTrackFinderConfig cfg;
  cfg.applyMaterial = false;
  const rz::Stepper stepper(cfg.radialField, 0.);
  const rz::Propagator propagator(cfg, layout, stepper);
  auto state = start();
  state.v[eRzDir0] = state.v[eRzDir2] = std::sqrt(0.5);
  state.anchor = state.v;
  rz::PropagationState propagation;
  propagator.initialize(propagation, state.v);
  RzTrackCandidate candidate;
  std::uint32_t count = 0;
  for (const auto& crossing :
       propagator.sensitiveCrossings(state, propagation, candidate)) {
    BOOST_REQUIRE_EQUAL(crossing.layer, count);
    BOOST_CHECK_EQUAL(crossing.stop, 2 * count + 1);
    // A filtered direction parallel to z must invalidate the cylinder solve.
    state.v[eRzDir0] = 0.;
    state.v[eRzDir2] = 1.;
    state.anchor = state.v;
    ++count;
  }
  BOOST_CHECK(propagation.status == rz::PropagationStatus::NoTarget);
  BOOST_CHECK_EQUAL(count, 2u);
  BOOST_CHECK_EQUAL(candidate.stops, 4u);
  BOOST_CHECK_CLOSE(state.v[eRzPos0], 20., 1e-10);
  BOOST_CHECK_EQUAL(state.v[eRzPos2], 20.);
}

BOOST_AUTO_TEST_CASE(PassiveMaterialPrecedesYieldAndBreakStopsTransport) {
  RzLayout layout = discLayout();
  for (auto& surface : layout.surfaces) {
    surface.materialEdges = {1., 200.};
    surface.materialBands.emplace_back(
        Material::fromMolarDensity(93.7f, 465.2f, 28.0855f, 14.f, 0.083f), 1.f);
  }
  const RzTrackFinderConfig cfg;
  const rz::Stepper stepper(cfg.radialField, 0.);
  const rz::Propagator propagator(cfg, layout, stepper);
  auto state = start();
  rz::PropagationState propagation;
  propagator.initialize(propagation, state.v);
  RzTrackCandidate candidate;
  std::uint32_t count = 0;
  for (const auto& crossing :
       propagator.sensitiveCrossings(state, propagation, candidate)) {
    BOOST_CHECK_EQUAL(crossing.layer, 0u);
    BOOST_CHECK_EQUAL(crossing.stop, 1u);
    BOOST_CHECK_GT(state.v[eRzQOverP], 1.);
    BOOST_CHECK_GT(state.c(eRzPos0, eRzPos0), 0.);
    BOOST_CHECK_GT(state.c(eRzQOverP, eRzQOverP), 0.);
    BOOST_CHECK(state.pending.empty());
    ++count;
    break;
  }
  BOOST_CHECK(propagation.status == rz::PropagationStatus::Active);
  BOOST_CHECK_EQUAL(count, 1u);
  BOOST_CHECK_EQUAL(candidate.stops, 2u);
  BOOST_CHECK_EQUAL(state.v[eRzPos2], 10.);
}

BOOST_AUTO_TEST_CASE(ExhaustedTurningBudgetYieldsNothing) {
  const RzLayout layout = discLayout();
  RzTrackFinderConfig cfg;
  cfg.maxTurningAngle = 0.;
  const rz::Stepper stepper(cfg.radialField, 2. * UnitConstants::T);
  const rz::Propagator propagator(cfg, layout, stepper);
  auto state = start();
  state.bz = state.anchorBz = 2. * UnitConstants::T;
  rz::PropagationState propagation;
  propagator.initialize(propagation, state.v);
  RzTrackCandidate candidate;
  for ([[maybe_unused]] const auto& crossing :
       propagator.sensitiveCrossings(state, propagation, candidate)) {
    BOOST_FAIL("Propagation exceeded its turning budget");
  }
  BOOST_CHECK(propagation.status == rz::PropagationStatus::TurningLimit);
  BOOST_CHECK_EQUAL(candidate.stops, 0u);
  BOOST_CHECK_EQUAL(state.v[eRzPos2], 0.);
}

BOOST_AUTO_TEST_CASE(MaterialFailureStopsPropagation) {
  RzLayout layout = discLayout();
  auto& surface = layout.surfaces.front();
  surface.materialEdges = {1., 200.};
  surface.materialBands.emplace_back(
      Material::fromMolarDensity(93.7f, 465.2f, 28.0855f, 14.f, 0.083f),
      10000.f);
  const RzTrackFinderConfig cfg;
  const rz::Stepper stepper(cfg.radialField, 0.);
  const rz::Propagator propagator(cfg, layout, stepper);
  auto state = start();
  rz::PropagationState propagation;
  propagator.initialize(propagation, state.v);
  RzTrackCandidate candidate;
  for ([[maybe_unused]] const auto& crossing :
       propagator.sensitiveCrossings(state, propagation, candidate)) {
    BOOST_FAIL("Propagation continued after exhausting the particle's energy");
  }
  BOOST_CHECK(propagation.status == rz::PropagationStatus::MaterialFailure);
  BOOST_CHECK_EQUAL(candidate.stops, 1u);
  BOOST_CHECK(!propagator.advance(state, propagation, candidate));
  BOOST_CHECK_EQUAL(candidate.stops, 1u);
}

BOOST_AUTO_TEST_SUITE_END()

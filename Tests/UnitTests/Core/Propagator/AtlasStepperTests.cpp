// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldError.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsTests/CommonHelpers/Assertions.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

using VectorHelpers::makeVector4;
using Covariance = BoundMatrix;
using Jacobian = BoundMatrix;
using Stepper = AtlasStepper;

// epsilon for floating point comparisons
static constexpr auto eps = 1024 * std::numeric_limits<double>::epsilon();

// propagation settings
static constexpr auto stepSize = 10_mm;
static constexpr Direction navDir = Direction::Backward();
static auto magneticField =
    std::make_shared<ConstantBField>(Vector3(0.1_T, -0.2_T, 2_T));

// initial parameter state
static const Vector4 pos4(1_mm, -1_mm, 2_mm, 2_ns);
static const Vector3 pos = pos4.segment<3>(ePos0);
static const auto time = pos4[eTime];
static const Vector3 unitDir = Vector3(-2, 2, 1).normalized();
static constexpr auto absMom = 1_GeV;
static constexpr auto charge = -1_e;
static const auto particleHypothesis = ParticleHypothesis::pion();
static const Covariance cov = Covariance::Identity();

// context objects
static const auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();
static const MagneticFieldContext magCtx;

BOOST_AUTO_TEST_SUITE(PropagatorSuite)

// test state construction from parameters w/o covariance
BOOST_AUTO_TEST_CASE(ConstructState) {
  BoundTrackParameters cp = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, std::nullopt, particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options);
  stepper.initialize(state, cp);

  BOOST_CHECK(!state.cov.has_value());
  BOOST_CHECK_EQUAL(state.pVector[0], pos.x());
  BOOST_CHECK_EQUAL(state.pVector[1], pos.y());
  BOOST_CHECK_EQUAL(state.pVector[2], pos.z());
  BOOST_CHECK_EQUAL(state.pVector[3], time);
  CHECK_CLOSE_ABS(state.pVector[4], unitDir.x(), eps);
  CHECK_CLOSE_ABS(state.pVector[5], unitDir.y(), eps);
  CHECK_CLOSE_ABS(state.pVector[6], unitDir.z(), eps);
  BOOST_CHECK_EQUAL(state.pVector[7], charge / absMom);
  BOOST_CHECK_EQUAL(state.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(state.stepSize.value(), stepSize);
  BOOST_CHECK_EQUAL(state.previousStepSize, 0.);
}

// test state construction from parameters w/ covariance
BOOST_AUTO_TEST_CASE(ConstructStateWithCovariance) {
  BoundTrackParameters cp = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options);
  stepper.initialize(state, cp);

  BOOST_CHECK(state.cov.has_value());
  BOOST_CHECK_EQUAL(*state.cov, cov);
  BOOST_CHECK_EQUAL(state.pVector[0], pos.x());
  BOOST_CHECK_EQUAL(state.pVector[1], pos.y());
  BOOST_CHECK_EQUAL(state.pVector[2], pos.z());
  BOOST_CHECK_EQUAL(state.pVector[3], time);
  CHECK_CLOSE_ABS(state.pVector[4], unitDir.x(), eps);
  CHECK_CLOSE_ABS(state.pVector[5], unitDir.y(), eps);
  CHECK_CLOSE_ABS(state.pVector[6], unitDir.z(), eps);
  BOOST_CHECK_EQUAL(state.pVector[7], charge / absMom);
  BOOST_CHECK_EQUAL(state.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(state.stepSize.value(), stepSize);
  BOOST_CHECK_EQUAL(state.previousStepSize, 0.);
}

// test stepper getters for particle state
BOOST_AUTO_TEST_CASE(Getters) {
  BoundTrackParameters cp = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options);
  stepper.initialize(state, cp);

  CHECK_CLOSE_ABS(stepper.position(state), pos, eps);
  CHECK_CLOSE_ABS(stepper.time(state), time, eps);
  CHECK_CLOSE_ABS(stepper.direction(state), unitDir, eps);
  CHECK_CLOSE_ABS(stepper.absoluteMomentum(state), absMom, eps);
  BOOST_CHECK_EQUAL(stepper.charge(state), charge);
}

// test stepper update methods with bound state as input
BOOST_AUTO_TEST_CASE(UpdateFromBound) {
  BoundTrackParameters cp = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options);
  stepper.initialize(state, cp);

  auto newPos4 = (pos4 + Vector4(1_mm, 2_mm, 3_mm, 20_ns)).eval();
  auto newPos = newPos4.segment<3>(ePos0);
  auto newTime = newPos4[eTime];
  auto newUnitDir = (unitDir + Vector3(1, -1, -1)).normalized();
  auto newAbsMom = 0.9 * absMom;

  // example surface and bound parameters at the updated position
  std::shared_ptr<PlaneSurface> plane =
      CurvilinearSurface(newPos, newUnitDir).planeSurface();
  auto params =
      BoundTrackParameters::create(geoCtx, plane, newPos4, newUnitDir,
                                   charge / absMom, cov, particleHypothesis)
          .value();
  FreeVector freeParams;
  freeParams[eFreePos0] = newPos4[ePos0];
  freeParams[eFreePos1] = newPos4[ePos1];
  freeParams[eFreePos2] = newPos4[ePos2];
  freeParams[eFreeTime] = newPos4[eTime];
  freeParams[eFreeDir0] = newUnitDir[eMom0];
  freeParams[eFreeDir1] = newUnitDir[eMom1];
  freeParams[eFreeDir2] = newUnitDir[eMom2];
  freeParams[eFreeQOverP] = charge / newAbsMom;

  BOOST_CHECK(params.covariance().has_value());
  stepper.update(state, freeParams, params.parameters(), *params.covariance(),
                 *plane);
  CHECK_CLOSE_ABS(stepper.position(state), newPos, eps);
  CHECK_CLOSE_ABS(stepper.time(state), newTime, eps);
  CHECK_CLOSE_ABS(stepper.direction(state), newUnitDir, eps);
  CHECK_CLOSE_ABS(stepper.absoluteMomentum(state), newAbsMom, eps);
  BOOST_CHECK_EQUAL(stepper.charge(state), charge);
}

// test stepper update methods with individual components as input
BOOST_AUTO_TEST_CASE(UpdateFromComponents) {
  BoundTrackParameters cp = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options);
  stepper.initialize(state, cp);

  auto newPos = (pos + Vector3(1_mm, 2_mm, 3_mm)).eval();
  auto newTime = time + 20_ns;
  auto newUnitDir = (unitDir + Vector3(1, -1, -1)).normalized();
  auto newAbsMom = 0.9 * absMom;

  stepper.update(state, newPos, newUnitDir, charge / newAbsMom, newTime);
  CHECK_CLOSE_ABS(stepper.position(state), newPos, eps);
  CHECK_CLOSE_ABS(stepper.time(state), newTime, eps);
  CHECK_CLOSE_ABS(stepper.direction(state), newUnitDir, eps);
  CHECK_CLOSE_ABS(stepper.absoluteMomentum(state), newAbsMom, eps);
  BOOST_CHECK_EQUAL(stepper.charge(state), charge);
}

// test building a bound state object from the stepper state
BOOST_AUTO_TEST_CASE(BuildBound) {
  BoundTrackParameters cp = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options);
  stepper.initialize(state, cp);

  // example surface at the current state position
  std::shared_ptr<PlaneSurface> plane =
      CurvilinearSurface(pos, unitDir).planeSurface();

  const Jacobian jac = stepper.transportToBound(state, *plane).value();
  const auto pars = stepper.boundParameters(state, *plane).value();
  const double pathLength = stepper.pathLength(state);
  // check parameters
  CHECK_CLOSE_ABS(pars.position(geoCtx), pos, eps);
  CHECK_CLOSE_ABS(pars.time(), time, eps);
  CHECK_CLOSE_ABS(pars.momentum(), absMom * unitDir, eps);
  BOOST_CHECK_EQUAL(pars.charge(), charge);
  BOOST_CHECK(pars.covariance().has_value());
  CHECK_CLOSE_ABS(*pars.covariance(), cov, eps);
  // check Jacobian. should be identity since there was no propagation yet
  CHECK_CLOSE_ABS(jac, Jacobian(Jacobian::Identity()), eps);
  // check propagation length
  CHECK_CLOSE_ABS(pathLength, 0., eps);
}

// test building a curvilinear state object from the stepper state
BOOST_AUTO_TEST_CASE(BuildCurvilinear) {
  BoundTrackParameters cp = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options);
  stepper.initialize(state, cp);

  const Jacobian jac = stepper.transportToCurvilinear(state);
  const auto pars = stepper.curvilinearParameters(state);
  const double pathLength = stepper.pathLength(state);
  // check parameters
  CHECK_CLOSE_ABS(pars.position(geoCtx), pos, eps);
  CHECK_CLOSE_ABS(pars.time(), time, eps);
  CHECK_CLOSE_ABS(pars.momentum(), absMom * unitDir, eps);
  BOOST_CHECK_EQUAL(pars.charge(), charge);
  BOOST_CHECK(pars.covariance().has_value());
  CHECK_CLOSE_ABS(*pars.covariance(), cov, eps);
  // check Jacobian. should be identity since there was no propagation yet
  CHECK_CLOSE_ABS(jac, Jacobian(Jacobian::Identity()), eps);
  // check propagation length
  CHECK_CLOSE_ABS(pathLength, 0., eps);
}

// test step method without covariance transport
BOOST_AUTO_TEST_CASE(Step) {
  BoundTrackParameters cp = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  auto state = stepper.makeState(options);
  stepper.initialize(state, cp);
  state.cov.reset();

  // ensure step does not result in an error
  auto res = stepper.step(state, Direction::Backward(), nullptr);
  BOOST_CHECK(res.ok());

  // extract the actual step size
  auto h = res.value();
  BOOST_CHECK_EQUAL(state.stepSize.value(), stepSize);
  BOOST_CHECK_EQUAL(state.stepSize.value(), h * navDir);

  // check that the position has moved
  auto deltaPos = (stepper.position(state) - pos).eval();
  BOOST_CHECK_LT(0, deltaPos.norm());
  // check that time has changed
  auto deltaTime = stepper.time(state) - time;
  BOOST_CHECK_LT(0, std::abs(deltaTime));
  // check that the direction was rotated
  auto projDir = stepper.direction(state).dot(unitDir);
  BOOST_CHECK_LT(projDir, 1);

  // momentum and charge should not change
  CHECK_CLOSE_ABS(stepper.absoluteMomentum(state), absMom, eps);
  BOOST_CHECK_EQUAL(stepper.charge(state), charge);
}

// test step method with covariance transport
BOOST_AUTO_TEST_CASE(StepWithCovariance) {
  BoundTrackParameters cp = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  auto state = stepper.makeState(options);
  stepper.initialize(state, cp);

  // ensure step does not result in an error
  auto res = stepper.step(state, Direction::Backward(), nullptr);
  BOOST_CHECK(res.ok());

  // extract the actual step size
  auto h = res.value();
  BOOST_CHECK_EQUAL(state.stepSize.value(), stepSize);
  BOOST_CHECK_EQUAL(state.stepSize.value(), h * navDir);

  // check that the position has moved
  auto deltaPos = (stepper.position(state) - pos).eval();
  BOOST_CHECK_LT(0, deltaPos.norm());
  // check that time has changed
  auto deltaTime = stepper.time(state) - time;
  BOOST_CHECK_LT(0, std::abs(deltaTime));
  // check that the direction was rotated
  auto projDir = stepper.direction(state).dot(unitDir);
  BOOST_CHECK_LT(projDir, 1);

  // momentum and charge should not change
  CHECK_CLOSE_ABS(stepper.absoluteMomentum(state), absMom, eps);
  BOOST_CHECK_EQUAL(stepper.charge(state), charge);

  stepper.transportToCurvilinear(state);
  BOOST_CHECK_NE(*state.cov, cov);
}

// test state reset method
BOOST_AUTO_TEST_CASE(Reset) {
  BoundTrackParameters cp = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  auto state = stepper.makeState(options);
  stepper.initialize(state, cp);

  // ensure step does not result in an error
  BOOST_CHECK(stepper.step(state, Direction::Backward(), nullptr).ok());

  // Construct the parameters
  Vector3 newPos(1.5, -2.5, 3.5);
  auto newAbsMom = 4.2 * absMom;
  double newTime = 7.5;
  double newCharge = 1.;
  BoundMatrix newCov = 8.5 * Covariance::Identity();
  cp = BoundTrackParameters::createCurvilinear(makeVector4(newPos, newTime),
                                               unitDir, newCharge / newAbsMom,
                                               newCov, particleHypothesis);
  FreeVector freeParams = transformBoundToFreeParameters(
      cp.referenceSurface(), geoCtx, cp.parameters());

  auto copyState = [&](auto& field, const auto& other) {
    using field_t = std::decay_t<decltype(field)>;
    std::decay_t<decltype(other)> copy = stepper.makeState(options);
    stepper.initialize(copy, cp);

    copy.newfield = other.newfield;
    copy.field = other.field;
    copy.pVector = other.pVector;
    copy.derivativeReady = other.derivativeReady;
    copy.cov = other.cov;
    copy.pathAccumulated = other.pathAccumulated;
    copy.stepSize = other.stepSize;
    copy.previousStepSize = other.previousStepSize;

    copy.fieldCache = MagneticFieldProvider::Cache(
        std::in_place_type<typename field_t::Cache>,
        other.fieldCache.template as<typename field_t::Cache>());

    return copy;
  };

  // Reset all possible parameters
  Stepper::State stateCopy = copyState(*magneticField, state);
  BOOST_CHECK(cp.covariance().has_value());
  stepper.initialize(stateCopy, cp.parameters(), *cp.covariance(),
                     cp.particleHypothesis(), cp.referenceSurface());
  // Test all components
  BOOST_CHECK(stateCopy.cov.has_value());
  BOOST_CHECK_EQUAL(*stateCopy.cov, newCov);
  BOOST_CHECK_EQUAL(stepper.position(stateCopy),
                    freeParams.template segment<3>(eFreePos0));
  CHECK_CLOSE_ABS(stepper.direction(stateCopy),
                  freeParams.template segment<3>(eFreeDir0).normalized(), eps);
  BOOST_CHECK_EQUAL(stepper.absoluteMomentum(stateCopy),
                    std::abs(1. / freeParams[eFreeQOverP]));
  BOOST_CHECK_EQUAL(stepper.charge(stateCopy), newCharge);
  BOOST_CHECK_EQUAL(stepper.time(stateCopy), freeParams[eFreeTime]);
  BOOST_CHECK_EQUAL(stateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(stateCopy.stepSize.value(), stepSize);
  BOOST_CHECK_EQUAL(stateCopy.previousStepSize, state.previousStepSize);

  // Reset using different surface shapes
  // 1) Disc surface
  // Setting some parameters
  newPos << 0.5, -1.5, 0.;
  newAbsMom *= 1.23;
  newTime = 8.4;
  newCharge = -1.;
  newCov = 10.9 * Covariance::Identity();
  Transform3 trafo = Transform3::Identity();
  auto disc = Surface::makeShared<DiscSurface>(trafo);
  auto boundDisc = BoundTrackParameters::create(
                       geoCtx, disc, makeVector4(newPos, newTime), unitDir,
                       newCharge / newAbsMom, newCov, particleHypothesis)
                       .value();

  // Reset the state and test
  Stepper::State stateDisc = copyState(*magneticField, state);
  BOOST_CHECK(boundDisc.covariance().has_value());
  stepper.initialize(stateDisc, boundDisc.parameters(), *boundDisc.covariance(),
                     cp.particleHypothesis(), boundDisc.referenceSurface());

  CHECK_NE_COLLECTIONS(stateDisc.pVector, stateCopy.pVector);
  CHECK_NE_COLLECTIONS(stateDisc.pVector, state.pVector);

  // 2) Perigee surface
  // Setting some parameters
  newPos << -2.06155, -2.06155, 3.5;
  newAbsMom *= 0.45;
  newTime = 2.3;
  newCharge = 1.;
  newCov = 8.7 * Covariance::Identity();
  auto perigee = Surface::makeShared<PerigeeSurface>(trafo);
  auto boundPerigee =
      BoundTrackParameters::create(
          geoCtx, perigee, makeVector4(newPos, newTime), unitDir,
          newCharge / newAbsMom, newCov, particleHypothesis)
          .value();

  // Reset the state and test
  Stepper::State statePerigee = copyState(*magneticField, state);
  BOOST_CHECK(boundPerigee.covariance().has_value());
  stepper.initialize(statePerigee, boundPerigee.parameters(),
                     *boundPerigee.covariance(), cp.particleHypothesis(),
                     boundPerigee.referenceSurface());
  CHECK_NE_COLLECTIONS(statePerigee.pVector, stateCopy.pVector);
  CHECK_NE_COLLECTIONS(statePerigee.pVector, state.pVector);
  CHECK_NE_COLLECTIONS(statePerigee.pVector, stateDisc.pVector);

  // 3) Straw surface
  // Use the same parameters as for previous Perigee surface
  auto straw = Surface::makeShared<StrawSurface>(trafo);
  auto boundStraw = BoundTrackParameters::create(
                        geoCtx, straw, makeVector4(newPos, newTime), unitDir,
                        newCharge / newAbsMom, newCov, particleHypothesis)
                        .value();

  // Reset the state and test
  Stepper::State stateStraw = copyState(*magneticField, state);
  BOOST_CHECK(boundStraw.covariance().has_value());
  stepper.initialize(stateStraw, boundStraw.parameters(),
                     *boundStraw.covariance(), cp.particleHypothesis(),
                     boundStraw.referenceSurface());
  CHECK_NE_COLLECTIONS(stateStraw.pVector, stateCopy.pVector);
  CHECK_NE_COLLECTIONS(stateStraw.pVector, state.pVector);
  CHECK_NE_COLLECTIONS(stateStraw.pVector, stateDisc.pVector);
  BOOST_CHECK_EQUAL_COLLECTIONS(
      stateStraw.pVector.begin(), stateStraw.pVector.end(),
      statePerigee.pVector.begin(), statePerigee.pVector.end());
}

BOOST_AUTO_TEST_CASE(StepSize) {
  BoundTrackParameters cp = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options);
  stepper.initialize(state, cp);

  stepper.updateStepSize(state, -5_cm, ConstrainedStep::Type::Navigator);
  BOOST_CHECK_EQUAL(state.previousStepSize, stepSize);
  BOOST_CHECK_EQUAL(state.stepSize.value(), -5_cm);

  stepper.releaseStepSize(state, ConstrainedStep::Type::Navigator);
  BOOST_CHECK_EQUAL(state.stepSize.value(), stepSize);
}

// test step size modification with target surfaces
BOOST_AUTO_TEST_CASE(StepSizeSurface) {
  BoundTrackParameters cp = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);

  Stepper stepper(magneticField);

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;

  Stepper::State state = stepper.makeState(options);
  stepper.initialize(state, cp);

  auto distance = 10_mm;
  std::shared_ptr<PlaneSurface> target =
      CurvilinearSurface(pos + navDir * distance * unitDir, unitDir)
          .planeSurface();

  stepper.updateSurfaceStatus(
      state, *target, 0, navDir, BoundaryTolerance::Infinite(),
      s_onSurfaceTolerance, ConstrainedStep::Type::Navigator);
  BOOST_CHECK_EQUAL(state.stepSize.value(ConstrainedStep::Type::Navigator),
                    distance);

  const auto getNavigationTarget = [&](const Surface& s,
                                       const BoundaryTolerance& bt) {
    auto [intersection, intersectionIndex] =
        s.intersect(geoCtx, stepper.position(state),
                    navDir * stepper.direction(state), bt)
            .closestWithIndex();
    return NavigationTarget(intersection, intersectionIndex, s, bt);
  };

  // test the step size modification in the context of a surface
  stepper.updateStepSize(
      state, getNavigationTarget(*target, BoundaryTolerance::Infinite()),
      navDir, ConstrainedStep::Type::Navigator);
  BOOST_CHECK_EQUAL(state.stepSize.value(), distance);

  // start with a different step size
  state.stepSize.setUser(navDir * stepSize);
  stepper.releaseStepSize(state, ConstrainedStep::Type::Navigator);
  stepper.updateStepSize(
      state, getNavigationTarget(*target, BoundaryTolerance::Infinite()),
      navDir, ConstrainedStep::Type::Navigator);
  BOOST_CHECK_EQUAL(state.stepSize.value(), navDir * stepSize);
}

// Record lookups in a non-uniform field so stale cached values are observable.
class RecordingField final : public MagneticFieldProvider {
 public:
  // A non-symmetric, divergence-free linear field exercises every derivative.
  SquareMatrix3 gradient = [] {
    SquareMatrix3 g;
    g << 0.003, 0.001, 0.002, 0.01, -0.002, 0.004, 0.005, 0.006, -0.001;
    return SquareMatrix3(g * (1_T / 1_mm));
  }();
  mutable std::vector<Vector3> lookups;
  mutable std::size_t failAt = std::numeric_limits<std::size_t>::max();

  Cache makeCache(const MagneticFieldContext& /*mctx*/) const override {
    return Cache(std::in_place_type<int>, 0);
  }

  Result<Vector3> getField(const Vector3& position,
                           Cache& /*cache*/) const override {
    lookups.push_back(position);
    if (lookups.size() == failAt) {
      return Result<Vector3>::failure(MagneticFieldError::OutOfBounds);
    }
    return Result<Vector3>::success(Vector3(0., 0., 1_T) + gradient * position);
  }
};

BOOST_AUTO_TEST_CASE(ReinitializeAndUpdateInvalidateField) {
  auto field = std::make_shared<RecordingField>();
  Stepper stepper(field);
  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;
  auto state = stepper.makeState(options);
  const auto initial = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);
  stepper.initialize(state, initial);
  BOOST_REQUIRE(stepper.step(state, Direction::Forward(), nullptr).ok());

  const auto updated = BoundTrackParameters::createCurvilinear(
      pos4 + Vector4(100_mm, 0., 0., 0.), unitDir, -charge / absMom,
      std::nullopt, particleHypothesis);
  stepper.initialize(state, updated);
  BOOST_CHECK(!stepper.hasCovariance(state));
  field->lookups.clear();
  BOOST_REQUIRE(stepper.step(state, Direction::Forward(), nullptr).ok());
  CHECK_CLOSE_ABS(field->lookups.front(), updated.position(geoCtx), eps);
  // No Jacobian should be propagated after removing the covariance.
  for (std::size_t i = eFreeSize; i < 56; ++i) {
    BOOST_CHECK_EQUAL(state.pVector[i], 0.);
  }
  auto fresh = stepper.makeState(options);
  stepper.initialize(fresh, updated);
  BOOST_REQUIRE(stepper.step(fresh, Direction::Forward(), nullptr).ok());
  CHECK_CLOSE_ABS(stepper.position(state), stepper.position(fresh), eps);
  CHECK_CLOSE_ABS(stepper.direction(state), stepper.direction(fresh), eps);

  // The component update must also invalidate a cached field.
  stepper.update(state, pos, unitDir, charge / absMom, time);
  field->lookups.clear();
  BOOST_REQUIRE(stepper.step(state, Direction::Forward(), nullptr).ok());
  CHECK_CLOSE_ABS(field->lookups.front(), pos, eps);

  // A bound update must accept the supplied charge sign and enable covariance.
  const FreeVector free = transformBoundToFreeParameters(
      updated.referenceSurface(), geoCtx, updated.parameters());
  stepper.update(state, free, updated.parameters(), cov,
                 updated.referenceSurface());
  BOOST_CHECK_EQUAL(stepper.qOverP(state), -charge / absMom);
  BOOST_CHECK(stepper.hasCovariance(state));
  field->lookups.clear();
  BOOST_REQUIRE(stepper.step(state, Direction::Forward(), nullptr).ok());
  CHECK_CLOSE_ABS(field->lookups.front(), updated.position(geoCtx), eps);
}

BOOST_AUTO_TEST_CASE(FieldQueriesAndFailuresDoNotPoisonCachedField) {
  auto field = std::make_shared<RecordingField>();
  Stepper stepper(field);
  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;
  auto state = stepper.makeState(options);
  const auto initial = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);
  stepper.initialize(state, initial);
  BOOST_REQUIRE(stepper.step(state, Direction::Forward(), nullptr).ok());
  const auto parameters = state.pVector;
  const Vector3 acceptedField = state.field;
  const double path = state.pathAccumulated;
  BOOST_REQUIRE(stepper.getField(state, Vector3(1000_mm, 0., 0.)).ok());
  CHECK_CLOSE_ABS(state.field, acceptedField, 0.);

  // Fail at the last stage, after a successful middle-stage lookup.
  field->failAt = field->lookups.size() + 2;
  const auto failed = stepper.step(state, Direction::Forward(), nullptr);
  BOOST_REQUIRE(!failed.ok());
  BOOST_CHECK(failed.error() == MagneticFieldError::OutOfBounds);
  CHECK_CLOSE_ABS(state.field, acceptedField, 0.);
  BOOST_CHECK_EQUAL_COLLECTIONS(state.pVector.begin(), state.pVector.end(),
                                parameters.begin(), parameters.end());
  BOOST_CHECK_EQUAL(state.pathAccumulated, path);
  BOOST_CHECK_EQUAL(state.nSteps, 1u);
  BOOST_CHECK_EQUAL(state.statistics.nSuccessfulSteps, 1u);
}

BOOST_AUTO_TEST_CASE(NeutralParticleDoesNotBend) {
  Stepper stepper(magneticField);
  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;
  auto state = stepper.makeState(options);
  const auto neutral =
      ParticleHypothesis::pion0().withMomentumHypothesis(absMom);
  const auto initial =
      BoundTrackParameters::createCurvilinear(pos4, unitDir, 0., cov, neutral);
  stepper.initialize(state, initial);
  const auto result = stepper.step(state, Direction::Forward(), nullptr);
  BOOST_REQUIRE(result.ok());
  CHECK_CLOSE_ABS(stepper.direction(state), unitDir, eps);
  CHECK_CLOSE_ABS(stepper.position(state), Vector3(pos + *result * unitDir),
                  eps);
  const double dtds = std::hypot(1., neutral.mass() / absMom);
  CHECK_CLOSE_ABS(stepper.time(state), time + *result * dtds, eps);
  const auto jac = stepper.transportToCurvilinear(state);
  CHECK_CLOSE_ABS(jac(eBoundPhi, eBoundQOverP), 0., eps);
  CHECK_CLOSE_ABS(jac(eBoundTheta, eBoundQOverP), 0., eps);
  CHECK_CLOSE_ABS(jac(eBoundTime, eBoundQOverP), 0., eps);
  BOOST_CHECK_EQUAL(stepper.qOverP(state), 0.);
}

BOOST_AUTO_TEST_CASE(TransportBeforeSteppingIncludesPathCorrections) {
  Stepper stepper(magneticField);
  Stepper::Options options(geoCtx, magCtx);
  auto state = stepper.makeState(options);
  const auto initial = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);
  stepper.initialize(state, initial);
  auto target = CurvilinearSurface(pos, Vector3::UnitZ()).planeSurface();
  const auto jac = stepper.transportToBound(state, *target).value();
  const auto frame =
      initial.referenceSurface().referenceFrame(geoCtx, pos, unitDir);
  const double dtds = std::hypot(1., particleHypothesis.mass() / absMom);
  for (int i = 0; i < 2; ++i) {
    CHECK_CLOSE_ABS(jac(eBoundTime, i), -dtds * frame(2, i) / unitDir.z(), eps);
  }
  // Reanchoring must allow another transport without an intervening step.
  const auto inverse = stepper.transportToCurvilinear(state);
  CHECK_CLOSE_ABS(BoundMatrix(inverse * jac), BoundMatrix::Identity().eval(),
                  1e-12);
  CHECK_CLOSE_ABS(*state.cov, cov, 1e-12);
}

BOOST_AUTO_TEST_CASE(CurvilinearFrameNearBeam) {
  Stepper stepper(magneticField);
  Stepper::Options options(geoCtx, magCtx);
  auto state = stepper.makeState(options);
  const Vector3 nearBeam = Vector3(1e-3, 2e-3, 1.).normalized();
  const auto initial = BoundTrackParameters::createCurvilinear(
      pos4, nearBeam, charge / absMom, cov, particleHypothesis);
  stepper.initialize(state, initial);
  const auto jac = stepper.transportToCurvilinear(state);
  CHECK_CLOSE_ABS(jac, BoundMatrix::Identity().eval(), 1e-10);
  CHECK_CLOSE_ABS(*state.cov, cov, 1e-10);
}

BOOST_AUTO_TEST_CASE(StepAdjustmentLimitsAndZeroStep) {
  Stepper stepper(magneticField);
  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = 1_m;
  options.stepTolerance = 1e-20;
  options.maxRungeKuttaStepTrials = 0;
  auto state = stepper.makeState(options);
  const auto initial = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);
  stepper.initialize(state, initial);
  const auto parameters = state.pVector;
  auto result = stepper.step(state, Direction::Forward(), nullptr);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() == EigenStepperError::StepSizeAdjustmentFailed);
  BOOST_CHECK_EQUAL(state.statistics.nRejectedSteps, 1u);
  BOOST_CHECK_EQUAL_COLLECTIONS(state.pVector.begin(), state.pVector.end(),
                                parameters.begin(), parameters.end());

  state.options.stepSizeCutOff = 1_m;
  stepper.initialize(state, initial);
  result = stepper.step(state, Direction::Forward(), nullptr);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() == EigenStepperError::StepSizeStalled);

  stepper.initialize(state, initial);
  stepper.updateStepSize(state, 0., ConstrainedStep::Type::Navigator);
  result = stepper.step(state, Direction::Forward(), nullptr);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(*result, 0.);
  BOOST_CHECK_EQUAL(state.nSteps, 0u);
  BOOST_CHECK_EQUAL(state.statistics.nAttemptedSteps, 0u);
  BOOST_CHECK_EQUAL(state.statistics.nSuccessfulSteps, 0u);
  BOOST_CHECK_EQUAL_COLLECTIONS(state.pVector.begin(), state.pVector.end(),
                                parameters.begin(), parameters.end());
}

BOOST_AUTO_TEST_CASE(StraightLineModeSkipsFieldQueries) {
  auto field = std::make_shared<RecordingField>();
  field->failAt = 1;
  Stepper stepper(field);
  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;
  options.useMagneticField = false;
  auto state = stepper.makeState(options);
  const auto initial = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);
  stepper.initialize(state, initial);
  BOOST_REQUIRE(stepper.prepareCurvilinearState(state));
  const auto result = stepper.step(state, Direction::Backward(), nullptr);
  BOOST_REQUIRE(result.ok());
  CHECK_CLOSE_ABS(stepper.position(state), Vector3(pos - stepSize * unitDir),
                  eps);
  CHECK_CLOSE_ABS(stepper.direction(state), unitDir, eps);
  BOOST_CHECK(field->lookups.empty());
  BOOST_CHECK_EQUAL(state.statistics.nAttemptedSteps, 1u);
  BOOST_CHECK_EQUAL(state.statistics.nSuccessfulSteps, 1u);
  BOOST_CHECK_EQUAL(state.statistics.nReverseSteps, 0u);
  BOOST_CHECK_EQUAL(state.statistics.pathLength, -stepSize);
  const auto jac = stepper.transportToCurvilinear(state);
  CHECK_CLOSE_ABS(jac(eBoundPhi, eBoundQOverP), 0., eps);
  CHECK_CLOSE_ABS(jac(eBoundTheta, eBoundQOverP), 0., eps);
  CHECK_CLOSE_ABS(jac(eBoundLoc0, eBoundPhi),
                  -stepSize * unitDir.head<2>().norm(), eps);
}

BOOST_AUTO_TEST_CASE(HelixShortcutMatchesConstantFieldAndRefreshes) {
  Stepper stepper(magneticField);
  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;
  auto normal = stepper.makeState(options);
  options.helixStep = 1_m;
  auto helix = stepper.makeState(options);
  const auto initial = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);
  stepper.initialize(normal, initial);
  stepper.initialize(helix, initial);
  for (int i = 0; i < 5; ++i) {
    BOOST_REQUIRE(stepper.step(normal, Direction::Forward(), nullptr).ok());
    BOOST_REQUIRE(stepper.step(helix, Direction::Forward(), nullptr).ok());
  }
  for (std::size_t i = 0; i < normal.pVector.size(); ++i) {
    CHECK_CLOSE_ABS(normal.pVector[i], helix.pVector[i], eps);
  }

  // Even consecutive short steps must sample at their own starting position.
  auto field = std::make_shared<RecordingField>();
  Stepper varying(field);
  auto state = varying.makeState(options);
  varying.initialize(state, initial);
  BOOST_REQUIRE(varying.step(state, Direction::Forward(), nullptr).ok());
  BOOST_REQUIRE_EQUAL(field->lookups.size(), 1u);
  const Vector3 nextStart = varying.position(state);
  BOOST_REQUIRE(varying.step(state, Direction::Forward(), nullptr).ok());
  BOOST_REQUIRE_EQUAL(field->lookups.size(), 2u);
  CHECK_CLOSE_ABS(field->lookups.back(), nextStart, eps);
}

BOOST_AUTO_TEST_CASE(ZeroFieldGradientMatchesFastJacobian) {
  Stepper::Config config{magneticField, {}};
  config.fieldGradient = [](const Vector3&, MagneticFieldProvider::Cache&) {
    return Result<SquareMatrix3>::success(SquareMatrix3::Zero());
  };
  Stepper stepper(config);
  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;
  auto fast = stepper.makeState(options);
  options.useFieldGradient = true;
  auto gradient = stepper.makeState(options);
  const auto initial = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);
  stepper.initialize(fast, initial);
  stepper.initialize(gradient, initial);
  for (int i = 0; i < 5; ++i) {
    BOOST_REQUIRE(stepper.step(fast, Direction::Forward(), nullptr).ok());
    BOOST_REQUIRE(stepper.step(gradient, Direction::Forward(), nullptr).ok());
  }
  for (std::size_t i = 0; i < fast.pVector.size(); ++i) {
    CHECK_CLOSE_ABS(fast.pVector[i], gradient.pVector[i], 1e-11);
  }
}

BOOST_AUTO_TEST_CASE(FieldGradientJacobianMatchesFiniteDifferences) {
  auto field = std::make_shared<RecordingField>();
  Stepper::Config config{field, {}};
  std::size_t gradientCalls = 0;
  config.fieldGradient = [&](const Vector3&, MagneticFieldProvider::Cache&) {
    ++gradientCalls;
    return Result<SquareMatrix3>::success(field->gradient);
  };
  Stepper analytic(config);
  Stepper finiteDifference(field);
  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = 5_mm;
  options.useFieldGradient = true;
  // Gradient transport takes precedence over the constant-field shortcut.
  options.helixStep = 1_m;
  const auto initial = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);
  const auto& surface = initial.referenceSurface();
  auto propagate = [&](const Stepper& stepper, const BoundVector& parameters) {
    auto state = stepper.makeState(options);
    stepper.initialize(state, parameters, cov, particleHypothesis, surface);
    for (int i = 0; i < 5; ++i) {
      const auto result = stepper.step(state, Direction::Forward(), nullptr);
      BOOST_REQUIRE(result.ok());
      BOOST_REQUIRE_EQUAL(*result, options.maxStepSize);
    }
    return state;
  };
  auto state = propagate(analytic, initial.parameters());
  BOOST_CHECK_EQUAL(gradientCalls, 15u);
  auto numericalGradient = propagate(finiteDifference, initial.parameters());
  for (std::size_t i = 0; i < state.pVector.size(); ++i) {
    CHECK_CLOSE_ABS(state.pVector[i], numericalGradient.pVector[i], 1e-10);
  }
  BoundToFreeMatrix jac =
      Eigen::Map<BoundToFreeMatrix>(state.pVector.data() + eFreeSize);
  jac.block<3, 1>(eFreePos0, eBoundQOverP) /= analytic.qOverP(state);
  jac.block<3, 1>(eFreeDir0, eBoundQOverP) /= analytic.qOverP(state);
  BoundToFreeMatrix numericalJac;
  for (std::uint32_t column = 0; column < eBoundSize; ++column) {
    const double delta = column < 2 ? 1e-3 : 1e-5;
    BoundVector plus = initial.parameters();
    BoundVector minus = initial.parameters();
    plus[column] += delta;
    minus[column] -= delta;
    auto plusState = propagate(analytic, plus);
    auto minusState = propagate(analytic, minus);
    numericalJac.col(column) =
        (Eigen::Map<const FreeVector>(plusState.pVector.data()) -
         Eigen::Map<const FreeVector>(minusState.pVector.data())) /
        (2. * delta);
  }
  CHECK_CLOSE_ABS(jac, numericalJac, 1e-7);
  // Check that the test resolves the terms omitted by the fast path.
  options.useFieldGradient = false;
  options.helixStep = 0.;
  auto withoutGradient = propagate(analytic, initial.parameters());
  const BoundToFreeMatrix fastJac =
      Eigen::Map<BoundToFreeMatrix>(withoutGradient.pVector.data() + eFreeSize);
  BOOST_CHECK_GT((fastJac.leftCols<2>() - numericalJac.leftCols<2>()).norm(),
                 1e-5);
}

BOOST_AUTO_TEST_CASE(GradientFailureLeavesParametersUnchanged) {
  auto field = std::make_shared<RecordingField>();
  Stepper::Config config{field, {}};
  config.fieldGradient = [](const Vector3&, MagneticFieldProvider::Cache&) {
    return Result<SquareMatrix3>::failure(MagneticFieldError::OutOfBounds);
  };
  Stepper stepper(config);
  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = stepSize;
  options.useFieldGradient = true;
  auto state = stepper.makeState(options);
  const auto initial = BoundTrackParameters::createCurvilinear(
      pos4, unitDir, charge / absMom, cov, particleHypothesis);
  stepper.initialize(state, initial);
  const auto parameters = state.pVector;
  const auto result = stepper.step(state, Direction::Forward(), nullptr);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() == MagneticFieldError::OutOfBounds);
  BOOST_CHECK_EQUAL_COLLECTIONS(state.pVector.begin(), state.pVector.end(),
                                parameters.begin(), parameters.end());
  BOOST_CHECK_EQUAL(state.statistics.nSuccessfulSteps, 0u);

  Stepper fallback(field);
  state.options.fieldGradientStep = 0.;
  const auto invalid = fallback.step(state, Direction::Forward(), nullptr);
  BOOST_REQUIRE(!invalid.ok());
  BOOST_CHECK(invalid.error() == EigenStepperError::StepInvalid);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

namespace bdata = boost::unit_test::data;

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;

Acts::GeometryContext gctx;
Acts::MagneticFieldContext mctx;

auto makeEigenPropagator(const Vector3& bField) {
  return Propagator{EigenStepper<>{std::make_shared<ConstantBField>(bField)},
                    VoidNavigator{}};
}

auto makeStraightLinePropagator() {
  return Propagator{StraightLineStepper{}, VoidNavigator{}};
}

auto makeDist(double a, double b, int seed) {
  return bdata::random(
      (bdata::engine = std::mt19937{}, bdata::seed = seed,
       bdata::distribution = std::uniform_real_distribution<double>(a, b)));
}

const auto bFieldDist = makeDist(0, 3_T, 11);
const auto directionDist = makeDist(0, 1, 12);
const auto distanceDist = makeDist(1_mm, 10_mm, 13);
const auto momentumDist = makeDist(0.1_GeV, 10_GeV, 14);
const auto chargeDist = makeDist(-1, 1, 15);

// Test the invariance of the angle covariance with and without a magnetic
// field.
//
// The test is performed by propagating a track with a given direction and
// momentum through a magnetic field and comparing the covariance of the
// direction at the end of the propagation with the covariance of the direction
// at the end of a straight line propagation.
BOOST_DATA_TEST_CASE(angle_cov_field_invariance,
                     (bFieldDist ^ bFieldDist ^ bFieldDist ^ directionDist ^
                      directionDist ^ directionDist ^ distanceDist ^
                      momentumDist ^ chargeDist) ^
                         bdata::xrange(100),
                     Bx, By, Bz, Dx, Dy, Dz, s, p, q_, index) {
  (void)index;
  const Vector3 bField{Bx, By, Bz};
  const Vector3 direction = Vector3{Dx, Dy, Dz}.normalized();
  const ParticleHypothesis particle = ParticleHypothesis::pion();
  const double q = q_ > 0 ? 1 : -1;

  const double qOverP = particle.qOverP(p, q);

  // Initial track parameters
  auto startSurface =
      Surface::makeShared<PlaneSurface>(Vector3{0, 0, 0}, direction);
  auto startParameters = BoundTrackParameters(
      startSurface,
      (BoundVector() << 0_mm, 0_mm, VectorHelpers::phi(direction),
       VectorHelpers::theta(direction), qOverP, 0_ns)
          .finished(),
      (BoundVector() << 1_um, 1_um, 1_mrad, 1_mrad, 1 / 1_MeV, 1_ps)
          .finished()
          .array()
          .square()
          .matrix()
          .asDiagonal(),
      particle);

  PropagatorOptions<> options(gctx, mctx);
  options.direction = Direction::Forward;
  options.pathLimit = s;

  // Reference propagation
  auto eigenPropagator = makeEigenPropagator(bField);
  auto eigenRes = eigenPropagator.propagate(startParameters, options);
  BOOST_CHECK(eigenRes.ok());
  BOOST_CHECK(eigenRes->endParameters.has_value());
  BOOST_CHECK(eigenRes->endParameters.value().covariance().has_value());
  auto eigenParameters = eigenRes->endParameters.value();

  PropagatorOptions<> straightLineBackwardOptions(gctx, mctx);
  options.direction = Direction::Backward;
  options.pathLimit = -s;

  // Now we propagate backwards with a straight line propagator so we can later
  // propagate forwards again to compare the covariances.
  auto straightLinePropagator = makeStraightLinePropagator();
  auto straightLineBackwardRes = straightLinePropagator.propagate(
      eigenParameters, straightLineBackwardOptions);
  BOOST_CHECK(straightLineBackwardRes.ok());
  BOOST_CHECK(straightLineBackwardRes->endParameters.has_value());
  BOOST_CHECK(
      straightLineBackwardRes->endParameters.value().covariance().has_value());
  auto straightLineBackwardParameters =
      straightLineBackwardRes->endParameters.value();

  // Propagate forward again
  auto straightLineForwardRes =
      straightLinePropagator.propagate(eigenParameters, options);
  BOOST_CHECK(straightLineForwardRes.ok());
  BOOST_CHECK(straightLineForwardRes->endParameters.has_value());
  BOOST_CHECK(
      straightLineForwardRes->endParameters.value().covariance().has_value());
  auto straightLineForwardParameters =
      straightLineForwardRes->endParameters.value();

  // Check if we get the same direction parameters and covariance

  {
    Vector2 exp = eigenParameters.parameters().segment<2>(eBoundPhi);
    Vector2 obs =
        straightLineForwardParameters.parameters().segment<2>(eBoundPhi);
    CHECK_CLOSE_ABS(exp, obs, 1e-7);
  }

  {
    SquareMatrix2 exp =
        eigenParameters.covariance().value().block<2, 2>(eBoundPhi, eBoundPhi);
    SquareMatrix2 obs =
        straightLineForwardParameters.covariance().value().block<2, 2>(
            eBoundPhi, eBoundPhi);
    CHECK_CLOSE_ABS(exp, obs, 1e-7);
  }
}

}  // namespace Test
}  // namespace Acts

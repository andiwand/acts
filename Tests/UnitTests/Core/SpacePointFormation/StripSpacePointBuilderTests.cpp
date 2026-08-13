// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/SpacePointFormation/SpacePointFormationError.hpp"
#include "Acts/SpacePointFormation/StripSpacePointBuilder.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <array>
#include <cmath>
#include <memory>
#include <numbers>

using namespace Acts;
using namespace Acts::StripSpacePointBuilder;
using namespace Acts::UnitLiterals;

namespace {

/// Stereo angle between the two strip layers, as in a barrel strip module.
constexpr double stereo = 0.04;

/// Inner strip: at x = 100, along z, centred on the x axis.
StripEnds innerStrip(double halfLength) {
  return StripEnds{Vector3(100, 0, halfLength), Vector3(100, 0, -halfLength)};
}

/// Outer strip: at x = 110, tilted by the stereo angle, centred on
/// (110, y, z). With a vertex at the origin this gives
///   m = (10 / (11 * halfLengthInner)) * (z - y / tan(stereo))
///   n = -y / (halfLengthOuter * sin(stereo))
StripEnds outerStrip(double y, double z, double halfLength) {
  const Vector3 centre(110, y, z);
  const Vector3 dir(0, std::sin(stereo), std::cos(stereo));
  return StripEnds{centre + halfLength * dir, centre - halfLength * dir};
}

/// Place the outer strip so that the space point parameters take the requested
/// values. Inverts the two relations above.
StripEnds outerStripFor(double m, double n, double halfLengthInner,
                        double halfLengthOuter) {
  const double y = -n * halfLengthOuter * std::sin(stereo);
  const double z = m * (11. * halfLengthInner) / 10. + y / std::tan(stereo);
  return outerStrip(y, z, halfLengthOuter);
}

const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

// two different strip pitches, so that the two directions cannot be confused
constexpr double var1 = 4.7e-4;
constexpr double var2 = 9.0e-4;

/// A plane whose local axes are given explicitly, centred at @p center
std::shared_ptr<PlaneSurface> makePlane(const Vector3& center,
                                        const Vector3& loc0,
                                        const Vector3& loc1) {
  Transform3 transform = Transform3::Identity();
  transform.matrix().block<3, 1>(0, 0) = loc0.normalized();
  transform.matrix().block<3, 1>(0, 1) = loc1.normalized();
  transform.matrix().block<3, 1>(0, 2) = loc0.cross(loc1).normalized();
  transform.matrix().block<3, 1>(0, 3) = center;
  return Surface::makeShared<PlaneSurface>(
      transform, std::make_shared<RectangleBounds>(50_mm, 50_mm));
}

/// The two local variances of a strip crossing, from the information matrix of
/// the two measurements built and inverted here
Vector2 referenceVariances(const double theta) {
  const Vector2 n1(1, 0);
  const Vector2 n2(std::cos(theta), std::sin(theta));
  const SquareMatrix2 information =
      n1 * n1.transpose() / var1 + n2 * n2.transpose() / var2;
  return information.inverse().diagonal();
}

/// A barrel module, whose only direction reaching `varianceZ` is the one
/// along z: pointing loc0 and then loc1 at it reads the two out one at a time
Vector2 varianceZR(const bool precisionAlongZ, const double theta) {
  const double phi = 0.4;
  const double r = 400_mm;
  const Vector3 position(r * std::cos(phi), r * std::sin(phi), 300_mm);
  const Vector3 rPhi(-std::sin(phi), std::cos(phi), 0);
  const Vector3 z = Vector3::UnitZ();

  auto surface = precisionAlongZ ? makePlane(position, z, rPhi)
                                 : makePlane(position, rPhi, z);
  return StripSpacePointBuilder::computeVarianceZR(gctx, *surface, position,
                                                   var1, var2, theta);
}

// orthogonal strips and the ITk strip stereo angle
constexpr std::array<double, 3> testAngles{std::numbers::pi / 2, 0.4, 40e-3};

}  // namespace

BOOST_AUTO_TEST_SUITE(StripSpacePointBuilderTests)

// A track from the origin through the centre of both strips puts the space
// point at the centre of the inner strip.
BOOST_AUTO_TEST_CASE(ConstrainedCentralHit) {
  const StripEnds first = innerStrip(30);
  const StripEnds second = outerStripFor(0., 0., 30, 30);

  ConstrainedOptions options;
  const Result<Vector3> sp =
      computeConstrainedSpacePoint(first, second, options);

  BOOST_REQUIRE(sp.ok());
  CHECK_CLOSE_ABS(sp->x(), 100., 1e-6);
  CHECK_CLOSE_ABS(sp->z(), 0., 1e-6);
}

// Well inside both strips: accepted, and the space point sits at m along the
// inner strip.
BOOST_AUTO_TEST_CASE(ConstrainedInsideBothStrips) {
  const StripEnds first = innerStrip(30);
  const StripEnds second = outerStripFor(0.5, -0.3, 30, 30);

  ConstrainedOptions options;
  const Result<Vector3> sp =
      computeConstrainedSpacePoint(first, second, options);

  BOOST_REQUIRE(sp.ok());
  CHECK_CLOSE_ABS(sp->z(), 0.5 * 30., 1e-6);
}

// Far outside stays rejected however generous the gap tolerance is.
BOOST_AUTO_TEST_CASE(ConstrainedFarOutsideIsRejected) {
  const StripEnds first = innerStrip(30);
  const StripEnds second = outerStripFor(10., 0.2, 30, 30);

  ConstrainedOptions options;
  options.stripLengthGapTolerance = 5.;

  const Result<Vector3> sp =
      computeConstrainedSpacePoint(first, second, options);

  BOOST_CHECK(!sp.ok());
  BOOST_CHECK_EQUAL(sp.error(), SpacePointFormationError::OutsideRelaxedLimits);
}

// Just beyond the end of the inner strip, with the outer one well inside.
BOOST_AUTO_TEST_CASE(ConstrainedRecoversOneSidedOvershoot) {
  const StripEnds first = innerStrip(30);
  const StripEnds second = outerStripFor(1.05, -0.5, 30, 30);

  ConstrainedOptions options;
  options.stripLengthTolerance = 0.01;

  // No gap tolerance: correctly outside the limits.
  options.stripLengthGapTolerance = 0.;
  const Result<Vector3> strict =
      computeConstrainedSpacePoint(first, second, options);
  BOOST_CHECK(!strict.ok());
  BOOST_CHECK_EQUAL(strict.error(),
                    SpacePointFormationError::OutsideRelaxedLimits);

  // With a gap tolerance the pair has to be recovered onto the strip end.
  options.stripLengthGapTolerance = 5.;
  const Result<Vector3> relaxed =
      computeConstrainedSpacePoint(first, second, options);
  BOOST_REQUIRE(relaxed.ok());
  CHECK_CLOSE_ABS(relaxed->z(), 30., 1e-6);
}

// The gap tolerance is a length in its own right: a crossing exactly that far
// beyond the end of a strip is the last one recovered. `m` reaches one at the
// end of the strip rather than at the far end of it, so reading the tolerance
// against the whole bottom-to-top length would recover only half of what it
// asks for -- and half a tolerance derived from a momentum is a detector that
// silently loses its soft tracks.
BOOST_AUTO_TEST_CASE(GapToleranceIsALengthPastTheStripEnd) {
  constexpr double halfLength = 30.;
  constexpr double overshoot = 3.;

  const StripEnds first = innerStrip(halfLength);
  const StripEnds second =
      outerStripFor(1. + overshoot / halfLength, -0.5, halfLength, halfLength);

  ConstrainedOptions options;
  options.stripLengthTolerance = 0.;

  options.stripLengthGapTolerance = overshoot * 1.001;
  BOOST_CHECK(computeConstrainedSpacePoint(first, second, options).ok());

  options.stripLengthGapTolerance = overshoot * 0.999;
  const Result<Vector3> tight =
      computeConstrainedSpacePoint(first, second, options);
  BOOST_CHECK(!tight.ok());
  BOOST_CHECK_EQUAL(tight.error(),
                    SpacePointFormationError::OutsideRelaxedLimits);
}

// The gap tolerance is a length, so a short strip gets a proportionally larger
// tolerance on its parameter than a long one.
BOOST_AUTO_TEST_CASE(GapToleranceUsesEachStripLength) {
  const StripEnds first = innerStrip(50);
  const StripEnds second = outerStripFor(0., -1.05, 50, 5);

  ConstrainedOptions options;
  options.stripLengthTolerance = 0.01;
  options.stripLengthGapTolerance = 1.;

  // 1 mm is 1% of the 100 mm inner strip but 10% of the 10 mm outer strip
  const Result<Vector3> sp =
      computeConstrainedSpacePoint(first, second, options);
  BOOST_REQUIRE(sp.ok());

  // Too small for either strip: rejected.
  options.stripLengthGapTolerance = 0.1;
  const Result<Vector3> tight =
      computeConstrainedSpacePoint(first, second, options);
  BOOST_CHECK(!tight.ok());
  BOOST_CHECK_EQUAL(tight.error(),
                    SpacePointFormationError::OutsideRelaxedLimits);
}

// After the shift both parameters still have to be on their strips.
BOOST_AUTO_TEST_CASE(ConstrainedRecoveryStillChecksLimits) {
  const StripEnds first = innerStrip(30);
  const StripEnds second = outerStripFor(1.05, -1.05, 30, 30);

  ConstrainedOptions options;
  options.stripLengthTolerance = 0.01;
  options.stripLengthGapTolerance = 5.;

  const Result<Vector3> sp =
      computeConstrainedSpacePoint(first, second, options);
  BOOST_CHECK(!sp.ok());
  BOOST_CHECK_EQUAL(sp.error(), SpacePointFormationError::OutsideLimits);
}

BOOST_AUTO_TEST_CASE(PrecisionDirection) {
  for (const double theta : testAngles) {
    BOOST_TEST_CONTEXT("theta = " << theta) {
      BOOST_CHECK_CLOSE(varianceZR(true, theta)[0],
                        referenceVariances(theta)[0], 1e-6);
      BOOST_CHECK_CLOSE(varianceZR(true, theta)[0], var1, 1e-6);
    }
  }
}

/// Along the strips the crossing is located to 1 / sin(theta) of a pitch, so
/// this direction degrades as the two become parallel
BOOST_AUTO_TEST_CASE(AlongStripDirection) {
  for (const double theta : testAngles) {
    BOOST_TEST_CONTEXT("theta = " << theta) {
      BOOST_CHECK_CLOSE(varianceZR(false, theta)[0],
                        referenceVariances(theta)[1], 1e-6);
    }
  }
}

/// Orthogonal strips measure the two local directions independently
BOOST_AUTO_TEST_CASE(OrthogonalStrips) {
  const double theta = std::numbers::pi / 2;
  BOOST_CHECK_CLOSE(varianceZR(true, theta)[0], var1, 1e-6);
  BOOST_CHECK_CLOSE(varianceZR(false, theta)[0], var2, 1e-6);
}

/// A barrel module carries no radial information either way
BOOST_AUTO_TEST_CASE(NoRadialVariance) {
  BOOST_CHECK_SMALL(varianceZR(true, 40e-3)[1], 1e-12);
  BOOST_CHECK_SMALL(varianceZR(false, 40e-3)[1], 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()

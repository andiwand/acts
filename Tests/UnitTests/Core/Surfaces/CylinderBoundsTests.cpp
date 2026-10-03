// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(SurfacesSuite)
/// Unit test for creating compliant/non-compliant CylinderBounds object

BOOST_AUTO_TEST_CASE(CylinderBoundsConstruction) {
  /// Test default construction
  // default construction is deleted

  const double radius = 0.5;
  const double halfZ = 10.;
  const double halfPhi = std::numbers::pi / 2.;
  const double averagePhi = std::numbers::pi / 2.;

  BOOST_CHECK_EQUAL(CylinderBounds(radius, halfZ).type(),
                    SurfaceBounds::eCylinder);
  BOOST_CHECK_EQUAL(CylinderBounds(radius, halfZ, halfPhi).type(),
                    SurfaceBounds::eCylinder);
  BOOST_CHECK_EQUAL(CylinderBounds(radius, halfZ, halfPhi, averagePhi).type(),
                    SurfaceBounds::eCylinder);

  /// Test copy construction;
  CylinderBounds cylinderBounds(radius, halfZ);
  CylinderBounds copyConstructedCylinderBounds(cylinderBounds);
  BOOST_CHECK_EQUAL(copyConstructedCylinderBounds, cylinderBounds);
}

BOOST_AUTO_TEST_CASE(CylinderBoundsRecreation) {
  const double radius = 0.5;
  const double halfZ = 10.;

  // Test construction with radii and default sector
  auto original = CylinderBounds(radius, halfZ);
  auto valvector = original.values();
  std::array<double, CylinderBounds::eSize> values{};
  std::copy_n(valvector.begin(), CylinderBounds::eSize, values.begin());
  CylinderBounds recreated(values);
  BOOST_CHECK_EQUAL(original, recreated);
}

BOOST_AUTO_TEST_CASE(CylinderBoundsException) {
  const double radius = 0.5;
  const double halfZ = 10.;
  const double halfPhi = std::numbers::pi / 2.;
  const double averagePhi = std::numbers::pi / 2.;

  /// Negative radius
  BOOST_CHECK_THROW(CylinderBounds(-radius, halfZ, halfPhi, averagePhi),
                    std::logic_error);

  /// Negative half length in z
  BOOST_CHECK_THROW(CylinderBounds(radius, -halfZ, halfPhi, averagePhi),
                    std::logic_error);

  /// Negative half sector in phi
  BOOST_CHECK_THROW(CylinderBounds(radius, halfZ, -halfPhi, averagePhi),
                    std::logic_error);

  /// Half sector in phi out of bounds
  BOOST_CHECK_THROW(CylinderBounds(radius, halfZ, 4., averagePhi),
                    std::logic_error);

  /// Phi position out of bounds
  BOOST_CHECK_THROW(CylinderBounds(radius, halfZ, halfPhi, 4.),
                    std::logic_error);
}

/// Unit tests for CylinderBounds properties
BOOST_AUTO_TEST_CASE(CylinderBoundsProperties) {
  // CylinderBounds object of radius 0.5 and halfZ 20
  const double radius = 0.5;
  const double halfZ = 20.;                      // != 10.
  const double halfPhi = std::numbers::pi / 4.;  // != pi/2
  const double averagePhi = 0.;                  // != pi/2

  CylinderBounds cylinderBoundsObject(radius, halfZ);
  CylinderBounds cylinderBoundsSegment(radius, halfZ, halfPhi, averagePhi);

  /// Test for type()
  BOOST_CHECK_EQUAL(cylinderBoundsObject.type(), SurfaceBounds::eCylinder);

  /// Test for inside(), 2D coords are r or phi ,z? : needs clarification
  const Vector2 origin{0., 0.};
  const Vector2 atPiBy2{std::numbers::pi / 2., 0.};
  const Vector2 atPi{std::numbers::pi, 0.};
  const Vector2 beyondEnd{0, 30.};
  const Vector2 unitZ{0., 1.};
  const Vector2 unitPhi{1., 0.};
  const BoundaryTolerance tolerance = BoundaryTolerance::AbsoluteEuclidean(0.1);

  BOOST_CHECK(cylinderBoundsObject.inside(atPiBy2, tolerance));
  BOOST_CHECK(!cylinderBoundsSegment.inside(unitPhi, tolerance));
  BOOST_CHECK(cylinderBoundsObject.inside(origin, tolerance));

  /// Test for r()
  CHECK_CLOSE_REL(cylinderBoundsObject.get(CylinderBounds::eR), radius, 1e-6);

  /// Test for averagePhi
  CHECK_CLOSE_OR_SMALL(cylinderBoundsObject.get(CylinderBounds::eAveragePhi),
                       averagePhi, 1e-6, 1e-6);

  /// Test for halfPhiSector
  CHECK_CLOSE_REL(cylinderBoundsSegment.get(CylinderBounds::eHalfPhiSector),
                  halfPhi,
                  1e-6);  // fail

  /// Test for halflengthZ (NOTE: Naming violation)
  CHECK_CLOSE_REL(cylinderBoundsObject.get(CylinderBounds::eHalfLengthZ), halfZ,
                  1e-6);

  /// Test for dump
  boost::test_tools::output_test_stream dumpOutput;
  cylinderBoundsObject.toStream(dumpOutput);
  BOOST_CHECK(dumpOutput.is_equal(
      "Acts::CylinderBounds: (radius, halfLengthZ, halfPhiSector, "
      "averagePhi) = (0.5000000, 20.0000000, 3.1415927, 0.0000000)"));
}

/// Unit test for testing CylinderBounds assignment
BOOST_AUTO_TEST_CASE(CylinderBoundsAssignment) {
  const double radius = 0.5;
  const double halfZ = 20.;  // != 10.

  CylinderBounds cylinderBoundsObject(radius, halfZ);
  CylinderBounds assignedCylinderBounds(10.5, 6.6);
  assignedCylinderBounds = cylinderBoundsObject;

  BOOST_CHECK_EQUAL(assignedCylinderBounds.get(CylinderBounds::eR),
                    cylinderBoundsObject.get(CylinderBounds::eR));
  BOOST_CHECK_EQUAL(assignedCylinderBounds, cylinderBoundsObject);
}

BOOST_AUTO_TEST_CASE(CylinderBoundsCenter) {
  const double radius = 5.0;
  const double halfZ = 10.0;

  // Test full cylinder
  CylinderBounds fullCylinder(radius, halfZ);
  Vector2 center = fullCylinder.center();
  CHECK_CLOSE_ABS(center, Vector2(0., 0.), 1e-6);

  // Test cylinder with average phi offset
  const double averagePhi = std::numbers::pi / 4.;
  CylinderBounds offsetCylinder(radius, halfZ, std::numbers::pi, averagePhi);
  Vector2 centerOffset = offsetCylinder.center();
  CHECK_CLOSE_ABS(centerOffset, Vector2(radius * averagePhi, 0.), 1e-6);
}

BOOST_AUTO_TEST_CASE(CylinderBoundsCoversFullAzimuth) {
  const double radius = 0.5;
  const double halfZ = 10.;

  BOOST_CHECK(CylinderBounds(radius, halfZ).coversFullAzimuth());
  BOOST_CHECK(
      CylinderBounds(radius, halfZ, std::numbers::pi).coversFullAzimuth());
  BOOST_CHECK(
      CylinderBounds(radius, halfZ, std::nextafter(std::numbers::pi, 0.))
          .coversFullAzimuth());
  BOOST_CHECK(!CylinderBounds(radius, halfZ, std::numbers::pi / 4.)
                   .coversFullAzimuth());
}

BOOST_AUTO_TEST_CASE(CylinderSectorBoundaryCoordinates) {
  const CylinderBounds bounds(10., 100., 0.5, 0.2);
  const Vector2 query(8., 0.);
  CHECK_CLOSE_ABS(bounds.closestPoint(query, SquareMatrix2::Identity()),
                  Vector2(7., 0.), 1e-12);
  CHECK_CLOSE_ABS(bounds.distance(query), 1., 1e-12);
  BOOST_CHECK(!bounds.inside(query, BoundaryTolerance::None()));
  BOOST_CHECK(bounds.inside(query, BoundaryTolerance::AbsoluteEuclidean(1.5)));
  BOOST_CHECK(!bounds.inside(query, BoundaryTolerance::AbsoluteEuclidean(0.5)));
  BOOST_CHECK(bounds.inside(Vector2(6., 0.),
                            BoundaryTolerance::AbsoluteEuclidean(-0.5)));
  BOOST_CHECK(!bounds.inside(Vector2(6., 0.),
                             BoundaryTolerance::AbsoluteEuclidean(-1.5)));
  BOOST_CHECK(bounds.inside(query, BoundaryTolerance::Infinite()));

  const double period = 20. * std::numbers::pi;
  for (int turn : {-3, 0, 2}) {
    const Vector2 shiftedQuery = query + Vector2(turn * period, 0.);
    CHECK_CLOSE_ABS(
        bounds.closestPoint(shiftedQuery, SquareMatrix2::Identity()) -
            shiftedQuery,
        Vector2(-1., 0.), 1e-12);
  }

  // A sector spanning the principal phi seam has a nearby boundary on either
  // representation of the query.
  const CylinderBounds seam(10., 100., 0.5, 3.);
  const Vector2 seamQuery(10. * (3.6 - 2. * std::numbers::pi), 0.);
  CHECK_CLOSE_ABS(seam.distance(seamQuery), 1., 1e-12);
  BOOST_CHECK(
      seam.inside(seamQuery, BoundaryTolerance::AbsoluteEuclidean(1.5)));
}

BOOST_AUTO_TEST_CASE(ClosedCylinderHasNoPhiBoundary) {
  const CylinderBounds bounds(10., 20.);
  const Vector2 query(10. * std::numbers::pi, 0.);
  BOOST_CHECK(bounds.inside(query));
  CHECK_CLOSE_ABS(bounds.distance(query), 20., 1e-12);
  BOOST_CHECK(bounds.inside(query, BoundaryTolerance::AbsoluteEuclidean(-19.)));
  BOOST_CHECK(
      !bounds.inside(query, BoundaryTolerance::AbsoluteEuclidean(-21.)));

  SquareMatrix2 metric;
  metric << 4., 1., 1., 2.;
  const Vector2 outside(query.x(), 22.);
  CHECK_CLOSE_ABS(bounds.closestPoint(outside, metric) - outside,
                  Vector2(0.5, -2.), 1e-12);
  BOOST_CHECK(bounds.inside(outside, BoundaryTolerance::Chi2Bound(metric, 8.)));
  BOOST_CHECK(
      !bounds.inside(outside, BoundaryTolerance::Chi2Bound(metric, 6.)));
}

BOOST_AUTO_TEST_CASE(CylinderPeriodicMetricProjection) {
  const double radius = 2.;
  const double halfZ = 3.;
  const double halfRPhi = radius * 0.7;
  const double center = radius * 2.9;
  const double period = 2. * std::numbers::pi * radius;
  const CylinderBounds bounds(radius, halfZ, 0.7, 2.9);
  for (double correlation : {-9., 0., 9.}) {
    SquareMatrix2 metric;
    metric << 1., correlation, correlation, 100.;
    for (const Vector2& query : {Vector2(0., 0.), Vector2(5., 8.),
                                 Vector2(-10., -12.), Vector2(30., 1.)}) {
      // Independent reference: exhaustively project onto many unrolled copies
      // of the rectangle, including copies beyond the adjacent phi periods.
      double reference = std::numeric_limits<double>::infinity();
      for (int turn = -30; turn <= 30; ++turn) {
        const double x = center + turn * period;
        const Vector2 candidate =
            detail::VerticesHelper::computeClosestPointOnAlignedBox(
                Vector2(x - halfRPhi, -halfZ), Vector2(x + halfRPhi, halfZ),
                query, metric);
        const Vector2 delta = candidate - query;
        reference = std::min(reference, delta.dot(metric * delta));
      }
      const Vector2 delta = bounds.closestPoint(query, metric) - query;
      CHECK_CLOSE_ABS(delta.dot(metric * delta), reference, 1e-9);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

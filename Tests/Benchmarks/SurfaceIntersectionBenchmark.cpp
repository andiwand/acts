// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/LineSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceImpl.hpp"
#include "Acts/Tests/CommonHelpers/BenchmarkTools.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <cmath>
#include <random>

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts::Test {

// Some randomness & number crunching
unsigned int ntests = 10;
unsigned int nrepts = 2000;
const BoundaryTolerance boundaryTolerance = BoundaryTolerance::None();
const bool testPlane = true;
const bool testDisc = false;
const bool testCylinder = false;
const bool testStraw = false;

// Create a test context
GeometryContext tgContext = GeometryContext();

// Create a test plane in 10 m distance
// Some random transform
Transform3 at = Transform3::Identity() * Translation3(0_m, 0_m, 10_m) *
                AngleAxis3(0.15, Vector3(1.2, 1.2, 0.12).normalized());

// Define the Plane surface
auto rb = std::make_shared<RectangleBounds>(1_m, 1_m);
auto aPlane = Surface::makeShared<PlaneSurface>(at, std::move(rb));

// Define the Disc surface
auto db = std::make_shared<RadialBounds>(0.2_m, 1.2_m);
auto aDisc = Surface::makeShared<DiscSurface>(at, std::move(db));

// Define a Cylinder surface
auto cb = std::make_shared<CylinderBounds>(10_m, 100_m);
auto aCylinder = Surface::makeShared<CylinderSurface>(at, std::move(cb));

// Define a Straw surface
auto aStraw = Surface::makeShared<StrawSurface>(at, 50_cm, 2_m);

// The number of repeat surfaces
unsigned int nRepeatSurface = 10;

// The origin of our attempts for plane, disc and cylinder
Vector3 origin(0., 0., 0.);

// The origin for straw/line attempts
Vector3 originStraw(0.3_m, -0.2_m, 11_m);

template <typename surface_t>
MicroBenchmarkResult intersectionTest(const surface_t& surface,
                                      const Vector3& direction) {
  return Acts::Test::microBenchmark(
      [&]() -> SurfaceMultiIntersection {
        return surface.intersect(tgContext, origin, direction,
                                 boundaryTolerance);
      },
      nrepts);
}

MicroBenchmarkResult intersectionTest2(const Surface& surface,
                                       const Vector3& direction) {
  return Acts::Test::microBenchmark(
      [&]() -> SurfaceMultiIntersection {
        return surface.intersectImpl(tgContext, origin, direction,
                                     boundaryTolerance, s_onSurfaceTolerance);
      },
      nrepts);
}

MicroBenchmarkResult multipleIntersectionTest(
    const std::vector<const Surface*>& surface, const Vector3& direction) {
  std::vector<SurfaceMultiIntersection> intersections;

  return Acts::Test::microBenchmark(
      [&]() -> const std::vector<SurfaceMultiIntersection>& {
        intersections.clear();
        for (const Surface* s : surface) {
          auto intersection =
              s->intersect(tgContext, origin, direction, boundaryTolerance);
          intersections.push_back(intersection);
        }
        return intersections;
      },
      nrepts / 10);
}

MicroBenchmarkResult multipleIntersectionTest2(
    const std::vector<const Surface*>& surface, const Vector3& direction) {
  std::vector<SurfaceMultiIntersection> intersections;

  return Acts::Test::microBenchmark(
      [&]() -> const std::vector<SurfaceMultiIntersection>& {
        intersections.clear();
        for (const Surface* s : surface) {
          auto intersection =
              s->intersectImpl(tgContext, origin, direction, boundaryTolerance,
                               s_onSurfaceTolerance);
          intersections.push_back(intersection);
        }
        return intersections;
      },
      nrepts / 10);
}

BOOST_DATA_TEST_CASE(
    benchmark_surface_intersections,
    bdata::random((bdata::engine = std::mt19937(), bdata::seed = 21,
                   bdata::distribution =
                       std::uniform_real_distribution<double>(-M_PI, M_PI))) ^
        bdata::random((bdata::engine = std::mt19937(), bdata::seed = 22,
                       bdata::distribution =
                           std::uniform_real_distribution<double>(-0.3, 0.3))) ^
        bdata::xrange(ntests),
    phi, theta, index) {
  (void)index;

  Vector3 direction = makeDirectionFromPhiTheta(phi, theta);

  std::cout << std::endl
            << "Single surface intersections benchmarking theta=" << theta
            << ", phi=" << phi << ", boundaryTolerance=" << boundaryTolerance
            << "..." << std::endl;
  if (testPlane) {
    std::cout << "- Plane: "
              << intersectionTest<PlaneSurface>(*aPlane, direction)
              << std::endl;
  }
  if (testDisc) {
    std::cout << "- Disc: " << intersectionTest<DiscSurface>(*aDisc, direction)
              << std::endl;
  }
  if (testCylinder) {
    std::cout << "- Cylinder: "
              << intersectionTest<CylinderSurface>(*aCylinder, direction)
              << std::endl;
  }
  if (testStraw) {
    std::cout << "- Straw: "
              << intersectionTest<StrawSurface>(*aStraw, direction)
              << std::endl;
  }
  std::cout << "Done." << std::endl;

  std::cout << std::endl
            << "Single virtual surface intersections benchmarking theta="
            << theta << ", phi=" << phi
            << ", boundaryTolerance=" << boundaryTolerance << "..."
            << std::endl;
  if (testPlane) {
    std::cout << "- Plane: "
              << intersectionTest(dynamic_cast<const Surface&>(*aPlane),
                                  direction)
              << std::endl;
  }
  if (testDisc) {
    std::cout << "- Disc: "
              << intersectionTest(dynamic_cast<const Surface&>(*aDisc),
                                  direction)
              << std::endl;
  }
  if (testCylinder) {
    std::cout << "- Cylinder: "
              << intersectionTest(dynamic_cast<const Surface&>(*aCylinder),
                                  direction)
              << std::endl;
  }
  if (testStraw) {
    std::cout << "- Straw: "
              << intersectionTest(dynamic_cast<const Surface&>(*aStraw),
                                  direction)
              << std::endl;
  }
  std::cout << "Done." << std::endl;

  std::cout << std::endl
            << "Single de-virtual surface intersections benchmarking theta="
            << theta << ", phi=" << phi
            << ", boundaryTolerance=" << boundaryTolerance << "..."
            << std::endl;
  if (testPlane) {
    std::cout << "- Plane: " << intersectionTest2(*aPlane, direction)
              << std::endl;
  }
  if (testDisc) {
    std::cout << "- Disc: " << intersectionTest2(*aDisc, direction)
              << std::endl;
  }
  if (testCylinder) {
    std::cout << "- Cylinder: " << intersectionTest2(*aCylinder, direction)
              << std::endl;
  }
  if (testStraw) {
    std::cout << "- Straw: " << intersectionTest2(*aStraw, direction)
              << std::endl;
  }
  std::cout << "Done." << std::endl;

  std::vector<const Surface*> surfaces;
  if (testPlane) {
    for (unsigned int i = 0; i < nRepeatSurface; i++) {
      surfaces.push_back(aPlane.get());
    }
  }
  if (testDisc) {
    for (unsigned int i = 0; i < nRepeatSurface; i++) {
      surfaces.push_back(aDisc.get());
    }
  }
  if (testCylinder) {
    for (unsigned int i = 0; i < nRepeatSurface; i++) {
      surfaces.push_back(aCylinder.get());
    }
  }
  if (testStraw) {
    for (unsigned int i = 0; i < nRepeatSurface; i++) {
      surfaces.push_back(aStraw.get());
    }
  }

  std::cout << std::endl
            << "Multiple virtual surface intersections benchmarking theta="
            << theta << ", phi=" << phi
            << ", boundaryTolerance=" << boundaryTolerance << "..."
            << std::endl;
  std::cout << "- All: " << multipleIntersectionTest(surfaces, direction)
            << std::endl;
  std::cout << "Done." << std::endl;

  std::cout << std::endl
            << "Multiple de-virtual surface intersections benchmarking theta="
            << theta << ", phi=" << phi
            << ", boundaryTolerance=" << boundaryTolerance << "..."
            << std::endl;
  std::cout << "- All: " << multipleIntersectionTest2(surfaces, direction)
            << std::endl;
  std::cout << "Done." << std::endl;
}

}  // namespace Acts::Test

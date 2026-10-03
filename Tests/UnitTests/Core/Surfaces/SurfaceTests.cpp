// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

#include <iomanip>
#include <memory>
#include <sstream>

#include "SurfaceStub.hpp"

namespace Acts {
/// Mock track object with minimal methods implemented for compilation
class MockTrack {
 public:
  MockTrack(const Vector3& mom, const Vector3& pos) : m_mom(mom), m_pos(pos) {
    // nop
  }

  Vector3 momentum() const { return m_mom; }

  Vector3 position() const { return m_pos; }

 private:
  Vector3 m_mom;
  Vector3 m_pos;
};
}  // namespace Acts

using namespace Acts;

namespace ActsTests {

// Create a test context
GeometryContext tgContext = GeometryContext::dangerouslyDefaultConstruct();

BOOST_AUTO_TEST_SUITE(SurfacesSuite)

BOOST_AUTO_TEST_CASE(SurfaceBoundaryHookOwnsOnSurfaceChecks) {
  struct BoundaryPlane : PlaneSurface {
    explicit BoundaryPlane(std::shared_ptr<const PlanarBounds> bounds)
        : GeometryObject(),
          PlaneSurface(Transform3::Identity(), std::move(bounds)) {}

    mutable unsigned int checks = 0;
    mutable Vector2 lastPosition = Vector2::Zero();

    bool insideBounds(const Vector2& position,
                      const BoundaryTolerance& tolerance) const override {
      ++checks;
      lastPosition = position;
      // Example surface-specific acceptance layered over existing bounds.
      return position.x() >= 0. && Surface::insideBounds(position, tolerance);
    }
  };
  const auto bounds = std::make_shared<const RectangleBounds>(5., 10.);
  const BoundaryPlane plane(bounds);
  const Surface& surface = plane;
  const RegularSurface& regular = plane;
  const Vector3 direction = Vector3::UnitZ();
  for (const double x : {-1., 1.}) {
    const Vector3 position{x, 2., 0.};
    const bool expected = x > 0.;
    BOOST_CHECK_EQUAL(surface.isOnSurface(tgContext, position, direction,
                                          BoundaryTolerance::None()),
                      expected);
    CHECK_CLOSE_ABS(plane.lastPosition, Vector2(x, 2.), 1e-12);
    BOOST_CHECK_EQUAL(
        regular.isOnSurface(tgContext, position, BoundaryTolerance::None()),
        expected);
    CHECK_CLOSE_ABS(plane.lastPosition, Vector2(x, 2.), 1e-12);
    BOOST_CHECK_EQUAL(surface
                          .intersect(tgContext, position - direction, direction,
                                     BoundaryTolerance::None())
                          .closest()
                          .isValid(),
                      expected);
  }
  BOOST_CHECK_EQUAL(plane.checks, 6u);
  // Points off the geometric surface must fail before region evaluation.
  BOOST_CHECK(!surface.isOnSurface(tgContext, {1., 2., 1.}, direction));
  BOOST_CHECK(!regular.isOnSurface(tgContext, {1., 2., 1.}));
  BOOST_CHECK_EQUAL(plane.checks, 6u);
  BOOST_CHECK_EQUAL(plane.boundsPtr().get(), bounds.get());
}

BOOST_AUTO_TEST_CASE(SurfaceBoundaryHookPreservesLegacyBoundsOverrides) {
  struct CustomBounds : RectangleBounds {
    CustomBounds() : RectangleBounds(5., 10.) {}
    using RectangleBounds::inside;
    mutable unsigned int checks = 0;
    bool inside(const Vector2&, const BoundaryTolerance&) const override {
      ++checks;
      return false;
    }
  };
  const auto bounds = std::make_shared<const CustomBounds>();
  const auto plane =
      Surface::makeShared<PlaneSurface>(Transform3::Identity(), bounds);
  const Surface& surface = *plane;
  const RegularSurface& regular = *plane;
  const Vector3 position{1., 2., 0.};
  BOOST_CHECK(bounds->inside(Vector2(1., 2.)));
  BOOST_CHECK(!surface.insideBounds({1., 2.}));
  BOOST_CHECK(!surface.isOnSurface(tgContext, position, Vector3::UnitZ()));
  BOOST_CHECK(!regular.isOnSurface(tgContext, position));
  BOOST_CHECK(!surface
                   .intersect(tgContext, {1., 2., -1.}, Vector3::UnitZ(),
                              BoundaryTolerance::None())
                   .closest()
                   .isValid());
  BOOST_CHECK_EQUAL(bounds->checks, 4u);
}

/// todo: make test fixture; separate out different cases

/// Unit test for creating compliant/non-compliant Surface object
BOOST_AUTO_TEST_CASE(SurfaceConstruction) {
  // SurfaceStub s;
  BOOST_CHECK_EQUAL(Surface::Other, SurfaceStub().type());
  SurfaceStub original;
  BOOST_CHECK_EQUAL(Surface::Other, SurfaceStub(original).type());
  Translation3 translation{0., 1., 2.};
  Transform3 transform(translation);
  BOOST_CHECK_EQUAL(Surface::Other,
                    SurfaceStub(tgContext, original, transform).type());
  // need some cruft to make the next one work
  auto pTransform = Transform3(translation);
  std::shared_ptr<const Acts::PlanarBounds> p =
      std::make_shared<const RectangleBounds>(5., 10.);
  DetectorElementStub detElement{pTransform, p, 0.2, nullptr};
  BOOST_CHECK_EQUAL(Surface::Other, SurfaceStub(detElement).type());
}

/// Unit test for testing Surface properties
BOOST_AUTO_TEST_CASE(SurfaceProperties) {
  // build a test object , 'surface'
  std::shared_ptr<const Acts::PlanarBounds> pPlanarBound =
      std::make_shared<const RectangleBounds>(5., 10.);
  Vector3 reference{0., 1., 2.};
  Translation3 translation{0., 1., 2.};
  auto pTransform = Transform3(translation);
  auto pLayer = PlaneLayer::create(pTransform, pPlanarBound, nullptr);
  auto pMaterial =
      std::make_shared<const HomogeneousSurfaceMaterial>(makePercentSlab());
  DetectorElementStub detElement{pTransform, pPlanarBound, 0.2, pMaterial};
  SurfaceStub surface(detElement);

  // associatedDetectorElement
  BOOST_CHECK_EQUAL(surface.surfacePlacement(), &detElement);

  // test associatelayer, associatedLayer
  surface.associateLayer(*pLayer);
  BOOST_CHECK_EQUAL(surface.associatedLayer(), pLayer.get());

  // associated Material is not set to the surface
  // it is set to the detector element surface though
  BOOST_CHECK_NE(surface.surfaceMaterial(), pMaterial.get());

  // center()
  CHECK_CLOSE_OR_SMALL(reference, surface.center(tgContext), 1e-6, 1e-9);

  // insideBounds
  Vector2 localPosition{0.1, 3.};
  BOOST_CHECK(surface.insideBounds(localPosition));
  Vector2 outside{20., 20.};
  BOOST_CHECK(surface.insideBounds(
      outside));  // should return false, but doesn't because SurfaceStub has
                  // "no bounds" hard-coded
  Vector3 mom{100., 200., 300.};

  // isOnSurface
  BOOST_CHECK(surface.isOnSurface(tgContext, reference, mom,
                                  BoundaryTolerance::Infinite()));
  BOOST_CHECK(surface.isOnSurface(
      tgContext, reference, mom,
      BoundaryTolerance::None()));  // need to improve bounds()

  // referenceFrame()
  RotationMatrix3 unitary;
  unitary << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  auto referenceFrame =
      surface.referenceFrame(tgContext, Vector3{1, 2, 3}.normalized(), mom);
  BOOST_CHECK_EQUAL(referenceFrame, unitary);

  // normal()
  auto normal = surface.normal(tgContext, Vector3{1, 2, 3}.normalized(),
                               Vector3::UnitZ());
  Vector3 zero{0., 0., 0.};
  BOOST_CHECK_EQUAL(zero, normal);

  // pathCorrection is pure virtual

  // surfaceMaterial()
  auto pNewMaterial =
      std::make_shared<const HomogeneousSurfaceMaterial>(makePercentSlab());
  surface.assignSurfaceMaterial(pNewMaterial);
  BOOST_CHECK_EQUAL(surface.surfaceMaterial(), pNewMaterial.get());

  CHECK_CLOSE_OR_SMALL(surface.localToGlobalTransform(tgContext), pTransform,
                       1e-6, 1e-9);

  // type() is pure virtual
}

BOOST_AUTO_TEST_CASE(SurfaceToStreamPreservesStreamState) {
  std::ostringstream stream;
  stream << std::scientific << std::showpos << std::setfill('#')
         << std::setprecision(3);
  stream.width(17);

  const auto flags = stream.flags();
  const auto precision = stream.precision();
  const auto width = stream.width();
  const auto fill = stream.fill();

  auto bounds = std::make_shared<const RectangleBounds>(5., 10.);
  auto surface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity(), bounds);

  stream << surface->toStream(tgContext);

  BOOST_CHECK(stream.flags() == flags);
  BOOST_CHECK_EQUAL(stream.precision(), precision);
  BOOST_CHECK_EQUAL(stream.width(), width);
  BOOST_CHECK_EQUAL(stream.fill(), fill);
}

BOOST_AUTO_TEST_CASE(EqualityOperators) {
  // build some test objects
  std::shared_ptr<const Acts::PlanarBounds> pPlanarBound =
      std::make_shared<const RectangleBounds>(5., 10.);
  Vector3 reference{0., 1., 2.};
  Translation3 translation1{0., 1., 2.};
  Translation3 translation2{1., 1., 2.};
  auto pTransform1 = Transform3(translation1);
  auto pTransform2 = Transform3(translation2);

  // build a planeSurface to be compared
  auto planeSurface =
      Surface::makeShared<PlaneSurface>(pTransform1, pPlanarBound);
  auto pLayer = PlaneLayer::create(pTransform1, pPlanarBound, nullptr);
  auto pMaterial =
      std::make_shared<const HomogeneousSurfaceMaterial>(makePercentSlab());
  DetectorElementStub detElement1{pTransform1, pPlanarBound, 0.2, pMaterial};
  DetectorElementStub detElement2{pTransform1, pPlanarBound, 0.3, pMaterial};
  DetectorElementStub detElement3{pTransform2, pPlanarBound, 0.3, pMaterial};

  SurfaceStub surface1(detElement1);
  SurfaceStub surface2(detElement1);  // 1 and 2 are the same
  SurfaceStub surface3(detElement2);  // 3 differs in thickness
  SurfaceStub surface4(detElement3);  // 4 has a different transform and id
  SurfaceStub surface5(detElement1);
  surface5.assignSurfaceMaterial(pMaterial);  // 5 has non-null surface material

  BOOST_CHECK(surface1 == surface2);
  BOOST_CHECK(surface1 != surface3);
  BOOST_CHECK(surface1 != surface4);
  BOOST_CHECK(surface1 != surface5);
  BOOST_CHECK(surface1 != *planeSurface);

  // Test the getSharedPtr
  const auto surfacePtr = Surface::makeShared<const SurfaceStub>(detElement1);
  const auto sharedSurfacePtr = surfacePtr->getSharedPtr();
  BOOST_CHECK(*surfacePtr == *sharedSurfacePtr);
}
BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

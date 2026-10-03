// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/BoundFactory.hpp"

#include <array>
#include <memory>
#include <span>

using namespace Acts;

namespace ActsTests {
BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(BoundFactoryReusesEquivalentConcreteBounds) {
  SurfaceBoundFactory factory;
  const auto first = factory.makeBounds<RectangleBounds>(2., 3.);
  const auto equal = factory.makeBounds<RectangleBounds>(2., 3.);
  const auto different = factory.makeBounds<RectangleBounds>(2., 4.);
  BOOST_REQUIRE(first);
  BOOST_CHECK_EQUAL(first, equal);
  BOOST_CHECK_NE(first, different);
  BOOST_CHECK_EQUAL(factory.size(), 2u);
  const std::shared_ptr<const PlanarBounds> base =
      std::make_shared<const RectangleBounds>(2., 3.);
  BOOST_CHECK_EQUAL(factory.insert(base), first);

  VolumeBoundFactory volumes;
  BOOST_CHECK_EQUAL(volumes.makeBounds<CuboidVolumeBounds>(1., 2., 3.),
                    volumes.makeBounds<CuboidVolumeBounds>(1., 2., 3.));
  BOOST_CHECK_EQUAL(volumes.size(), 1u);
}

BOOST_AUTO_TEST_CASE(BoundFactoryPreservesPolygonDynamicTypes) {
  const std::array<Vector2, 4> vertices{Vector2(0., 0.), Vector2(1., 0.),
                                        Vector2(1., 1.), Vector2(0., 1.)};
  for (bool dynamicFirst : {false, true}) {
    SurfaceBoundFactory factory;
    const auto fixed = std::make_shared<const ConvexPolygonBounds<4>>(vertices);
    const auto dynamic =
        std::make_shared<const ConvexPolygonBounds<PolygonDynamic>>(
            std::span<const Vector2>(vertices));
    if (dynamicFirst) {
      factory.insert(dynamic);
    } else {
      factory.insert(fixed);
    }
    const auto fixedResult = factory.insert(fixed);
    const auto dynamicResult = factory.insert(dynamic);
    BOOST_REQUIRE(fixedResult);
    BOOST_REQUIRE(dynamicResult);
    BOOST_CHECK_EQUAL(factory.size(), 2u);
    BOOST_CHECK_EQUAL(factory.makeBounds<ConvexPolygonBounds<4>>(vertices),
                      fixedResult);
    BOOST_CHECK_EQUAL(factory.makeBounds<ConvexPolygonBounds<PolygonDynamic>>(
                          std::span<const Vector2>(vertices)),
                      dynamicResult);
    const std::shared_ptr<const PlanarBounds> base = dynamic;
    BOOST_CHECK_EQUAL(factory.insert(base), dynamicResult);
  }
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests

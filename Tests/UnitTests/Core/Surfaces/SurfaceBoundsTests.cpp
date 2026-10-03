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
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Surfaces/detail/BoundaryCheck.hpp"
#include "Acts/Surfaces/detail/PolarChart.hpp"

#include <cstddef>
#include <iomanip>
#include <memory>
#include <numbers>
#include <numeric>
#include <ostream>
#include <sstream>
#include <vector>

namespace Acts {

/// Class to implement pure virtual method of SurfaceBounds for testing only
class SurfaceBoundsStub : public SurfaceBounds {
 public:
  /// Implement ctor and pure virtual methods of SurfaceBounds
  explicit SurfaceBoundsStub(std::size_t nValues = 0) : m_values(nValues, 0) {
    std::iota(m_values.begin(), m_values.end(), 0);
  }

#if defined(__GNUC__) && (__GNUC__ == 13 || __GNUC__ == 14) && \
    !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
  SurfaceBoundsStub(const SurfaceBoundsStub& other) = default;
  SurfaceBoundsStub& operator=(const SurfaceBoundsStub& other) = default;
#if defined(__GNUC__) && (__GNUC__ == 13 || __GNUC__ == 14) && \
    !defined(__clang__)
#pragma GCC diagnostic pop
#endif

  BoundsType type() const final { return eOther; }

  bool isCartesian() const final { return true; }

  SquareMatrix2 boundToCartesianJacobian(const Vector2& lposition) const final {
    static_cast<void>(lposition);
    return SquareMatrix2::Identity();
  }

  SquareMatrix2 boundToCartesianMetric(const Vector2& lposition) const final {
    static_cast<void>(lposition);
    return SquareMatrix2::Identity();
  }

  std::vector<double> values() const final { return m_values; }

  bool inside(const Vector2& lposition) const final {
    static_cast<void>(lposition);
    return true;
  }

  Vector2 closestPoint(const Vector2& lposition,
                       const SquareMatrix2& metric) const final {
    static_cast<void>(metric);
    return lposition;
  }

  Vector2 center() const final { return Vector2(0.0, 0.0); }

  bool inside(const Vector2& lposition,
              const BoundaryTolerance& boundaryTolerance) const final {
    static_cast<void>(lposition);
    static_cast<void>(boundaryTolerance);
    return true;
  }

  std::ostream& toStream(std::ostream& sl) const final {
    sl << "SurfaceBoundsStub";
    return sl;
  }

 private:
  std::vector<double> m_values;
};

}  // namespace Acts

using namespace Acts;

namespace ActsTests {

namespace {

// A region with no inheritance, chart API or ownership requirements. Its
// boundary y=0 has an analytic projection for a general positive metric.
struct HalfPlaneRegion {
  bool inside(const Vector2& p) const { return p.y() <= 0.; }

  Vector2 closestPoint(const Vector2& p, const SquareMatrix2& metric) const {
    return {p.x() + metric(0, 1) / metric(0, 0) * p.y(), 0.};
  }
};

struct StreamState {
  std::ios_base::fmtflags flags;
  std::streamsize precision;
  std::streamsize width;
  char fill;
};

StreamState setNonDefaultStreamState(std::ostringstream& stream) {
  stream << std::scientific << std::showpos << std::setfill('#')
         << std::setprecision(3);
  stream.width(17);
  return {stream.flags(), stream.precision(), stream.width(), stream.fill()};
}

void checkStreamState(const std::ostringstream& stream,
                      const StreamState& state) {
  BOOST_CHECK(stream.flags() == state.flags);
  BOOST_CHECK_EQUAL(stream.precision(), state.precision);
  BOOST_CHECK_EQUAL(stream.width(), state.width);
  BOOST_CHECK_EQUAL(stream.fill(), state.fill);
}

void checkStreamStatePreserved(const SurfaceBounds& bounds) {
  std::ostringstream stream;
  const auto state = setNonDefaultStreamState(stream);

  bounds.toStream(stream);

  checkStreamState(stream, state);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(SurfacesSuite)

BOOST_AUTO_TEST_CASE(LegacyBoundsOverridesRemainEffective) {
  // Downstream bounds can override the old tolerance/distance interface, even
  // when all geometric primitives of their base bounds are final.
  struct CustomBounds : RadialBounds {
    CustomBounds() : RadialBounds(1., 4.) {}
    using RadialBounds::inside;
    bool inside(const Vector2&, const BoundaryTolerance&) const override {
      return false;
    }
    double distance(const Vector2&) const override { return 123.; }
  };
  const auto bounds = std::make_shared<const CustomBounds>();
  const auto surface =
      Surface::makeShared<DiscSurface>(Transform3::Identity(), bounds);
  const Surface& base = *surface;
  const GeometryContext context =
      GeometryContext::dangerouslyDefaultConstruct();
  const Vector2 position{2., 0.};
  BOOST_CHECK(bounds->inside(position));
  BOOST_CHECK(!base.insideBounds(position, BoundaryTolerance::None()));
  BOOST_CHECK(!base.isOnSurface(context, Vector3{2., 0., 0.}, Vector3::UnitZ(),
                                BoundaryTolerance::None()));
  BOOST_CHECK_EQUAL(base.distanceToBoundary(position), 123.);
  BOOST_CHECK_EQUAL(surface->boundsPtr().get(), bounds.get());
}

BOOST_AUTO_TEST_CASE(ExplicitChartBoundaryTolerance) {
  const HalfPlaneRegion region;
  const auto chart = detail::PolarChart::toCartesianJacobian;
  SquareMatrix2 weight;
  weight << 4., 1., 1., 2.;

  for (const double phi : {0.3, -0.3}) {
    const Vector2 position{2., phi};
    const SquareMatrix2 jacobian = chart(position);
    const SquareMatrix2 inverse = jacobian.inverse();
    const SquareMatrix2 cartesianWeight =
        inverse.transpose() * weight * inverse;
    // Minimize delta^T W delta subject to delta_y = -phi.
    const double minimumChi2 = phi * phi * (2. - 1. / 4.);
    const double sign = phi > 0. ? 1. : -1.;
    for (const double factor : {0.99, 1.01}) {
      const double threshold = sign * factor * minimumChi2;
      const bool expected = phi > 0. ? factor > 1. : factor < 1.;
      // Identity charts retain the analytic bound-coordinate result. The
      // legacy polar Chi2Bound convention is left for a numerical follow-up.
      BOOST_CHECK_EQUAL(
          detail::insideWithTolerance(
              region, position, BoundaryTolerance::Chi2Bound(weight, threshold),
              [](const Vector2&) -> SquareMatrix2 {
                return SquareMatrix2::Identity();
              }),
          expected);
      BOOST_CHECK_EQUAL(
          detail::insideWithTolerance(
              region, position,
              BoundaryTolerance::Chi2Cartesian(cartesianWeight, threshold),
              chart),
          expected);

      // The polar Euclidean metric is diag(1, r^2); distance = r*abs(phi).
      BOOST_CHECK_EQUAL(
          detail::insideWithTolerance(
              region, position,
              BoundaryTolerance::AbsoluteEuclidean(sign * factor * 0.6), chart),
          expected);
    }
    BOOST_CHECK_CLOSE(
        detail::boundaryDistance(region, position,
                                 detail::PolarChart::metric(position)),
        0.6, 1e-10);
  }
}

BOOST_AUTO_TEST_CASE(ExplicitChartBoundaryShortcuts) {
  const HalfPlaneRegion region;
  const auto unavailableChart = [](const Vector2&) -> SquareMatrix2 {
    throw std::runtime_error("This check must not request a chart");
  };
  BOOST_CHECK(detail::insideWithTolerance(
      region, {0., 1.}, BoundaryTolerance::Infinite(), unavailableChart));
  BOOST_CHECK(!detail::insideWithTolerance(
      region, {0., 1.}, BoundaryTolerance::None(), unavailableChart));
  BOOST_CHECK(detail::insideWithTolerance(
      region, {0., -1.}, BoundaryTolerance::None(), unavailableChart));
  BOOST_CHECK(detail::insideWithTolerance(
      region, {0., -1.}, BoundaryTolerance::AbsoluteEuclidean(1.),
      unavailableChart));
  // A singular forward chart does not require an inverse. The local polar
  // metric cannot distinguish angles at r=0; it is not a finite Cartesian map.
  const Vector2 origin{0., 1.};
  BOOST_CHECK_SMALL(detail::boundaryDistance(
                        region, origin, detail::PolarChart::metric(origin)),
                    1e-15);
  BOOST_CHECK(detail::insideWithTolerance(
      region, origin, BoundaryTolerance::AbsoluteEuclidean(0.1),
      detail::PolarChart::toCartesianJacobian));
  BOOST_CHECK(!detail::insideWithTolerance(
      region, {2., 0.3}, BoundaryTolerance::AbsoluteEuclidean(-0.1),
      detail::PolarChart::toCartesianJacobian));
}

/// Unit test for creating compliant/non-compliant SurfaceBounds object
BOOST_AUTO_TEST_CASE(SurfaceBoundsConstruction) {
  SurfaceBoundsStub u;
  SurfaceBoundsStub s(1);  // would act as std::size_t cast to SurfaceBounds
  SurfaceBoundsStub t(s);
  SurfaceBoundsStub v(u);
}

BOOST_AUTO_TEST_CASE(SurfaceBoundsProperties) {
  SurfaceBoundsStub surface(5);
  std::vector<double> reference{0, 1, 2, 3, 4};
  const auto& boundValues = surface.values();
  BOOST_CHECK_EQUAL_COLLECTIONS(reference.cbegin(), reference.cend(),
                                boundValues.cbegin(), boundValues.cend());
}

/// Unit test for testing SurfaceBounds properties
BOOST_AUTO_TEST_CASE(SurfaceBoundsEquality) {
  SurfaceBoundsStub surface(1);
  SurfaceBoundsStub copiedSurface(surface);
  SurfaceBoundsStub differentSurface(2);
  BOOST_CHECK_EQUAL(surface, copiedSurface);
  BOOST_CHECK_NE(surface, differentSurface);

  SurfaceBoundsStub assignedSurface;
  assignedSurface = surface;
  BOOST_CHECK_EQUAL(surface, assignedSurface);

  const auto& surfaceboundValues = surface.values();
  const auto& assignedboundValues = assignedSurface.values();
  BOOST_CHECK_EQUAL_COLLECTIONS(
      surfaceboundValues.cbegin(), surfaceboundValues.cend(),
      assignedboundValues.cbegin(), assignedboundValues.cend());
}

BOOST_AUTO_TEST_CASE(SurfaceBoundsToStreamPreservesStreamState) {
  std::vector<std::unique_ptr<SurfaceBounds>> bounds;
  bounds.push_back(
      std::make_unique<AnnulusBounds>(7.2, 12., 0.7, 1.3, Vector2{-2., 2.}));
  bounds.push_back(std::make_unique<ConeBounds>(std::numbers::pi / 8., 3., 6.));
  bounds.push_back(std::make_unique<CylinderBounds>(0.5, 10.));
  bounds.push_back(std::make_unique<DiamondBounds>(10., 20., 15., 5., 7.));
  bounds.push_back(std::make_unique<DiscTrapezoidBounds>(1., 5., 2., 6., 0.));
  bounds.push_back(std::make_unique<EllipseBounds>(1., 2., 3., 4.,
                                                   std::numbers::pi / 2., 0.));
  bounds.push_back(std::make_unique<LineBounds>(0.5, 20.));
  bounds.push_back(std::make_unique<RadialBounds>(1., 5.));
  bounds.push_back(std::make_unique<RectangleBounds>(10., 5.));
  bounds.push_back(std::make_unique<TrapezoidBounds>(1., 6., 2.));

  for (const auto& bound : bounds) {
    checkStreamStatePreserved(*bound);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

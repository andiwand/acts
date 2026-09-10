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
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/TrackFinding/Rz/RzBound.hpp"
#include "Acts/TrackFinding/Rz/RzTransport.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <memory>

using namespace Acts;
using namespace Acts::Experimental;
using namespace Acts::UnitLiterals;
namespace bdata = boost::unit_test::data;

namespace {

const double kBz = 2_T;
const RzHelix helix{kBz};

RzVector makeState(double pt, double phi, double theta, double charge) {
  RzVector v;
  v[eRzPos0] = 1_mm;
  v[eRzPos1] = -2_mm;
  v[eRzPos2] = 5_mm;
  v[eRzDir0] = std::sin(theta) * std::cos(phi);
  v[eRzDir1] = std::sin(theta) * std::sin(phi);
  v[eRzDir2] = std::cos(theta);
  v[eRzQOverP] = charge * std::sin(theta) / pt;
  return v;
}

const auto ptSamples = bdata::make({0.2_GeV, 0.5_GeV, 1_GeV, 10_GeV, 100_GeV});
const auto phiSamples = bdata::make({-2.9, -1.0, 0.3, 2.5});
const auto thetaSamples = bdata::make({0.3, 1.2, 1.6, 2.7});
const auto chargeSamples = bdata::make({-1., 1.});

}  // namespace

BOOST_AUTO_TEST_SUITE(RzTransportSuite)

BOOST_DATA_TEST_CASE(StepMatchesRungeKutta,
                     ptSamples * phiSamples * thetaSamples * chargeSamples, pt,
                     phi, theta, charge) {
  const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  const MagneticFieldContext mctx;
  using Stepper = EigenStepper<>;
  using Prop = Propagator<Stepper, VoidNavigator>;
  Prop propagator(Stepper(std::make_shared<ConstantBField>(Vector3(0, 0, kBz))),
                  VoidNavigator());
  auto cylinder = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(), std::make_shared<CylinderBounds>(120_mm, 3_m));

  RzVector v = makeState(pt, phi, theta, charge);
  const std::optional<double> s = helix.pathToCylinder(v, 120_mm);
  BOOST_REQUIRE(s.has_value());

  const Vector3 pos = v.segment<3>(eRzPos0);
  const Vector3 dir = v.segment<3>(eRzDir0);
  const BoundTrackParameters start = BoundTrackParameters::createCurvilinear(
      Vector4(pos.x(), pos.y(), pos.z(), 0.), dir, v[eRzQOverP], std::nullopt,
      ParticleHypothesis::pion());
  Prop::Options<> options(gctx, mctx);
  options.stepping.stepTolerance = 1e-8;
  const auto result = propagator.propagate(start, *cylinder, options);
  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE(result->endParameters.has_value());

  helix.step(v, *s);
  CHECK_CLOSE_ABS(*s, result->pathLength, 1e-4_mm);
  CHECK_CLOSE_ABS(v.segment<3>(eRzPos0), result->endParameters->position(gctx),
                  1e-4_mm);
  CHECK_CLOSE_ABS(v.segment<3>(eRzDir0), result->endParameters->direction(),
                  1e-6);
  CHECK_CLOSE_ABS(std::hypot(v[eRzPos0], v[eRzPos1]), 120_mm, 1e-9_mm);
}

BOOST_DATA_TEST_CASE(FixedPathJacobian,
                     ptSamples * phiSamples * thetaSamples * chargeSamples, pt,
                     phi, theta, charge) {
  const RzVector v0 = makeState(pt, phi, theta, charge);
  const double s = 87_mm;
  const RzMatrix j = helix.stepJacobian(v0, s);

  // central differences, each parameter stepped by its own scale
  RzVector eps;
  eps << 1e-3, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6, 1e-5 * std::abs(v0[eRzQOverP]);
  for (unsigned int col = 0; col < eRzSize; ++col) {
    RzVector up = v0;
    RzVector dn = v0;
    up[col] += eps[col];
    dn[col] -= eps[col];
    helix.step(up, s);
    helix.step(dn, s);
    const RzVector num = (up - dn) / (2. * eps[col]);
    for (unsigned int row = 0; row < eRzSize; ++row) {
      const double scale = std::max(1., std::abs(num[row]));
      BOOST_CHECK_MESSAGE(std::abs(j(row, col) - num[row]) < 1e-5 * scale,
                          "J(" << row << "," << col << ") = " << j(row, col)
                               << " vs " << num[row]);
    }
  }
}

BOOST_DATA_TEST_CASE(ConstrainedJacobian,
                     ptSamples * phiSamples * thetaSamples * chargeSamples, pt,
                     phi, theta, charge) {
  const double radius = 120_mm;
  const RzVector v0 = makeState(pt, phi, theta, charge);
  auto land = [&](const RzVector& start) {
    RzVector w = start;
    const std::optional<double> s = helix.pathToCylinder(w, radius);
    BOOST_REQUIRE(s.has_value());
    helix.step(w, *s);
    return std::pair{w, *s};
  };
  const auto [v1, s] = land(v0);

  RzMatrix j = helix.stepJacobian(v0, s);
  const Vector3 normal(v1[eRzPos0] / radius, v1[eRzPos1] / radius, 0.);
  RzHelix::constrainToSurface(j, helix.derivative(v1), normal);

  RzVector eps;
  eps << 1e-3, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6, 1e-5 * std::abs(v0[eRzQOverP]);
  for (unsigned int col = 0; col < eRzSize; ++col) {
    RzVector up = v0;
    RzVector dn = v0;
    up[col] += eps[col];
    dn[col] -= eps[col];
    const RzVector num = (land(up).first - land(dn).first) / (2. * eps[col]);
    for (unsigned int row = 0; row < eRzSize; ++row) {
      const double scale = std::max(1., std::abs(num[row]));
      BOOST_CHECK_MESSAGE(std::abs(j(row, col) - num[row]) < 1e-5 * scale,
                          "J(" << row << "," << col << ") = " << j(row, col)
                               << " vs " << num[row]);
    }
  }
  // the landing point stays on the cylinder whatever the start varies by
  const Eigen::Matrix<double, 1, eRzSize> onSurface =
      normal.transpose() * j.block<3, eRzSize>(eRzPos0, 0);
  BOOST_CHECK_LT(onSurface.cwiseAbs().maxCoeff(), 1e-9);
}

BOOST_DATA_TEST_CASE(DiscPlanePerigee, ptSamples * phiSamples * chargeSamples,
                     pt, phi, charge) {
  const RzVector v0 = makeState(pt, phi, 0.4, charge);

  const std::optional<double> sDisc = helix.pathToDisc(v0, 300_mm);
  BOOST_REQUIRE(sDisc.has_value());
  RzVector w = v0;
  helix.step(w, *sDisc);
  CHECK_CLOSE_ABS(w[eRzPos2], 300_mm, 1e-12_mm);

  const Vector3 point(40_mm, 25_mm, 200_mm);
  const Vector3 normal = Vector3(0.6, 0.3, 0.74).normalized();
  const std::optional<double> sPlane = helix.pathToPlane(v0, point, normal);
  BOOST_REQUIRE(sPlane.has_value());
  w = v0;
  helix.step(w, *sPlane);
  CHECK_CLOSE_ABS(normal.dot(w.segment<3>(eRzPos0) - point), 0., 1e-9_mm);

  const double sPerigee = helix.pathToPerigee(v0);
  w = v0;
  helix.step(w, sPerigee);
  // at the perigee the transverse momentum is tangent to the circle around
  // the axis
  CHECK_CLOSE_ABS(w[eRzPos0] * w[eRzDir0] + w[eRzPos1] * w[eRzDir1], 0.,
                  1e-9_mm);
  BOOST_CHECK_LT(std::hypot(w[eRzPos0], w[eRzPos1]),
                 std::hypot(v0[eRzPos0], v0[eRzPos1]) + 1e-9);
}

BOOST_DATA_TEST_CASE(NewtonMatchesClosedForm,
                     ptSamples * phiSamples * thetaSamples * chargeSamples, pt,
                     phi, theta, charge) {
  const RzVector v = makeState(pt, phi, theta, charge);
  for (const double radius : {30_mm, 120_mm, 500_mm, 1000_mm}) {
    const std::optional<double> fast = helix.pathToCylinder(v, radius);
    const std::optional<double> exact =
        helix.pathToCylinderClosedForm(v, radius);
    BOOST_CHECK_EQUAL(fast.has_value(), exact.has_value());
    if (fast.has_value() && exact.has_value()) {
      CHECK_CLOSE_ABS(*fast, *exact, 1e-6_mm);
    }
  }
}

BOOST_AUTO_TEST_CASE(StiffTrackIsStraight) {
  RzVector v = makeState(1_TeV, 0.7, 1.3, 1.);
  const std::optional<double> curved = helix.pathToCylinder(v, 30_mm);
  BOOST_REQUIRE(curved.has_value());
  v[eRzQOverP] = 0.;
  const std::optional<double> straight = helix.pathToCylinder(v, 30_mm);
  BOOST_REQUIRE(straight.has_value());
  CHECK_CLOSE_REL(*curved, *straight, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(StepJacobianOntoMatchesDense) {
  for (const double pt : {0.3_GeV, 1_GeV, 20_GeV}) {
    for (const double theta : {0.4, 1.2, 2.5}) {
      const RzVector v = makeState(pt, 0.7, theta, 1.);
      const double s = 37_mm;
      RzVector end = v;
      helix.step(end, s);
      for (const Vector3& normal :
           {Vector3(std::cos(0.9), std::sin(0.9), 0.), Vector3(0., 0., 1.),
            Vector3(0.6, 0., 0.8)}) {
        RzMatrix dense = helix.stepJacobian(v, s);
        RzHelix::constrainToSurface(dense, helix.derivative(end), normal);
        const RzHelix::StepJacobian sparse =
            helix.stepJacobianOnto(v, s, end, normal);
        CHECK_CLOSE_ABS(sparse.dense(), dense, 1e-12);
        RzMatrix c = RzMatrix::Random();
        c = (c * c.transpose()).eval();
        CHECK_CLOSE_ABS(sparse.transport(c),
                        RzMatrix(dense * c * dense.transpose()),
                        1e-9 * c.norm());
        const Eigen::Matrix<double, 3, eRzSize> rows = sparse.positionRows();
        const Eigen::Matrix<double, 3, eRzSize> denseRows =
            dense.block<3, eRzSize>(eRzPos0, 0);
        CHECK_CLOSE_ABS(rows, denseRows, 1e-12);
      }
    }
  }
}

// The closed-form bound conversion against what the surface says, on a
// tilted plane and on an annulus disc, for a state on the surface

namespace {

/// A state on the surface at a local point, with a direction, and a
/// covariance that is full rank on the seven components
std::pair<RzVector, RzMatrix> stateOn(const Surface& surface,
                                      const GeometryContext& gctx,
                                      const Vector2& local,
                                      const Vector3& direction) {
  RzVector v;
  v.segment<3>(eRzPos0) = surface.localToGlobal(gctx, local, direction);
  v.segment<3>(eRzDir0) = direction.normalized();
  v[eRzQOverP] = -0.4;
  RzMatrix c;
  for (unsigned int a = 0; a < eRzSize; ++a) {
    for (unsigned int b = 0; b < eRzSize; ++b) {
      c(a, b) = 0.01 * (1. + a) * (1. + b) + (a == b ? 0.3 : 0.);
    }
  }
  return {v, c};
}

RzBoundState throughSurface(const Surface& surface, const GeometryContext& gctx,
                            const RzVector& v, const RzMatrix& c) {
  constexpr std::array<unsigned int, eRzSize> freeOf = {
      eFreePos0, eFreePos1, eFreePos2,  eFreeDir0,
      eFreeDir1, eFreeDir2, eFreeQOverP};
  const Vector3 position = v.segment<3>(eRzPos0);
  const Vector3 direction = v.segment<3>(eRzDir0);
  const auto bound = transformFreeToBoundParameters(
      position, 0., direction, v[eRzQOverP], surface, gctx);
  BOOST_REQUIRE(bound.ok());
  FreeMatrix freeCov = FreeMatrix::Zero();
  for (unsigned int a = 0; a < eRzSize; ++a) {
    for (unsigned int b = 0; b < eRzSize; ++b) {
      freeCov(freeOf[a], freeOf[b]) = c(a, b);
    }
  }
  const FreeToBoundMatrix j =
      surface.freeToBoundJacobian(gctx, position, direction);
  BoundMatrix cov = j * freeCov * j.transpose();
  cov(eBoundTime, eBoundTime) = 1.;
  return {*bound, cov};
}

void expectSame(const RzBoundState& a, const RzBoundState& b) {
  for (unsigned int i = 0; i < eBoundSize; ++i) {
    BOOST_CHECK_SMALL(a.parameters[i] - b.parameters[i], 1e-12);
  }
  for (unsigned int r = 0; r < eBoundSize; ++r) {
    for (unsigned int s = 0; s < eBoundSize; ++s) {
      BOOST_CHECK_SMALL(a.covariance(r, s) - b.covariance(r, s),
                        1e-10 * std::max(1., std::abs(b.covariance(r, s))));
    }
  }
}

}  // namespace

BOOST_AUTO_TEST_CASE(BoundOnCartesianModuleMatchesSurface) {
  const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  Transform3 transform = Transform3::Identity();
  transform.translation() = Vector3(120., -35., 410.);
  transform.linear() =
      (AngleAxis3(0.4, Vector3::UnitZ()) * AngleAxis3(1.1, Vector3::UnitX()) *
       AngleAxis3(-0.3, Vector3::UnitY()))
          .toRotationMatrix();
  auto surface = Surface::makeShared<PlaneSurface>(
      transform, std::make_shared<RectangleBounds>(20., 40.));
  RzModule module;
  module.center = transform.translation();
  module.u = transform.rotation().col(0);
  module.v = transform.rotation().col(1);
  module.normal = transform.rotation().col(2);
  module.halfU = 20.;
  module.halfV = 40.;
  for (const Vector2& local : {Vector2(3., -7.), Vector2(-15., 31.)}) {
    for (const Vector3& direction :
         {Vector3(0.3, 0.2, 0.9), Vector3(-0.8, 0.1, 0.4)}) {
      const auto [v, c] = stateOn(*surface, gctx, local, direction);
      const std::optional<RzBoundState> closed = rzBoundOnModule(module, v, c);
      BOOST_REQUIRE(closed.has_value());
      expectSame(*closed, throughSurface(*surface, gctx, v, c));
    }
  }
}

BOOST_AUTO_TEST_CASE(BoundOnPolarModuleMatchesSurface) {
  const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  Transform3 transform = Transform3::Identity();
  transform.translation() = Vector3(0., 0., 2700.);
  transform.linear() = AngleAxis3(0.7, Vector3::UnitZ()).toRotationMatrix();
  // an annulus with its own origin at the surface origin, so that the bound
  // coordinates are plain polar there
  auto surface = Surface::makeShared<DiscSurface>(
      transform, std::make_shared<AnnulusBounds>(400., 500., -0.1, 0.1,
                                                 Vector2(0., 0.), 0.));
  RzModule module;
  // the box centre is not the surface origin; the closed form finds the
  // origin back from it
  module.localCenter = Vector2(450., 3.);
  module.u = transform.rotation().col(0);
  module.v = transform.rotation().col(1);
  module.normal = transform.rotation().col(2);
  module.center = transform.translation() + module.localCenter.x() * module.u +
                  module.localCenter.y() * module.v;
  module.polar = true;
  module.polarIsPlain = true;
  for (const Vector2& local : {Vector2(430., 0.05), Vector2(495., -0.08)}) {
    for (const Vector3& direction :
         {Vector3(0.15, 0.05, 0.98), Vector3(0.2, -0.1, 0.9)}) {
      const auto [v, c] = stateOn(*surface, gctx, local, direction);
      const std::optional<RzBoundState> closed = rzBoundOnModule(module, v, c);
      BOOST_REQUIRE(closed.has_value());
      expectSame(*closed, throughSurface(*surface, gctx, v, c));
    }
  }
  module.polarIsPlain = false;
  const auto [v, c] =
      stateOn(*surface, gctx, Vector2(430., 0.), Vector3(0., 0., 1.));
  BOOST_CHECK(!rzBoundOnModule(module, v, c).has_value());
}

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <cstdint>

#include "RzMaterial.hpp"
#include "RzMeasurement.hpp"

using namespace Acts;
using namespace Acts::Experimental;
namespace rz = Acts::Experimental::detail::rz;

BOOST_AUTO_TEST_SUITE(RzFilteringSuite)

BOOST_AUTO_TEST_CASE(TwoScatteringPlanesGiveSignedProjectedNoise) {
  const Vector3 direction = Vector3::UnitZ();
  const Vector3 normal(0., 0.6, 0.8);
  for (double sign : {-1., 1.}) {
    rz::ProcessNoise noise;
    noise.varAngle = 4.;
    noise.varQOverP = 0.25;
    noise.advance(sign * 2.);
    noise.varAngle += 9.;
    noise.varQOverP += 0.5;
    noise.advance(sign * 3.);

    // Two independent kicks, five and three length units before this plane.
    RzMatrix expected = RzMatrix::Zero();
    for (std::uint32_t axis = 0; axis < 2; ++axis) {
      expected(eRzPos0 + axis, eRzPos0 + axis) = 4. * 25. + 9. * 9.;
      expected(eRzDir0 + axis, eRzDir0 + axis) = 4. + 9.;
      expected(eRzPos0 + axis, eRzDir0 + axis) = sign * (4. * 5. + 9. * 3.);
      expected(eRzDir0 + axis, eRzPos0 + axis) =
          expected(eRzPos0 + axis, eRzDir0 + axis);
    }
    expected(eRzQOverP, eRzQOverP) = 0.75;
    RzMatrix projection = RzMatrix::Identity();
    projection.block<3, 3>(eRzPos0, eRzPos0) -=
        direction * normal.transpose() / normal.dot(direction);
    expected = (projection * expected * projection.transpose()).eval();

    RzMatrix actual = RzMatrix::Zero();
    noise.apply(actual, direction, normal);
    BOOST_CHECK_SMALL((actual - expected).norm(), 1e-12);
    BOOST_CHECK(noise.empty());
    noise.apply(actual, direction, normal);
    BOOST_CHECK_SMALL((actual - expected).norm(), 1e-12);
  }
}

BOOST_AUTO_TEST_CASE(MeasurementUpdateMatchesDensePosterior) {
  const auto check = []<std::int32_t N>() {
    rz::State state;
    state.v = RzVector::Zero();
    state.v[eRzPos0] = 10.;
    state.v[eRzPos1] = 20.;
    state.v[eRzDir2] = 1.;
    state.v[eRzQOverP] = 1.;
    state.anchor = state.v;
    state.c = RzMatrix::Zero();
    state.c.block<2, 2>(eRzPos0, eRzPos0) << 4., 1., 1., 9.;
    state.c(eRzQOverP, eRzQOverP) = 0.25;

    rz::Placed measurement;
    measurement.position = Vector3(12., 17., 0.);
    measurement.u = Vector3::UnitX();
    measurement.v = Vector3::UnitY();
    measurement.normal = Vector3::UnitZ();
    measurement.cov00 = 1.;
    measurement.cov01 = 0.2;
    measurement.cov11 = 4.;
    measurement.halfV = 100.;
    measurement.maxDistance = 50.;
    measurement.pixel = N == 2;

    Eigen::Matrix<double, N, eRzSize> h =
        Eigen::Matrix<double, N, eRzSize>::Zero();
    h.template leftCols<N>().setIdentity();
    SquareMatrix2 measuredCovariance;
    measuredCovariance << 1., 0.2, 0.2, 4.;
    const Vector2 residual(2., -3.);
    const auto inverse = (h * state.c * h.transpose() +
                          measuredCovariance.template topLeftCorner<N, N>())
                             .inverse()
                             .eval();
    const auto gain = (state.c * h.transpose() * inverse).eval();
    const RzVector expectedParameters =
        state.v + gain * residual.template head<N>();
    const RzMatrix expectedCovariance = state.c - gain * h * state.c;
    const double expectedChi2 =
        residual.template head<N>().dot(inverse * residual.template head<N>());

    const rz::MeasurementEvaluator evaluator(50., 60., 2.);
    const auto evaluation =
        evaluator.evaluate(state, measurement, false, false);
    BOOST_REQUIRE(evaluation.has_value());
    BOOST_CHECK_CLOSE(evaluation->chi2, expectedChi2, 1e-10);
    rz::kalmanUpdate(state, *evaluation);
    BOOST_CHECK_SMALL((state.v - expectedParameters).norm(), 1e-12);
    BOOST_CHECK_SMALL((state.c - expectedCovariance).norm(), 1e-12);
    BOOST_CHECK(state.anchor == state.v);
  };
  check.template operator()<1>();
  check.template operator()<2>();
}

BOOST_AUTO_TEST_SUITE_END()

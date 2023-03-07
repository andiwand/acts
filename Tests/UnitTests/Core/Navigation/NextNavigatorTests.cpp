// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/Portal.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/NextNavigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <memory>

Acts::GeometryContext tgContext;
Acts::MagneticFieldContext mfContext;

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(NextNavigator) {
  // TODO get some geometry

  using ActionListType = Acts::ActionList<>;
  using AbortListType = Acts::AbortList<>;

  auto bField = std::make_shared<Acts::ConstantBField>(Acts::Vector3{0, 0, 2 * Acts::UnitConstants::T});

  auto stepper = Acts::EigenStepper<>(bField);
  auto navigator = Acts::Experimental::NextNavigator({});
  auto options = Acts::PropagatorOptions<ActionListType, AbortListType>(tgContext, mfContext);
  auto propagator = Acts::Propagator<Acts::EigenStepper<>, Acts::Experimental::NextNavigator>(stepper, navigator);

  // define start parameters
  double pT = 10 * Acts::UnitConstants::GeV;
  double phi = 0 * Acts::UnitConstants::rad;
  double theta = 0 * Acts::UnitConstants::rad;
  Acts::Vector4 pos(0, 0, 0, 0);
  Acts::Vector3 mom(pT * cos(phi), pT * sin(phi), pT / tan(theta));
  Acts::CurvilinearTrackParameters start(pos, mom, mom.norm(), +1);
  // propagate to the cylinder surface
  const auto& result = propagator.propagate(start, options).value();
}

BOOST_AUTO_TEST_SUITE_END()

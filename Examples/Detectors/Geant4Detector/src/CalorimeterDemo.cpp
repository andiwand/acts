// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Portal.hpp"
#include "Acts/Detector/PortalHelper.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/NextNavigator.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"

#include <memory>

Acts::GeometryContext tgContext;
Acts::MagneticFieldContext mfContext;

using namespace Acts;
using namespace ActsExamples;

int main() {
  auto geoGdml = ActsExamples::GdmlDetectorConstruction(
      "/home/andreas/cern/scripts/calo_mockup/lar_cylinder/lar_cylinder.gdml");
  auto g4World = geoGdml.Construct();

  auto g4WorldConfig = ActsExamples::Geant4::Geant4Detector::Config();
  g4WorldConfig.name = "Calo Demo";
  g4WorldConfig.g4World = g4World;

  auto g4detector = ActsExamples::Geant4::Geant4Detector();

  auto [detector, surfaces, detectorElements] =
      g4detector.constructDetector(g4WorldConfig, Acts::getDummyLogger());

  auto bField = std::make_shared<Acts::ConstantBField>(
      Acts::Vector3(0, 0, 2 * Acts::UnitConstants::T));

  using ActionListType = Acts::ActionList<>;
  using AbortListType = Acts::AbortList<>;

  Acts::Experimental::NextNavigator::Config navCfg;
  navCfg.detector = detector.get();

  auto stepper = Acts::EigenStepper<>(bField);
  auto navigator = Acts::Experimental::NextNavigator(
      navCfg,
      Acts::getDefaultLogger("NextNavigator", Acts::Logging::Level::VERBOSE));
  auto options = Acts::PropagatorOptions<ActionListType, AbortListType>(
      tgContext, mfContext);
  auto propagator =
      Acts::Propagator<Acts::EigenStepper<>, Acts::Experimental::NextNavigator>(
          stepper, navigator,
          Acts::getDefaultLogger("Propagator", Acts::Logging::Level::VERBOSE));

  // define start parameters
  Acts::Vector4 pos(0, 0, 0, 0);
  Acts::Vector3 mom(0, 0, 10);
  Acts::CurvilinearTrackParameters start(pos, mom, mom.norm(), +1);
  // propagate to the cylinder surface
  propagator.propagate(start, options);
}

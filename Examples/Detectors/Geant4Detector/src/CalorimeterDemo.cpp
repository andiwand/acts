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
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/NextNavigator.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Plugins/Geant4/Geant4Converters.hpp"
#include "Acts/Propagator/DefaultExtension.hpp"
#include "Acts/Propagator/DenseEnvironmentExtension.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"

#include <memory>

#include "G4LogicalVolume.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"

Acts::GeometryContext tgContext;
Acts::MagneticFieldContext mfContext;

namespace {

void visitG4PhysicalVolumes(
    G4VPhysicalVolume *physVol,
    std::function<void(G4VPhysicalVolume *physVol)> visiter) {
  auto logVol = physVol->GetLogicalVolume();
  for (std::size_t d = 0; d < logVol->GetNoDaughters(); ++d) {
    auto daughter = logVol->GetDaughter(d);
    visiter(daughter);
    visitG4PhysicalVolumes(daughter, visiter);
  }
}

std::shared_ptr<Acts::Experimental::Detector> buildDetector(
    G4VPhysicalVolume *g4World) {
  G4VPhysicalVolume *g4CaloPhys = nullptr;

  visitG4PhysicalVolumes(g4World, [&](G4VPhysicalVolume *physVol) {
    if (G4StrUtil::contains(physVol->GetName(), "LAr_phys")) {
      g4CaloPhys = physVol;
    }
  });

  if (g4CaloPhys == nullptr) {
    std::cout << "calo volume not found" << std::endl;
    return {};
  }

  auto g4CaloLog = g4CaloPhys->GetLogicalVolume();
  auto g4CaloMaterial = g4CaloLog->GetMaterial();
  auto g4CaloSolid = dynamic_cast<G4Tubs *>(g4CaloLog->GetSolid());

  auto caloBounds = Acts::Geant4VolumeConverter{}.cylinderBounds(*g4CaloSolid);
  auto caloTransform = Acts::Geant4AlgebraConverter{}.transform(*g4CaloPhys);
  auto caloMaterial = Acts::Geant4MaterialConverter{}.material(*g4CaloMaterial);

  std::cout << "calo bounds " << *caloBounds << std::endl;
  std::cout << "calo transform " << caloTransform.matrix() << std::endl;

  for (auto surface : caloBounds->orientedSurfaces()) {
    std::cout << "calo bound surface " << std::tie(*surface.first, tgContext)
              << std::endl;
  }

  auto caloVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), tgContext,
      "Calo Volume", caloTransform, std::move(caloBounds),
      std::vector<std::shared_ptr<Acts::Surface>>(),
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>(),
      Acts::Experimental::tryAllPortalsAndSurfaces());
  caloVolume->assignVolumeMaterial(
      std::make_shared<Acts::HomogeneousVolumeMaterial>(
          std::move(caloMaterial)));

  auto detectorVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      Acts::Experimental::defaultPortalAndSubPortalGenerator(), tgContext,
      "Detector Volume", Acts::Transform3::Identity(),
      std::make_unique<Acts::CuboidVolumeBounds>(50000, 50000, 50000),
      std::vector<std::shared_ptr<Acts::Surface>>(),
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>(
          {caloVolume}),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  auto detector = Acts::Experimental::Detector::makeShared(
      "Detector",
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>(
          {detectorVolume}),
      Acts::Experimental::tryAllVolumes());

  return detector;
}

}  // namespace

int main() {
  auto geoGdml = ActsExamples::GdmlDetectorConstruction(
      "/home/andreas/cern/scripts/calo_mockup/lar_cylinder/lar_cylinder.gdml");
  auto g4World = geoGdml.Construct();

  auto detector = buildDetector(g4World);
  if (!detector) {
    std::cout << "detector not built" << std::endl;
    return 1;
  }

  auto bField = std::make_shared<Acts::ConstantBField>(
      Acts::Vector3(0, 0, 2 * Acts::UnitConstants::T));

  Acts::Experimental::NextNavigator::Config navCfg;
  navCfg.detector = detector.get();

  using Stepper = Acts::EigenStepper<Acts::StepperExtensionList<
      Acts::DefaultExtension, Acts::DenseEnvironmentExtension>>;
  using Navigator = Acts::Experimental::NextNavigator;
  using ActionListType = Acts::ActionList<>;
  using AbortListType = Acts::AbortList<>;
  using PropagatorOptions =
      Acts::DenseStepperPropagatorOptions<ActionListType, AbortListType>;
  using Propagator = Acts::Propagator<Stepper, Navigator>;

  auto stepper = Stepper(bField);
  auto navigator = Navigator(
      navCfg,
      Acts::getDefaultLogger("NextNavigator", Acts::Logging::Level::VERBOSE));
  auto options = PropagatorOptions(tgContext, mfContext);
  auto propagator = Propagator(
      stepper, navigator,
      Acts::getDefaultLogger("Propagator", Acts::Logging::Level::VERBOSE));

  // define start parameters
  Acts::Vector4 pos(0, 0, 0, 0);
  Acts::Vector3 mom(10, 0, 0);
  Acts::CurvilinearTrackParameters start(pos, mom, mom.norm(), +1);
  // propagate to the cylinder surface
  propagator.propagate(start, options);
}

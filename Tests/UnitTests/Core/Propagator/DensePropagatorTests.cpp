// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/EigenStepperDenseExtension.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <fstream>
#include <memory>
#include <utility>

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts::Test {

const GeometryContext geoCtx;
const MagneticFieldContext magCtx;

inline Material makeLiquidArgon() {
  return Material::fromMassDensity(140.0343868497_mm, 857.0639538668_mm,
                                   39.9476933511, 18, 1.396 * 1_g / 1_cm3);
}

inline Material makeIron() {
  return Material::fromMassDensity(17.57493465097_mm, 169.9030027586_mm,
                                   55.845110798, 26, 7.874 * 1_g / 1_cm3);
}

inline std::tuple<std::shared_ptr<const TrackingGeometry>,
                  std::vector<const Surface*>>
makeDetector(Material material, double thickness) {
  CuboidVolumeBuilder::Config conf;
  conf.position = {0., 0., 0.};
  conf.length = {4_m, 2_m, 2_m};

  {
    CuboidVolumeBuilder::VolumeConfig start;
    start.position = {-1_m, 0, 0};
    start.length = {1_m, 1_m, 1_m};
    start.name = "start";

    conf.volumeCfg.push_back(start);
  }

  if (thickness < 1_m) {
    CuboidVolumeBuilder::VolumeConfig gap;
    gap.position = {-0.25 * (1_m + thickness), 0, 0};
    gap.length = {0.5 * (1_m - thickness), 1_m, 1_m};
    gap.name = "gap1";

    conf.volumeCfg.push_back(gap);
  }

  {
    CuboidVolumeBuilder::VolumeConfig dense;
    dense.position = {0, 0, 0};
    dense.length = {thickness, 1_m, 1_m};
    dense.volumeMaterial =
        std::make_shared<const HomogeneousVolumeMaterial>(material);
    dense.name = "dense";

    conf.volumeCfg.push_back(dense);
  }

  if (thickness < 1_m) {
    CuboidVolumeBuilder::VolumeConfig gap;
    gap.position = {0.25 * (1_m + thickness), 0, 0};
    gap.length = {0.5 * (1_m - thickness), 1_m, 1_m};
    gap.name = "gap2";

    conf.volumeCfg.push_back(gap);
  }

  {
    CuboidVolumeBuilder::SurfaceConfig surface;
    surface.position = {1.5_m, 0, 0};
    surface.rotation =
        Eigen::AngleAxisd(std::numbers::pi / 2, Eigen::Vector3d::UnitY())
            .matrix();
    surface.rBounds = std::make_shared<RectangleBounds>(1_m, 1_m);

    CuboidVolumeBuilder::LayerConfig layer;
    layer.surfaceCfg.push_back(surface);

    CuboidVolumeBuilder::VolumeConfig end;
    end.position = {1_m, 0, 0};
    end.length = {1_m, 1_m, 1_m};
    end.layerCfg.push_back(layer);
    end.name = "end";

    conf.volumeCfg.push_back(end);
  }

  CuboidVolumeBuilder cvb(conf);

  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto&) {
        return cvb.trackingVolume(context, inner, nullptr);
      });
  auto detector = TrackingGeometryBuilder(tgbCfg).trackingGeometry(geoCtx);

  std::vector<const Surface*> surfaces;
  detector->visitSurfaces(
      [&](const Surface* surface) { surfaces.push_back(surface); });

  return {std::move(detector), std::move(surfaces)};
}

auto makePropagator(std::shared_ptr<const TrackingGeometry> detector,
                    std::shared_ptr<const MagneticFieldProvider> bfield) {
  // using Stepper = EigenStepper<EigenStepperDenseExtension>;
  using Stepper = SympyStepper;
  using Navigator = Navigator;
  using Propagator = Propagator<Stepper, Navigator>;

  Navigator navigator({std::move(detector), true, true, false},
                      getDefaultLogger("nav", Logging::INFO));
  Stepper stepper(std::move(bfield));
  Propagator propagator(std::move(stepper), std::move(navigator),
                        getDefaultLogger("prop", Logging::INFO));

  return propagator;
}

BOOST_DATA_TEST_CASE(dense_propagator_test,
                     bdata::make({1_GeV, 10_GeV, 100_GeV}), p) {
  const double q = 1;

  auto bfield = std::make_shared<ConstantBField>(Vector3{0., 0., 0.});
  auto material = makeLiquidArgon();
  auto [detector, surfaces] = makeDetector(material, 1000_mm);

  auto propagator = makePropagator(detector, bfield);

  auto particleHypothesis = ParticleHypothesis::muon();
  double qOverP = particleHypothesis.qOverP(p, q);
  CurvilinearTrackParameters startParams(
      Vector4(-1.5_m, 0, 0, 0), Vector3(1, 0, 0), qOverP,
      BoundVector::Constant(1e-16).asDiagonal(), particleHypothesis);

  decltype(propagator)::Options<> options(geoCtx, magCtx);
  options.maxSteps = 10000;
  options.stepping.maxStepSize = 1_m;
  options.stepping.dense.meanEnergyLoss = true;

  const Surface& target = *surfaces.back();

  auto result = propagator.propagate(startParams, target, options);

  BOOST_CHECK(result.ok());
  CHECK_CLOSE_REL(3_m, result->pathLength, 1e-6);
  BOOST_CHECK(result->endParameters);

  BoundTrackParameters endParams = result->endParameters.value();

  BOOST_CHECK(endParams.covariance());
  CHECK_CLOSE_ABS(startParams.position(geoCtx) + Vector3(3_m, 0, 0),
                  endParams.position(geoCtx), 1e-6);
  CHECK_CLOSE_ABS(startParams.direction(), endParams.direction(), 1e-6);

  const auto& cov = endParams.covariance().value();

  double endP = endParams.absoluteMomentum();
  double endVarX = cov(eBoundLoc0, eBoundLoc0);
  double endVarY = cov(eBoundLoc1, eBoundLoc1);
  double endVarQOverP = cov(eBoundQOverP, eBoundQOverP);
  double endVarP =
      std::pow(q / std::pow(endParams.qOverP(), 2), 2) * endVarQOverP;

  std::cout << "input p = " << p << std::endl;
  std::cout << "output p = " << endP << std::endl;
  std::cout << "output t = " << endParams.time() << std::endl;
  std::cout << "output std x = " << std::sqrt(endVarX) << std::endl;
  std::cout << "output std y = " << std::sqrt(endVarY) << std::endl;
  std::cout << "output std q/p = " << std::sqrt(endVarQOverP) << std::endl;
  std::cout << "output std p = " << std::sqrt(endVarP) << std::endl;

  float theta0 = computeMultipleScatteringTheta0(
      MaterialSlab(material, 1_m), particleHypothesis.absolutePdg(),
      particleHypothesis.mass(), qOverP, particleHypothesis.absoluteCharge());

  std::cout << "theta0 = " << theta0 << std::endl;

  double eloss = computeEnergyLossMean(
      MaterialSlab(material, 1_m), particleHypothesis.absolutePdg(),
      particleHypothesis.mass(), qOverP, particleHypothesis.absoluteCharge());

  std::cout << "eloss = " << eloss << std::endl;

  std::cout << std::endl;
}

void chart_eloss(double q, double p_min, double p_max, int n, Material material,
                 std::ostream& out) {
  out << "p_initial,total_mean,bethe,total_mode,landau_sigma" << std::endl;

  MaterialSlab slab(material, 1);

  auto particleHypothesis = ParticleHypothesis::muon();

  for (int i = 0; i < n; ++i) {
    const double p = std::pow(
        10, std::log10(p_min) +
                (std::log10(p_max) - std::log10(p_min)) * (1.0 * i / n));

    double qOverP = particleHypothesis.qOverP(p, q);

    double bethe =
        computeEnergyLossBethe(slab, particleHypothesis.mass(), qOverP,
                               particleHypothesis.absoluteCharge());

    double total_mean = computeEnergyLossMean(
        slab, particleHypothesis.absolutePdg(), particleHypothesis.mass(),
        qOverP, particleHypothesis.absoluteCharge());

    double total_mode = computeEnergyLossMode(
        slab, particleHypothesis.absolutePdg(), particleHypothesis.mass(),
        qOverP, particleHypothesis.absoluteCharge());

    double landau_sigma =
        computeEnergyLossLandauSigma(slab, particleHypothesis.mass(), qOverP,
                                     particleHypothesis.absoluteCharge());

    out << p << "," << total_mean << "," << bethe << "," << total_mode << ","
        << landau_sigma << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(chart_eloss_lar) {
  const double q = 1;
  const double p_min = 50_MeV;
  const double p_max = 1_TeV;
  const int n = 1000;
  const auto material = makeLiquidArgon();

  std::ofstream file("eloss_lar.csv");
  chart_eloss(q, p_min, p_max, n, material, file);
}

BOOST_AUTO_TEST_CASE(chart_eloss_fe) {
  const double q = 1;
  const double p_min = 50_MeV;
  const double p_max = 1_TeV;
  const int n = 1000;
  const auto material = makeIron();

  std::ofstream file("eloss_fe.csv");
  chart_eloss(q, p_min, p_max, n, material, file);
}

void chart_msc(double q, double p_min, double p_max, int n, Material material,
               double thickness, std::ostream& out) {
  out << "p_initial,p_final,sigma" << std::endl;

  auto bfield = std::make_shared<ConstantBField>(Vector3{0., 0., 0.});
  auto [detector, surfaces] = makeDetector(material, thickness);

  auto propagator = makePropagator(detector, bfield);

  auto particleHypothesis = ParticleHypothesis::muon();

  decltype(propagator)::Options<> options(geoCtx, magCtx);
  options.maxSteps = 10000;
  options.stepping.maxStepSize = 1_m;
  options.stepping.dense.meanEnergyLoss = true;

  const Surface& target = *surfaces.back();

  for (int i = 0; i < n; ++i) {
    const double p = std::pow(
        10, std::log10(p_min) +
                (std::log10(p_max) - std::log10(p_min)) * (1.0 * i / n));

    double qOverP = particleHypothesis.qOverP(p, q);

    CurvilinearTrackParameters startParams(
        Vector4(-1.5_m, 0, 0, 0), Vector3(1, 0, 0), qOverP,
        BoundVector::Constant(1e-16).asDiagonal(), particleHypothesis);

    auto result = propagator.propagate(startParams, target, options);

    BOOST_CHECK(result.ok());
    CHECK_CLOSE_REL(3_m, result->pathLength, 1e-6);
    BOOST_CHECK(result->endParameters);

    BoundTrackParameters endParams = result->endParameters.value();

    BOOST_CHECK(endParams.covariance());
    CHECK_CLOSE_ABS(startParams.position(geoCtx) + Vector3(3_m, 0, 0),
                    endParams.position(geoCtx), 1e-6);
    CHECK_CLOSE_ABS(startParams.direction(), endParams.direction(), 1e-6);

    const auto& cov = endParams.covariance().value();

    double endVarX = cov(eBoundLoc0, eBoundLoc0);

    out << p << "," << endParams.absoluteMomentum() << "," << std::sqrt(endVarX)
        << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(chart_msc_lar) {
  const double q = 1;
  const double p_min = 1_GeV;
  const double p_max = 1_TeV;
  const int n = 1000;
  const auto material = makeLiquidArgon();
  const double thickness = 100_mm;

  std::ofstream file("msc_lar.csv");
  chart_msc(q, p_min, p_max, n, material, thickness, file);
}

BOOST_AUTO_TEST_CASE(chart_msc_fe) {
  const double q = 1;
  const double p_min = 1_GeV;
  const double p_max = 1_TeV;
  const int n = 1000;
  const auto material = makeIron();
  const double thickness = 100_mm;

  std::ofstream file("msc_fe.csv");
  chart_msc(q, p_min, p_max, n, material, thickness, file);
}

}  // namespace Acts::Test

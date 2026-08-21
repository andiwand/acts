// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Isolated navigation benchmark.
///
/// The full-chain A/B has a noise floor around 1% on the CKF, and navigation is
/// only a sixth of it, so anything below a few percent *of navigation* is
/// invisible there. This drives the Gen1 navigator directly over a geometry
/// with real surface arrays, so that navigation is essentially all of what is
/// timed.
///
/// Two measurements:
///   - `propagate`  : a straight-line propagation through the whole geometry,
///                    i.e. navigation plus the minimal stepping needed to move
///                    along the ray.
///   - `resolve`    : `Layer::compatibleSurfaces` on its own, the single
///                    largest navigation component.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsTests/CommonHelpers/BenchmarkTools.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"

#include <numbers>
#include <random>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(NavigationBenchmark)

BOOST_AUTO_TEST_CASE(NavigatorThroughput) {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  MagneticFieldContext mctx;

  CylindricalTrackingGeometry cGeometry(gctx);
  std::shared_ptr<const TrackingGeometry> geometry = cGeometry();

  // Rays from the origin, flat in eta and phi, fixed for reproducibility
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> phiDist(-std::numbers::pi,
                                                 std::numbers::pi);
  std::uniform_real_distribution<double> etaDist(-2.5, 2.5);
  constexpr std::size_t nRays = 500;
  std::vector<Vector3> directions;
  directions.reserve(nRays);
  for (std::size_t i = 0; i < nRays; ++i) {
    const double phi = phiDist(rng);
    const double theta = 2 * std::atan(std::exp(-etaDist(rng)));
    directions.push_back({std::sin(theta) * std::cos(phi),
                          std::sin(theta) * std::sin(phi), std::cos(theta)});
  }

  // ---- navigation driven by a straight-line propagation ----
  using Stepper = StraightLineStepper;
  using Aborters = ActorList<EndOfWorldReached>;
  using Prop = Propagator<Stepper, Navigator>;
  Navigator::Config navCfg;
  navCfg.trackingGeometry = geometry;
  Navigator navigator(navCfg);
  Prop propagator(Stepper(), std::move(navigator));

  std::size_t sink = 0;
  std::size_t nProp = 0, nOk = 0, maxSteps = 0;
  auto propResult = microBenchmark(
      [&](const Vector3& direction) {
        Prop::Options<Aborters> options(gctx, mctx);
        options.maxSteps = 1000;
        auto start = BoundTrackParameters::createCurvilinear(
            Vector4::Zero(), direction, 1_e / 1_GeV, std::nullopt,
            ParticleHypothesis::pion());
        auto result = propagator.propagate(start, options);
        ++nProp;
        if (result.ok()) {
          ++nOk;
          sink += result->steps;
          maxSteps = std::max(maxSteps, result->steps);
        }
        return sink;
      },
      directions, 40);
  std::cout << "navigation.propagate: " << propResult << std::endl;
  std::cout << "propagations=" << nProp << " ok=" << nOk
            << " meanSteps=" << (nOk > 0 ? sink / nOk : 0)
            << " maxSteps=" << maxSteps << std::endl;

  // ---- Layer::compatibleSurfaces on its own ----
  std::vector<const Layer*> layers;
  geometry->visitVolumes([&](const TrackingVolume* volume) {
    if (volume->confinedLayers() == nullptr) {
      return;
    }
    for (const auto& layer : volume->confinedLayers()->arrayObjects()) {
      if (layer->surfaceArray() != nullptr) {
        layers.push_back(layer.get());
      }
    }
  });
  std::cout << "layers with a surface array: " << layers.size() << std::endl;

  NavigationOptions<Surface> navOpts;
  navOpts.resolveSensitive = true;
  navOpts.resolveMaterial = true;
  navOpts.resolvePassive = false;
  navOpts.nearLimit = 0;
  navOpts.farLimit = std::numeric_limits<double>::max();

  std::size_t sink2 = 0;
  auto resolveResult = microBenchmark(
      [&](const Vector3& direction) {
        for (const Layer* layer : layers) {
          auto candidates = layer->compatibleSurfaces(gctx, Vector3::Zero(),
                                                      direction, navOpts);
          sink2 += candidates.size();
        }
        return sink2;
      },
      directions, 2000);
  std::cout << "navigation.compatibleSurfaces: " << resolveResult << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

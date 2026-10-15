// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/SpacePointFormation/detail/StripSpacePointCalibrationImpl.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/EventConfig.hpp"
#include "ActsFatras/Synthetic/EventGenerator.hpp"
#include "ActsFatras/Synthetic/SyntheticEvent.hpp"
#include "ActsFatras/Synthetic/detail/StripReadout.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <random>
#include <ranges>
#include <stdexcept>
#include <string>
#include <vector>

using namespace ActsFatras::Synthetic;
using ActsFatras::Synthetic::detail::readStrip;
using ActsFatras::Synthetic::detail::StripHit;
using ActsFatras::Synthetic::detail::StripLayer;
using ActsFatras::Synthetic::detail::stripLayers;
using namespace Acts::UnitLiterals;

namespace ActsTests {

namespace {

/// `MeasurementConfig::stripGapParameter`, at the value Athena's
/// `ITkSiSpacePointMakerToolCfg` sets.
constexpr float kGapParameter = 0.0015f;

/// The ITk's strip barrel, near enough: 26 mrad between the sensors of a
/// module, 2 mm between them and 75.5 um strips.
StripSensor makeSensor() {
  return StripSensor{.stereoAngle = 26e-3f,
                     .pitch = 75.5_um,
                     .moduleGap = 2._mm,
                     .halfLength = 24.1_mm};
}

/// A detector of one pixel barrel and one strip barrel, which is the least that
/// tells the two readouts apart. `strips` names the module the strip layers are
/// built of, the way a shipped description does.
DetectorDescription makeTestDescription(bool strips = true) {
  const Acts::MaterialSlab sensor =
      materialSlab(93.7f, 465.2f, 28.0855f, 14.f, 8.2925e-5f, 0.015f * 93.7f);

  DetectorDescription description;

  SubsystemDescription pixels;
  pixels.name = "pixel";
  BarrelDescription pixelBarrel;
  pixelBarrel.cylinders.push_back(
      CylinderDescription{.radius = 100.f,
                          .halfLengthZ = 400.f,
                          .material = SurfaceMaterial{sensor}});
  pixels.barrels.push_back(std::move(pixelBarrel));

  SubsystemDescription stripSystem;
  stripSystem.name = "strip";
  BarrelDescription stripBarrel;
  stripBarrel.cylinders.push_back(
      CylinderDescription{.radius = 400.f,
                          .halfLengthZ = 1000.f,
                          .material = SurfaceMaterial{sensor}});
  stripSystem.barrels.push_back(std::move(stripBarrel));
  EndcapDescription stripEndcap;
  stripEndcap.discs.push_back(
      DiscDescription{.absZ = 1500.f,
                      .rings = {RingBounds{380.f, 950.f}},
                      .material = SurfaceMaterial{sensor}});
  stripSystem.endcaps.push_back(std::move(stripEndcap));

  description.subsystems = {std::move(pixels), std::move(stripSystem)};

  if (strips) {
    description.sensors["short"] = makeSensor();
    description.subsystems[1].barrels[0].cylinders[0].sensor = "short";
    description.subsystems[1].endcaps[0].discs[0].sensor = "short";
  }
  return description;
}

DetectorLayout makeTestLayout(bool strips = true) {
  return makeLayout(makeTestDescription(strips));
}

/// @return the unit vector at these angles
std::array<float, 3> directionOf(const float phi, const float theta) {
  return {std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
          std::cos(theta)};
}

}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasSyntheticStripReadoutSuite)

/// A layer is a strip layer because it names a module type the description
/// holds. Nothing else marks one, and a name the description does not hold is a
/// mistake rather than a layer with no readout.
BOOST_AUTO_TEST_CASE(ANamedSensorIsWhatMakesALayerAStripLayer) {
  BOOST_CHECK(stripLayers(makeTestLayout(/*strips=*/false)).empty());

  const DetectorLayout layout = makeTestLayout();
  const std::vector<StripLayer> strips = stripLayers(layout);
  BOOST_REQUIRE_EQUAL(strips.size(), layout.layers.size());
  std::size_t counted = 0;
  for (std::size_t index = 0; index < layout.layers.size(); ++index) {
    const bool isStrip =
        layout.subsystems[layout.layers[index].subsystem] == "strip";
    BOOST_CHECK_EQUAL(strips[index].strip, isStrip);
    BOOST_CHECK_EQUAL(layout.layers[index].sensor.has_value(), isStrip);
    counted += isStrip ? 1 : 0;
  }
  // one barrel cylinder and the disc on either side of the interaction point
  BOOST_CHECK_EQUAL(counted, 3u);

  DetectorDescription unknown = makeTestDescription();
  unknown.subsystems[1].barrels[0].cylinders[0].sensor = "scintillator";
  BOOST_CHECK_THROW(makeLayout(unknown), std::invalid_argument);
}

/// Naming the module rather than restating it is the point: two layers of one
/// name are one module, and a ring may still name its own.
BOOST_AUTO_TEST_CASE(SensorsAreSharedByName) {
  DetectorDescription description = makeTestDescription();
  BOOST_REQUIRE_EQUAL(description.sensors.size(), 1u);

  // the barrel and the endcap both name "short", so both resolve to it
  const DetectorLayout shared = makeLayout(description);
  for (const DetectorLayer& layer : shared.layers) {
    if (layer.sensor.has_value()) {
      BOOST_CHECK(*layer.sensor == makeSensor());
    }
  }

  // a ring overrides its disc, which is how an endcap of unequal rings is said
  StripSensor longer = makeSensor();
  longer.halfLength = 2.f * makeSensor().halfLength;
  description.sensors["long"] = longer;
  description.subsystems[1].endcaps[0].discs[0].rings[0].sensor = "long";

  const DetectorLayout mixed = makeLayout(description);
  for (const DetectorLayer& layer : mixed.layers) {
    if (!layer.sensor.has_value()) {
      continue;
    }
    const bool disc = layer.shape == SurfaceShape::Disc;
    BOOST_CHECK_EQUAL(layer.sensor->halfLength,
                      disc ? longer.halfLength : makeSensor().halfLength);
  }

  // and a selection still resolves what it kept
  const std::vector<std::string> stripOnly{"strip"};
  const DetectorLayout selected =
      makeLayout(selectSubsystems(description, stripOnly));
  BOOST_CHECK_EQUAL(std::ranges::count_if(selected.layers,
                                          [](const DetectorLayer& l) {
                                            return l.sensor.has_value();
                                          }),
                    3);
}

/// The stereo angle sets the resolution: sharp across the strips, an order of
/// magnitude worse along them, and both follow the pitch.
BOOST_AUTO_TEST_CASE(ResolutionFollowsTheStereoAngle) {
  const auto barrelOf = [](const DetectorLayout& layout) {
    const std::vector<StripLayer> strips = stripLayers(layout);
    const auto found = std::ranges::find_if(
        strips, [](const StripLayer& s) { return s.strip; });
    BOOST_REQUIRE(found != strips.end());
    return *found;
  };

  const StripLayer barrel = barrelOf(makeTestLayout());
  // 75.5 um over the square root of twelve, then over sqrt(2) sin(theta/2)
  BOOST_CHECK_CLOSE(std::sqrt(barrel.varianceAlong), 1.185_mm, 1.);
  // and over sqrt(2) cos(theta/2), two strips measuring one coordinate
  BOOST_CHECK_CLOSE(std::sqrt(barrel.varianceAcross), 15.4_um, 1.);

  // halving the angle doubles the error along the strip and leaves the other
  DetectorDescription narrower = makeTestDescription();
  narrower.sensors["short"].stereoAngle = 13e-3f;
  const StripLayer narrow = barrelOf(makeLayout(narrower));
  BOOST_CHECK_CLOSE(std::sqrt(narrow.varianceAlong),
                    2.f * std::sqrt(barrel.varianceAlong), 0.1);
  BOOST_CHECK_CLOSE(std::sqrt(narrow.varianceAcross),
                    std::sqrt(barrel.varianceAcross), 0.1);
}

/// The pair the generator emits is the module's real geometry: given the
/// direction and perfect strips, resolving it gives the crossing back exactly.
/// This is what says the two strip lines are where they should be.
BOOST_AUTO_TEST_CASE(ThePairIsTheModuleGeometry) {
  DetectorDescription description = makeTestDescription();
  // strips of no width, so nothing but the geometry is left
  description.sensors["short"].pitch = 0.f;
  const std::vector<StripLayer> strips = stripLayers(makeLayout(description));
  const StripLayer& strip = *std::ranges::find_if(
      strips, [](const StripLayer& s) { return s.strip; });

  std::mt19937 rng{12345};
  constexpr float radius = 400.f;
  for (std::size_t trial = 0; trial < 200; ++trial) {
    const float phi = 0.4f * static_cast<float>(trial % 7) - 1.2f;
    const float theta =
        std::numbers::pi_v<float> * (0.30f + 0.02f * (trial % 11));
    // a direction with a transverse slope, i.e. not one the beam spot would
    // have guessed
    const std::array<float, 3> direction = directionOf(phi + 0.12f, theta);
    const std::array<float, 3> position{radius * std::cos(phi),
                                        radius * std::sin(phi), 120.f};

    const std::optional<StripHit> read = readStrip(
        rng, strip, kGapParameter, /*cylinder=*/true, position, direction);
    BOOST_REQUIRE(read.has_value());

    std::array<float, 3> moved{};
    BOOST_REQUIRE(Acts::detail::calibrateOuterStripSpacePoint(
        direction,
        Acts::detail::deriveOuterStripSpacePointCalibrationDetails(
            read->strips),
        moved, 2.f));
    // the calibrated point sits on the outer sensor, half a gap out
    const float incidence =
        direction[0] * std::cos(phi) + direction[1] * std::sin(phi);
    const float trueZ = position[2] + strip.halfGap / incidence * direction[2];
    BOOST_CHECK_SMALL(moved[2] - trueZ, 1e-2f);
  }
}

/// What the whole exercise rests on. Resolving the pair from the beam spot is
/// *biased* for a track that did not come straight out of it, by millimetres
/// along the strip; the same pair given the true direction is not. The random
/// part is the same either way -- two strips at 26 mrad cannot do better than a
/// millimetre along themselves, and no amount of knowing the direction changes
/// that -- so what a calibration buys is the bias and not the resolution.
BOOST_AUTO_TEST_CASE(CalibrationRemovesTheVertexBias) {
  const DetectorLayout layout = makeTestLayout();
  const std::vector<StripLayer> strips = stripLayers(layout);
  const StripLayer& strip = *std::ranges::find_if(
      strips, [](const StripLayer& s) { return s.strip; });

  const auto measure = [&](const float slope) {
    std::mt19937 rng{12345};
    constexpr float radius = 400.f;
    double nominal = 0., nominal2 = 0., calibrated = 0., calibrated2 = 0.;
    std::size_t resolved = 0;
    for (std::size_t trial = 0; trial < 4000; ++trial) {
      const float phi = 0.013f * static_cast<float>(trial % 479) - 3.1f;
      const float theta =
          std::numbers::pi_v<float> * (0.30f + 0.02f * (trial % 11));
      const std::array<float, 3> direction = directionOf(phi + slope, theta);
      const std::array<float, 3> position{radius * std::cos(phi),
                                          radius * std::sin(phi),
                                          2.f * static_cast<float>(trial % 97)};

      const std::optional<StripHit> read = readStrip(
          rng, strip, kGapParameter, /*cylinder=*/true, position, direction);
      if (!read.has_value()) {
        continue;
      }
      ++resolved;

      const float incidence =
          direction[0] * std::cos(phi) + direction[1] * std::sin(phi);
      // the nominal point sits on the inner sensor and the calibrated one on
      // the outer, so each is compared against the track where it is
      const float innerZ =
          position[2] - strip.halfGap / incidence * direction[2];
      const float outerZ =
          position[2] + strip.halfGap / incidence * direction[2];
      nominal += read->hit.z - innerZ;
      nominal2 += (read->hit.z - innerZ) * (read->hit.z - innerZ);

      std::array<float, 3> moved{};
      BOOST_REQUIRE(Acts::detail::calibrateOuterStripSpacePoint(
          direction,
          Acts::detail::deriveOuterStripSpacePointCalibrationDetails(
              read->strips),
          moved, 2.f));
      calibrated += moved[2] - outerZ;
      calibrated2 += (moved[2] - outerZ) * (moved[2] - outerZ);
    }
    BOOST_REQUIRE_GT(resolved, 100u);
    const double n = static_cast<double>(resolved);
    struct Result {
      double nominalBias, nominalRms, calibratedBias, calibratedRms;
      double resolvedFraction;
    };
    return Result{
        nominal / n, std::sqrt(nominal2 / n - (nominal / n) * (nominal / n)),
        calibrated / n,
        std::sqrt(calibrated2 / n - (calibrated / n) * (calibrated / n)),
        n / 4000.};
  };

  // A track straight out of the beam spot: the assumption reconstruction makes
  // is the truth, so neither is biased.
  const auto radial = measure(0.f);
  BOOST_CHECK_SMALL(radial.nominalBias, 0.2);
  BOOST_CHECK_SMALL(radial.calibratedBias, 0.2);
  // and almost all of them resolve, the pair landing on both sensors
  BOOST_CHECK_GT(radial.resolvedFraction, 0.9);

  // A track with a transverse slope of 0.12 rad, which is a GeV at this radius
  // in a two tesla field. The beam spot now guesses the direction wrong and the
  // nominal point walks along the strip by centimetres.
  const auto curved = measure(0.12f);
  BOOST_CHECK_GT(std::abs(curved.nominalBias), 5.);
  BOOST_CHECK_SMALL(curved.calibratedBias, 0.2);
  // It costs no acceptance, though, because the tolerance the layer derives is
  // this same walk evaluated at `stripGapParameter`: the slope it admits is
  // `kGapParameter * radius`, 0.6 rad here, and a GeV is a fifth of it.
  BOOST_CHECK_GT(curved.resolvedFraction, 0.9);

  // And the random part is the stereo-limited millimetre either way, so the
  // generator has to carry the pair: an inflated Gaussian would have the same
  // width and no bias to remove.
  for (const double rms : {radial.nominalRms, radial.calibratedRms,
                           curved.nominalRms, curved.calibratedRms}) {
    BOOST_CHECK_GT(rms, 0.5);
    BOOST_CHECK_LT(rms, 3.);
  }

  // Softer than the tolerance was derived for and the pair stops resolving,
  // which is the acceptance the parameter buys and what it costs.
  const auto accepted = [&](const float slope) {
    std::mt19937 rng{12345};
    constexpr float radius = 400.f;
    std::size_t resolved = 0;
    for (std::size_t trial = 0; trial < 4000; ++trial) {
      const float phi = 0.013f * static_cast<float>(trial % 479) - 3.1f;
      const float theta =
          std::numbers::pi_v<float> * (0.30f + 0.02f * (trial % 11));
      const std::array<float, 3> position{radius * std::cos(phi),
                                          radius * std::sin(phi),
                                          2.f * static_cast<float>(trial % 97)};
      resolved += readStrip(rng, strip, kGapParameter, /*cylinder=*/true,
                            position, directionOf(phi + slope, theta))
                      .has_value();
    }
    return static_cast<double>(resolved) / 4000.;
  };
  BOOST_CHECK_LT(accepted(1.5f), 0.05);
}

/// A track along the module resolves nothing, which is a strip measuring
/// nothing rather than a point at infinity.
BOOST_AUTO_TEST_CASE(GrazingTrackResolvesNothing) {
  const DetectorLayout layout = makeTestLayout();
  const std::vector<StripLayer> strips = stripLayers(layout);
  const StripLayer& strip = *std::ranges::find_if(
      strips, [](const StripLayer& s) { return s.strip; });

  std::mt19937 rng{1};
  const std::array<float, 3> position{400.f, 0.f, 0.f};
  // straight along the cylinder, so it never leaves the module plane
  BOOST_CHECK(!readStrip(rng, strip, kGapParameter, /*cylinder=*/true, position,
                         std::array<float, 3>{0.f, 0.f, 1.f})
                   .has_value());
}

/// The generator emits the two kinds of space point in two collections, and
/// only the strip one carries the pair.
BOOST_AUTO_TEST_CASE(TheTwoCollectionsAreSeparate) {
  const DetectorLayout layout = makeTestLayout();
  EventConfig config;
  config.generation.pileup = 20;
  config.simulation.propagation.maxTurns = 0.5f;

  const Event event = generateEvent(layout, config);
  BOOST_CHECK_GT(event.spacePoints.size(), 0u);
  BOOST_CHECK_GT(event.stripSpacePoints.size(), 0u);
  BOOST_CHECK(event.stripSpacePoints.hasColumns(
      Acts::SpacePointColumns::StripCalibrationDetails));
  BOOST_CHECK(!event.spacePoints.hasColumns(
      Acts::SpacePointColumns::StripCalibrationDetails));

  // each collection holds only its own kind of layer, in one index space
  const auto check = [&layout](const Acts::SpacePointContainer& points,
                               const std::string& subsystem) {
    const auto layers = points.column<std::uint32_t>("layerId");
    for (const auto sp : points) {
      const std::uint32_t index = sp.extra(layers);
      BOOST_REQUIRE_LT(index, layout.layers.size());
      BOOST_CHECK_EQUAL(layout.subsystems[layout.layers[index].subsystem],
                        subsystem);
    }
  };
  check(event.spacePoints, "pixel");
  check(event.stripSpacePoints, "strip");

  // with the sensors taken off, everything is a pixel again
  const Event pixelOnly =
      generateEvent(makeTestLayout(/*strips=*/false), config);
  BOOST_CHECK_GT(pixelOnly.spacePoints.size(), event.spacePoints.size());
  BOOST_CHECK_EQUAL(pixelOnly.stripSpacePoints.size(), 0u);
}

/// A seeder takes one container, so the two are gathered into one. What the
/// gathering has to preserve is what a caller looks a point up by: its layer,
/// its particle, and the pair a strip point carries.
BOOST_AUTO_TEST_CASE(TheTwoCollectionsGatherIntoOne) {
  const DetectorLayout layout = makeTestLayout();
  EventConfig config;
  config.generation.pileup = 20;
  config.simulation.propagation.maxTurns = 0.5f;

  const Event event = generateEvent(layout, config);
  const std::size_t pixels = event.spacePoints.size();
  const std::size_t strips = event.stripSpacePoints.size();
  BOOST_REQUIRE_GT(pixels, 0u);
  BOOST_REQUIRE_GT(strips, 0u);

  BOOST_CHECK_EQUAL(selectSpacePoints(event, SpacePointSelection::Pixel).size(),
                    pixels);
  BOOST_CHECK_EQUAL(selectSpacePoints(event, SpacePointSelection::Strip).size(),
                    strips);

  const Acts::SpacePointContainer both =
      selectSpacePoints(event, SpacePointSelection::Combined);
  BOOST_CHECK_EQUAL(both.size(), pixels + strips);
  // the pair column comes with the strips and is dead weight on a pixel point
  BOOST_CHECK(
      both.hasColumns(Acts::SpacePointColumns::StripCalibrationDetails));
  BOOST_CHECK(
      !selectSpacePoints(event, SpacePointSelection::Pixel)
           .hasColumns(Acts::SpacePointColumns::StripCalibrationDetails));

  const auto layers = both.column<std::uint32_t>("layerId");
  const auto particles = both.column<std::uint32_t>("particleId");
  const auto sourceLayers =
      event.stripSpacePoints.column<std::uint32_t>("layerId");
  const auto sourceParticles =
      event.stripSpacePoints.column<std::uint32_t>("particleId");
  for (const auto sp : both) {
    const bool strip = sp.index() >= pixels;
    BOOST_CHECK_EQUAL(layout.layers[sp.extra(layers)].sensor.has_value(),
                      strip);
    if (!strip) {
      continue;
    }
    // and a strip point still says which of the event's it was made from
    const auto source = event.stripSpacePoints[sp.copiedFromIndex()];
    BOOST_CHECK_EQUAL(sp.extra(layers), source.extra(sourceLayers));
    BOOST_CHECK_EQUAL(sp.extra(particles), source.extra(sourceParticles));
    BOOST_CHECK_EQUAL(sp.z(), source.z());
    BOOST_CHECK_EQUAL(sp.outerStripCalibrationDetails().outerCenter[2],
                      source.outerStripCalibrationDetails().outerCenter[2]);
  }
}

/// How a real track fares: from the beam spot with a beam-spot z spread, bent
/// by the field, met at the innermost strip barrel. The two things that decide
/// whether strips are usable at all are how often a pair resolves and how well
/// it lands, and both are functions of momentum.
BOOST_AUTO_TEST_CASE(AcceptanceAndResolutionAgainstMomentum) {
  const DetectorLayout layout = makeTestLayout();
  const std::vector<StripLayer> strips = stripLayers(layout);
  const StripLayer& strip = *std::ranges::find_if(
      strips, [](const StripLayer& s) { return s.strip; });
  constexpr float radius = 400.f;

  std::mt19937 rng{1};
  std::normal_distribution<float> z0Of(0.f, 50.f);
  std::uniform_real_distribution<float> phiOf(-3.14f, 3.14f);
  std::uniform_real_distribution<float> etaOf(-1.f, 1.f);

  for (const float pt : {0.2f, 0.5f, 1.f, 2.f, 5.f, 100.f}) {
    std::size_t tried = 0, resolved = 0;
    double sum = 0., sum2 = 0.;
    for (std::size_t n = 0; n < 20000; ++n) {
      const float z0 = z0Of(rng);
      const float phi0 = phiOf(rng);
      const float theta = 2.f * std::atan(std::exp(-etaOf(rng)));
      const float cotTheta = std::cos(theta) / std::sin(theta);
      // a helix out of (0, 0, z0) in a two tesla field
      const float bendRadius = pt * 1000.f / (0.3f * 2.f);
      const float half = radius / (2.f * bendRadius);
      if (half >= 1.f) {
        continue;
      }
      const float gamma = 2.f * std::asin(half);
      const float phiPos = phi0 + 0.5f * gamma;
      const std::array<float, 3> position{radius * std::cos(phiPos),
                                          radius * std::sin(phiPos),
                                          z0 + bendRadius * gamma * cotTheta};
      if (std::abs(position[2]) > 900.f) {
        continue;
      }
      const std::array<float, 3> direction = directionOf(phi0 + gamma, theta);
      ++tried;

      const std::optional<StripHit> read = readStrip(
          rng, strip, kGapParameter, /*cylinder=*/true, position, direction);
      if (!read.has_value()) {
        continue;
      }
      ++resolved;
      const float incidence =
          direction[0] * std::cos(phiPos) + direction[1] * std::sin(phiPos);
      const float trueZ =
          position[2] - strip.halfGap / incidence * direction[2];
      sum += read->hit.z - trueZ;
      sum2 += (read->hit.z - trueZ) * (read->hit.z - trueZ);
    }
    BOOST_REQUIRE_GT(tried, 100u);
    const double n = std::max<std::size_t>(1, resolved);
    BOOST_TEST_MESSAGE("pt " << pt << " GeV: resolved "
                             << 100. * static_cast<double>(resolved) /
                                    static_cast<double>(tried)
                             << "%  bias " << sum / n << "  rms "
                             << std::sqrt(sum2 / n - (sum / n) * (sum / n))
                             << " mm");
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

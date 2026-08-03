// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/TrackFitting/GlobalChiSquareFitter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

#include <cmath>
#include <numbers>
#include <unordered_map>
#include <vector>

using namespace Acts;
using namespace Acts::Experimental;
using namespace Acts::UnitLiterals;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(TrackFittingSuite)

// The extended parameter vector packs the parameters of each material surface
// after the bound parameters. Check that the offsets of the four possible
// configurations are consistent, non-overlapping and inside the system.
BOOST_AUTO_TEST_CASE(Gx2fParameterLayoutOffsets) {
  constexpr std::size_t nMat = 4;

  // Nothing fitted
  {
    const Gx2fParameterLayout layout{false, false, nMat};
    BOOST_CHECK(!layout.fitMaterial());
    BOOST_CHECK_EQUAL(layout.stride(), 0u);
    // Without any fitted parameters the surfaces do not count
    BOOST_CHECK_EQUAL(layout.nMaterialSurfaces(), 0u);
    BOOST_CHECK_EQUAL(layout.nDims(), eBoundSize);
  }

  // Scattering only, the historical layout
  {
    const Gx2fParameterLayout layout{true, false, nMat};
    BOOST_CHECK(layout.fitMaterial());
    BOOST_CHECK_EQUAL(layout.stride(), 2u);
    BOOST_CHECK_EQUAL(layout.nDims(), eBoundSize + 2 * nMat);
    for (std::size_t k = 0; k < nMat; k++) {
      BOOST_CHECK_EQUAL(layout.scatteringOffset(k), eBoundSize + 2 * k);
    }
  }

  // Energy loss only
  {
    const Gx2fParameterLayout layout{false, true, nMat};
    BOOST_CHECK(layout.fitMaterial());
    BOOST_CHECK_EQUAL(layout.stride(), 1u);
    BOOST_CHECK_EQUAL(layout.nDims(), eBoundSize + nMat);
    for (std::size_t k = 0; k < nMat; k++) {
      BOOST_CHECK_EQUAL(layout.energyLossOffset(k), eBoundSize + k);
    }
  }

  // Both
  {
    const Gx2fParameterLayout layout{true, true, nMat};
    BOOST_CHECK_EQUAL(layout.stride(), 3u);
    BOOST_CHECK_EQUAL(layout.nDims(), eBoundSize + 3 * nMat);

    std::vector<std::size_t> seen;
    for (std::size_t k = 0; k < nMat; k++) {
      seen.push_back(layout.scatteringOffset(k));
      seen.push_back(layout.scatteringOffset(k) + 1);
      seen.push_back(layout.energyLossOffset(k));
    }

    // Strictly increasing, hence non-overlapping, and all inside the system
    for (std::size_t i = 0; i < seen.size(); i++) {
      BOOST_CHECK_LT(seen[i], layout.nDims());
      BOOST_CHECK_GE(seen[i], eBoundSize);
      if (i > 0) {
        BOOST_CHECK_LT(seen[i - 1], seen[i]);
      }
    }
  }
}

// The deterministic q/p correction must have the right sign and magnitude, and
// must degrade gracefully in the corner cases.
BOOST_AUTO_TEST_CASE(Gx2fQOverPOffsetSign) {
  const MaterialSlab slab{makeSilicon(), 5_mm};
  const auto muon = ParticleHypothesis::muon();
  const double qOverP = 1. / 1_GeV;

  const double offsetForward = computeGx2fQOverPOffset(
      slab, muon, qOverP, Direction::Forward(), Gx2fEnergyLossMode::Mode);

  // Forward propagation loses energy, so |p| decreases and |q/p| grows. The
  // offset therefore carries the sign of q/p.
  BOOST_CHECK_GT(offsetForward, 0.);

  const double offsetBackward = computeGx2fQOverPOffset(
      slab, muon, qOverP, Direction::Backward(), Gx2fEnergyLossMode::Mode);
  BOOST_CHECK_LT(offsetBackward, 0.);
  // To first order the two are mirror images
  CHECK_CLOSE_REL(offsetForward, -offsetBackward, 1e-2);

  // A negatively charged particle has a negative q/p, and the offset follows it
  const double offsetNegative = computeGx2fQOverPOffset(
      slab, muon, -qOverP, Direction::Forward(), Gx2fEnergyLossMode::Mode);
  BOOST_CHECK_LT(offsetNegative, 0.);
  CHECK_CLOSE_REL(offsetForward, -offsetNegative, 1e-6);

  // The mean includes the full radiative term, the mode only 15% of it, so the
  // mean loss is the larger one
  const double offsetMean = computeGx2fQOverPOffset(
      slab, muon, qOverP, Direction::Forward(), Gx2fEnergyLossMode::Mean);
  BOOST_CHECK_GT(offsetMean, offsetForward);

  // Vacuum does not change the momentum
  const MaterialSlab vacuum = MaterialSlab::Nothing();
  BOOST_CHECK_EQUAL(
      computeGx2fQOverPOffset(vacuum, muon, qOverP, Direction::Forward(),
                              Gx2fEnergyLossMode::Mode),
      0.);

  // A neutral particle does not lose energy by ionisation, and the material
  // functions are not defined for it
  BOOST_CHECK_EQUAL(
      computeGx2fQOverPOffset(slab, ParticleHypothesis::pion0(), qOverP,
                              Direction::Forward(), Gx2fEnergyLossMode::Mode),
      0.);
}

// A very thick slab drives the particle into the momentum floor. This must stay
// finite rather than produce a NaN or an infinite q/p.
BOOST_AUTO_TEST_CASE(Gx2fQOverPOffsetMomentumFloor) {
  const MaterialSlab slab{makeIron(), 10_m};
  const auto muon = ParticleHypothesis::muon();
  const double qOverP = 1. / 1_GeV;

  const double offset = computeGx2fQOverPOffset(
      slab, muon, qOverP, Direction::Forward(), Gx2fEnergyLossMode::Mode);

  BOOST_CHECK(std::isfinite(offset));
  // Floored at 10 MeV, so q/p saturates at 1/(10 MeV)
  CHECK_CLOSE_REL(qOverP + offset, 1. / 10_MeV, 1e-6);
}

// At high momentum the q/p offset is many orders of magnitude smaller than q/p
// itself. Computing the difference in single precision would destroy it, so
// check that enough significant digits survive.
BOOST_AUTO_TEST_CASE(Gx2fQOverPOffsetHighMomentum) {
  const MaterialSlab slab{makeSilicon(), 5_mm};
  const auto muon = ParticleHypothesis::muon();
  const double qOverP = 1. / 100_GeV;

  const double offset = computeGx2fQOverPOffset(
      slab, muon, qOverP, Direction::Forward(), Gx2fEnergyLossMode::Mean);

  BOOST_CHECK_GT(offset, 0.);
  // The offset is tiny compared to q/p, but must not be lost to rounding
  BOOST_CHECK_LT(offset, 1e-3 * qOverP);
  BOOST_CHECK_GT(offset, 1e-9 * qOverP);
}

// The energy loss penalty must land on the energy loss column, and must not
// disturb the scattering entries.
BOOST_AUTO_TEST_CASE(AddMaterialToGx2fSumsEnergyLoss) {
  // A minimal stand-in for a track state. addMaterialToGx2fSums only needs the
  // reference surface geometry id and the smoothed parameters.
  struct SurfaceStub {
    GeometryIdentifier m_geoId;
    GeometryIdentifier geometryId() const { return m_geoId; }
  };
  struct TrackStateStub {
    SurfaceStub m_surface;
    BoundVector m_smoothed;
    const SurfaceStub& referenceSurface() const { return m_surface; }
    const BoundVector& smoothed() const { return m_smoothed; }
  };

  const GeometryIdentifier geoId =
      GeometryIdentifier().withVolume(1).withLayer(2);

  BoundVector smoothed = BoundVector::Zero();
  smoothed[eBoundTheta] = std::numbers::pi / 2.;  // sin(theta) == 1
  const TrackStateStub trackState{SurfaceStub{geoId}, smoothed};

  const double invCovScattering = 400.;
  const double invCovQOverP = 1e6;
  const double deltaQOverP = 3e-4;
  const double scatteringPhi = 1e-3;
  const double scatteringTheta = 2e-3;

  BoundVector angles = BoundVector::Zero();
  angles[eBoundPhi] = scatteringPhi;
  angles[eBoundTheta] = scatteringTheta;

  Gx2fMaterialProperties properties{angles, invCovScattering, true};
  properties.qOverPOffset() = deltaQOverP;
  properties.invCovarianceQOverP() = invCovQOverP;

  std::unordered_map<GeometryIdentifier, Gx2fMaterialProperties> materialMap;
  materialMap.emplace(geoId, properties);

  const Gx2fParameterLayout layout{true, true, 2};
  Gx2fSystem system{layout};

  // Handle the second of the two material surfaces, to catch offset mistakes
  constexpr std::size_t nMaterialsHandled = 1;
  addMaterialToGx2fSums(system, nMaterialsHandled, materialMap, trackState,
                        *getDefaultLogger("Gx2fComponentTests", Logging::INFO));

  const std::size_t scatteringPos = layout.scatteringOffset(nMaterialsHandled);
  const std::size_t energyLossPos = layout.energyLossOffset(nMaterialsHandled);

  // The energy loss column
  CHECK_CLOSE_REL(system.aMatrix()(energyLossPos, energyLossPos), invCovQOverP,
                  1e-12);
  CHECK_CLOSE_REL(system.bVector()(energyLossPos), -invCovQOverP * deltaQOverP,
                  1e-12);

  // The scattering columns are untouched by the energy loss contribution
  CHECK_CLOSE_REL(system.aMatrix()(scatteringPos, scatteringPos),
                  invCovScattering, 1e-12);
  CHECK_CLOSE_REL(system.aMatrix()(scatteringPos + 1, scatteringPos + 1),
                  invCovScattering, 1e-12);

  // chi2 collects all three penalties
  const double expectedChi2 =
      invCovQOverP * deltaQOverP * deltaQOverP +
      invCovScattering * scatteringPhi * scatteringPhi +
      invCovScattering * scatteringTheta * scatteringTheta;
  CHECK_CLOSE_REL(system.chi2(), expectedChi2, 1e-12);

  // Nothing leaked into the first material surface
  BOOST_CHECK_EQUAL(
      system.aMatrix()(layout.energyLossOffset(0), layout.energyLossOffset(0)),
      0.);
}

// The solved deltas must be routed to the right material surface and the right
// parameter within it.
BOOST_AUTO_TEST_CASE(UpdateGx2fParamsEnergyLoss) {
  const GeometryIdentifier geoId0 =
      GeometryIdentifier().withVolume(1).withLayer(2);
  const GeometryIdentifier geoId1 =
      GeometryIdentifier().withVolume(1).withLayer(4);

  const Gx2fParameterLayout layout{true, true, 2};

  std::unordered_map<GeometryIdentifier, Gx2fMaterialProperties> materialMap;
  materialMap.emplace(geoId0,
                      Gx2fMaterialProperties{BoundVector::Zero(), 1., true});
  materialMap.emplace(geoId1,
                      Gx2fMaterialProperties{BoundVector::Zero(), 1., true});
  const std::vector<GeometryIdentifier> geoIdVector{geoId0, geoId1};

  Eigen::VectorXd delta = Eigen::VectorXd::Zero(layout.nDims());
  delta[eBoundLoc0] = 0.5;
  delta[layout.scatteringOffset(0)] = 1e-3;
  delta[layout.scatteringOffset(0) + 1] = 2e-3;
  delta[layout.energyLossOffset(0)] = 3e-4;
  delta[layout.scatteringOffset(1)] = 4e-3;
  delta[layout.scatteringOffset(1) + 1] = 5e-3;
  delta[layout.energyLossOffset(1)] = 6e-4;

  BoundTrackParameters params = BoundTrackParameters::createCurvilinear(
      Vector4::Zero(), Vector3::UnitX(), 1. / 1_GeV, std::nullopt,
      ParticleHypothesis::muon());
  const double loc0Before = params.parameters()[eBoundLoc0];

  updateGx2fParams(params, delta, layout, materialMap, geoIdVector);

  CHECK_CLOSE_REL(params.parameters()[eBoundLoc0], loc0Before + 0.5, 1e-12);

  CHECK_CLOSE_REL(materialMap.at(geoId0).scatteringAngles()[eBoundPhi], 1e-3,
                  1e-12);
  CHECK_CLOSE_REL(materialMap.at(geoId0).scatteringAngles()[eBoundTheta], 2e-3,
                  1e-12);
  CHECK_CLOSE_REL(materialMap.at(geoId0).qOverPOffset(), 3e-4, 1e-12);

  CHECK_CLOSE_REL(materialMap.at(geoId1).scatteringAngles()[eBoundPhi], 4e-3,
                  1e-12);
  CHECK_CLOSE_REL(materialMap.at(geoId1).scatteringAngles()[eBoundTheta], 5e-3,
                  1e-12);
  CHECK_CLOSE_REL(materialMap.at(geoId1).qOverPOffset(), 6e-4, 1e-12);

  // The expectation is untouched by the fit update; only the deviation moves
  BOOST_CHECK_EQUAL(materialMap.at(geoId0).expectedQOverPOffset(), 0.);
  CHECK_CLOSE_REL(materialMap.at(geoId0).totalQOverPOffset(), 3e-4, 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

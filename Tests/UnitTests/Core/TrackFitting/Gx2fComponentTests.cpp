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
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
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

// Material parameter blocks must fit within the system without overlap.
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

// Check the signs for both charges and propagation directions.
BOOST_AUTO_TEST_CASE(Gx2fQOverPOffsetSign) {
  const MaterialSlab slab{makeSilicon(), 5_mm};
  const auto muon = ParticleHypothesis::muon();
  const double qOverP = 1. / 1_GeV;

  const double offsetForward = computeGx2fQOverPOffset(
      slab, muon, qOverP, Direction::Forward(), Gx2fEnergyLossMode::Mode);

  // Forward energy loss increases |q/p|.
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

  // The mean loss exceeds the mode.
  const double offsetMean = computeGx2fQOverPOffset(
      slab, muon, qOverP, Direction::Forward(), Gx2fEnergyLossMode::Mean);
  BOOST_CHECK_GT(offsetMean, offsetForward);

  // Vacuum does not change the momentum
  const MaterialSlab vacuum = MaterialSlab::Nothing();
  BOOST_CHECK_EQUAL(
      computeGx2fQOverPOffset(vacuum, muon, qOverP, Direction::Forward(),
                              Gx2fEnergyLossMode::Mode),
      0.);

  // Skip ionisation for neutral particles.
  BOOST_CHECK_EQUAL(
      computeGx2fQOverPOffset(slab, ParticleHypothesis::pion0(), qOverP,
                              Direction::Forward(), Gx2fEnergyLossMode::Mode),
      0.);
}

// Thick material must saturate at the momentum floor.
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

// Preserve small q/p offsets at high momentum.
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

// Route each material penalty to its parameter column.
BOOST_AUTO_TEST_CASE(AddMaterialToGx2fSumsEnergyLoss) {
  // Only the surface ID and smoothed parameters are needed.
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

// Route fitted updates to their material surfaces.
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

// Mixed units must not hide measured parameters or alter unmeasured covariance.
BOOST_AUTO_TEST_CASE(Gx2fEquilibratedSystem) {
  Gx2fSystem system{Gx2fParameterLayout{false, true, 1}};
  Eigen::VectorXd diagonal = Eigen::VectorXd::Ones(system.nDims());
  diagonal << 1e-12, 1e12, 4., 9., 0., 16., 1e24;
  system.aMatrix() = diagonal.asDiagonal();
  system.aMatrix()(eBoundLoc0, eBoundSize) = 5e5;
  system.aMatrix()(eBoundSize, eBoundLoc0) = 5e5;
  Eigen::VectorXd expected = Eigen::VectorXd::Ones(system.nDims());
  expected[eBoundQOverP] = 0.;
  expected[eBoundSize] = 1e-18;
  system.bVector() = system.aMatrix() * expected;
  CHECK_CLOSE_ABS(computeGx2fDeltaParams(system), expected, 1e-10);

  BOOST_CHECK_EQUAL(system.findRequiredNdf(), 5u);
  const Eigen::MatrixXd before = system.aMatrix();
  BoundMatrix covariance = 7. * BoundMatrix::Identity();
  updateGx2fCovarianceParams(covariance, system);
  CHECK_CLOSE_REL(covariance(eBoundLoc0, eBoundLoc0), 4. / 3. * 1e12, 1e-12);
  CHECK_CLOSE_REL(covariance(eBoundLoc1, eBoundLoc1), 1e-12, 1e-12);
  BOOST_CHECK_EQUAL(covariance(eBoundQOverP, eBoundQOverP), 7.);
  BOOST_CHECK_EQUAL(covariance(eBoundTime, eBoundTime), 1. / 16.);
  CHECK_CLOSE_ABS(system.aMatrix(), before, 0.);
}

// Include skipped transport segments and material columns on the measured
// state.
BOOST_AUTO_TEST_CASE(Gx2fSystemTransportAndLocalMaterial) {
  TrackContainer tracks{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto track = tracks.makeTrack();
  const auto geoId = GeometryIdentifier().withVolume(1).withLayer(2);
  auto surface = Surface::makeShared<PlaneSurface>(Transform3::Identity());
  surface->assignGeometryId(geoId);

  auto skipped = track.appendTrackState(Gx2fConstants::trackStateMask);
  skipped.setReferenceSurface(surface);
  skipped.jacobian() = 2. * BoundMatrix::Identity();

  auto measured = track.appendTrackState(Gx2fConstants::trackStateMask);
  measured.setReferenceSurface(surface);
  measured.jacobian() = 3. * BoundMatrix::Identity();
  measured.smoothed().setZero();
  measured.smoothed()[eBoundTheta] = std::numbers::pi / 2.;
  measured.typeFlags().setIsMeasurement();
  measured.typeFlags().setHasMaterial();
  measured.allocateCalibrated(1);
  measured.calibrated<1>()[0] = 2.;
  measured.calibratedCovariance<1>()(0, 0) = 1.;
  measured.setProjectorSubspaceIndices(std::array{eBoundQOverP});
  track.linkForward();

  Gx2fMaterialProperties material{BoundVector::Zero(), 0., true};
  material.invCovarianceQOverP() = 4.;
  const std::unordered_map<GeometryIdentifier, Gx2fMaterialProperties>
      materials{{geoId, material}};
  const auto logger = getDefaultLogger("Gx2fComponentTests", Logging::INFO);
  for (const bool fitEnergyLoss : {false, true}) {
    const Gx2fParameterLayout layout{false, fitEnergyLoss, 1};
    Gx2fSystem system{layout};
    std::vector<GeometryIdentifier> ids;
    fillGx2fSystem(track, system, materials, ids, *logger);
    BOOST_CHECK_EQUAL(system.aMatrix()(eBoundQOverP, eBoundQOverP), 36.);
    BOOST_CHECK_EQUAL(system.bVector()[eBoundQOverP], 12.);
    BOOST_CHECK_EQUAL(system.chi2(), 4.);
    if (fitEnergyLoss) {
      const auto offset = layout.energyLossOffset(0);
      BOOST_REQUIRE_EQUAL(ids.size(), 1u);
      BOOST_CHECK(ids.front() == geoId);
      BOOST_CHECK_EQUAL(system.aMatrix()(eBoundQOverP, offset), 6.);
      BOOST_CHECK_EQUAL(system.aMatrix()(offset, offset), 5.);
      BOOST_CHECK_EQUAL(system.bVector()[offset], 2.);
    } else {
      BOOST_CHECK(ids.empty());
    }
  }
}

BOOST_AUTO_TEST_CASE(Gx2fQOverPOffsetUnsupportedHypotheses) {
  const MaterialSlab slab{makeSilicon(), 5_mm};
  for (const auto& particle :
       {ParticleHypothesis::pion0(), ParticleHypothesis::chargedGeantino(),
        ParticleHypothesis::muon().withMomentumHypothesis(1_GeV)}) {
    BOOST_CHECK_EQUAL(
        computeGx2fQOverPOffset(slab, particle, 1. / 1_GeV,
                                Direction::Forward(), Gx2fEnergyLossMode::Mean),
        0.);
  }
  BOOST_CHECK_EQUAL(
      computeGx2fQOverPOffset(slab, ParticleHypothesis::muon(), 0.,
                              Direction::Forward(), Gx2fEnergyLossMode::Mean),
      0.);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

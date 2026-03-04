// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Fatras/FatrasSimulation.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsFatras/Kernel/InteractionList.hpp"
#include "ActsFatras/Kernel/Simulation.hpp"
#include "ActsFatras/Physics/Decay/NoDecay.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/PhotonConversion.hpp"
#include "ActsFatras/Physics/StandardInteractions.hpp"
#include "ActsFatras/Selectors/SelectorHelpers.hpp"
#include "ActsFatras/Selectors/SurfaceSelectors.hpp"

#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ActsExamples {

namespace {

/// Simple struct to select surfaces where hits should be generated.
struct HitSurfaceSelector {
  bool sensitive = false;
  bool material = false;
  bool passive = false;

  /// Check if the surface should be used.
  bool operator()(const Acts::Surface &surface) const {
    // sensitive/material are not mutually exclusive
    bool isSensitive = surface.isSensitive();
    bool isMaterial = surface.surfaceMaterial() != nullptr;
    // passive should be an orthogonal category
    bool isPassive = !(isSensitive || isMaterial);
    return (isSensitive && sensitive) || (isMaterial && material) ||
           (isPassive && passive);
  }
};

}  // namespace

// Same interface as `ActsFatras::Simulation` but with concrete types.
struct detail::FatrasSimulation {
  virtual ~FatrasSimulation() = default;
  virtual Acts::Result<std::vector<ActsFatras::FailedParticle>> simulate(
      const Acts::GeometryContext &geoCtx,
      const Acts::MagneticFieldContext &magCtx, RandomEngine &rng,
      const SimParticleContainer &inputParticles,
      SimParticleContainer &simulatedParticles,
      SimHitContainer &hits) const = 0;
};

namespace {

// Magnetic-field specific PIMPL implementation.
//
// This always uses the SympyStepper for charged particle propagation and is
// thus limited to propagation in vacuum at the moment.
struct FatrasSimulationT final : detail::FatrasSimulation {
  using CutPMin = ActsFatras::Min<ActsFatras::Casts::P>;

  // typedefs for charge particle simulation
  // propagate charged particles numerically in the given magnetic field
  using ChargedStepper = Acts::SympyStepper;
  using ChargedPropagator = Acts::Propagator<ChargedStepper, Acts::Navigator>;
  // charged particles w/ standard em physics list and selectable hits
  using ChargedSelector = CutPMin;
  using ChargedSimulation = ActsFatras::SingleParticleSimulation<
      ChargedPropagator, ActsFatras::StandardChargedElectroMagneticInteractions,
      HitSurfaceSelector, ActsFatras::NoDecay>;

  // typedefs for neutral particle simulation
  // propagate neutral particles with just straight lines
  using NeutralStepper = Acts::StraightLineStepper;
  using NeutralPropagator = Acts::Propagator<NeutralStepper, Acts::Navigator>;
  // neutral particles w/ photon conversion and no hits
  using NeutralSelector = CutPMin;
  using NeutralInteractions =
      ActsFatras::InteractionList<ActsFatras::PhotonConversion>;
  using NeutralSimulation = ActsFatras::SingleParticleSimulation<
      NeutralPropagator, NeutralInteractions, ActsFatras::NoSurface,
      ActsFatras::NoDecay>;

  // combined simulation type
  using Simulation = ActsFatras::Simulation<ChargedSelector, ChargedSimulation,
                                            NeutralSelector, NeutralSimulation>;

  Simulation simulation;

  FatrasSimulationT(const ActsExamples::FatrasSimulation::Config &cfg,
                    const Acts::Logger &logger)
      : simulation(
            ChargedSimulation(
                ChargedPropagator(ChargedStepper(cfg.magneticField),
                                  Acts::Navigator({cfg.trackingGeometry},
                                                  logger.clone("SimNav")),
                                  logger.clone("SimProp")),
                logger.clone("Simulation")),
            NeutralSimulation(
                NeutralPropagator(NeutralStepper(),
                                  Acts::Navigator({cfg.trackingGeometry},
                                                  logger.clone("SimNav")),
                                  logger.clone("SimProp")),
                logger.clone("Simulation"))) {
    using namespace ActsFatras;
    using namespace ActsFatras::detail;
    // apply the configuration

    // minimal p cut on input particles and as is-alive check for interactions
    simulation.selectCharged.valMin = cfg.pMin;
    simulation.selectNeutral.valMin = cfg.pMin;
    simulation.charged.interactions =
        makeStandardChargedElectroMagneticInteractions(cfg.pMin);

    // processes are enabled by default
    if (!cfg.emScattering) {
      simulation.charged.interactions.template disable<StandardScattering>();
    }
    if (!cfg.emEnergyLossIonisation) {
      simulation.charged.interactions.template disable<StandardBetheBloch>();
    }
    if (!cfg.emEnergyLossRadiation) {
      simulation.charged.interactions.template disable<StandardBetheHeitler>();
    }
    if (!cfg.emPhotonConversion) {
      simulation.neutral.interactions.template disable<PhotonConversion>();
    }

    // configure hit surfaces for charged particles
    simulation.charged.selectHitSurface.sensitive = cfg.generateHitsOnSensitive;
    simulation.charged.selectHitSurface.material = cfg.generateHitsOnMaterial;
    simulation.charged.selectHitSurface.passive = cfg.generateHitsOnPassive;

    simulation.charged.maxStepSize = cfg.maxStepSize;
    simulation.charged.pathLimit = cfg.pathLimit;
    simulation.neutral.maxStepSize = cfg.maxStepSize;
    simulation.neutral.pathLimit = cfg.pathLimit;
  }
  ~FatrasSimulationT() final = default;

  Acts::Result<std::vector<ActsFatras::FailedParticle>> simulate(
      const Acts::GeometryContext &geoCtx,
      const Acts::MagneticFieldContext &magCtx, RandomEngine &rng,
      const SimParticleContainer &inputParticles,
      SimParticleContainer &simulatedParticles,
      SimHitContainer &hits) const final {
    return simulation.simulate(geoCtx, magCtx, rng, inputParticles,
                               simulatedParticles, hits);
  }
};

}  // namespace

FatrasSimulation::FatrasSimulation(Config cfg,
                                   std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("FatrasSimulation", std::move(logger)), m_cfg(std::move(cfg)) {
  ACTS_LOG_WITH_LOGGER(
      this->logger(), Acts::Logging::DEBUG,
      "hits on sensitive surfaces: " << m_cfg.generateHitsOnSensitive);
  ACTS_LOG_WITH_LOGGER(
      this->logger(), Acts::Logging::DEBUG,
      "hits on material surfaces: " << m_cfg.generateHitsOnMaterial);
  ACTS_LOG_WITH_LOGGER(
      this->logger(), Acts::Logging::DEBUG,
      "hits on passive surfaces: " << m_cfg.generateHitsOnPassive);

  if (!m_cfg.generateHitsOnSensitive && !m_cfg.generateHitsOnMaterial &&
      !m_cfg.generateHitsOnPassive) {
    ACTS_LOG_WITH_LOGGER(
        this->logger(), Acts::Logging::WARNING,
        "FatrasSimulation not configured to generate any hits!");
  }

  if (!m_cfg.trackingGeometry) {
    throw std::invalid_argument{"Missing tracking geometry"};
  }
  if (!m_cfg.magneticField) {
    throw std::invalid_argument{"Missing magnetic field"};
  }
  if (!m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers tool");
  }

  // construct the simulation for the specific magnetic field
  m_sim = std::make_unique<FatrasSimulationT>(m_cfg, this->logger());

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_outputParticles.initialize(m_cfg.outputParticles);
  m_outputSimHits.initialize(m_cfg.outputSimHits);
}

// explicit destructor needed for the PIMPL implementation to work
FatrasSimulation::~FatrasSimulation() = default;

ProcessCode FatrasSimulation::execute(const AlgorithmContext &ctx) const {
  // read input containers
  const auto &inputParticles = m_inputParticles(ctx);

  ACTS_DEBUG(inputParticles.size() << " input particles");

  // prepare output containers
  SimParticleContainer simulatedParticles;
  SimHitContainer hits;
  // reserve appropriate resources
  simulatedParticles.reserve(inputParticles.size());
  hits.reserve(inputParticles.size() * m_cfg.averageHitsPerParticle);

  // run the simulation w/ a local random generator
  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);
  auto ret = m_sim->simulate(ctx.geoContext, ctx.magFieldContext, rng,
                             inputParticles, simulatedParticles, hits);
  // fatal error leads to panic
  if (!ret.ok()) {
    ACTS_FATAL("event " << ctx.eventNumber << " simulation failed with error "
                        << ret.error().message());
    return ProcessCode::ABORT;
  }
  // failed particles are just logged. assumes that failed particles are due
  // to edge-cases representing a tiny fraction of the event; not due to a
  // fundamental issue.
  for (const auto &failed : ret.value()) {
    ACTS_ERROR("event " << ctx.eventNumber << " particle " << failed.particleId
                        << " failed to simulate with error " << failed.error
                        << ": " << failed.error.message());
  }

  ACTS_DEBUG(simulatedParticles.size() << " simulated particles");
  ACTS_DEBUG(hits.size() << " simulated hits");

  // store ordered output containers
  m_outputParticles(ctx, std::move(simulatedParticles));
  m_outputSimHits(ctx, std::move(hits));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

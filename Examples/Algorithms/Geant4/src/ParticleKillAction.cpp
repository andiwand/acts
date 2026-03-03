// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/ParticleKillAction.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Geant4/AlgebraConverters.hpp"

#include <ostream>
#include <utility>

#include <G4RunManager.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4Track.hh>
#include <G4UnitsTable.hh>
#include <G4VPhysicalVolume.hh>

namespace ActsExamples::Geant4 {

constexpr double convertTime = Acts::UnitConstants::ns / CLHEP::ns;

ParticleKillAction::ParticleKillAction(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4UserSteppingAction(), m_cfg(cfg), m_logger(std::move(logger)) {}

void ParticleKillAction::UserSteppingAction(const G4Step* stepPtr) {
  assert(stepPtr != nullptr);
  const G4Step& step = *stepPtr;

  assert(step.GetTrack() != nullptr);

  G4Track& track = *step.GetTrack();

  const double time = convertTime * track.GetGlobalTime();
  const bool isSecondary =
      track.GetDynamicParticle()->GetPrimaryParticle() == nullptr;

  const bool outOfVolume =
      m_cfg.volume &&
      !m_cfg.volume->inside(eventStore().geoContext,
                            convertPosition(track.GetPosition()));
  const bool outOfTime = time > m_cfg.maxTime;
  const bool invalidSecondary = m_cfg.secondaries && isSecondary;

  if (outOfVolume || outOfTime || invalidSecondary) {
    ACTS_DEBUG("Kill track with internal track ID "
               << track.GetTrackID() << " at "
               << convertPosition(track.GetPosition()) << " and global time "
               << time / Acts::UnitConstants::ns << "ns and isSecondary "
               << isSecondary);
    track.SetTrackStatus(G4TrackStatus::fStopAndKill);
  }

  // store the outcome of the particle
  // check if we have a particle assigned to track
  if (auto trackIt = eventStore().trackIdMapping.find(track.GetTrackID());
      trackIt != eventStore().trackIdMapping.end()) {
    const SimParticleIndex simParticleIndex = trackIt->second;
    MutableSimParticle simParticle =
        eventStore().particles.at(simParticleIndex);

    using enum ActsFatras::SimulationOutcome;

    // set the outcome of the particle
    if (outOfVolume) {
      simParticle.simulationOutcome() = KilledVolumeExit;
    } else if (outOfTime) {
      simParticle.simulationOutcome() = KilledTime;
    } else if (invalidSecondary) {
      simParticle.simulationOutcome() = KilledSecondaryParticle;
    } else if (track.GetTrackStatus() == G4TrackStatus::fStopAndKill) {
      simParticle.simulationOutcome() = KilledInteraction;
    }
  } else {
    ACTS_DEBUG("Track ID " << track.GetTrackID() << " not found");
    ++eventStore().particleIdNotFound;
  }
}

}  // namespace ActsExamples::Geant4

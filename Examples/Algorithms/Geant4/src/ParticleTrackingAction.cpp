// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/ParticleTrackingAction.hpp"

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"
#include "ActsFatras/EventData/ParticleContainer.hpp"

#include <cassert>
#include <ostream>
#include <unordered_map>
#include <utility>

#include <G4ParticleDefinition.hh>
#include <G4RunManager.hh>
#include <G4Track.hh>
#include <G4UnitsTable.hh>

namespace ActsExamples::Geant4 {

namespace {

// Unit conversions G4->::ACTS
constexpr double convertTime = Acts::UnitConstants::ns / CLHEP::ns;
constexpr double convertLength = Acts::UnitConstants::mm / CLHEP::mm;
constexpr double convertEnergy = Acts::UnitConstants::GeV / CLHEP::GeV;

template <typename StateAccessor>
void convertState(const G4Track& track,
                  ActsFatras::ParticleStateProxy<StateAccessor, false> state) {
  const G4ThreeVector pPosition = convertLength * track.GetPosition();
  const G4double pTime = convertTime * track.GetGlobalTime();
  const G4ThreeVector pDirection = track.GetMomentumDirection();
  const G4double p = convertEnergy * track.GetKineticEnergy();

  state.fourPosition() =
      Acts::Vector4(pPosition.x(), pPosition.y(), pPosition.z(), pTime);
  state.absoluteMomentum() = p;
  state.direction() =
      Acts::Vector3(pDirection.x(), pDirection.y(), pDirection.z());
}

}  // namespace

ParticleTrackingAction::ParticleTrackingAction(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4UserTrackingAction(), m_cfg(cfg), m_logger(std::move(logger)) {}

void ParticleTrackingAction::PreUserTrackingAction(const G4Track* trackPtr) {
  assert(trackPtr != nullptr);
  const G4Track& track = *trackPtr;

  // If this is not the case, there are unhandled cases of particle stopping in
  // the SensitiveSteppingAction
  // TODO We could also merge the remaining hits to a hit here, but it would be
  // nicer to investigate, if we can handle all particle stop conditions in the
  // SensitiveSteppingAction... This seems to happen O(1) times in a ttbar
  // event, so seems not to be too problematic
  if (!eventStore().hitBuffer.empty()) {
    ACTS_WARNING("Hit buffer not empty after track. Clearing ...");
    eventStore().hitBuffer.clear();
  }

  const std::optional<SimBarcode> barcode =
      makeBarcode(track.GetTrackID(), track.GetParentID());

  // There is already a warning printed in the makeBarcode function if this
  // indicates a failure
  if (!barcode) {
    return;
  }

  MutableSimParticle simParticle = convertGenerated(track, *barcode);
  eventStore().trackIdMapping[track.GetTrackID()] = simParticle.index();
}

void ParticleTrackingAction::PostUserTrackingAction(const G4Track* trackPtr) {
  assert(trackPtr != nullptr);
  const G4Track& track = *trackPtr;

  // The initial particle maybe was not registered because of a particle ID
  // collision
  if (!eventStore().trackIdMapping.contains(track.GetTrackID())) {
    ACTS_WARNING("Particle ID for track ID " << track.GetTrackID()
                                             << " not registered. Skip");
    return;
  }

  const SimParticleIndex simParticleId =
      eventStore().trackIdMapping.at(track.GetTrackID());
  MutableSimParticle simParticle = eventStore().particles.at(simParticleId);
  convertFinal(track, simParticle);
}

std::optional<SimBarcode> ParticleTrackingAction::makeBarcode(
    G4int trackId, G4int parentId) const {
  // We already have this particle registered (it is one of the input particles
  // or we are making a final particle state)
  if (eventStore().trackIdMapping.contains(trackId)) {
    return std::nullopt;
  }

  if (!eventStore().trackIdMapping.contains(parentId)) {
    ACTS_DEBUG("Parent track ID " << parentId
                                  << " not found, cannot build barcode");
    ++eventStore().parentIdNotFound;
    return std::nullopt;
  }

  const SimParticleIndex simParticleId =
      eventStore().trackIdMapping.at(parentId);
  const SimParticle parentParticle =
      eventStore().particles.at(simParticleId).asConst();
  SimBarcode barcode = parentParticle.barcode().makeDescendant();
  const SimBarcode key = parentParticle.barcode().withoutSubparticle();
  ++eventStore().subparticleMap[key];
  barcode = barcode.withSubParticle(eventStore().subparticleMap[key]);

  return barcode;
}

MutableSimParticle ParticleTrackingAction::convertGenerated(
    const G4Track& track, const SimBarcode& barcode) const {
  // Get all the information from the Track
  const G4ParticleDefinition* particleDef = track.GetParticleDefinition();
  const G4int pdg = particleDef->GetPDGEncoding();
  const G4double charge = particleDef->GetPDGCharge();
  const G4double mass = convertEnergy * particleDef->GetPDGMass();

  // Now create the Particle
  MutableSimParticle simParticle = eventStore().particles.createParticle();
  simParticle.mass() = mass;
  simParticle.charge() = charge;
  simParticle.barcode() = barcode;
  simParticle.pdg() = Acts::PdgParticle{pdg};

  convertState(track, simParticle.generationState());

  return simParticle;
}

void ParticleTrackingAction::convertFinal(
    const G4Track& track, MutableSimParticle simParticle) const {
  convertState(track, simParticle.simulationState());
}

}  // namespace ActsExamples::Geant4

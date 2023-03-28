// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <memory>
#include <mutex>
#include <string>

class TFile;
class TTree;

namespace ActsExamples {

/// Write ambiguity resolution performance measures.
class AmbiguityResolutionPerformanceWriter final
    : public WriterT<ConstTrackContainer> {
 public:
  struct Config {
    /// Input reconstructed proto tracks collection.
    std::string inputTracks;
    /// Input particles collection.
    std::string inputParticles;
    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;
    /// Output filename.
    std::string filePath = "performance_ambiguity_resolution.root";
    /// Output file mode
    std::string fileMode = "RECREATE";
    /// Output tree name
    std::string treeName = "ambiguity_resolution";
  };

  /// Constructor
  /// @param config the configuration
  /// @param level The log level
  AmbiguityResolutionPerformanceWriter(Config config,
                                       Acts::Logging::Level level);
  ~AmbiguityResolutionPerformanceWriter();

  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  using SimParticleContainer = ActsExamples::SimParticleContainer;
  using HitParticlesMap = ActsExamples::IndexMultimap<ActsFatras::Barcode>;

  ProcessCode writeT(const AlgorithmContext& ctx,
                     const ConstTrackContainer& tracks) override;

  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<HitParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};

  std::mutex m_writeMutex;
  TFile* m_outputFile = nullptr;
  TTree* m_outputTree = nullptr;

  /// Event identifier.
  uint32_t m_eventId = 0;
  /// Event-unique particle identifier a.k.a barcode.
  std::vector<uint64_t> m_particleId;
  /// Particle type a.k.a. PDG particle number
  std::vector<int32_t> m_particleType;
  /// Production process type, i.e. what generated the particle.
  std::vector<uint32_t> m_process;
  /// Production position components in mm.
  std::vector<float> m_vx;
  std::vector<float> m_vy;
  std::vector<float> m_vz;
  // Production time in ns.
  std::vector<float> m_vt;
  /// Total momentum in GeV
  std::vector<float> m_p;
  /// Momentum components in GeV.
  std::vector<float> m_px;
  std::vector<float> m_py;
  std::vector<float> m_pz;
  /// Mass in GeV.
  std::vector<float> m_m;
  /// Charge in e.
  std::vector<float> m_q;
  // Derived kinematic quantities
  /// Direction pseudo-rapidity.
  std::vector<float> m_eta;
  /// Direction angle in the transverse plane.
  std::vector<float> m_phi;
  /// Transverse momentum in GeV.
  std::vector<float> m_pt;
  // Decoded particle identifier; see Barcode definition for details.
  std::vector<uint32_t> m_vertexPrimary;
  std::vector<uint32_t> m_vertexSecondary;
  std::vector<uint32_t> m_particle;
  std::vector<uint32_t> m_generation;
  std::vector<uint32_t> m_subParticle;
};

}  // namespace ActsExamples

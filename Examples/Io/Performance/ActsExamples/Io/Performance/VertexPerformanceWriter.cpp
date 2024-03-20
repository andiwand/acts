// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Performance/VertexPerformanceWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/MultiIndex.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ios>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <ostream>
#include <set>
#include <stdexcept>
#include <utility>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

namespace ActsExamples {

VertexPerformanceWriter::VertexPerformanceWriter(
    const VertexPerformanceWriter::Config& config, Acts::Logging::Level level)
    : WriterT(config.inputVertices, "VertexPerformanceWriter", level),
      m_cfg(config) {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }
  if (m_cfg.inputTruthVertices.empty()) {
    throw std::invalid_argument("Collection with truth vertices missing");
  }
  if (m_cfg.inputAllTruthParticles.empty()) {
    throw std::invalid_argument("Collection with all truth particles missing");
  }
  if (m_cfg.inputSelectedTruthParticles.empty()) {
    throw std::invalid_argument(
        "Collection with selected truth particles missing");
  }

  m_inputTruthVertices.initialize(m_cfg.inputTruthVertices);
  m_inputAllTruthParticles.initialize(m_cfg.inputAllTruthParticles);
  m_inputSelectedTruthParticles.initialize(m_cfg.inputSelectedTruthParticles);

  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_inputTracks.initialize(m_cfg.inputTracks);

  // Setup ROOT I/O
  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + path + "'");
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  m_outputTree->Branch("event_nr", &m_eventNr);

  m_outputTree->Branch("nRecoVtx", &m_nRecoVtx);
  m_outputTree->Branch("nTrueVtx", &m_nTrueVtx);
  m_outputTree->Branch("nMergedVtx", &m_nMergedVtx);
  m_outputTree->Branch("nSplitVtx", &m_nSplitVtx);
  m_outputTree->Branch("nVtxDetectorAcceptance", &m_nVtxDetAcceptance);
  m_outputTree->Branch("nVtxReconstructable", &m_nVtxReconstructable);

  m_outputTree->Branch("recoX", &m_recoX);
  m_outputTree->Branch("recoY", &m_recoY);
  m_outputTree->Branch("recoZ", &m_recoZ);
  m_outputTree->Branch("recoT", &m_recoT);

  m_outputTree->Branch("covXX", &m_covXX);
  m_outputTree->Branch("covYY", &m_covYY);
  m_outputTree->Branch("covZZ", &m_covZZ);
  m_outputTree->Branch("covTT", &m_covTT);
  m_outputTree->Branch("covXY", &m_covXY);
  m_outputTree->Branch("covXZ", &m_covXZ);
  m_outputTree->Branch("covXT", &m_covXT);
  m_outputTree->Branch("covYZ", &m_covYZ);
  m_outputTree->Branch("covYT", &m_covYT);
  m_outputTree->Branch("covZT", &m_covZT);

  m_outputTree->Branch("nTracksRecoVtx", &m_nTracksOnRecoVertex);
  m_outputTree->Branch("recoVertexTotalTrackWeight",
                       &m_recoVertexTotalTrackWeight);
  m_outputTree->Branch("sumPt2", &m_sumPt2);

  m_outputTree->Branch("vertex_primary", &m_vertexPrimary);
  m_outputTree->Branch("vertex_secondary", &m_vertexSecondary);

  m_outputTree->Branch("truthVertexMatchRatio", &m_truthVertexMatchRatio);

  m_outputTree->Branch("nTracksTruthVtx", &m_nTracksOnTruthVertex);

  m_outputTree->Branch("recoVertexClassification", &m_recoVertexClassification);

  m_outputTree->Branch("truthX", &m_truthX);
  m_outputTree->Branch("truthY", &m_truthY);
  m_outputTree->Branch("truthZ", &m_truthZ);
  m_outputTree->Branch("truthT", &m_truthT);

  m_outputTree->Branch("resX", &m_resX);
  m_outputTree->Branch("resY", &m_resY);
  m_outputTree->Branch("resZ", &m_resZ);
  m_outputTree->Branch("resT", &m_resT);

  m_outputTree->Branch("pullX", &m_pullX);
  m_outputTree->Branch("pullY", &m_pullY);
  m_outputTree->Branch("pullZ", &m_pullZ);
  m_outputTree->Branch("pullT", &m_pullT);

  m_outputTree->Branch("trk_weight", &m_trkWeight);

  m_outputTree->Branch("trk_recoPhi", &m_recoPhi);
  m_outputTree->Branch("trk_recoTheta", &m_recoTheta);
  m_outputTree->Branch("trk_recoQOverP", &m_recoQOverP);
  m_outputTree->Branch("trk_recoPhiFitted", &m_recoPhiFitted);
  m_outputTree->Branch("trk_recoThetaFitted", &m_recoThetaFitted);
  m_outputTree->Branch("trk_recoQOverPFitted", &m_recoQOverPFitted);

  m_outputTree->Branch("trk_particleId", &m_particleId);

  m_outputTree->Branch("trk_truthPhi", &m_truthPhi);
  m_outputTree->Branch("trk_truthTheta", &m_truthTheta);
  m_outputTree->Branch("trk_truthQOverP", &m_truthQOverP);

  m_outputTree->Branch("trk_resPhi", &m_resPhi);
  m_outputTree->Branch("trk_resTheta", &m_resTheta);
  m_outputTree->Branch("trk_resQOverP", &m_resQOverP);
  m_outputTree->Branch("trk_momOverlap", &m_momOverlap);
  m_outputTree->Branch("trk_resPhiFitted", &m_resPhiFitted);
  m_outputTree->Branch("trk_resThetaFitted", &m_resThetaFitted);
  m_outputTree->Branch("trk_resQOverPFitted", &m_resQOverPFitted);
  m_outputTree->Branch("trk_momOverlapFitted", &m_momOverlapFitted);

  m_outputTree->Branch("trk_pullPhi", &m_pullPhi);
  m_outputTree->Branch("trk_pullTheta", &m_pullTheta);
  m_outputTree->Branch("trk_pullQOverP", &m_pullQOverP);
  m_outputTree->Branch("trk_pullPhiFitted", &m_pullPhiFitted);
  m_outputTree->Branch("trk_pullThetaFitted", &m_pullThetaFitted);
  m_outputTree->Branch("trk_pullQOverPFitted", &m_pullQOverPFitted);
}

VertexPerformanceWriter::~VertexPerformanceWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode VertexPerformanceWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  return ProcessCode::SUCCESS;
}

int VertexPerformanceWriter::getNumberOfReconstructableVertices(
    const SimParticleContainer& collection) const {
  // map for finding frequency
  std::map<int, int> fmap;

  std::vector<int> reconstructableTruthVertices;

  // traverse the array for frequency
  for (const auto& p : collection) {
    int secVtxId = p.particleId().vertexSecondary();
    if (secVtxId != 0) {
      // truthparticle from secondary vtx
      continue;
    }
    int priVtxId = p.particleId().vertexPrimary();
    fmap[priVtxId]++;
  }

  // iterate over the map
  for (auto it : fmap) {
    // Require at least 2 tracks
    if (it.second > 1) {
      reconstructableTruthVertices.push_back(it.first);
    }
  }

  return reconstructableTruthVertices.size();
}

int VertexPerformanceWriter::getNumberOfTruePriVertices(
    const SimParticleContainer& collection) const {
  // Vector to store indices of all primary vertices
  std::set<int> allPriVtxIds;
  for (const auto& p : collection) {
    int priVtxId = p.particleId().vertexPrimary();
    int secVtxId = p.particleId().vertexSecondary();
    if (secVtxId != 0) {
      // truthparticle from secondary vtx
      continue;
    }
    // Insert to set, removing duplicates
    allPriVtxIds.insert(priVtxId);
  }
  // Size of set corresponds to total number of primary vertices
  return allPriVtxIds.size();
}

ProcessCode VertexPerformanceWriter::writeT(const AlgorithmContext& ctx,
                                            const VertexContainer& vertices) {
  const double nan = std::numeric_limits<double>::quiet_NaN();

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Read truth particle input collection
  const auto& allTruthParticles = m_inputAllTruthParticles(ctx);

  // Read truth vertex input collection
  const auto& truthVertices = m_inputTruthVertices(ctx);
  // Get number of generated true primary vertices
  m_nTrueVtx = truthVertices.size();

  m_nRecoVtx = vertices.size();
  m_nMergedVtx = 0;
  m_nSplitVtx = 0;

  ACTS_VERBOSE("Total number of generated truth particles in event : "
               << allTruthParticles.size());
  ACTS_VERBOSE(
      "Total number of generated truth primary vertices : " << m_nTrueVtx);
  ACTS_DEBUG("Number of reco vertices in event: " << m_nRecoVtx);

  // Read selected truth particle input collection
  const auto& selectedTruthParticles = m_inputSelectedTruthParticles(ctx);
  // Get number of detector-accepted true primary vertices
  m_nVtxDetAcceptance = getNumberOfTruePriVertices(selectedTruthParticles);

  ACTS_VERBOSE("Total number of selected truth particles in event : "
               << selectedTruthParticles.size());
  ACTS_VERBOSE("Total number of detector-accepted truth primary vertices : "
               << m_nVtxDetAcceptance);

  TrackParametersContainer trackParameters;
  std::vector<SimParticle> associatedTruthParticles;

  // Get the event number
  m_eventNr = ctx.eventNumber;

  // If we don't know which truth particle corresponds to which track a
  // priori, we check how many hits particles and tracks share. We match the
  // particle to the track if a fraction of more than truthMatchProbMin of
  // hits that contribute to the track come from the particle. Note that not
  // all tracksatVertex have matching parameters in trackParameters in this
  // case. Equivalently, one could say that not all tracksAtVertex will be
  // assigned to a truth particle.
  {
    std::vector<ParticleHitCount> particleHitCounts;

    const auto& hitParticlesMap = m_inputMeasurementParticlesMap(ctx);

    for (const auto& track : m_inputTracks(ctx)) {
      if (!track.hasReferenceSurface()) {
        continue;
      }

      identifyContributingParticles(hitParticlesMap, track, particleHitCounts);
      ActsFatras::Barcode majorityParticleId =
          particleHitCounts.front().particleId;
      std::size_t nMajorityHits = particleHitCounts.front().hitCount;

      if ((double)nMajorityHits / track.nMeasurements() <
          m_cfg.trackMatchThreshold) {
        continue;
      }

      auto it = std::find_if(allTruthParticles.begin(), allTruthParticles.end(),
                             [&](const auto& tp) {
                               return tp.particleId() == majorityParticleId;
                             });

      if (it == allTruthParticles.end()) {
        continue;
      }

      trackParameters.push_back(track.createParametersAtReference());
      const auto& majorityParticle = *it;
      associatedTruthParticles.push_back(majorityParticle);
    }
  }

  if (trackParameters.size() != associatedTruthParticles.size()) {
    ACTS_ERROR(
        "Number of track parameters and associated truth particles do not "
        "match. ("
        << trackParameters.size() << " != " << associatedTruthParticles.size()
        << ").");
  }

  // Get number of track-associated true primary vertices
  m_nVtxReconstructable =
      getNumberOfReconstructableVertices(SimParticleContainer(
          associatedTruthParticles.begin(), associatedTruthParticles.end()));

  ACTS_INFO(
      "Total number of reconstructed tracks : " << trackParameters.size());
  ACTS_INFO("Total number of reco track-associated truth particles in event : "
            << associatedTruthParticles.size());
  ACTS_INFO("Total number of reco track-associated truth primary vertices : "
            << m_nVtxReconstructable);

  // We compare the reconstructed momenta to the true momenta at the vertex. For
  // this, we propagate the reconstructed tracks to the PCA of the true vertex
  // position. Setting up propagator:
  Acts::EigenStepper<> stepper(m_cfg.bField);
  using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
  auto propagator = std::make_shared<Propagator>(stepper);

  auto findTruthParticle =
      [&](const auto& tracksAtVtx) -> std::optional<SimParticle> {
    // Track parameters before the vertex fit
    Acts::BoundTrackParameters origTrack =
        *(tracksAtVtx.originalParams.template as<Acts::BoundTrackParameters>());

    // Finding the matching parameters in the container of all track
    // parameters. This allows us to identify the corresponding particle,
    // since we expect trackParameters and associatedTruthParticles to
    // align.
    for (std::size_t i = 0; i < trackParameters.size(); ++i) {
      const auto& params = trackParameters[i].parameters();
      if (origTrack.parameters() == params) {
        return associatedTruthParticles[i];
      }
    }

    return std::nullopt;
  };

  // Get number of contributing tracks (i.e., tracks with a weight above
  // threshold)
  auto weightHighEnough = [this](const Acts::TrackAtVertex& trkAtVtx) {
    return trkAtVtx.trackWeight > m_cfg.minTrkWeight;
  };

  auto calcSumPt2 = [this](const Acts::Vertex& vtx) {
    double sumPt2 = 0;
    for (const auto& trk : vtx.tracks()) {
      if (trk.trackWeight > m_cfg.minTrkWeight) {
        double pt = trk.originalParams.as<Acts::BoundTrackParameters>()
                        ->transverseMomentum();
        sumPt2 += pt * pt;
      }
    }
    return sumPt2;
  };

  struct ToTruthMatching {
    std::optional<SimVertexBarcode> vertexId;
    double totalTrackWeight{};
    double matchFraction{};

    RecoVertexClassification classification{RecoVertexClassification::Unknown};
  };
  struct ToRecoMatching {
    std::size_t recoIndex{};

    double recoSumPt2{};
  };

  std::vector<ToTruthMatching> recoToTruthMatching;
  std::map<SimVertexBarcode, ToRecoMatching> truthToRecoMatching;

  // Do truth matching for each reconstructed vertex
  for (const auto& [vtxIndex, vtx] : Acts::enumerate(vertices)) {
    // Reconstructed tracks that contribute to the reconstructed vertex
    const auto& tracksAtVtx = vtx.tracks();

    // Containers for storing truth particles and truth vertices that
    // contribute
    // to the reconstructed vertex
    std::vector<std::pair<SimVertexBarcode, double>> contributingTruthVertices;

    double totalTrackWeight = 0;
    for (const auto& trk : tracksAtVtx) {
      if (trk.trackWeight < m_cfg.minTrkWeight) {
        continue;
      }

      totalTrackWeight += trk.trackWeight;

      std::optional<SimParticle> particleOpt = findTruthParticle(trk);
      if (!particleOpt.has_value()) {
        ACTS_VERBOSE("Track has no matching truth particle.");
      } else {
        contributingTruthVertices.emplace_back(
            SimBarcode(particleOpt->particleId()).vertexId(), trk.trackWeight);
      }
    }

    // Find true vertex that contributes most to the reconstructed vertex
    std::map<SimVertexBarcode, std::pair<int, double>> fmap;
    for (const auto& [vtxId, weight] : contributingTruthVertices) {
      ++fmap[vtxId].first;
      fmap[vtxId].second += weight;
    }
    double maxOccurrence = -1;
    SimVertexBarcode maxOccurrenceId = -1;
    for (const auto& [vtxId, counter] : fmap) {
      if (counter.second > maxOccurrence) {
        maxOccurrenceId = vtxId;
        maxOccurrence = counter.second;
      }
    }

    double sumPt2 = calcSumPt2(vtx);

    double vertexMatchFraction =
        fmap[maxOccurrenceId].second / totalTrackWeight;
    RecoVertexClassification recoVertexClassification =
        RecoVertexClassification::Unknown;

    if (vertexMatchFraction >= m_cfg.vertexMatchThreshold) {
      recoVertexClassification = RecoVertexClassification::Clean;
    } else {
      recoVertexClassification = RecoVertexClassification::Merged;
    }

    recoToTruthMatching.push_back({maxOccurrenceId, totalTrackWeight,
                                   vertexMatchFraction,
                                   recoVertexClassification});

    auto& recoToTruth = recoToTruthMatching.back();

    // We have to decide if this reco vertex is a split vertex.
    if (auto it = truthToRecoMatching.find(maxOccurrenceId);
        it != truthToRecoMatching.end()) {
      // This truth vertex is already matched to a reconstructed vertex so we
      // are dealing with a split vertex.

      // We have to decide which of the two reconstructed vertices is the
      // split vertex.
      if (sumPt2 <= it->second.recoSumPt2) {
        // Since the sumPt2 is smaller we can simply call this a split vertex

        recoToTruth.classification = RecoVertexClassification::Split;

        // Keep the existing truth to reco matching
      } else {
        // The sumPt2 is larger, so we call the other vertex a split vertex.

        auto& otherRecoToTruth = recoToTruthMatching.at(it->second.recoIndex);
        // Swap the classification
        recoToTruth.classification = otherRecoToTruth.classification;
        otherRecoToTruth.classification = RecoVertexClassification::Split;

        // Overwrite the truth to reco matching
        it->second = {vtxIndex, sumPt2};
      }
    } else {
      truthToRecoMatching[maxOccurrenceId] = {vtxIndex, sumPt2};
    }
  }

  // Loop over reconstructed vertices and see if they can be matched to a true
  // vertex.
  for (const auto& [vtxIndex, vtx] : Acts::enumerate(vertices)) {
    // Reconstructed tracks that contribute to the reconstructed vertex
    const auto& tracksAtVtx = vtx.tracks();

    const auto& toTruthMatching = recoToTruthMatching[vtxIndex];

    // Helper function for computing the pull
    auto pull =
        [&](const Acts::ActsScalar& diff, const Acts::ActsScalar& variance,
            const std::string& variableStr, const bool& afterFit = true) {
          if (variance <= 0) {
            std::string tempStr;
            if (afterFit) {
              tempStr = "after";
            } else {
              tempStr = "before";
            }
            ACTS_WARNING("Nonpositive variance "
                         << tempStr << " vertex fit: Var(" << variableStr
                         << ") = " << variance << " <= 0.");
            return std::numeric_limits<double>::quiet_NaN();
          }
          double std = std::sqrt(variance);
          return diff / std;
        };

    m_recoX.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreePos0]);
    m_recoY.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreePos1]);
    m_recoZ.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreePos2]);
    m_recoT.push_back(vtx.fullPosition()[Acts::FreeIndices::eFreeTime]);

    Acts::ActsScalar varX = vtx.fullCovariance()(Acts::FreeIndices::eFreePos0,
                                                 Acts::FreeIndices::eFreePos0);
    Acts::ActsScalar varY = vtx.fullCovariance()(Acts::FreeIndices::eFreePos1,
                                                 Acts::FreeIndices::eFreePos1);
    Acts::ActsScalar varZ = vtx.fullCovariance()(Acts::FreeIndices::eFreePos2,
                                                 Acts::FreeIndices::eFreePos2);
    Acts::ActsScalar varTime = vtx.fullCovariance()(
        Acts::FreeIndices::eFreeTime, Acts::FreeIndices::eFreeTime);

    m_covXX.push_back(varX);
    m_covYY.push_back(varY);
    m_covZZ.push_back(varZ);
    m_covTT.push_back(varTime);
    m_covXY.push_back(vtx.fullCovariance()(Acts::FreeIndices::eFreePos0,
                                           Acts::FreeIndices::eFreePos1));
    m_covXZ.push_back(vtx.fullCovariance()(Acts::FreeIndices::eFreePos0,
                                           Acts::FreeIndices::eFreePos2));
    m_covXT.push_back(vtx.fullCovariance()(Acts::FreeIndices::eFreePos0,
                                           Acts::FreeIndices::eFreeTime));
    m_covYZ.push_back(vtx.fullCovariance()(Acts::FreeIndices::eFreePos1,
                                           Acts::FreeIndices::eFreePos2));
    m_covYT.push_back(vtx.fullCovariance()(Acts::FreeIndices::eFreePos1,
                                           Acts::FreeIndices::eFreeTime));
    m_covZT.push_back(vtx.fullCovariance()(Acts::FreeIndices::eFreePos2,
                                           Acts::FreeIndices::eFreeTime));

    double sumPt2 = calcSumPt2(vtx);
    m_sumPt2.push_back(sumPt2);

    double totalTrackWeight = toTruthMatching.totalTrackWeight;
    m_recoVertexTotalTrackWeight.push_back(totalTrackWeight);

    RecoVertexClassification recoVertexClassification =
        toTruthMatching.classification;
    m_recoVertexClassification.push_back(
        static_cast<int>(recoVertexClassification));

    if (recoVertexClassification == RecoVertexClassification::Merged) {
      ++m_nMergedVtx;
    } else if (recoVertexClassification == RecoVertexClassification::Split) {
      ++m_nSplitVtx;
    }

    unsigned int nTracksOnRecoVertex =
        std::count_if(tracksAtVtx.begin(), tracksAtVtx.end(), weightHighEnough);
    m_nTracksOnRecoVertex.push_back(nTracksOnRecoVertex);

    // Saving truth information for the reconstructed vertex
    bool truthInfoWritten = false;
    std::optional<Acts::Vector4> truthPos;
    if (toTruthMatching.vertexId.has_value()) {
      auto iTruthVertex = truthVertices.find(toTruthMatching.vertexId.value());
      if (iTruthVertex == truthVertices.end()) {
        ACTS_ERROR("Truth vertex not found.");
      } else {
        const SimVertex& truthVertex = *iTruthVertex;

        // Count number of reconstructible tracks on truth vertex
        int nTracksOnTruthVertex = 0;
        for (const auto& particle : selectedTruthParticles) {
          if (particle.particleId().vertexId() == truthVertex.vertexId()) {
            ++nTracksOnTruthVertex;
          }
        }
        m_nTracksOnTruthVertex.push_back(nTracksOnTruthVertex);

        m_truthVertexMatchRatio.push_back(toTruthMatching.matchFraction);

        m_vertexPrimary.push_back(truthVertex.vertexId().vertexPrimary());
        m_vertexSecondary.push_back(truthVertex.vertexId().vertexSecondary());

        const Acts::Vector4& truePos = truthVertex.position4;
        truthPos = truePos;
        m_truthX.push_back(truePos[Acts::FreeIndices::eFreePos0]);
        m_truthY.push_back(truePos[Acts::FreeIndices::eFreePos1]);
        m_truthZ.push_back(truePos[Acts::FreeIndices::eFreePos2]);
        m_truthT.push_back(truePos[Acts::FreeIndices::eFreeTime]);

        const Acts::Vector4 diffPos = vtx.fullPosition() - truePos;
        m_resX.push_back(diffPos[Acts::FreeIndices::eFreePos0]);
        m_resY.push_back(diffPos[Acts::FreeIndices::eFreePos1]);
        m_resZ.push_back(diffPos[Acts::FreeIndices::eFreePos2]);
        m_resT.push_back(diffPos[Acts::FreeIndices::eFreeTime]);

        m_pullX.push_back(
            pull(diffPos[Acts::FreeIndices::eFreePos0], varX, "X"));
        m_pullY.push_back(
            pull(diffPos[Acts::FreeIndices::eFreePos1], varY, "Y"));
        m_pullZ.push_back(
            pull(diffPos[Acts::FreeIndices::eFreePos2], varZ, "Z"));
        m_pullT.push_back(
            pull(diffPos[Acts::FreeIndices::eFreeTime], varTime, "T"));

        truthInfoWritten = true;
      }
    }
    if (!truthInfoWritten) {
      m_nTracksOnTruthVertex.push_back(-1);

      m_truthVertexMatchRatio.push_back(-1);

      m_vertexPrimary.push_back(-1);
      m_vertexSecondary.push_back(-1);

      m_truthX.push_back(nan);
      m_truthY.push_back(nan);
      m_truthZ.push_back(nan);
      m_truthT.push_back(nan);

      m_resX.push_back(nan);
      m_resY.push_back(nan);
      m_resZ.push_back(nan);
      m_resT.push_back(nan);

      m_pullX.push_back(nan);
      m_pullY.push_back(nan);
      m_pullZ.push_back(nan);
      m_pullT.push_back(nan);
    }

    // Saving the reconstructed/truth momenta. The reconstructed momenta
    // are taken at the PCA to the truth vertex position -> we need to
    // perform a propagation.

    // Get references to inner vectors where all track variables corresponding
    // to the current vertex will be saved
    auto& innerTrkWeight = m_trkWeight.emplace_back();

    auto& innerRecoPhi = m_recoPhi.emplace_back();
    auto& innerRecoTheta = m_recoTheta.emplace_back();
    auto& innerRecoQOverP = m_recoQOverP.emplace_back();

    auto& innerRecoPhiFitted = m_recoPhiFitted.emplace_back();
    auto& innerRecoThetaFitted = m_recoThetaFitted.emplace_back();
    auto& innerRecoQOverPFitted = m_recoQOverPFitted.emplace_back();

    auto& innerParticleId = m_particleId.emplace_back();

    auto& innerTruthPhi = m_truthPhi.emplace_back();
    auto& innerTruthTheta = m_truthTheta.emplace_back();
    auto& innerTruthQOverP = m_truthQOverP.emplace_back();

    auto& innerResPhi = m_resPhi.emplace_back();
    auto& innerResTheta = m_resTheta.emplace_back();
    auto& innerResQOverP = m_resQOverP.emplace_back();

    auto& innerResPhiFitted = m_resPhiFitted.emplace_back();
    auto& innerResThetaFitted = m_resThetaFitted.emplace_back();
    auto& innerResQOverPFitted = m_resQOverPFitted.emplace_back();

    auto& innerMomOverlap = m_momOverlap.emplace_back();
    auto& innerMomOverlapFitted = m_momOverlapFitted.emplace_back();

    auto& innerPullPhi = m_pullPhi.emplace_back();
    auto& innerPullTheta = m_pullTheta.emplace_back();
    auto& innerPullQOverP = m_pullQOverP.emplace_back();

    auto& innerPullPhiFitted = m_pullPhiFitted.emplace_back();
    auto& innerPullThetaFitted = m_pullThetaFitted.emplace_back();
    auto& innerPullQOverPFitted = m_pullQOverPFitted.emplace_back();

    // Perigee at the true vertex position
    std::shared_ptr<Acts::PerigeeSurface> perigeeSurface;
    if (truthPos.has_value()) {
      perigeeSurface =
          Acts::Surface::makeShared<Acts::PerigeeSurface>(truthPos->head<3>());
    }
    // Lambda for propagating the tracks to the PCA
    auto propagateToVtx =
        [&](const auto& params) -> std::optional<Acts::BoundTrackParameters> {
      if (!perigeeSurface) {
        return std::nullopt;
      }

      auto intersection =
          perigeeSurface
              ->intersect(ctx.geoContext, params.position(ctx.geoContext),
                          params.direction(), Acts::BoundaryCheck(false))
              .closest();

      // Setting the geometry/magnetic field context for the event
      Acts::PropagatorOptions pOptions(ctx.geoContext, ctx.magFieldContext);
      pOptions.direction =
          Acts::Direction::fromScalarZeroAsPositive(intersection.pathLength());

      auto result = propagator->propagate(params, *perigeeSurface, pOptions);
      if (!result.ok()) {
        ACTS_ERROR("Propagation to true vertex position failed.");
        return std::nullopt;
      }
      auto& paramsAtVtx = *result->endParameters;
      return std::make_optional(paramsAtVtx);
    };

    for (const auto& trk : tracksAtVtx) {
      if (trk.trackWeight < m_cfg.minTrkWeight) {
        continue;
      }

      innerTrkWeight.push_back(trk.trackWeight);

      Acts::Vector3 trueUnitDir = Acts::Vector3::Zero();
      Acts::Vector3 trueMom = Acts::Vector3::Zero();

      const auto& particleOpt = findTruthParticle(trk);
      if (particleOpt.has_value()) {
        const auto& particle = *particleOpt;

        innerParticleId.push_back(particle.particleId().value());

        trueUnitDir = particle.direction();
        trueMom.head<2>() = Acts::makePhiThetaFromDirection(trueUnitDir);
        trueMom[2] = particle.qOverP();

        innerTruthPhi.push_back(trueMom[0]);
        innerTruthTheta.push_back(trueMom[1]);
        innerTruthQOverP.push_back(trueMom[2]);
      } else {
        ACTS_VERBOSE("Track has no matching truth particle.");

        innerParticleId.push_back(-1);

        innerTruthPhi.push_back(nan);
        innerTruthTheta.push_back(nan);
        innerTruthQOverP.push_back(nan);
      }

      // Save track parameters before the vertex fit
      const auto paramsAtVtx = propagateToVtx(
          *(trk.originalParams.as<Acts::BoundTrackParameters>()));
      if (paramsAtVtx.has_value()) {
        Acts::Vector3 recoMom =
            paramsAtVtx->parameters().segment(Acts::eBoundPhi, 3);
        const Acts::ActsMatrix<3, 3>& momCov =
            paramsAtVtx->covariance()->template block<3, 3>(Acts::eBoundPhi,
                                                            Acts::eBoundPhi);
        innerRecoPhi.push_back(recoMom[0]);
        innerRecoTheta.push_back(recoMom[1]);
        innerRecoQOverP.push_back(recoMom[2]);

        if (particleOpt.has_value()) {
          Acts::Vector3 diffMom = recoMom - trueMom;
          // Accounting for the periodicity of phi. We overwrite the
          // previously computed value for better readability.
          diffMom[0] = Acts::detail::difference_periodic(recoMom(0), trueMom(0),
                                                         2 * M_PI);
          innerResPhi.push_back(diffMom[0]);
          innerResTheta.push_back(diffMom[1]);
          innerResQOverP.push_back(diffMom[2]);

          innerPullPhi.push_back(pull(diffMom[0], momCov(0, 0), "phi", false));
          innerPullTheta.push_back(
              pull(diffMom[1], momCov(1, 1), "theta", false));
          innerPullQOverP.push_back(
              pull(diffMom[2], momCov(2, 2), "q/p", false));

          const auto& recoUnitDir = paramsAtVtx->direction();
          double overlap = trueUnitDir.dot(recoUnitDir);
          innerMomOverlap.push_back(overlap);
        } else {
          innerResPhi.push_back(nan);
          innerResTheta.push_back(nan);
          innerResQOverP.push_back(nan);
          innerPullPhi.push_back(nan);
          innerPullTheta.push_back(nan);
          innerPullQOverP.push_back(nan);
          innerMomOverlap.push_back(nan);
        }
      } else {
        innerRecoPhi.push_back(nan);
        innerRecoTheta.push_back(nan);
        innerRecoQOverP.push_back(nan);
        innerResPhi.push_back(nan);
        innerResTheta.push_back(nan);
        innerResQOverP.push_back(nan);
        innerPullPhi.push_back(nan);
        innerPullTheta.push_back(nan);
        innerPullQOverP.push_back(nan);
        innerMomOverlap.push_back(nan);
      }

      // Save track parameters after the vertex fit
      const auto paramsAtVtxFitted = propagateToVtx(trk.fittedParams);
      if (paramsAtVtxFitted.has_value()) {
        Acts::Vector3 recoMomFitted =
            paramsAtVtxFitted->parameters().segment(Acts::eBoundPhi, 3);
        const Acts::ActsMatrix<3, 3>& momCovFitted =
            paramsAtVtxFitted->covariance()->block<3, 3>(Acts::eBoundPhi,
                                                         Acts::eBoundPhi);
        innerRecoPhiFitted.push_back(recoMomFitted[0]);
        innerRecoThetaFitted.push_back(recoMomFitted[1]);
        innerRecoQOverPFitted.push_back(recoMomFitted[2]);

        if (particleOpt.has_value()) {
          Acts::Vector3 diffMomFitted = recoMomFitted - trueMom;
          // Accounting for the periodicity of phi. We overwrite the
          // previously computed value for better readability.
          diffMomFitted[0] = Acts::detail::difference_periodic(
              recoMomFitted(0), trueMom(0), 2 * M_PI);
          innerResPhiFitted.push_back(diffMomFitted[0]);
          innerResThetaFitted.push_back(diffMomFitted[1]);
          innerResQOverPFitted.push_back(diffMomFitted[2]);

          innerPullPhiFitted.push_back(
              pull(diffMomFitted[0], momCovFitted(0, 0), "phi"));
          innerPullThetaFitted.push_back(
              pull(diffMomFitted[1], momCovFitted(1, 1), "theta"));
          innerPullQOverPFitted.push_back(
              pull(diffMomFitted[2], momCovFitted(2, 2), "q/p"));

          const auto& recoUnitDirFitted = paramsAtVtxFitted->direction();
          double overlapFitted = trueUnitDir.dot(recoUnitDirFitted);
          innerMomOverlapFitted.push_back(overlapFitted);
        } else {
          innerResPhiFitted.push_back(nan);
          innerResThetaFitted.push_back(nan);
          innerResQOverPFitted.push_back(nan);
          innerPullPhiFitted.push_back(nan);
          innerPullThetaFitted.push_back(nan);
          innerPullQOverPFitted.push_back(nan);
          innerMomOverlapFitted.push_back(nan);
        }
      } else {
        innerRecoPhiFitted.push_back(nan);
        innerRecoThetaFitted.push_back(nan);
        innerRecoQOverPFitted.push_back(nan);
        innerResPhiFitted.push_back(nan);
        innerResThetaFitted.push_back(nan);
        innerResQOverPFitted.push_back(nan);
        innerPullPhiFitted.push_back(nan);
        innerPullThetaFitted.push_back(nan);
        innerPullQOverPFitted.push_back(nan);
        innerMomOverlapFitted.push_back(nan);
      }
    }
  }

  // fill the variables
  m_outputTree->Fill();

  m_recoX.clear();
  m_recoY.clear();
  m_recoZ.clear();
  m_recoT.clear();
  m_covXX.clear();
  m_covYY.clear();
  m_covZZ.clear();
  m_covTT.clear();
  m_covXY.clear();
  m_covXZ.clear();
  m_covXT.clear();
  m_covYZ.clear();
  m_covYT.clear();
  m_covZT.clear();
  m_nTracksOnRecoVertex.clear();
  m_recoVertexTotalTrackWeight.clear();
  m_sumPt2.clear();
  m_vertexPrimary.clear();
  m_vertexSecondary.clear();
  m_truthVertexMatchRatio.clear();
  m_nTracksOnTruthVertex.clear();
  m_recoVertexClassification.clear();
  m_truthX.clear();
  m_truthY.clear();
  m_truthZ.clear();
  m_truthT.clear();
  m_resX.clear();
  m_resY.clear();
  m_resZ.clear();
  m_resT.clear();
  m_pullX.clear();
  m_pullY.clear();
  m_pullZ.clear();
  m_pullT.clear();

  m_trkWeight.clear();
  m_recoPhi.clear();
  m_recoTheta.clear();
  m_recoQOverP.clear();
  m_recoPhiFitted.clear();
  m_recoThetaFitted.clear();
  m_recoQOverPFitted.clear();
  m_particleId.clear();
  m_truthPhi.clear();
  m_truthTheta.clear();
  m_truthQOverP.clear();
  m_resPhi.clear();
  m_resTheta.clear();
  m_resQOverP.clear();
  m_resPhiFitted.clear();
  m_resThetaFitted.clear();
  m_resQOverPFitted.clear();
  m_momOverlap.clear();
  m_momOverlapFitted.clear();
  m_pullPhi.clear();
  m_pullTheta.clear();
  m_pullQOverP.clear();
  m_pullPhiFitted.clear();
  m_pullThetaFitted.clear();
  m_pullQOverPFitted.clear();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

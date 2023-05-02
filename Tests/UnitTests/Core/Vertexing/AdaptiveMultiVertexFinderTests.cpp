// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Vertexing/AdaptiveGridDensityVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/GridDensityVertexFinder.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"

#include <chrono>

#include "VertexingDataHelper.hpp"

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;

using Covariance = BoundSymMatrix;
using Propagator = Acts::Propagator<EigenStepper<>>;
using Linearizer = HelicalTrackLinearizer<Propagator>;

// Create a test context
GeometryContext geoContext = GeometryContext();
MagneticFieldContext magFieldContext = MagneticFieldContext();

const std::string toolString = "AMVF";

/// @brief AMVF test with Gaussian seed finder
BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_finder_test) {
  // Set debug mode
  bool debugMode = false;
  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));

  // Set up EigenStepper
  // EigenStepper<> stepper(bField);
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // IP 3D Estimator
  using IPEstimator = ImpactPointEstimator<BoundTrackParameters, Propagator>;

  IPEstimator::Config ipEstimatorCfg(bField, propagator);
  IPEstimator ipEstimator(ipEstimatorCfg);

  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer>;

  Fitter::Config fitterCfg(ipEstimator);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  Fitter fitter(fitterCfg);

  using SeedFinder =
      TrackDensityVertexFinder<Fitter,
                               GaussianTrackDensity<BoundTrackParameters>>;

  SeedFinder seedFinder;

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEstimator,
                              std::move(linearizer), bField);

  // TODO: test this as well!
  // finderConfig.useBeamSpotConstraint = false;

  Finder finder(finderConfig);
  Finder::State state;

  auto csvData = readTracksAndVertexCSV(toolString);
  auto tracks = std::get<TracksData>(csvData);

  if (debugMode) {
    std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
    int maxCout = 10;
    int count = 0;
    for (const auto& trk : tracks) {
      std::cout << count << ". track: " << std::endl;
      std::cout << "params: " << trk << std::endl;
      count++;
      if (count == maxCout) {
        break;
      }
    }
  }

  std::vector<const BoundTrackParameters*> tracksPtr;
  for (const auto& trk : tracks) {
    tracksPtr.push_back(&trk);
  }

  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);

  vertexingOptions.vertexConstraint = std::get<BeamSpotData>(csvData);

  auto t1 = std::chrono::system_clock::now();
  auto findResult = finder.find(tracksPtr, vertexingOptions, state);
  auto t2 = std::chrono::system_clock::now();

  auto timediff =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex<BoundTrackParameters>> allVertices = *findResult;

  if (debugMode) {
    std::cout << "Time needed: " << timediff << " ms." << std::endl;
    std::cout << "Number of vertices reconstructed: " << allVertices.size()
              << std::endl;

    int count = 0;
    for (const auto& vtx : allVertices) {
      count++;
      std::cout << count << ". Vertex at position: " << vtx.position()[0]
                << ", " << vtx.position()[1] << ", " << vtx.position()[2]
                << std::endl;
      std::cout << count << ". Vertex with cov: " << vtx.covariance()
                << std::endl;
      std::cout << "\t with n tracks: " << vtx.tracks().size() << std::endl;
    }
  }

  // Test expected outcomes from athena implementation
  // Number of reconstructed vertices
  auto verticesInfo = std::get<VerticesData>(csvData);
  const int expNRecoVertices = verticesInfo.size();

  BOOST_CHECK_EQUAL(allVertices.size(), expNRecoVertices);

  for (int i = 0; i < expNRecoVertices; i++) {
    auto recoVtx = allVertices[i];
    auto expVtx = verticesInfo[i];
    CHECK_CLOSE_ABS(recoVtx.position(), expVtx.position, 0.001_mm);
    CHECK_CLOSE_ABS(recoVtx.covariance(), expVtx.covariance, 0.001_mm);
    BOOST_CHECK_EQUAL(recoVtx.tracks().size(), expVtx.nTracks);
    CHECK_CLOSE_ABS(recoVtx.tracks()[0].trackWeight, expVtx.trk1Weight, 0.003);
    CHECK_CLOSE_ABS(recoVtx.tracks()[0].vertexCompatibility, expVtx.trk1Comp,
                    0.003);
  }
}

// Dummy user-defined InputTrack type
struct InputTrack {
  InputTrack(const BoundTrackParameters& params, int id)
      : m_parameters(params), m_id(id) {}

  const BoundTrackParameters& parameters() const { return m_parameters; }
  // store e.g. link to original objects here

  int id() const { return m_id; }

 private:
  BoundTrackParameters m_parameters;

  // Some test track ID
  int m_id;
};

/// @brief AMVF test with user-defined input track type
BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_finder_usertype_test) {
  // Set debug mode
  bool debugMode = false;
  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));

  // Set up EigenStepper
  // EigenStepper<> stepper(bField);
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // Create a custom std::function to extract BoundTrackParameters from
  // user-defined InputTrack
  std::function<BoundTrackParameters(InputTrack)> extractParameters =
      [](const InputTrack& params) { return params.parameters(); };

  // IP 3D Estimator
  using IPEstimator = ImpactPointEstimator<InputTrack, Propagator>;

  IPEstimator::Config ipEstimatorCfg(bField, propagator);
  IPEstimator ipEstimator(ipEstimatorCfg);

  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<InputTrack, Linearizer>;

  Fitter::Config fitterCfg(ipEstimator);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  Fitter fitter(fitterCfg, extractParameters);

  using SeedFinder =
      TrackDensityVertexFinder<Fitter, GaussianTrackDensity<InputTrack>>;

  SeedFinder seedFinder(extractParameters);

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEstimator,
                              std::move(linearizer), bField);
  Finder::State state;

  Finder finder(finderConfig, extractParameters);

  auto csvData = readTracksAndVertexCSV(toolString);
  auto tracks = std::get<TracksData>(csvData);

  std::vector<InputTrack> userTracks;
  int idCount = 0;
  for (const auto& trk : tracks) {
    userTracks.push_back(InputTrack(trk, idCount));
    idCount++;
  }

  if (debugMode) {
    std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
    int maxCout = 10;
    int count = 0;
    for (const auto& trk : tracks) {
      std::cout << count << ". track: " << std::endl;
      std::cout << "params: " << trk << std::endl;
      count++;
      if (count == maxCout) {
        break;
      }
    }
  }

  std::vector<const InputTrack*> userTracksPtr;
  for (const auto& trk : userTracks) {
    userTracksPtr.push_back(&trk);
  }

  VertexingOptions<InputTrack> vertexingOptions(geoContext, magFieldContext);

  Vertex<InputTrack> constraintVtx;
  constraintVtx.setPosition(std::get<BeamSpotData>(csvData).position());
  constraintVtx.setCovariance(std::get<BeamSpotData>(csvData).covariance());

  vertexingOptions.vertexConstraint = constraintVtx;

  auto findResult = finder.find(userTracksPtr, vertexingOptions, state);

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex<InputTrack>> allVertices = *findResult;

  if (debugMode) {
    std::cout << "Number of vertices reconstructed: " << allVertices.size()
              << std::endl;

    int count = 0;
    for (const auto& vtx : allVertices) {
      count++;
      std::cout << count << ". Vertex at position: " << vtx.position()[0]
                << ", " << vtx.position()[1] << ", " << vtx.position()[2]
                << std::endl;
      std::cout << count << ". Vertex with cov: " << vtx.covariance()
                << std::endl;
      std::cout << "\t with n tracks: " << vtx.tracks().size() << std::endl;
    }
    for (auto& trk : allVertices[0].tracks()) {
      std::cout << "Track ID at first vertex: " << trk.originalParams->id()
                << std::endl;
    }
  }

  auto verticesInfo = std::get<VerticesData>(csvData);
  const int expNRecoVertices = verticesInfo.size();

  BOOST_CHECK_EQUAL(allVertices.size(), expNRecoVertices);

  for (int i = 0; i < expNRecoVertices; i++) {
    auto recoVtx = allVertices[i];
    auto expVtx = verticesInfo[i];
    CHECK_CLOSE_ABS(recoVtx.position(), expVtx.position, 0.001_mm);
    CHECK_CLOSE_ABS(recoVtx.covariance(), expVtx.covariance, 0.001_mm);
    BOOST_CHECK_EQUAL(recoVtx.tracks().size(), expVtx.nTracks);
    CHECK_CLOSE_ABS(recoVtx.tracks()[0].trackWeight, expVtx.trk1Weight, 0.003);
    CHECK_CLOSE_ABS(recoVtx.tracks()[0].vertexCompatibility, expVtx.trk1Comp,
                    0.003);
  }
}

/// @brief AMVF test with grid seed finder
BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_finder_grid_seed_finder_test) {
  // Set debug mode
  bool debugMode = false;
  if (debugMode) {
    std::cout << "Starting AMVF test with grid seed finder..." << std::endl;
  }
  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));

  // Set up EigenStepper
  // EigenStepper<> stepper(bField);
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // IP Estimator
  using IPEstimator = ImpactPointEstimator<BoundTrackParameters, Propagator>;

  IPEstimator::Config ipEstCfg(bField, propagator);
  IPEstimator ipEst(ipEstCfg);

  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer>;

  Fitter::Config fitterCfg(ipEst);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  Fitter fitter(fitterCfg);

  using SeedFinder = GridDensityVertexFinder<4000, 55>;
  SeedFinder::Config seedFinderCfg(250);
  seedFinderCfg.cacheGridStateForTrackRemoval = true;

  SeedFinder seedFinder(seedFinderCfg);

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEst,
                              std::move(linearizer), bField);

  // TODO: test this as well!
  // finderConfig.useBeamSpotConstraint = false;

  Finder finder(finderConfig);
  Finder::State state;

  auto csvData = readTracksAndVertexCSV(toolString);
  auto tracks = std::get<TracksData>(csvData);

  if (debugMode) {
    std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
    int maxCout = 10;
    int count = 0;
    for (const auto& trk : tracks) {
      std::cout << count << ". track: " << std::endl;
      std::cout << "params: " << trk << std::endl;
      count++;
      if (count == maxCout) {
        break;
      }
    }
  }

  std::vector<const BoundTrackParameters*> tracksPtr;
  for (const auto& trk : tracks) {
    tracksPtr.push_back(&trk);
  }

  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);

  vertexingOptions.vertexConstraint = std::get<BeamSpotData>(csvData);

  auto t1 = std::chrono::system_clock::now();
  auto findResult = finder.find(tracksPtr, vertexingOptions, state);
  auto t2 = std::chrono::system_clock::now();

  auto timediff =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex<BoundTrackParameters>> allVertices = *findResult;

  if (debugMode) {
    std::cout << "Time needed: " << timediff << " ms." << std::endl;
    std::cout << "Number of vertices reconstructed: " << allVertices.size()
              << std::endl;

    int count = 0;
    for (const auto& vtx : allVertices) {
      count++;
      std::cout << count << ". Vertex at position: " << vtx.position()[0]
                << ", " << vtx.position()[1] << ", " << vtx.position()[2]
                << std::endl;
      std::cout << count << ". Vertex with cov: " << vtx.covariance()
                << std::endl;
      std::cout << "\t with n tracks: " << vtx.tracks().size() << std::endl;
    }
  }
  // Test expected outcomes from athena implementation
  // Number of reconstructed vertices
  auto verticesInfo = std::get<VerticesData>(csvData);
  const int expNRecoVertices = verticesInfo.size();

  BOOST_CHECK_EQUAL(allVertices.size(), expNRecoVertices);
  std::vector<bool> vtxFound(expNRecoVertices, false);

  for (const auto& vtx : allVertices) {
    double vtxZ = vtx.position()[2];
    double diffZ = 1e5;
    int foundVtxIdx = -1;
    for (int i = 0; i < expNRecoVertices; i++) {
      if (not vtxFound[i]) {
        if (std::abs(vtxZ - verticesInfo[i].position[2]) < diffZ) {
          diffZ = std::abs(vtxZ - verticesInfo[i].position[2]);
          foundVtxIdx = i;
        }
      }
    }
    if (diffZ < 0.5_mm) {
      vtxFound[foundVtxIdx] = true;
      CHECK_CLOSE_ABS(vtx.tracks().size(), verticesInfo[foundVtxIdx].nTracks,
                      1);
    }
  }
  for (bool found : vtxFound) {
    BOOST_CHECK_EQUAL(found, true);
  }
}

/// @brief AMVF test with adaptive grid seed finder
BOOST_AUTO_TEST_CASE(
    adaptive_multi_vertex_finder_adaptive_grid_seed_finder_test) {
  // Set debug mode
  bool debugMode = false;
  if (debugMode) {
    std::cout << "Starting AMVF test with adaptive grid seed finder..."
              << std::endl;
  }
  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));

  // Set up EigenStepper
  // EigenStepper<> stepper(bField);
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // IP Estimator
  using IPEstimator = ImpactPointEstimator<BoundTrackParameters, Propagator>;

  IPEstimator::Config ipEstCfg(bField, propagator);
  IPEstimator ipEst(ipEstCfg);

  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer>;

  Fitter::Config fitterCfg(ipEst);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  Fitter fitter(fitterCfg);

  using SeedFinder = AdaptiveGridDensityVertexFinder<55>;
  SeedFinder::Config seedFinderCfg(0.05);
  seedFinderCfg.cacheGridStateForTrackRemoval = true;

  SeedFinder seedFinder(seedFinderCfg);

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEst,
                              std::move(linearizer), bField);

  Finder finder(finderConfig);
  Finder::State state;

  auto csvData = readTracksAndVertexCSV(toolString);
  auto tracks = std::get<TracksData>(csvData);

  if (debugMode) {
    std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
    int maxCout = 10;
    int count = 0;
    for (const auto& trk : tracks) {
      std::cout << count << ". track: " << std::endl;
      std::cout << "params: " << trk << std::endl;
      count++;
      if (count == maxCout) {
        break;
      }
    }
  }

  std::vector<const BoundTrackParameters*> tracksPtr;
  for (const auto& trk : tracks) {
    tracksPtr.push_back(&trk);
  }

  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);

  vertexingOptions.vertexConstraint = std::get<BeamSpotData>(csvData);

  auto t1 = std::chrono::system_clock::now();
  auto findResult = finder.find(tracksPtr, vertexingOptions, state);
  auto t2 = std::chrono::system_clock::now();

  auto timediff =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex<BoundTrackParameters>> allVertices = *findResult;

  if (debugMode) {
    std::cout << "Time needed: " << timediff << " ms." << std::endl;
    std::cout << "Number of vertices reconstructed: " << allVertices.size()
              << std::endl;

    int count = 0;
    for (const auto& vtx : allVertices) {
      count++;
      std::cout << count << ". Vertex at position: " << vtx.position()[0]
                << ", " << vtx.position()[1] << ", " << vtx.position()[2]
                << std::endl;
      std::cout << count << ". Vertex with cov: " << vtx.covariance()
                << std::endl;
      std::cout << "\t with n tracks: " << vtx.tracks().size() << std::endl;
    }
  }
  // Test expected outcomes from athena implementation
  // Number of reconstructed vertices
  auto verticesInfo = std::get<VerticesData>(csvData);
  const int expNRecoVertices = verticesInfo.size();

  BOOST_CHECK_EQUAL(allVertices.size(), expNRecoVertices);
  std::vector<bool> vtxFound(expNRecoVertices, false);

  for (const auto& vtx : allVertices) {
    double vtxZ = vtx.position()[2];
    double diffZ = 1e5;
    int foundVtxIdx = -1;
    for (int i = 0; i < expNRecoVertices; i++) {
      if (not vtxFound[i]) {
        if (std::abs(vtxZ - verticesInfo[i].position[2]) < diffZ) {
          diffZ = std::abs(vtxZ - verticesInfo[i].position[2]);
          foundVtxIdx = i;
        }
      }
    }
    if (diffZ < 0.5_mm) {
      vtxFound[foundVtxIdx] = true;
      CHECK_CLOSE_ABS(vtx.tracks().size(), verticesInfo[foundVtxIdx].nTracks,
                      2);
    }
  }
  for (bool found : vtxFound) {
    BOOST_CHECK_EQUAL(found, true);
  }
}

/// @brief AMVF test with Gaussian seed finder
BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_finder_test_with_timing) {
  using IPEstimator = ImpactPointEstimator<BoundTrackParameters, Propagator>;
  using Fitter = AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer>;
  using SeedFinder = TrackDensityVertexFinder<Fitter,
                              GaussianTrackDensity<BoundTrackParameters>>;
  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  bool printLoadedData = false;
  bool printSeed = false;
  bool printAssignedTracks = false;
  bool printMultimap = false;
  bool printImpactParams = false;
  bool printCompatibilities = false;
  bool printWeights = false;
  bool printFittedVtx = false;
  bool printIfGood = false;

  //load data
  std::string toolStringTiming = "truth";
  auto csvData = readTracksAndVertexCSV(toolStringTiming, "vertexing_event_mu3");

  //retrieve track and vertex data
  auto tracks = std::get<TracksData>(csvData);
  auto trueVtxInfo = std::get<VerticesData>(csvData);

  //make vector of pointers to tracks
  std::vector<const Fitter::InputTrack_t*> tracksPtr;
  for (const auto& trk : tracks) {
    tracksPtr.push_back(&trk);
  }
  
  //print track parameters and true vertex positions:
  if (printLoadedData) {
    std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
    int maxCount = 12;
    int trkCount = 0;
    for (const auto& trk : tracks) {
      std::cout << trkCount+1 << ". track: " << std::endl;
      std::cout << "params: \n" << trk << std::endl;
      //std::cout << "covariance: \n" << trk.covariance().value() << std::endl;
      trkCount++;
      if (trkCount == maxCount) {
        break;
      }
    }
    int vertCount = 0;
    for (const auto& vtx : trueVtxInfo) {
      vertCount++;
      std::cout << vertCount+1 << ". Vertex at position:\n" << vtx.position
                << std::endl;
    }
  }

  //define fitter and finder (two instances of classes, which we use for vertexing)
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));
  
  //what is this?
  EigenStepper<> stepper(bField); 

  // Set up propagator with void navigator - what is this?
  auto propagator = std::make_shared<Propagator>(stepper);

  //estimates the impact point, i.e., the point of closest approach of the track to a vertex (?)
  IPEstimator::Config ipEstimatorCfg(bField, propagator);
  IPEstimator ipEstimator(ipEstimatorCfg);

  //weight function depends on a "temperature", which gets gradually decreased 
  //this procedure called "annealing" reduces effects of outliers
  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);
  //vertexingOptions.vertexConstraint = std::get<BeamSpotData>(csvData); //if uncomment need to set useBeamSpotConstraint = true

  Fitter::Config fitterCfg(ipEstimator);
  fitterCfg.annealingTool = annealingUtility;
  fitterCfg.doSmoothing = true; //should this be true?
  Fitter fitter(fitterCfg);
  Fitter::State fitterState(*bField, vertexingOptions.magFieldContext);

  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  SeedFinder seedFinder;
  SeedFinder::State seedFinderState;
  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEstimator,
                              std::move(linearizer), bField);
  finderConfig.useBeamSpotConstraint = false; // TODO: test this as well!
  Finder finder(finderConfig);
  Finder::State state;

  //do a first estimate for the vertex (so-called "seed")
  std::vector<const Fitter::InputTrack_t*> seedTracks = tracksPtr; //tracks from which we compute the seed
  Vertex<Fitter::InputTrack_t> currentConstraint = vertexingOptions.vertexConstraint; 
  std::vector<const Fitter::InputTrack_t*> removedSeedTracks; //TODO: all compatible tracks with last vertex that need to be removed from seed tracks
  auto seedResult = finder.doSeeding(seedTracks, currentConstraint, vertexingOptions,
                                     seedFinderState, removedSeedTracks);
  
  auto seed = *std::make_unique<Vertex<Fitter::InputTrack_t>>(*seedResult);
  if (printSeed){
    std::cout << "\n position of seed:\n" << seed.fullPosition() << "\n";
  }

  //tracks that can be assigned to vertex
  std::vector<const Fitter::InputTrack_t*> searchTracks = tracksPtr;
  //assign tracks to vertex
  auto prepResult = finder.canPrepareVertexForFit(searchTracks, seedTracks,
                                           seed, currentConstraint,
                                           fitterState, vertexingOptions);
  if(printAssignedTracks){
    bool printParams = true;
    if (printParams){
      std::cout << "\nparameters of assigned tracks:\n\n";
    }
    else {
      std::cout << "\npositions and time of assigned tracks:\n\n";
    }
    for (const auto& trk : fitterState.vtxInfoMap[&seed].trackLinks){
      if (printParams){
        std::cout << *trk << "\n\n";
      }
      else{
        auto params = finder.m_extractParameters(*trk);
        auto pos = params.position(vertexingOptions.geoContext);
        std::cout << pos << " time: "<< params.time() << "\n\n";
      }
    }
  }

  //map every track to all of its connected vertices
  fitterState.addVertexToMultiMap(seed); 
  if (printMultimap){
    std::cout << "\nmultimap of tracks to verts:\n\n";
    for (auto itr = fitterState.trackToVerticesMultiMap.begin(); 
        itr != fitterState.trackToVerticesMultiMap.end(); ++itr){
            auto params = finder.m_extractParameters(*(itr->first));
            auto pos = params.position(vertexingOptions.geoContext);
            std::cout << "track position: \n" << pos << "\nvertex position:\n" << (itr->second->fullPosition()) << "\n\n"; 
        }
  }

  //estimate the track parameters in the point of closest approach to the vertex (i.e., the impact point)
  // The seed position (head<3> returns only the first 3 elements of the vector)
  //the Perigee parametrization is done wrt the vertex seed position rather than wrt the origin
  const Vector3& seedPos = fitterState.vtxInfoMap[&seed].seedPosition.template head<3>();
  if (printImpactParams) std::cout << "\nimpact parameters:\n\n";
  for (const auto& trk : fitterState.vtxInfoMap[&seed].trackLinks) {
    auto res = finder.config().vertexFitter.m_cfg.ipEst.estimate3DImpactParameters(
        vertexingOptions.geoContext, vertexingOptions.magFieldContext,
        finder.m_extractParameters(*trk), seedPos, fitterState.ipState);
    // Set ip3dParams
    fitterState.vtxInfoMap[&seed].ip3dParams.emplace(trk, res.value());
    if (printImpactParams) std::cout << "\n" << res.value() << "\n";
  }

  //we only consider one vertex here
  std::vector<Vertex<Fitter::InputTrack_t>*> verticesToFit;
  verticesToFit.push_back(&seed);
  fitterState.vertexCollection = verticesToFit;

  //compute compatibility of all tracks with vertex (low compatibility means that tracks are compatible)
  //timing info does not enter!
  if (printCompatibilities){
    std::cout << "\n\nVertex Compatibilities:\n";
  }

  if (fitterState.vtxInfoMap[&seed].constraintVertex.fullCovariance() !=
          SymMatrix4::Zero()) {
        //std::cout << "\nSetting position to:\n" << state.vtxInfoMap[&seed].constraintVertex.fullPosition() << "\n";
        seed.setFullPosition(
            fitterState.vtxInfoMap[&seed].constraintVertex.fullPosition()); //only set z position!
        seed.setFitQuality(
            fitterState.vtxInfoMap[&seed].constraintVertex.fitQuality());
        seed.setFullCovariance(
            fitterState.vtxInfoMap[&seed].constraintVertex.fullCovariance());
        } 
  std::cout << "\n seed vertex covariance:\n" << seed.fullCovariance();
  double weight =
          1. / finder.config().vertexFitter.m_cfg.annealingTool.getWeight(fitterState.annealingState, 1.); //why do we set chi2=1?
  seed.setFullCovariance(seed.fullCovariance() * weight); //why do we do this? has no effect atm
  
  for (const auto& trk : fitterState.vtxInfoMap[&seed].trackLinks) {
    auto& trkAtVtx =
        fitterState.tracksAtVerticesMap.at(std::make_pair(trk, &seed));
    // Set compatibility with current vertex
    auto compRes = finder.config().vertexFitter.m_cfg.ipEst.get3dVertexCompatibility(
        vertexingOptions.geoContext, &(fitterState.vtxInfoMap[&seed].ip3dParams.at(trk)),
        VectorHelpers::position(seed.fullPosition()));
    trkAtVtx.vertexCompatibility = *compRes;
    if (printCompatibilities){
      std::cout << *compRes << "\n";
    }
  }
  finder.config().vertexFitter.setWeightsAndUpdate(fitterState, finder.config().linearizer, vertexingOptions);

  //compute track weights
  /*
  if (printWeights) {
    std::cout << "\n\nTrack Weights (minWeight = " << finder.config().vertexFitter.m_cfg.minWeight << ")\n";
  }
  else {
    std::cout << "\n";
  }
  if (printFittedVtx){
    std::cout << "Vertex positions when adding tracks to the fit:\n";
  }
  for (const auto& trk : fitterState.vtxInfoMap[&seed].trackLinks) {
    auto& trkAtVtx = fitterState.tracksAtVerticesMap.at(std::make_pair(trk, &seed));
    double currentTrkWeight = finder.config().vertexFitter.m_cfg.annealingTool.getWeight(
          fitterState.annealingState, trkAtVtx.vertexCompatibility,
          finder.config().vertexFitter.collectTrackToVertexCompatibilities(fitterState, trk));
    trkAtVtx.trackWeight = currentTrkWeight;
    if (printWeights) {
      std::cout << "weight: " << currentTrkWeight << "\n";
    }
    if (trkAtVtx.trackWeight > finder.config().vertexFitter.m_cfg.minWeight) {
        // Check if linearization state exists or needs to be relinearized
        if (not trkAtVtx.isLinearized || fitterState.vtxInfoMap[&seed].relinearize) {
          auto result = finder.config().linearizer.linearizeTrack(
              finder.m_extractParameters(*trk), fitterState.vtxInfoMap[&seed].oldPosition,
              vertexingOptions.geoContext, vertexingOptions.magFieldContext,
              fitterState.linearizerState);
          if (trkAtVtx.isLinearized) {
            fitterState.vtxInfoMap[&seed].linPoint = fitterState.vtxInfoMap[&seed].oldPosition;
          }
          trkAtVtx.linearizedState = *result;
          trkAtVtx.isLinearized = true;
          std::cout << "\n position jacobian:\n" << result->positionJacobian << "\n\n";
        }
        // Update the vertex with the new track
        KalmanVertexUpdater::updateVertexWithTrack<Fitter::InputTrack_t>(seed,
                                                                        trkAtVtx);
        if (printFittedVtx){
          std::cout << "position:\n" << seed.fullPosition() << "\n\n";
        }
    }
  }

  //check if vertex is good
  std::vector<Vertex<Fitter::InputTrack_t>*> allVerticesPtr {&seed};
  auto [nCompatibleTracks, isGoodVertex] =
        finder.checkVertexAndCompatibleTracks(seed, seedTracks, fitterState);
  
  bool knV = finder.keepNewVertex(seed, allVerticesPtr, fitterState);

  if (printIfGood) {
    std::cout << "\n\n Good Vertex: " << isGoodVertex << " Keep Vertex: " << knV <<"\n";
  }
  */
/*

  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));

  // Set up EigenStepper
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // IP 3D Estimator
  using IPEstimator = ImpactPointEstimator<BoundTrackParameters, Propagator>;

  IPEstimator::Config ipEstimatorCfg(bField, propagator);
  IPEstimator ipEstimator(ipEstimatorCfg);

  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer>;

  Fitter::Config fitterCfg(ipEstimator);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  Fitter fitter(fitterCfg);

  using SeedFinder =
      TrackDensityVertexFinder<Fitter,
                               GaussianTrackDensity<BoundTrackParameters>>;

  SeedFinder seedFinder;

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEstimator,
                              std::move(linearizer), bField);

  // TODO: test this as well!
  // finderConfig.useBeamSpotConstraint = false;

  Finder finder(finderConfig);
  Finder::State state;



  if (debugMode) {
    std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
    int maxCout = 12;
    int count = 0;
    for (const auto& trk : tracks) {
      std::cout << count << ". track: " << std::endl;
      std::cout << "params: \n" << trk << std::endl;
      std::cout << "covariance: \n" << trk.covariance().value() << std::endl;
      count++;
      if (count == maxCout) {
        break;
      }
    }
  }
  


  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);
  
  //work without vertexConstraint?
  vertexingOptions.vertexConstraint = std::get<BeamSpotData>(csvData); 

  const std::vector<const Fitter::InputTrack_t*>& origTracks = tracksPtr;

  Fitter::State fitterState(bField, magFieldContext);
  SeedFinder::State seedFinderState;

  std::vector<std::unique_ptr<Vertex<Fitter::InputTrack_t>>> allVertices;
  std::vector<Vertex<Fitter::InputTrack_t>*> allVerticesPtr;
  std::vector<const Fitter::InputTrack_t*> removedSeedTracks;
  */
  /*
  auto t1 = std::chrono::system_clock::now();
  auto findResult = finder.find(tracksPtr, vertexingOptions, state);
  auto t2 = std::chrono::system_clock::now();

  auto timediff =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex<BoundTrackParameters>> allVertices = *findResult;

  if (debugMode) {
    std::cout << "Time needed: " << timediff << " ms." << std::endl;
    std::cout << "Number of vertices reconstructed: " << allVertices.size()
              << std::endl;

    int count = 0;
    for (const auto& vtx : allVertices) {
      count++;
      std::cout << count << ". Vertex at position: " << vtx.position()[0]
                << ", " << vtx.position()[1] << ", " << vtx.position()[2]
                << std::endl;
      std::cout << count << ". Vertex with cov: " << vtx.covariance()
                << std::endl;
      std::cout << "\t with n tracks: " << vtx.tracks().size() << std::endl;
    }
  }

  // Test expected outcomes from athena implementation
  // Number of reconstructed vertices
  auto verticesInfo = std::get<VerticesData>(csvData);
  const int expNRecoVertices = verticesInfo.size();

  std::cout << "\n" << allVertices.size() << " " << expNRecoVertices << "\n";
  */
  /*
  BOOST_CHECK_EQUAL(allVertices.size(), expNRecoVertices);
  for (int i = 0; i < expNRecoVertices; i++) {
    auto recoVtx = allVertices[i];
    auto expVtx = verticesInfo[i];
    CHECK_CLOSE_ABS(recoVtx.position(), expVtx.position, 0.001_mm);
    CHECK_CLOSE_ABS(recoVtx.covariance(), expVtx.covariance, 0.001_mm);
    BOOST_CHECK_EQUAL(recoVtx.tracks().size(), expVtx.nTracks);
    CHECK_CLOSE_ABS(recoVtx.tracks()[0].trackWeight, expVtx.trk1Weight, 0.003);
    CHECK_CLOSE_ABS(recoVtx.tracks()[0].vertexCompatibility, expVtx.trk1Comp,
                    0.003);
  }
  */
}

}  // namespace Test
}  // namespace Acts

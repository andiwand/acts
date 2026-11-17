// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Throughput benchmark for the GBTS seeder on an ITk-like synthetic event.
///
/// The detector and the event come from `ActsFatras::Synthetic` through
/// SyntheticEventOptions.hpp, shared with the other seeding benchmarks. Only
/// the GBTS-specific geometry, namely the logical layer numbering and the layer
/// connection table, is built here.
///
/// Two timings are reported. The one-shot `createSeeds` is what an experiment
/// calls, and it includes filling the node storage. The second times that fill
/// alone, so that the graph stage is the difference between the two. A filled
/// storage cannot be reused instead: the graph builder writes its edge state
/// into the storage and only partly resets it.
///
/// The cuts are the ACTS defaults, which are what ATLAS runs on the ITk pixel
/// detector. The connection table is `<layout>-gbts-connections.txt` from the
/// Fatras data when it exists, read with the Athena layer ids (see
/// `athenaLayerIds`). The ITk one is learned from the generator truth with
/// `--train-connection-table` at its defaults, as Athena learned its own; an
/// Athena table can be passed with `--connection-table`. Without a file, or
/// with `--connection-table none`, a hand-written table is used, which misses
/// endcap links and loses a few points of efficiency forward.
///
/// @note A connection table belongs to the layout it was written for. A layout
///       describing the rings of an endcap has many more discs per side than
///       one outlining it, and a table whose reach along z does not cover them
///       costs half the efficiency. See `kDiscReach`.
///
/// @note GBTS returns seeds of four to eleven space points, so it is scored on
///       sharing `evaluateSeeds`' minimum with one primary rather than on every
///       space point matching. Scoring it the strict way costs four points of
///       efficiency that say nothing about the seeder.

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsGraphBuilder.hpp"
#include "Acts/Seeding/GbtsLayerConnection.hpp"
#include "Acts/Seeding/GbtsNodeStorage.hpp"
#include "Acts/Seeding/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/EventGenerator.hpp"
#include "ActsFatras/Synthetic/SeedingTruth.hpp"
#include "ActsTests/CommonHelpers/BenchmarkTools.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <unordered_set>
#include <stdexcept>
#include <iostream>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "SeedEtaReport.hpp"
#include "SyntheticEventOptions.hpp"

namespace po = boost::program_options;
namespace Exp = Acts::Experimental;

using namespace Acts::UnitLiterals;
using namespace ActsTests;
using namespace ActsFatras::Synthetic;

namespace {

/// Rings a disc may carry before the endcap encoding below runs out of room.
constexpr int kRingsPerDisc = 8;

/// How many discs back along z the connection table lets a disc reach. One per
/// ring set of the ITk pixel endcap, see `makeConnectionTable`.
constexpr int kDiscReach = 9;

/// GBTS numbers its logical layers as `gbtsId * 1000 + subIndex` and derives
/// the barrel flag from the leading digit as `layerId / 10000 == 8`, so the
/// barrel takes the eighties, one group per cylinder, with the eta module as
/// the sub-index.
///
/// The endcaps cannot do the same: a layout describing the ring structure of a
/// real endcap has far more than ten discs per side, so a group per disc would
/// run out of the nineties and collide with the barrel. Each endcap is instead
/// one group, 90 and 70, with the disc and the ring packed into the sub-index
/// together.
Exp::GbtsExperimentLayerId gbtsLayerId(const DetectorLayer& layer) {
  const Exp::GbtsExperimentLayerId layerIdx = layer.layer;
  const Exp::GbtsExperimentLayerId moduleIdx = layer.moduleIndex;
  if (layer.shape == SurfaceShape::Cylinder) {
    return (80 + layerIdx) * 1000 + moduleIdx;
  }
  const Exp::GbtsExperimentLayerId group =
      layer.side == SurfaceSide::Positive ? 90 : 70;
  return group * 1000 + layerIdx * kRingsPerDisc + moduleIdx;
}

/// The `gbtsLayerId` encoding above only holds within limits: the barrel group
/// numbers have to stay in the eighties for `layerId / 10000 == 8` to mean
/// barrel, and both sub-indices have to fit below the thousands digit.
/// `GbtsLayerIndex` is also 16 bits, which bounds the number of layers.
/// @param layout the layout to check
void checkEncodingFits(const DetectorLayout& layout) {
  int numBarrel = 0;
  int maxBarrelModule = 0;
  int maxDiscSubIndex = 0;
  for (const DetectorLayer& layer : layout.layers) {
    const auto layerIdx = static_cast<int>(layer.layer);
    const auto moduleIdx = static_cast<int>(layer.moduleIndex);
    if (layer.shape == SurfaceShape::Cylinder) {
      numBarrel = std::max(numBarrel, layerIdx + 1);
      maxBarrelModule = std::max(maxBarrelModule, moduleIdx);
    } else {
      if (moduleIdx >= kRingsPerDisc) {
        throw std::invalid_argument(
            "the GBTS layer id encoding takes at most eight rings per disc");
      }
      maxDiscSubIndex =
          std::max(maxDiscSubIndex, layerIdx * kRingsPerDisc + moduleIdx);
    }
  }
  if (numBarrel > 10) {
    throw std::invalid_argument(
        "the GBTS layer id encoding takes at most ten barrel layers");
  }
  if (maxBarrelModule >= 1000) {
    throw std::invalid_argument(
        "the GBTS layer id encoding takes at most a thousand modules per "
        "cylinder");
  }
  if (maxDiscSubIndex >= 1000) {
    throw std::invalid_argument(
        "the GBTS layer id encoding takes at most a hundred and twenty-five "
        "discs per side");
  }
  if (layout.layers.size() > 65535) {
    throw std::invalid_argument("GbtsLayerIndex is 16 bits wide");
  }
}

std::vector<Exp::GbtsLayerDescription> makeLayerDescriptions(
    const DetectorLayout& layout) {
  std::vector<Exp::GbtsLayerDescription> descriptions;
  descriptions.reserve(layout.layers.size());
  for (const DetectorLayer& layer : layout.layers) {
    Exp::GbtsLayerDescription description;
    description.id = gbtsLayerId(layer);
    description.type = layer.shape == SurfaceShape::Cylinder
                           ? Exp::GbtsLayerType::Barrel
                           : Exp::GbtsLayerType::Endcap;
    description.refCoord = layer.refCoord;
    description.minBound = layer.minBound;
    description.maxBound = layer.maxBound;
    descriptions.push_back(description);
  }
  return descriptions;
}

/// The connection table lists which layer may feed which. `src` is the outer
/// layer of a doublet, `dst` the inner one. The per-eta-bin compatibility is
/// worked out by GbtsGeometry from the beam spot constraint, so listing a
/// layer pair here only says that the pair is allowed at all.
std::vector<Exp::GbtsLayerConnection> makeConnections(
    const DetectorLayout& layout) {
  // group layers by (side, layer) so that all modules of a surface are linked
  auto modulesOf = [&](SurfaceSide side, int layerIdx) {
    std::vector<Exp::GbtsExperimentLayerId> ids;
    for (const DetectorLayer& layer : layout.layers) {
      if (layer.side == side &&
          layer.layer == static_cast<std::uint32_t>(layerIdx)) {
        ids.push_back(gbtsLayerId(layer));
      }
    }
    return ids;
  };

  std::vector<Exp::GbtsLayerConnection> connections;
  auto connect = [&](SurfaceSide outerSide, int outerLayer,
                     SurfaceSide innerSide, int innerLayer) {
    for (const Exp::GbtsExperimentLayerId src :
         modulesOf(outerSide, outerLayer)) {
      for (const Exp::GbtsExperimentLayerId dst :
           modulesOf(innerSide, innerLayer)) {
        connections.push_back({src, dst});
      }
    }
  };

  int numBarrel = 0;
  int numDiscs = 0;
  for (const DetectorLayer& layer : layout.layers) {
    if (layer.shape == SurfaceShape::Cylinder) {
      numBarrel = std::max(numBarrel, static_cast<int>(layer.layer) + 1);
    } else if (layer.side == SurfaceSide::Positive) {
      numDiscs = std::max(numDiscs, static_cast<int>(layer.layer) + 1);
    }
  }

  // barrel to barrel, adjacent and one skipped layer
  for (int i = 0; i < numBarrel; ++i) {
    for (int step : {1, 2}) {
      if (i + step < numBarrel) {
        connect(SurfaceSide::Barrel, i + step, SurfaceSide::Barrel, i);
      }
    }
  }
  for (const SurfaceSide side :
       {SurfaceSide::Positive, SurfaceSide::Negative}) {
    // Disc to disc. Discs are numbered outwards along z but consecutive ones do
    // not sit at consecutive radii, a layout describing the rings of an endcap
    // interleaving the ring sets, so `kDiscReach` has to cover the number of
    // ring sets rather than just the neighbours.
    for (int j = 0; j < numDiscs; ++j) {
      for (int step = 1; step <= kDiscReach; ++step) {
        if (j + step < numDiscs) {
          connect(side, j + step, side, j);
        }
      }
    }
    // forward transition: the innermost barrel layers are fed by the first
    // discs of each side, and by as many of them as a disc reaches back over
    for (int j = 0; j < kDiscReach; ++j) {
      for (int i = 0; i < numBarrel; ++i) {
        connect(side, j, SurfaceSide::Barrel, i);
      }
    }
  }

  return connections;
}


/// Athena's ITk pixel layer ids: barrel `(80 + layer) * 1000`, endcap
/// `(70|90 + ringLayer) * 1000 + ring`. A ring layer is one radius band of
/// rings, ordered by radius and then by the |z| it starts at (the inclined-like
/// low-z set before the endcap set); a ring is its rank in |z| within the band.
std::vector<Exp::GbtsExperimentLayerId> athenaLayerIds(
    const DetectorLayout& layout) {
  struct Band {
    long rMin;
    long rMax;
    std::vector<float> absZ;
  };
  std::vector<Band> bands;
  const auto key = [](float v) { return std::lround(v * 10.f); };
  for (const DetectorLayer& layer : layout.layers) {
    if (layer.shape == SurfaceShape::Cylinder ||
        layer.side != SurfaceSide::Positive) {
      continue;
    }
    auto it = std::ranges::find_if(bands, [&](const Band& b) {
      return b.rMin == key(layer.minBound) && b.rMax == key(layer.maxBound);
    });
    if (it == bands.end()) {
      bands.push_back({key(layer.minBound), key(layer.maxBound), {}});
      it = bands.end() - 1;
    }
    it->absZ.push_back(std::abs(layer.refCoord));
  }
  for (Band& b : bands) {
    std::ranges::sort(b.absZ);
  }
  std::ranges::sort(bands, [](const Band& a, const Band& b) {
    const long ga = std::lround(a.rMin / 300.);
    const long gb = std::lround(b.rMin / 300.);
    return ga != gb ? ga < gb : a.absZ.front() < b.absZ.front();
  });
  std::vector<Exp::GbtsExperimentLayerId> ids;
  for (const DetectorLayer& layer : layout.layers) {
    if (layer.shape == SurfaceShape::Cylinder) {
      ids.push_back((80 + layer.layer) * 1000);
      continue;
    }
    const auto band = std::ranges::find_if(bands, [&](const Band& b) {
      return b.rMin == key(layer.minBound) && b.rMax == key(layer.maxBound);
    });
    if (band == bands.end()) {
      throw std::invalid_argument("endcap ring without a positive-side twin");
    }
    const auto ring = std::ranges::find_if(band->absZ, [&](float z) {
      return std::abs(z - std::abs(layer.refCoord)) < 0.05f;
    });
    const auto ringLayer = static_cast<Exp::GbtsExperimentLayerId>(
        std::distance(bands.begin(), band));
    const auto base = layer.side == SurfaceSide::Positive ? 90 : 70;
    ids.push_back((base + ringLayer) * 1000 +
                  static_cast<Exp::GbtsExperimentLayerId>(
                      std::distance(band->absZ.begin(), ring)));
  }
  return ids;
}

/// Read the (src, dst) pairs and the eta bin width of an Athena GBTS
/// connection table, skipping the bin tables as Athena does.
std::vector<Exp::GbtsLayerConnection> readConnectionTable(
    const std::string& path, float& etaBinWidth) {
  std::ifstream in(path);
  std::uint32_t n = 0;
  if (!(in >> n >> etaBinWidth)) {
    throw std::runtime_error("cannot read " + path);
  }
  std::vector<std::pair<std::uint32_t, Exp::GbtsLayerConnection>> staged;
  for (std::uint32_t i = 0; i < n; ++i) {
    std::uint32_t index{}, stage{}, src{}, dst{}, height{}, width{}, entries{};
    if (!(in >> index >> stage >> src >> dst >> height >> width >> entries)) {
      throw std::runtime_error("bad connection header in " + path);
    }
    for (std::uint64_t b = 0; b < std::uint64_t{height} * width; ++b) {
      std::uint32_t unused{};
      in >> unused;
    }
    staged.push_back({stage, {src, dst}});
  }
  std::ranges::stable_sort(staged, {}, [](const auto& e) { return e.first; });
  std::vector<Exp::GbtsLayerConnection> out;
  for (const auto& e : staged) {
    out.push_back(e.second);
  }
  return out;
}


/// Learn a connection table from the generator truth, the way Athena trained
/// its own: the layer sequence of every primary above `minPt`, ordered away
/// from the origin with overlap hits merged, gives one link per adjacent pair
/// and one per pair skipping a layer. Links seen fewer than `minCount` times
/// are dropped.
std::map<std::pair<Exp::GbtsExperimentLayerId, Exp::GbtsExperimentLayerId>,
         std::size_t>
countTruthLinks(const Event& event,
                const std::vector<Exp::GbtsExperimentLayerId>& ids,
                float minPt, bool skipOne) {
  const auto particleColumn =
      event.spacePoints.column<std::uint32_t>("particleId");
  const auto layerColumn = event.spacePoints.column<std::uint32_t>("layerId");
  std::vector<std::vector<std::pair<float, std::uint32_t>>> hits(
      event.particles.size());
  for (const auto sp : event.spacePoints) {
    const std::uint32_t p = sp.extra(particleColumn);
    const auto& info = event.particles[p];
    if (!info.primary() || info.pt < minPt) {
      continue;
    }
    hits[p].emplace_back(std::hypot(sp.r(), sp.z()), sp.extra(layerColumn));
  }
  std::map<std::pair<Exp::GbtsExperimentLayerId, Exp::GbtsExperimentLayerId>,
           std::size_t>
      counts;
  for (auto& h : hits) {
    std::ranges::sort(h);
    std::vector<Exp::GbtsExperimentLayerId> seq;
    for (const auto& [d, layer] : h) {
      if (seq.empty() || seq.back() != ids[layer]) {
        seq.push_back(ids[layer]);
      }
    }
    for (std::size_t i = 0; i + 1 < seq.size(); ++i) {
      ++counts[{seq[i + 1], seq[i]}];
      if (skipOne && i + 2 < seq.size()) {
        ++counts[{seq[i + 2], seq[i]}];
      }
    }
  }
  return counts;
}

/// Write links in the Athena connection table format, one 1x1 bin table each.
void writeConnectionTable(
    const std::string& path,
    const std::vector<std::pair<Exp::GbtsExperimentLayerId,
                                Exp::GbtsExperimentLayerId>>& links,
    float etaBinWidth) {
  std::ofstream out(path);
  out << links.size() << " " << etaBinWidth << "\n";
  std::size_t index = 0;
  for (const auto& [src, dst] : links) {
    out << index++ << " 1 " << src << " " << dst << " 1 1 100\n100\n";
  }
}

/// GBTS reads the layer of a space point from a 16 bit column of its own name,
/// while the generator writes a 32 bit one. Copy it over once per event,
/// outside the timed region; `checkEncodingFits` makes sure that nothing is cut
/// off.
/// @param spacePoints the space points of the event
void addGbtsLayerIndex(Acts::SpacePointContainer& spacePoints) {
  auto layerIds = spacePoints.column<std::uint32_t>("layerId");
  auto gbtsLayers =
      spacePoints.createColumn<Exp::GbtsLayerIndex>("gbtsLayerIndex");
  for (auto sp : spacePoints) {
    sp.extra(gbtsLayers) = static_cast<Exp::GbtsLayerIndex>(sp.extra(layerIds));
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  SyntheticEventOptions::Values shared;
  std::size_t numRuns = 10;
  float minPt = 900.f;
  // Efficiency is counted over a harder threshold than the seeder is cut at, so
  // that the turn-on stays out of it.
  float truthPt = 1000.f;
  bool verbose = false;
  float maxSeedSplitEta = Exp::GraphBasedTrackSeeder::Config{}.maxSeedSplitEta;
  float maxInvRadDiffPerM = 0.7e-2f;
  int minSeedLevel = 3;
  bool addTriplets = false;
  float maxEtaAddTriplets = 1.5f;
  bool validateTriplets = true;
  std::string connectionTable;
  bool adaptiveCuts = true;
  std::string trainTable;
  std::size_t trainEvents = 100;
  std::size_t trainMinCount = 1000;
  float trainMinPt = 900.f;
  bool trainSkipOne = true;
  // The graph is capped at this many edges per space point. The synthetic event
  // is noise free, so more geometric doublet candidates survive the cuts than
  // in a real event; the default still leaves room for all of them on the ITk
  // at pile-up 200. Benchmarking against a saturated ceiling would hide any
  // work done past it, so lower it only to measure the truncated regime.
  float maxEdgesPerSP = Exp::GbtsGraphBuilder::Config{}.maxEdgesPerSP;

  try {
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message");
    SyntheticEventOptions::add(desc, shared);
    desc.add_options()("runs",
                       po::value<std::size_t>(&numRuns)->default_value(numRuns),
                       "number of benchmark runs")(
        "min-pt", po::value<float>(&minPt)->default_value(minPt),
        "seed momentum threshold in MeV")(
        "truth-pt", po::value<float>(&truthPt)->default_value(truthPt),
        "momentum in MeV a primary has to reach to be counted in the "
        "efficiency")(
        "max-edges-per-sp",
        po::value<float>(&maxEdgesPerSP)->default_value(maxEdgesPerSP),
        "ceiling on the number of graph edges per space point")(
        "max-seed-split-eta",
        po::value<float>(&maxSeedSplitEta)->default_value(maxSeedSplitEta),
        "seed split |eta| reach, 0 disables")(
        "max-inv-rad-diff",
        po::value<float>(&maxInvRadDiffPerM)->default_value(maxInvRadDiffPerM),
        "seed split curvature tolerance in 1/m")(
        "min-seed-level",
        po::value<int>(&minSeedLevel)->default_value(minSeedLevel),
        "minimum chain length in edges")(
        "add-triplets", po::value<bool>(&addTriplets)->default_value(addTriplets),
        "accept 3 SP seeds")(
        "max-eta-add-triplets",
        po::value<float>(&maxEtaAddTriplets)->default_value(maxEtaAddTriplets),
        "|eta| reach of add-triplets")(
        "validate-triplets",
        po::value<bool>(&validateTriplets)->default_value(validateTriplets),
        "barrel triplet pT/d0 validation")(
        "adaptive-cuts", po::value<bool>(&adaptiveCuts)->default_value(adaptiveCuts),
        "widen tau ratio cuts for skipped layers as Athena does")(
        "train-connection-table", po::value<std::string>(&trainTable),
        "learn a connection table from the truth of --train-events events "
        "(random seeds 1001...), write it here and exit")(
        "train-events", po::value<std::size_t>(&trainEvents)->default_value(trainEvents),
        "events to learn the connection table on")(
        "train-min-count",
        po::value<std::size_t>(&trainMinCount)->default_value(trainMinCount),
        "times a link has to be seen to be kept")(
        "train-min-pt", po::value<float>(&trainMinPt)->default_value(trainMinPt),
        "primaries to learn from, in MeV")(
        "train-skip-one",
        po::value<bool>(&trainSkipOne)->default_value(trainSkipOne),
        "also learn links that skip one layer")(
        "connection-table",
        po::value<std::string>(&connectionTable),
        "connection table in the Athena format, read with Athena layer ids; "
        "defaults to <layout>-gbts-connections.txt when it exists, none for "
        "the hand-written one")(
        "verbose", po::bool_switch(&verbose),
        "log the seeder's own statistics");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help") > 0) {
      std::cout << desc << std::endl;
      return 0;
    }
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }

  DetectorLayout layout;
  EventConfig eventConfig;
  try {
    layout = SyntheticEventOptions::makeLayout(shared);
    eventConfig = SyntheticEventOptions::makeConfig(shared);
    checkEncodingFits(layout);
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }
  if (!trainTable.empty()) {
    const auto ids = athenaLayerIds(layout);
    std::map<std::pair<Exp::GbtsExperimentLayerId, Exp::GbtsExperimentLayerId>,
             std::size_t>
        total;
    for (std::size_t k = 0; k < trainEvents; ++k) {
      EventConfig config = eventConfig;
      config.seed = static_cast<std::uint32_t>(1001 + k);
      for (const auto& [link, n] :
           countTruthLinks(generateEvent(layout, config), ids,
                           trainMinPt * 1_MeV / 1_GeV, trainSkipOne)) {
        total[link] += n;
      }
    }
    std::vector<std::pair<Exp::GbtsExperimentLayerId,
                          Exp::GbtsExperimentLayerId>>
        links;
    for (const auto& [link, n] : total) {
      if (n >= trainMinCount) {
        links.push_back(link);
      }
    }
    writeConnectionTable(trainTable, links, 0.2f);
    std::cout << "wrote " << links.size() << " of " << total.size()
              << " links to " << trainTable << std::endl;
    return 0;
  }
  if (connectionTable.empty()) {
    if (const auto shipped = layoutPath(shared, "-gbts-connections.txt");
        std::filesystem::exists(shipped)) {
      connectionTable = shipped.string();
    }
  }
  if (connectionTable == "none") {
    connectionTable.clear();
  }
  Event event = generateEvent(layout, eventConfig);
  applyHitInefficiency(event, shared.seed);
  addGbtsLayerIndex(event.spacePoints);
  const Acts::SpacePointContainer& spacePoints = event.spacePoints;

  float etaBinWidth = 0.2f;
  std::vector<Exp::GbtsLayerDescription> layers = makeLayerDescriptions(layout);
  std::vector<Exp::GbtsLayerConnection> connections;
  if (connectionTable.empty()) {
    connections = makeConnections(layout);
  } else {
    const auto ids = athenaLayerIds(layout);
    for (std::size_t i = 0; i < layers.size(); ++i) {
      layers[i].id = ids[i];
    }
    const auto all = readConnectionTable(connectionTable, etaBinWidth);
    std::unordered_set<Exp::GbtsExperimentLayerId> known(ids.begin(),
                                                         ids.end());
    for (const auto& c : all) {
      if (known.contains(c.src) && known.contains(c.dst)) {
        connections.push_back(c);
      }
    }
    std::cout << "connection table: " << connections.size() << " of "
              << all.size() << " links, eta bin width " << etaBinWidth
              << std::endl;
  }
  if (adaptiveCuts) {
    Exp::assignAdaptiveTauRatioCuts(connections, layers);
  }
  const auto geometry =
      std::make_shared<Exp::GbtsGeometry>(layers, connections, etaBinWidth);

  const auto makeLogger = [&](const std::string& name) {
    // Warnings are let through even when quiet: the graph builder stops once
    // `maxEdgesPerSP` is reached and says so at that level, and a truncated
    // graph costs efficiency without failing.
    return Acts::getDefaultLogger(name, verbose
                                            ? Acts::Logging::Level::DEBUG
                                            : Acts::Logging::Level::WARNING);
  };

  // The cuts are the ACTS defaults, which ATLAS runs the ITk pixel detector
  // with: `ActsTrk::GbtsSeedingTool` overrides only the connector file, the ML
  // lookup table and minPt. Two of its defaults differ and are set below.
  Exp::GraphBasedTrackSeeder::Config cfg;
  // ATLAS runs with the cluster width and its trained lookup table, which the
  // synthetic event has no cluster shapes to offer. It only adjusts the barrel
  // cot(theta) cuts, which is not where this event loses anything.
  cfg.useClusterWidthCuts = false;
  cfg.maxSeedSplitEta = maxSeedSplitEta;
  cfg.maxInvRadDiff = maxInvRadDiffPerM / 1_m;
  cfg.minSeedLevel = static_cast<std::uint8_t>(minSeedLevel);
  cfg.addTriplets = addTriplets;
  cfg.maxAbsEtaAddTriplets = maxEtaAddTriplets;

  Exp::GbtsGraphBuilder::Config graphCfg;
  graphCfg.minPt = minPt * 1_MeV;
  graphCfg.maxOuterRadius = 550.f;
  graphCfg.matchBeforeCreate = true;
  graphCfg.maxEdgesPerSP = maxEdgesPerSP;
  graphCfg.validateTriplets = validateTriplets;

  const Exp::GraphBasedTrackSeeder seeder(
      Exp::GraphBasedTrackSeeder::DerivedConfig(cfg), geometry,
      makeLogger("Gbts"));
  const Exp::GbtsGraphBuilder graphBuilder(graphCfg, geometry,
                                           makeLogger("GbtsGraphBuilder"));
  const Exp::GbtsTrackingFilter filter(Exp::GbtsTrackingFilter::Config{},
                                       geometry, makeLogger("GbtsFilter"));
  const Exp::GbtsRoiDescriptor roi(-5., 5., -200.f, 200.f);
  const Exp::GraphBasedTrackSeeder::Options options{.bFieldInZ =
                                                        eventConfig.bFieldZ};

  const EventSummary summary = summarize(event, truthPt * 1_MeV / 1_GeV);
  std::cout << "layers=" << layout.layers.size()
            << " etaBins=" << geometry->numBins()
            << " binGroups=" << geometry->binGroups().size() << "\n"
            << "spacePoints=" << summary.spacePoints
            << " primaryHits=" << summary.primaryHits
            << " secondaryHits=" << summary.secondaryHits
            << " primaries=" << summary.primaries
            << " secondaries=" << summary.secondaries
            << " seedable=" << summary.seedablePrimaries << std::endl;

  // matched outside the timed region, so that the truth lookup does not show
  // up in the measurement
  Acts::SeedContainer reference;
  reference.assignSpacePointContainer(spacePoints);
  seeder.createSeeds(spacePoints, roi, graphBuilder, filter, options,
                     reference);
  const SeedingSummary seedSummary =
      evaluateSeeds(event, reference, truthPt * 1_MeV / 1_GeV);

  std::cout << "seeds=" << seedSummary.seeds
            << " trueSeeds=" << seedSummary.trueSeeds << " efficiency="
            << static_cast<float>(seedSummary.matchedPrimaries) /
                   static_cast<float>(
                       std::max<std::size_t>(1, summary.seedablePrimaries))
            << std::endl;

  printSeedEtaReport(event, reference, truthPt * 1_MeV / 1_GeV, "gbts");
  if (numRuns == 0) {
    return 0;
  }

  // What an experiment calls: node building, graph and seed extraction in one.
  const auto full = microBenchmark(
      [&] {
        Acts::SeedContainer seeds;
        seeds.assignSpacePointContainer(spacePoints);
        seeder.createSeeds(spacePoints, roi, graphBuilder, filter, options,
                           seeds);
        assumeRead(seeds);
      },
      1, numRuns);
  std::cout << "full: " << full << std::endl;

  // The node storage fill alone, which copies every space point into the
  // storage and bins it. The graph stage is what the full call adds to it.
  const auto layerColumn =
      spacePoints.column<Exp::GbtsLayerIndex>("gbtsLayerIndex");
  const auto clusterWidthColumn = spacePoints.column<float>("clusterWidth");
  const auto localPositionColumn = spacePoints.column<float>("localPositionY");
  const auto nodes = microBenchmark(
      [&] {
        Exp::GbtsNodeStorage nodeStorage = seeder.makeNodeStorage();
        nodeStorage.extend(spacePoints, layerColumn, clusterWidthColumn,
                           localPositionColumn);
        nodeStorage.finalize();
        assumeRead(nodeStorage.numberOfNodes());
      },
      1, numRuns);
  std::cout << "nodes only: " << nodes << std::endl;

  return 0;
}

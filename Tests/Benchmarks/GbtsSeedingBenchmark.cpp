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
/// calls, and it includes building the graph nodes, which allocates a node
/// vector per layer inside the timed region. The second hoists the node
/// building out, so that a change to the graph stage cannot hide behind that
/// allocation.
///
/// The cuts are the ACTS defaults, which are what ATLAS runs on the ITk pixel
/// detector. The layer connection table below is not: it is written for these
/// layouts rather than trained on them, so the efficiency here is not a
/// statement about GBTS as ATLAS runs it. What it loses sits at the
/// barrel-endcap transition, where a hand-written table is weakest.
///
/// `--space-points` picks which of the two collections the event holds to seed
/// on. GBTS is a pixel seeder; running it on the strips of a detector, or on
/// both at once, is what the strip work is measuring against.
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

#include <cstdint>
#include <iostream>
#include <memory>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "SyntheticEventOptions.hpp"

namespace po = boost::program_options;
namespace Exp = Acts::Experimental;

using namespace Acts::UnitLiterals;
using namespace ActsTests;
using namespace ActsFatras::Synthetic;

namespace {

/// Rings a disc may carry before the endcap encoding below runs out of room.
constexpr std::uint32_t kRingsPerDisc = 8;

/// How many discs back along z the connection table lets a disc reach. One per
/// ring set of the ITk pixel endcap, see `makeConnectionTable`.
constexpr std::size_t kDiscReach = 9;

/// One logical layer of the detector: a barrel cylinder or an endcap disc with
/// every module of it. This is what a connection links; the `DetectorLayer`s a
/// layout holds are the modules.
struct LayerGroup {
  /// Which endcap it belongs to, or the barrel
  SurfaceSide side{};
  /// Which subsystem it belongs to
  std::uint16_t subsystem{};
  /// Radius of a cylinder and `|z|` of a disc: what the groups are ordered
  /// outwards by
  float refCoord{};
  /// How far along z a cylinder reaches, which is where the discs past it
  /// begin; unused for a disc
  float halfLengthZ{};
  /// GBTS ids of its modules
  std::vector<std::int32_t> ids;
};

/// The GBTS numbering of a layout.
///
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
///
/// A layer index belongs to its subsystem, so the pixels and the strips of one
/// detector both have a layer zero and neither the numbering nor the connection
/// table can be written in terms of it. Both are written in terms of the layout
/// as a whole instead, cylinders ordered outwards by radius and the discs of
/// each endcap outwards by `|z|`.
struct LayerNumbering {
  /// GBTS id of every layer of the layout, by its index
  std::vector<std::int32_t> ids;
  /// Its logical layers, each side ordered outwards
  std::vector<LayerGroup> groups;
};

/// Number the layers of a layout, and check that the encoding holds: the
/// barrel group numbers have to stay in the eighties for `layerId / 10000 == 8`
/// to mean barrel, and both sub-indices have to fit below the thousands digit.
/// `GbtsNode::layer` is also 16 bits, which bounds the number of layers.
///
/// @param layout the layout to number
/// @return the numbering
/// @throws std::invalid_argument if the layout does not fit the encoding
LayerNumbering numberLayers(const DetectorLayout& layout) {
  if (layout.layers.size() > 65535) {
    throw std::invalid_argument("GbtsNode::layer is 16 bits wide");
  }

  // one group per (subsystem, side, layer), which the modules of a surface
  // share
  using Key = std::tuple<SurfaceSide, std::uint16_t, std::uint32_t>;
  std::map<Key, LayerGroup> byKey;
  for (const DetectorLayer& layer : layout.layers) {
    if (layer.moduleIndex >= kRingsPerDisc &&
        layer.shape == SurfaceShape::Disc) {
      throw std::invalid_argument(
          "the GBTS layer id encoding takes at most eight rings per disc");
    }
    LayerGroup& group = byKey[Key{layer.side, layer.subsystem, layer.layer}];
    group.side = layer.side;
    group.subsystem = layer.subsystem;
    group.refCoord = std::abs(layer.refCoord);
    if (layer.shape == SurfaceShape::Cylinder) {
      group.halfLengthZ = std::max(
          group.halfLengthZ,
          std::max(std::abs(layer.minBound), std::abs(layer.maxBound)));
    }
  }

  // each side outwards, which is what the ids and the table are written in
  std::vector<LayerGroup*> ordered;
  ordered.reserve(byKey.size());
  for (auto& [key, group] : byKey) {
    ordered.push_back(&group);
  }
  std::ranges::stable_sort(
      ordered, [](const LayerGroup* a, const LayerGroup* b) {
        return std::tie(a->side, a->refCoord) < std::tie(b->side, b->refCoord);
      });
  std::map<const LayerGroup*, int> position;
  std::map<SurfaceSide, int> counts;
  for (LayerGroup* group : ordered) {
    position[group] = counts[group->side]++;
  }
  if (counts[SurfaceSide::Barrel] > 10) {
    throw std::invalid_argument(
        "the GBTS layer id encoding takes at most ten barrel layers");
  }

  LayerNumbering numbering;
  numbering.ids.reserve(layout.layers.size());
  for (const DetectorLayer& layer : layout.layers) {
    LayerGroup& group = byKey.at(Key{layer.side, layer.subsystem, layer.layer});
    const auto moduleIdx = static_cast<std::int32_t>(layer.moduleIndex);
    const int index = position.at(&group);
    std::int32_t id{};
    if (layer.shape == SurfaceShape::Cylinder) {
      if (moduleIdx >= 1000) {
        throw std::invalid_argument(
            "the GBTS layer id encoding takes at most a thousand modules per "
            "cylinder");
      }
      id = (80 + index) * 1000 + moduleIdx;
    } else {
      const std::int32_t subIndex =
          index * static_cast<std::int32_t>(kRingsPerDisc) + moduleIdx;
      if (subIndex >= 1000) {
        throw std::invalid_argument(
            "the GBTS layer id encoding takes at most a hundred and "
            "twenty-five discs per side");
      }
      id = (layer.side == SurfaceSide::Positive ? 90 : 70) * 1000 + subIndex;
    }
    numbering.ids.push_back(id);
    group.ids.push_back(id);
  }

  numbering.groups.reserve(ordered.size());
  for (const LayerGroup* group : ordered) {
    numbering.groups.push_back(*group);
  }
  return numbering;
}

std::vector<Exp::GbtsLayerDescription> makeLayerDescriptions(
    const DetectorLayout& layout, const LayerNumbering& numbering) {
  std::vector<Exp::GbtsLayerDescription> descriptions;
  descriptions.reserve(layout.layers.size());
  for (std::size_t index = 0; index < layout.layers.size(); ++index) {
    const DetectorLayer& layer = layout.layers[index];
    Exp::GbtsLayerDescription description;
    description.id = numbering.ids[index];
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
/// layer pair here only says that the pair is allowed at all, and a pair the
/// geometry cannot serve costs nothing beyond the time to reject it.
std::string makeConnectionTable(const LayerNumbering& numbering,
                                float etaBinWidth) {
  // The endcaps are taken one subsystem at a time, because `kDiscReach` counts
  // the ring sets of an endcap: measuring it over two endcaps interleaved in z
  // would give each of them less reach than it needs.
  std::vector<const LayerGroup*> barrel;
  std::map<SurfaceSide, std::map<std::uint16_t, std::vector<const LayerGroup*>>>
      endcaps;
  for (const LayerGroup& group : numbering.groups) {
    if (group.side == SurfaceSide::Barrel) {
      barrel.push_back(&group);
    } else {
      endcaps[group.side][group.subsystem].push_back(&group);
    }
  }

  std::vector<std::pair<std::int32_t, std::int32_t>> connections;
  auto connect = [&](const LayerGroup& outer, const LayerGroup& inner) {
    for (const std::int32_t src : outer.ids) {
      for (const std::int32_t dst : inner.ids) {
        connections.emplace_back(src, dst);
      }
    }
  };

  // barrel to barrel, adjacent and one skipped layer. Ordered by radius across
  // the subsystems, so the pixel to strip transition falls out of it.
  for (std::size_t i = 0; i < barrel.size(); ++i) {
    for (std::size_t step : {1u, 2u}) {
      if (i + step < barrel.size()) {
        connect(*barrel[i + step], *barrel[i]);
      }
    }
  }
  for (const auto& [side, sequences] : endcaps) {
    for (const auto& [subsystem, discs] : sequences) {
      // Disc to disc. Discs are ordered outwards along z but consecutive ones
      // do not sit at consecutive radii, a layout describing the rings of an
      // endcap interleaving the ring sets, so `kDiscReach` has to cover the
      // number of ring sets rather than just the neighbours.
      for (std::size_t j = 0; j < discs.size(); ++j) {
        for (std::size_t step = 1; step <= kDiscReach; ++step) {
          if (j + step < discs.size()) {
            connect(*discs[j + step], *discs[j]);
          }
        }
      }
      // Forward transition: a cylinder feeds the discs just past the end of it,
      // as many of them as the reach covers. A disc short of that end sits at a
      // smaller radius than the cylinder does and is inward of it, not outward.
      for (const LayerGroup* cylinder : barrel) {
        std::size_t taken = 0;
        for (const LayerGroup* disc : discs) {
          if (disc->refCoord > cylinder->halfLengthZ && taken++ < kDiscReach) {
            connect(*disc, *cylinder);
          }
        }
      }
    }
    // One endcap to the next along z, the strips of a detector picking up where
    // its pixels leave off. The same reach, counted into the other endcap.
    for (const auto& [outerSubsystem, outer] : sequences) {
      for (const auto& [innerSubsystem, inner] : sequences) {
        if (outerSubsystem == innerSubsystem) {
          continue;
        }
        for (const LayerGroup* disc : outer) {
          std::size_t taken = 0;
          for (const LayerGroup* below : std::ranges::reverse_view(inner)) {
            if (below->refCoord < disc->refCoord && taken++ < kDiscReach) {
              connect(*disc, *below);
            }
          }
        }
      }
    }
  }

  std::ostringstream os;
  os << connections.size() << " " << etaBinWidth << "\n";
  for (std::size_t k = 0; k < connections.size(); ++k) {
    // index stage src dst height width entries
    os << k << " " << k << " " << connections[k].first << " "
       << connections[k].second << " 0 0 0\n";
  }
  return os.str();
}

/// @param name what `--space-points` was given
/// @return the selection it names
/// @throws std::invalid_argument if it names none
SpacePointSelection parseSelection(const std::string& name) {
  if (name == "pixel") {
    return SpacePointSelection::Pixel;
  }
  if (name == "strip") {
    return SpacePointSelection::Strip;
  }
  if (name == "combined") {
    return SpacePointSelection::Combined;
  }
  throw std::invalid_argument("--space-points takes pixel, strip or combined");
}

/// How far out the graph has to reach, i.e. the radius of the outermost layer.
/// @param layout the layout
/// @return the radius
float outerRadius(const DetectorLayout& layout) {
  float radius = 0.f;
  for (const DetectorLayer& layer : layout.layers) {
    radius = std::max(radius, layer.shape == SurfaceShape::Cylinder
                                  ? layer.refCoord
                                  : layer.maxBound);
  }
  return radius;
}

}  // namespace

int main(int argc, char* argv[]) {
  SyntheticEventOptions::Values shared;
  std::size_t numRuns = 10;
  float minPt = 900.f;
  // Efficiency is counted over a harder threshold than the seeder is cut at, so
  // that the turn-on stays out of it.
  float truthPt = 1000.f;
  std::string selectionName = "pixel";
  bool calibrateStrips = true;
  float stripTolerance = 1.1f;
  bool verbose = false;
  // The synthetic event is noise free, so more geometric doublet candidates
  // survive the cuts than in a real event and the graph outgrows the 2000000
  // edges allowed by default. Benchmarking against a saturated ceiling would
  // hide any work done past it; pass the production value to measure the
  // truncated regime instead.
  std::uint32_t maxEdges = 6000000;

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
        "max-edges",
        po::value<std::uint32_t>(&maxEdges)->default_value(maxEdges),
        "ceiling on the number of graph edges")(
        "space-points",
        po::value<std::string>(&selectionName)->default_value(selectionName),
        "which of the event's space points to seed on: pixel, strip or "
        "combined")(
        "calibrate-strips",
        po::value<bool>(&calibrateStrips)->default_value(calibrateStrips),
        "resolve a strip endpoint against the direction of the doublet it is "
        "part of before cutting on it")(
        "strip-tolerance",
        po::value<float>(&stripTolerance)->default_value(stripTolerance),
        "how far outside its strips a crossing may be recovered")(
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
  LayerNumbering numbering;
  SpacePointSelection selection{};
  try {
    layout = SyntheticEventOptions::makeLayout(shared);
    eventConfig = SyntheticEventOptions::makeConfig(shared);
    numbering = numberLayers(layout);
    selection = parseSelection(selectionName);
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }
  const Event event = generateEvent(layout, eventConfig);
  // The seeder takes one container, so the two collections the generator emits
  // are gathered into the one the selection asks for.
  const Acts::SpacePointContainer spacePoints =
      selectSpacePoints(event, selection);

  constexpr float etaBinWidth = 0.2f;
  const std::string table = makeConnectionTable(numbering, etaBinWidth);
  std::istringstream tableStream(table);
  const Exp::GbtsLayerConnectionMap connections =
      Exp::GbtsLayerConnectionMap::fromStream(tableStream, false);
  const auto geometry = std::make_shared<Exp::GbtsGeometry>(
      makeLayerDescriptions(layout, numbering), connections);

  // The cuts are the ACTS defaults, which ATLAS runs the ITk pixel detector
  // with: `ActsTrk::GbtsSeedingTool` overrides only the connector file, the ML
  // lookup table and minPt. Two of its defaults differ and are set below; the
  // third, beamSpotCorrection, is read by nothing in the ACTS seeder.
  Exp::GraphBasedTrackSeeder::Config cfg;
  cfg.minPt = minPt * 1_MeV;
  cfg.minZ0 = -200.f;
  cfg.maxZ0 = 200.f;
  // How far out the graph reaches, which the ITk pixel default of 550 mm is the
  // measure of; the strips reach twice that.
  cfg.maxOuterRadius = outerRadius(layout);
  // ATLAS runs with the cluster width and its trained lookup table, which the
  // synthetic event has no cluster shapes to offer. It only adjusts the barrel
  // cot(theta) cuts, which is not where this event loses anything.
  cfg.useClusterWidthCuts = false;
  cfg.matchBeforeCreate = true;
  cfg.nMaxEdges = maxEdges;
  cfg.calibrateStrips = calibrateStrips;
  cfg.stripLengthTolerance = stripTolerance;

  const Exp::GraphBasedTrackSeeder::DerivedConfig derived(cfg);
  // Warnings are let through even when quiet: the seeder stops building the
  // graph once `nMaxEdges` is reached and says so at that level, and a
  // truncated graph costs efficiency without failing.
  const Exp::GraphBasedTrackSeeder seeder(
      derived, geometry,
      Acts::getDefaultLogger("Gbts", verbose ? Acts::Logging::Level::DEBUG
                                             : Acts::Logging::Level::WARNING));
  const Exp::GbtsTrackingFilter filter(Exp::GbtsTrackingFilter::Config{},
                                       geometry);
  const Exp::GbtsRoiDescriptor roi(-5., 5., cfg.minZ0, cfg.maxZ0);
  const Exp::GraphBasedTrackSeeder::Options options(eventConfig.bFieldZ);
  // A layer the description gives a strip sensor is read out as a stereo pair;
  // the rest are pixels.
  std::vector<bool> isPixelLayer(layout.layers.size());
  for (std::size_t index = 0; index < layout.layers.size(); ++index) {
    isPixelLayer[index] = !layout.layers[index].sensor.has_value();
  }

  const EventSummary summary =
      summarize(event, truthPt * 1_MeV / 1_GeV, selection);
  std::cout << "layers=" << layout.layers.size()
            << " spacePoints=" << spacePoints.size() << "\n"
            << "generated=" << summary.spacePoints
            << " strip=" << summary.stripSpacePoints
            << " primaryHits=" << summary.primaryHits
            << " secondaryHits=" << summary.secondaryHits
            << " primaries=" << summary.primaries
            << " secondaries=" << summary.secondaries
            << " seedable=" << summary.seedablePrimaries << std::endl;

  // matched outside the timed region, so that the truth lookup does not show
  // up in the measurement
  Acts::SeedContainer reference;
  reference.assignSpacePointContainer(spacePoints);
  seeder.createSeeds(spacePoints, roi, isPixelLayer, filter, options,
                     reference);
  const SeedingSummary seedSummary =
      evaluateSeeds(event, reference, truthPt * 1_MeV / 1_GeV);

  std::cout << "seeds=" << seedSummary.seeds
            << " trueSeeds=" << seedSummary.trueSeeds << " efficiency="
            << static_cast<float>(seedSummary.matchedPrimaries) /
                   static_cast<float>(
                       std::max<std::size_t>(1, summary.seedablePrimaries))
            << std::endl;

  // What an experiment calls: node building, graph and seed extraction in one.
  const auto full = microBenchmark(
      [&] {
        Acts::SeedContainer seeds;
        seeds.assignSpacePointContainer(spacePoints);
        seeder.createSeeds(spacePoints, roi, isPixelLayer, filter, options,
                           seeds);
        assumeRead(seeds);
      },
      1, numRuns);
  std::cout << "full: " << full << std::endl;

  // The same without the node building, which sorts every space point into its
  // eta and phi bin. Hoisting it out is what makes a change to the graph stage
  // visible.
  Exp::GbtsNodeStorage nodeStorage = seeder.makeNodeStorage(isPixelLayer);
  nodeStorage.extend(spacePoints, spacePoints.column<std::uint32_t>("layerId"),
                     spacePoints.column<float>("clusterWidth"),
                     spacePoints.column<float>("localPositionY"));
  nodeStorage.finalize();
  const auto graph = microBenchmark(
      [&] {
        Acts::SeedContainer seeds;
        seeds.assignSpacePointContainer(spacePoints);
        seeder.createSeeds(nodeStorage, roi, filter, options, seeds);
        assumeRead(seeds);
      },
      1, numRuns);
  std::cout << "graph only: " << graph << std::endl;

  return 0;
}

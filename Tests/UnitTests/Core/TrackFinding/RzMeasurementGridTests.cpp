// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/TrackFinding/Rz/RzLayout.hpp"
#include "Acts/TrackFinding/Rz/RzMeasurementGrid.hpp"

#include <cstdint>
#include <map>
#include <vector>

using namespace Acts;
using namespace Acts::Experimental;

namespace {

/// Five modules of which the odd ones are polar, which is all the grid reads
/// of a layout: how many modules there are and which of them carry a frame
/// per measurement.
RzLayout makeLayout() {
  RzLayout layout;
  layout.layers.resize(1);
  for (std::uint32_t i = 0; i < 5; ++i) {
    RzModule module;
    module.center = Vector3(i, 0., 0.);
    module.u = Vector3::UnitX();
    module.v = Vector3::UnitY();
    module.normal = Vector3::UnitZ();
    module.polar = (i % 2) == 1;
    module.layer = 0;
    layout.modules.push_back(module);
  }
  return layout;
}

bool isPolar(std::uint32_t module) {
  return (module % 2) == 1;
}

/// What the grid should hold: the sources added to each module, in order
using Added = std::map<std::uint32_t, std::vector<std::uint32_t>>;

/// Add one measurement, tagging both it and its frame with the source so that
/// a permutation that loses their pairing shows up
void add(RzMeasurementGrid& grid, std::uint32_t module, std::uint32_t source,
         Added& added) {
  RzMeasurement measurement;
  measurement.loc0 = source;
  measurement.source = source;
  RzMeasurementFrame frame;
  frame.u = Vector3(source, 0., 0.);
  frame.v = Vector3::UnitY();
  frame.normal = Vector3::UnitZ();
  BOOST_CHECK_EQUAL(grid.add(module, measurement, frame),
                    added[module].size());
  added[module].push_back(source);
}

void checkHolds(const RzMeasurementGrid& grid, const Added& added) {
  for (const auto& [module, sources] : added) {
    const RzModuleMeasurements on = grid.moduleRange(module);
    BOOST_TEST_CONTEXT("module " << module) {
      BOOST_REQUIRE_EQUAL(on.entries.size(), sources.size());
      BOOST_CHECK_EQUAL(!on.frames.empty(), isPolar(module));
      for (std::size_t i = 0; i < sources.size(); ++i) {
        BOOST_CHECK_EQUAL(on.entries[i].source, sources[i]);
        if (isPolar(module)) {
          BOOST_CHECK_EQUAL(on.frames[i].u.x(),
                            static_cast<double>(sources[i]));
        }
      }
    }
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(RzMeasurementGridSuite)

// A caller that adds each module's measurements in one run leaves nothing for
// `finalize` to do, whatever order it visits the modules in.
BOOST_AUTO_TEST_CASE(RunsInAnyModuleOrder) {
  const RzLayout layout = makeLayout();
  RzMeasurementGrid grid(layout);
  Added added;
  std::uint32_t source = 0;
  for (const std::uint32_t module : {3u, 0u, 4u, 1u}) {
    for (int i = 0; i < 3; ++i) {
      add(grid, module, source++, added);
    }
  }
  grid.finalize();

  BOOST_CHECK_EQUAL(grid.size(), 12u);
  checkHolds(grid, added);
  // a module nothing was added to is empty, not missing
  BOOST_CHECK(grid.moduleRange(2).entries.empty());
}

// The same measurements with every module split across the fill: `finalize`
// has to group them, and must keep each module's own order and carry the
// frames along with the entries they belong to.
BOOST_AUTO_TEST_CASE(InterleavedModulesAreGrouped) {
  const RzLayout layout = makeLayout();
  RzMeasurementGrid grid(layout);
  Added added;
  std::uint32_t source = 0;
  for (int i = 0; i < 3; ++i) {
    for (const std::uint32_t module : {3u, 0u, 4u, 1u}) {
      add(grid, module, source++, added);
    }
  }
  grid.finalize();

  BOOST_CHECK_EQUAL(grid.size(), 12u);
  checkHolds(grid, added);
  BOOST_CHECK(grid.moduleRange(2).entries.empty());
}

// A whole module at once, which is one run by construction
BOOST_AUTO_TEST_CASE(AddRange) {
  const RzLayout layout = makeLayout();
  RzMeasurementGrid grid(layout);

  std::vector<RzMeasurement> entries(2);
  entries[0].source = 7;
  entries[1].source = 8;
  std::vector<RzMeasurementFrame> frames(2);
  frames[0].u = Vector3(7., 0., 0.);
  frames[1].u = Vector3(8., 0., 0.);
  grid.addRange(1, entries, frames);
  grid.addRange(0, std::span(entries).first(1), {});
  grid.finalize();

  const RzModuleMeasurements polar = grid.moduleRange(1);
  BOOST_REQUIRE_EQUAL(polar.entries.size(), 2u);
  BOOST_REQUIRE_EQUAL(polar.frames.size(), 2u);
  BOOST_CHECK_EQUAL(polar.entries[1].source, 8u);
  BOOST_CHECK_EQUAL(polar.frames[1].u.x(), 8.);
  const RzModuleMeasurements cartesian = grid.moduleRange(0);
  BOOST_REQUIRE_EQUAL(cartesian.entries.size(), 1u);
  BOOST_CHECK(cartesian.frames.empty());
}

// The grid is filled once per event, so what one event left must not reach
// the next
BOOST_AUTO_TEST_CASE(ClearForgetsTheEvent) {
  const RzLayout layout = makeLayout();
  RzMeasurementGrid grid(layout);
  Added added;
  std::uint32_t source = 0;
  for (const std::uint32_t module : {1u, 3u, 1u}) {
    add(grid, module, source++, added);
  }
  grid.finalize();
  checkHolds(grid, added);

  grid.clear();
  BOOST_CHECK_EQUAL(grid.size(), 0u);
  BOOST_CHECK(grid.moduleRange(1).entries.empty());

  Added again;
  source = 100;
  for (const std::uint32_t module : {2u, 1u}) {
    add(grid, module, source++, again);
    add(grid, module, source++, again);
  }
  grid.finalize();
  BOOST_CHECK_EQUAL(grid.size(), 4u);
  checkHolds(grid, again);
}

// The accessor is what the finder is handed, and it has to see the same
// grouping the grid does
BOOST_AUTO_TEST_CASE(AccessorSeesTheSameGrid) {
  const RzLayout layout = makeLayout();
  RzMeasurementGrid grid(layout);
  Added added;
  std::uint32_t source = 0;
  for (int i = 0; i < 2; ++i) {
    for (const std::uint32_t module : {4u, 1u}) {
      add(grid, module, source++, added);
    }
  }
  grid.finalize();

  const RzMeasurementAccessor accessor = grid.accessor();
  for (const auto& [module, sources] : added) {
    const RzModuleMeasurements on = accessor(module);
    BOOST_REQUIRE_EQUAL(on.entries.size(), sources.size());
    for (std::size_t i = 0; i < sources.size(); ++i) {
      BOOST_CHECK_EQUAL(on.entries[i].source, sources[i]);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

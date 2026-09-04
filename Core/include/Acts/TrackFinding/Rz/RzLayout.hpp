// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// The RZ skeleton of a tracker the finder navigates: cylinders and discs,
/// each with the material it carries and, if sensitive, the modules that hang
/// off it. Built once from a Gen1 tracking geometry, whose layers already are
/// cylinders and discs.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/MaterialSlab.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

namespace Acts {
class Surface;
class TrackingGeometry;
class Logger;
}  // namespace Acts

namespace Acts::Experimental {

/// Stands for "none" wherever an index is optional.
constexpr std::uint32_t kRzNone = std::numeric_limits<std::uint32_t>::max();

/// The two shapes the layout is made of.
enum class RzShape : std::uint8_t {
  /// At a fixed radius, extending along z
  Cylinder,
  /// At a fixed z, extending in r
  Disc,
};

/// The effects of one material band on one particle species, tabulated
/// against the momentum so that a stop costs a lookup rather than the Bethe
/// and Highland formulas. Values are for the band's own thickness; a crossing
/// at an angle scales them by the path factor, the scattering with its
/// logarithmic correction kept exact.
struct RzMaterialTable {
  static constexpr std::uint32_t kBins = 48;
  /// Logarithmic momentum grid
  static constexpr double kMinP = 0.1 * UnitConstants::GeV;
  static constexpr double kMaxP = 200. * UnitConstants::GeV;

  std::array<float, kBins> theta0Sq{};
  std::array<float, kBins> energyLoss{};
  std::array<float, kBins> sigmaQOverPSq{};
  /// `ln(t / X0)` of the band, for the Highland factor at a path factor
  float logThicknessInX0{};

  static double logMinP() { return std::log(kMinP); }
  static double logStep() { return std::log(kMaxP / kMinP) / (kBins - 1); }
};

/// One RZ surface: a navigation stop, carrying material and possibly a layer.
struct RzSurface {
  RzShape shape{};
  /// Radius of a cylinder, signed z of a disc
  double refCoord{};
  /// Extent along the coordinate the surface extends in: signed z for a
  /// cylinder, r for a disc
  double minBound{};
  double maxBound{};
  /// Band edges along the extended coordinate, one more than there are bands;
  /// empty for a surface without material
  std::vector<double> materialEdges;
  /// What each band is made of
  std::vector<MaterialSlab> materialBands;
  /// The bands' effects tabulated, one per band, or empty
  std::vector<RzMaterialTable> materialTables;
  /// `Bz` along the surface, one value per `fieldBinWidth` from `minBound`,
  /// averaged over azimuth; empty for a constant field
  std::vector<double> bzTable;
  double fieldBinWidth{};

  /// The field at a crossing, or nothing if the surface carries no table
  /// @param along z on a cylinder, r on a disc
  /// @return `Bz`
  std::optional<double> bzAt(double along) const {
    if (bzTable.empty()) {
      return std::nullopt;
    }
    const auto bin =
        static_cast<std::int64_t>((along - minBound) / fieldBinWidth);
    const std::size_t i = static_cast<std::size_t>(std::clamp<std::int64_t>(
        bin, 0, static_cast<std::int64_t>(bzTable.size()) - 1));
    return bzTable[i];
  }
  /// Index into `RzLayout::layers` if sensitive
  std::uint32_t layer{kRzNone};
  /// Where it came from
  GeometryIdentifier geometryId;
  /// The surface it came from, for track states that need one
  std::shared_ptr<const Surface> surface;

  /// Whether a crossing at a position along the surface lands on it
  /// @param along z on a cylinder, r on a disc
  /// @return true if within the bounds
  bool contains(double along) const {
    return minBound <= along && along <= maxBound;
  }

  /// The material a crossing meets, or nothing
  /// @param along z on a cylinder, r on a disc
  /// @return the slab, or nullptr for none
  const MaterialSlab* materialAt(double along) const;

  /// The band a crossing meets
  /// @param along z on a cylinder, r on a disc
  /// @return the band index, or -1 for none
  int materialBandAt(double along) const;
};

/// A sensitive module: a plane with its frame, which is all the finder needs
/// to project a free state onto a measurement.
struct RzModule {
  Vector3 center{Vector3::Zero()};
  /// Local axes and normal in the global frame
  Vector3 u{Vector3::Zero()};
  Vector3 v{Vector3::Zero()};
  Vector3 normal{Vector3::Zero()};
  /// Half extents of the bounding box in the local frame; `halfV` is also
  /// the coordinate a strip does not measure
  double halfU{};
  double halfV{};
  /// Local coordinates are polar (r, phi) in the surface frame: an annulus
  bool polar{false};
  std::uint32_t layer{kRzNone};
  GeometryIdentifier geometryId;
  std::shared_ptr<const Surface> surface;
};

/// A sensitive layer and the binning its measurements are stored in.
struct RzLayer {
  std::uint32_t surface{kRzNone};
  std::uint32_t phiBins{};
  std::uint32_t alongBins{};
  double alongMin{};
  double alongMax{};
  /// Largest `RzModule::halfV` of the layer's modules, the window a strip
  /// search has to open along the strip
  double maxHalfV{};
  /// Largest offset of a module centre from the RZ surface: in r for a
  /// cylinder, in z for a disc. Staggered, inclined or double-sided modules
  /// all show up here, and the search window opens by it times the slope.
  double halfThickness{};
  /// Largest half diagonal of a module, the room a module lookup by centre
  /// has to leave
  double maxHalfExtent{};
  /// How far from the RZ surface a module of this layer may be met: the
  /// finder's own limit, or three half thicknesses for a layer whose
  /// modules spread (the ITk inclined section)
  double moduleDistance{};
  /// First global bin of this layer
  std::uint32_t binOffset{};

  double phiBinWidth() const;
  double alongBinWidth() const { return (alongMax - alongMin) / alongBins; }
};

struct RzLayout {
  std::vector<RzSurface> surfaces;
  /// Indices of the cylinders by increasing radius
  std::vector<std::uint32_t> cylinders;
  /// Indices of the discs by increasing z
  std::vector<std::uint32_t> discs;
  std::vector<RzLayer> layers;
  std::vector<RzModule> modules;
  std::unordered_map<GeometryIdentifier, std::uint32_t> moduleIndex;
  /// Total number of measurement bins over all layers
  std::uint32_t totalBins{};
  /// The modules binned by their centres in the layers' measurement binning:
  /// bin `b` holds `moduleOrder[moduleBinStart[b] .. moduleBinStart[b + 1])`
  std::vector<std::uint32_t> moduleBinStart;
  std::vector<std::uint32_t> moduleOrder;

  /// Bin of a point on a layer, in the layer's binning
  /// @param layer the layer
  /// @param phi azimuth
  /// @param along z on a cylinder, r on a disc
  /// @return the global bin
  std::uint32_t bin(std::uint32_t layer, double phi, double along) const;
  /// Beyond either of these a track has left the tracker
  double escapeRadius{};
  double escapeHalfZ{};
};

struct RzLayoutOptions {
  /// Keep only sensitive surfaces for which this returns true; empty keeps
  /// everything. Passive material is kept regardless.
  std::function<bool(const Surface&)> surfaceSelector;
  /// Samples along the extended coordinate the material is read at
  std::uint32_t materialSamples = 50;
  /// Samples in azimuth averaged into each band
  std::uint32_t phiSamples = 8;
  /// Relative difference in x/X0 below which neighbouring samples are one band
  double materialBandTolerance = 0.2;
  /// Measurement bins in azimuth per layer
  std::uint32_t phiBins = 64;
  /// Measurement bin width along the extended coordinate
  double alongBinWidth = 20 * UnitConstants::mm;
  /// Tabulate each band's effects for this particle, see `RzMaterialTable`.
  /// Measured to buy nothing once the rest is in, so off.
  bool materialTables = false;
  ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();
  /// The field, sampled once per surface into `RzSurface::bzTable`; empty
  /// for a constant field
  std::function<Vector3(const Vector3&)> fieldSampler;
  double fieldBinWidth = 100 * UnitConstants::mm;
  /// Floor of `RzLayer::moduleDistance`
  double moduleDistance = 50 * UnitConstants::mm;
};

/// Reduce a Gen1 tracking geometry to its RZ skeleton.
/// @param trackingGeometry the geometry
/// @param gctx the geometry context
/// @param options steering
/// @param logger a logger
/// @return the layout
RzLayout makeRzLayout(const TrackingGeometry& trackingGeometry,
                      const GeometryContext& gctx,
                      const RzLayoutOptions& options, const Logger& logger);

}  // namespace Acts::Experimental

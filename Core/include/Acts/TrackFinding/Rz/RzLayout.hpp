// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// RZ navigation surfaces and modules built from tracking geometry.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/TrackFinding/Rz/RzTypes.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
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

/// The two shapes the layout is made of.
enum class RzShape : std::uint8_t {
  /// At a fixed radius, extending along z
  Cylinder,
  /// At a fixed z, extending in r
  Disc,
};

/// Momentum-tabulated material effects for one band and particle species.
struct RzMaterialTable {
  static constexpr std::uint32_t kBins = 48;
  /// Logarithmic momentum grid
  static constexpr double kMinP = 0.1 * UnitConstants::GeV;
  static constexpr double kMaxP = 200. * UnitConstants::GeV;

  std::array<float, kBins> theta0Sq{};
  std::array<float, kBins> energyLoss{};
  std::array<float, kBins> sigmaQOverPSq{};
  /// `ln(t / X0)` for the Highland correction.
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
  /// The radial field along the surface, `B . r_hat`, binned and averaged as
  /// `bzTable`; empty for a constant field
  std::vector<double> brTable;
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
  /// The radial field at a crossing, zero without a table
  /// @param along z on a cylinder, r on a disc
  /// @return `Br`
  double brAt(double along) const {
    if (brTable.empty()) {
      return 0.;
    }
    const auto bin =
        static_cast<std::int64_t>((along - minBound) / fieldBinWidth);
    const std::size_t i = static_cast<std::size_t>(std::clamp<std::int64_t>(
        bin, 0, static_cast<std::int64_t>(brTable.size()) - 1));
    return brTable[i];
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

  /// The band a crossing meets
  /// @param along z on a cylinder, r on a disc
  /// @return the band index, or -1 for none
  std::int32_t materialBandAt(double along) const;
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
  /// The bound coordinates are not the module's own axes: they are polar
  /// (r, phi) in the surface frame, which is an annulus disc
  bool polar{false};
  /// Plain polar coordinates, verified against the surface during layout.
  bool polarIsPlain{false};
  /// The module's centre in the surface's own cartesian frame, which is what
  /// a plain polar measurement is placed against
  Vector2 localCenter{Vector2::Zero()};
  /// The centre in the surface's own bound coordinates, where a measurement
  /// leaves a coordinate it does not measure
  Vector2 boundCenter{Vector2::Zero()};
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
  /// Contiguous disc coordinates and bounds for navigation probes.
  std::vector<double> discCoord;
  std::vector<double> discMin;
  std::vector<double> discMax;
  /// Contiguous cylinder radii for navigation probes.
  std::vector<double> cylCoord;
  std::vector<RzLayer> layers;
  std::vector<RzModule> modules;
  std::unordered_map<GeometryIdentifier, std::uint32_t> moduleIndex;
  /// Total number of measurement bins over all layers
  std::uint32_t totalBins{};
  /// The modules binned by their centres in the layers' measurement binning:
  /// bin `b` holds `moduleOrder[moduleBinStart[b] .. moduleBinStart[b + 1])`
  std::vector<std::uint32_t> moduleBinStart;
  std::vector<std::uint32_t> moduleOrder;

  /// Visit every global bin of a layer a window touches
  /// @param layer the layer
  /// @param phi azimuth of the point
  /// @param along z on a cylinder, r on a disc
  /// @param halfPhi half width of the window in azimuth
  /// @param halfAlong half width along
  /// @param visitor called with each global bin
  template <typename visitor_t>
  void visitBins(std::uint32_t layer, double phi, double along, double halfPhi,
                 double halfAlong, visitor_t&& visitor) const {
    const RzLayer& l = layers[layer];
    const double phiWidth = l.phiBinWidth();
    const double alongWidth = l.alongBinWidth();
    const std::int32_t phiLo =
        static_cast<std::int32_t>(std::floor((phi - halfPhi) / phiWidth));
    const std::int32_t phiHi =
        static_cast<std::int32_t>(std::floor((phi + halfPhi) / phiWidth));
    const std::int32_t nPhi = static_cast<std::int32_t>(l.phiBins);
    const std::int32_t alongLoRaw = static_cast<std::int32_t>(
        std::floor((along - halfAlong - l.alongMin) / alongWidth));
    const std::int32_t alongHiRaw = static_cast<std::int32_t>(
        std::floor((along + halfAlong - l.alongMin) / alongWidth));
    const std::int32_t alongLo = std::max(alongLoRaw, 0);
    const std::int32_t alongHi =
        std::min(alongHiRaw, static_cast<std::int32_t>(l.alongBins) - 1);
    if (alongLo > alongHi) {
      return;
    }
    // a window wider than the full circle visits each bin once
    const std::int32_t phiCount = std::min(phiHi - phiLo + 1, nPhi);
    for (std::int32_t i = 0; i < phiCount; ++i) {
      std::int32_t p = (phiLo + i) % nPhi;
      if (p < 0) {
        p += nPhi;
      }
      for (std::int32_t a = alongLo; a <= alongHi; ++a) {
        visitor(l.binOffset + static_cast<std::uint32_t>(p) * l.alongBins +
                static_cast<std::uint32_t>(a));
      }
    }
  }

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
  /// Tabulate material effects by momentum; useful for material-rich layouts.
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

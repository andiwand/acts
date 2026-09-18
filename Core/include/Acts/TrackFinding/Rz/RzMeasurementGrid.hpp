// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Event measurements indexed by module, with an optional custom accessor.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/TrackFinding/Rz/RzLayout.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <cstdint>
#include <span>
#include <vector>

namespace Acts::Experimental {

/// Which of the module's two local coordinates a measurement holds.
enum class RzProjector : std::uint8_t {
  /// A strip along `v`, measuring `u`
  Loc0,
  /// A strip along `u`, measuring `v`
  Loc1,
  /// A pixel
  Both,
};

/// Local offset and covariance on a module. Polar entries carry their own
/// frame; an unmeasured coordinate is zero at the module centre.
struct RzMeasurement {
  /// Offset from the module centre along the first axis
  double loc0{};
  /// Offset from the module centre along the second axis
  double loc1{};
  /// Variance of `loc0`
  double cov00{};
  /// Covariance, used only by a pixel
  double cov01{};
  /// Variance of `loc1`
  double cov11{};
  /// Inverse lever arm for a polar angular variance; zero otherwise.
  double invLever{};
  /// Optional measured time and variance; zero variance means no time value.
  float time{};
  float timeVariance{};
  RzProjector projector{RzProjector::Both};
  /// The caller's index of the measurement
  std::uint32_t source{kRzNone};
};

/// Per-measurement axes for polar modules, where each strip has its own frame.
struct RzMeasurementFrame {
  /// The direction `loc0` is measured along
  Vector3 u{Vector3::Zero()};
  /// The direction `loc1` is measured along
  Vector3 v{Vector3::Zero()};
  /// `u x v`
  Vector3 normal{Vector3::Zero()};
  /// Coordinates of these axes in the module frame.
  double uU{1.};
  double uV{0.};
  double vU{0.};
  double vV{1.};
};

/// The measurements on one module, and their frames if it is polar.
struct RzModuleMeasurements {
  std::span<const RzMeasurement> entries;
  /// One frame per entry for a polar module, empty for a cartesian one, whose
  /// entries all share the module's own axes
  std::span<const RzMeasurementFrame> frames;
};

/// Fetch measurements for one crossed module.
/// @param module index into `RzLayout::modules`
/// @return the measurements on that module, contiguous
using RzMeasurementAccessor =
    Delegate<RzModuleMeasurements(std::uint32_t module)>;

/// The measurements of an event, binned by module.
class RzMeasurementGrid {
 public:
  explicit RzMeasurementGrid(const RzLayout& layout);

  const RzLayout& layout() const { return *m_layout; }

  /// Drop the measurements of the last event
  void clear();

  /// Make room for a known number of measurements, so that filling does not
  /// grow and copy. A caller knows how many it is about to hand over.
  /// @param n the number of measurements about to be added
  void reserve(std::size_t n);

  /// Add one measurement already on the module's own axes. This is the
  /// cheapest form and does no arithmetic: a caller whose measurements are
  /// cartesian offsets on the module needs nothing else.
  /// @param module index into `RzLayout::modules`
  /// @param measurement the entry
  /// @param frame its own axes, for a polar module only
  /// @return its index within the module
  std::uint32_t add(std::uint32_t module, const RzMeasurement& measurement,
                    const RzMeasurementFrame& frame = {});

  /// Add a whole module's measurements at once. A caller whose container is
  /// grouped by module hands over one group here, which is one run by
  /// construction, and the grid never sorts.
  /// @param module index into `RzLayout::modules`
  /// @param measurements the entries
  /// @param frames their axes, for a polar module only
  void addRange(std::uint32_t module,
                std::span<const RzMeasurement> measurements,
                std::span<const RzMeasurementFrame> frames = {});

  /// Convert and add a measurement in the surface's bound coordinates.
  /// @param module index into `RzLayout::modules`
  /// @param surface the module's surface, whose bounds give the local frame
  /// @param gctx the geometry context
  /// @param dim 1 or 2
  /// @param boundIndices which bound coordinate each measured value is, 0 for
  ///        `eBoundLoc0` and 1 for `eBoundLoc1`; the ones not measured are
  ///        taken at the module centre
  /// @param boundParams the measured values
  /// @param boundCov the covariance in the bound coordinates, row major,
  ///        `dim` by `dim`
  /// @param source the caller's index of the measurement
  /// @return its index within the module
  std::uint32_t addBound(std::uint32_t module, const Surface& surface,
                         const GeometryContext& gctx, std::uint8_t dim,
                         std::span<const std::uint8_t> boundIndices,
                         std::span<const double> boundParams,
                         std::span<const double> boundCov, std::uint32_t source,
                         double time = 0., double timeVariance = 0.);

  /// Convert bound coordinates to a module-frame entry without adding it.
  /// @param layout the layout the module belongs to
  /// @param module index into `RzLayout::modules`
  /// @param surface the module's surface, whose bounds give the local frame
  /// @param gctx the geometry context
  /// @param dim 1 or 2
  /// @param boundIndices which bound coordinate each measured value is
  /// @param boundParams the measured values
  /// @param boundCov the covariance in the bound coordinates, row major
  /// @param source the caller's index of the measurement
  /// @param frame filled with the measurement's own axes if the module is
  ///        polar, left alone if it is not
  /// @return the entry
  static RzMeasurement fromBound(const RzLayout& layout, std::uint32_t module,
                                 const Surface& surface,
                                 const GeometryContext& gctx, std::uint8_t dim,
                                 std::span<const std::uint8_t> boundIndices,
                                 std::span<const double> boundParams,
                                 std::span<const double> boundCov,
                                 std::uint32_t source,
                                 RzMeasurementFrame& frame);

  /// Group entries by module if they were added out of order.
  void finalize();

  std::size_t size() const { return m_entries.size(); }

  /// The measurements on one module
  /// @param module index into `RzLayout::modules`
  /// @return the entries, contiguous, with their frames if the module is polar
  RzModuleMeasurements moduleRange(std::uint32_t module) const {
    const ModuleBlock& block = m_blocks[module];
    RzModuleMeasurements r;
    r.entries = {m_entries.data() + block.begin, block.size};
    if (block.frame != kRzNone) {
      r.frames = {m_frames.data() + block.frame, block.size};
    }
    return r;
  }

  /// An accessor onto this grid, to hand to the finder
  /// @return the accessor
  RzMeasurementAccessor accessor() const {
    RzMeasurementAccessor a;
    a.connect<&RzMeasurementGrid::moduleRange>(this);
    return a;
  }

  /// Visit every global bin of a layer a window touches
  /// @param layout the layout
  /// @param layer the layer
  /// @param phi azimuth of the point
  /// @param along z on a cylinder, r on a disc
  /// @param halfPhi half width of the window in azimuth
  /// @param halfAlong half width along
  /// @param visitor called with each global bin
  template <typename visitor_t>
  static void visitBins(const RzLayout& layout, std::uint32_t layer, double phi,
                        double along, double halfPhi, double halfAlong,
                        visitor_t&& visitor) {
    const RzLayer& l = layout.layers[layer];
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

 private:
  /// Where one module's measurements sit in the grid. The search asks for a
  /// module and gets a span, so all it needs is where the module's run starts
  /// and how long it is — not that the modules themselves are in any order.
  struct ModuleBlock {
    /// Index of the module's first measurement
    std::uint32_t begin{};
    /// How many it holds
    std::uint32_t size{};
    /// Index of its first frame, `kRzNone` for a cartesian module, whose
    /// measurements all share the module's own axes
    std::uint32_t frame{kRzNone};
  };

  const RzLayout* m_layout{};
  std::vector<RzMeasurement> m_entries;
  /// One per measurement on a polar module, in the order they were added
  std::vector<RzMeasurementFrame> m_frames;
  std::vector<std::uint32_t> m_moduleOf;
  std::vector<ModuleBlock> m_blocks;
  /// Whether every module's measurements have arrived in one run, so that the
  /// blocks already describe them and nothing has to move
  bool m_contiguous{true};
};

}  // namespace Acts::Experimental

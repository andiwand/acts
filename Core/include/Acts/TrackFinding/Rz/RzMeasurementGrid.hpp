// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// The measurements of one event, held per module of an `RzLayout` and in the
/// module's own frame. The finder walks to a stop, asks the layout which
/// modules the track crosses there, and asks for those modules' measurements —
/// so the index is by module, not by a window in the layer.
///
/// A caller reaches the search either by filling `RzMeasurementGrid`, one
/// measurement or one module at a time, or by handing the finder its own
/// `RzMeasurementAccessor`, which is asked once per module the search crosses.

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

/// One measurement as the finder sees it: an offset from the module centre
/// along a pair of axes in the module's plane, and the variance of that
/// offset. A module is a plane, so two numbers place the measurement on it and
/// the module supplies the rest — the normal, the extent and the distance the
/// search may meet it at.
///
/// The axes are the module's own, `RzModule::u` and `RzModule::v`, unless the
/// module is polar, in which case each measurement carries its own
/// `RzMeasurementFrame`. A coordinate the measurement does not hold sits at
/// zero, the module centre, which is where a strip's extent is measured from.
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
  /// For a measurement whose bound frame is polar, the reciprocal of the lever
  /// arm its variance was converted with: what the angle stands for in length
  /// grows with the distance from the frame's origin, and the entry sits at
  /// its own radius while the track crosses at another. Zero for a cartesian
  /// frame, which has no such dependence.
  double invLever{};
  RzProjector projector{RzProjector::Both};
  /// The caller's index of the measurement
  std::uint32_t source{kRzNone};
};

/// The axes of a measurement whose bound coordinates are not the module's own.
/// An annulus strip measures an azimuth, so its frame turns with the strip and
/// no two strips on the module share one.
struct RzMeasurementFrame {
  /// The direction `loc0` is measured along
  Vector3 u{Vector3::Zero()};
  /// The direction `loc1` is measured along
  Vector3 v{Vector3::Zero()};
  /// `u x v`
  Vector3 normal{Vector3::Zero()};
};

/// The measurements on one module, and their frames if it is polar.
struct RzModuleMeasurements {
  std::span<const RzMeasurement> entries;
  /// One frame per entry for a polar module, empty for a cartesian one, whose
  /// entries all share the module's own axes
  std::span<const RzMeasurementFrame> frames;
};

/// Where the search gets a module's measurements. It asks once per module it
/// crosses, so an accessor that converts from a framework's own container pays
/// that conversion only for the modules a track reaches, and never for the
/// rest of the event.
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
  /// grouped by module hands over one group here, and the grid never sorts.
  /// @param module index into `RzLayout::modules`
  /// @param measurements the entries
  /// @param frames their axes, for a polar module only
  void addRange(std::uint32_t module,
                std::span<const RzMeasurement> measurements,
                std::span<const RzMeasurementFrame> frames = {});

  /// Add a measurement given in the surface's own bound coordinates, whatever
  /// local frame that is. A cartesian frame needs only the module centre,
  /// which the layout read off the surface once. A polar frame (an annulus
  /// strip measures an azimuth) needs the bounds' map to the cartesian frame,
  /// and this is the only path that touches the surface at all.
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
                         std::span<const double> boundCov,
                         std::uint32_t source);

  /// Convert a measurement in the surface's own bound coordinates into the
  /// module's frame, without adding it. A caller that keeps its own container
  /// and serves the search through an `RzMeasurementAccessor` uses this to
  /// build its entries.
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

  /// Group what was added by module. A grid filled in module order, or filled
  /// only through `addRange`, is already grouped and this does nothing.
  void finalize();

  std::size_t size() const { return m_entries.size(); }

  /// The measurements on one module
  /// @param module index into `RzLayout::modules`
  /// @return the entries, contiguous, with their frames if the module is polar
  RzModuleMeasurements moduleRange(std::uint32_t module) const {
    const std::uint32_t begin = m_moduleStart[module];
    const std::uint32_t end = m_moduleStart[module + 1];
    RzModuleMeasurements r;
    r.entries = {m_entries.data() + begin, m_entries.data() + end};
    if (!m_frames.empty() && m_layout->modules[module].polar) {
      r.frames = {m_frames.data() + begin, m_frames.data() + end};
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
    const int phiLo = static_cast<int>(std::floor((phi - halfPhi) / phiWidth));
    const int phiHi = static_cast<int>(std::floor((phi + halfPhi) / phiWidth));
    const int nPhi = static_cast<int>(l.phiBins);
    const int alongLoRaw = static_cast<int>(
        std::floor((along - halfAlong - l.alongMin) / alongWidth));
    const int alongHiRaw = static_cast<int>(
        std::floor((along + halfAlong - l.alongMin) / alongWidth));
    const int alongLo = std::max(alongLoRaw, 0);
    const int alongHi = std::min(alongHiRaw, static_cast<int>(l.alongBins) - 1);
    if (alongLo > alongHi) {
      return;
    }
    // a window wider than the full circle visits each bin once
    const int phiCount = std::min(phiHi - phiLo + 1, nPhi);
    for (int i = 0; i < phiCount; ++i) {
      int p = (phiLo + i) % nPhi;
      if (p < 0) {
        p += nPhi;
      }
      for (int a = alongLo; a <= alongHi; ++a) {
        visitor(l.binOffset + static_cast<std::uint32_t>(p) * l.alongBins +
                static_cast<std::uint32_t>(a));
      }
    }
  }

 private:
  void padFrames();

  const RzLayout* m_layout{};
  std::vector<RzMeasurement> m_entries;
  /// Parallel to `m_entries`, filled only once a polar module is added
  std::vector<RzMeasurementFrame> m_frames;
  std::vector<std::uint32_t> m_moduleOf;
  std::vector<std::uint32_t> m_moduleStart;
  /// Whether every add so far named a module at least as large as the one
  /// before, so that the entries are already grouped
  bool m_grouped{true};
  std::uint32_t m_lastModule{0};
};

}  // namespace Acts::Experimental

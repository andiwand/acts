// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// The measurements of one event in global coordinates, held per module of an
/// `RzLayout`. The finder walks to a stop, asks the layout which modules the
/// track crosses there, and looks at those modules' measurements — so the
/// index is by module, not by a window in the layer.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/TrackFinding/Rz/RzLayout.hpp"

#include <cstdint>
#include <span>
#include <vector>

namespace Acts::Experimental {

/// One measurement as the finder sees it: a point on a module plane and the
/// local frame to take a residual in.
struct RzMeasurement {
  /// The measured point in the global frame; an unmeasured coordinate sits at
  /// the module centre
  Vector3 position{Vector3::Zero()};
  /// Local axes of the module, the residual is taken along `u` and `v`
  Vector3 u{Vector3::Zero()};
  Vector3 v{Vector3::Zero()};
  Vector3 normal{Vector3::Zero()};
  /// Local covariance, `cov01` and `cov11` unused for a strip
  double cov00{};
  double cov01{};
  double cov11{};
  /// Half extent of a strip along `v`
  double halfV{};
  /// How far from the RZ stop the module may be met, from its layer
  double maxDistance{};
  /// 1 for a strip measuring `u`, 2 for a pixel
  std::uint8_t dim{};
  std::uint32_t module{kRzNone};
  /// The caller's index of the measurement
  std::uint32_t source{kRzNone};
};

class RzMeasurementGrid {
 public:
  explicit RzMeasurementGrid(const RzLayout& layout);

  const RzLayout& layout() const { return *m_layout; }

  /// Drop the measurements of the last event
  void clear();

  /// Add a measurement on a module.
  /// @param module index into `RzLayout::modules`
  /// @param dim 1 or 2
  /// @param localIndices which local coordinate each measured value is, 0 for
  ///        `u` and 1 for `v`
  /// @param localParams the measured values
  /// @param localCov the covariance, row major, `dim` by `dim`
  /// @param source the caller's index of the measurement
  void add(std::uint32_t module, std::uint8_t dim,
           std::span<const std::uint8_t> localIndices,
           std::span<const double> localParams,
           std::span<const double> localCov, std::uint32_t source);

  /// Add a measurement given in the global frame: the point, the unit
  /// directions the residual is taken along and the covariance in length
  /// units along them. This is the form for a module whose local
  /// coordinates are not cartesian (an annulus strip measures an angle) and
  /// what an experiment framework with its own local-to-global would call.
  /// @param module index into `RzLayout::modules`
  /// @param dim 1 or 2
  /// @param position the measured point, the unmeasured coordinate at the
  ///        module centre
  /// @param u the direction of the (first) measured coordinate
  /// @param v the second one, or for a strip the unmeasured one in the plane
  /// @param cov00 variance along `u`
  /// @param cov01 covariance, unused for a strip
  /// @param cov11 variance along `v`, unused for a strip
  /// @param source the caller's index of the measurement
  void add(std::uint32_t module, std::uint8_t dim, const Vector3& position,
           const Vector3& u, const Vector3& v, double cov00, double cov01,
           double cov11, std::uint32_t source);

  /// Add a measurement given in the surface's own bound coordinates, whatever
  /// local frame that is: the surface maps the point to the global frame and
  /// its bounds supply the map from the bound coordinates to the cartesian
  /// ones, so a polar frame (an annulus strip) needs nothing of the caller.
  /// This is the form for an experiment framework whose measurements are
  /// bound parameters on an ACTS surface.
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
  void addBound(std::uint32_t module, const Surface& surface,
                const GeometryContext& gctx, std::uint8_t dim,
                std::span<const std::uint8_t> boundIndices,
                std::span<const double> boundParams,
                std::span<const double> boundCov, std::uint32_t source);

  /// Sort what was added into the bins. Nothing can be searched before.
  void finalize();

  std::size_t size() const { return m_entries.size(); }

  const RzMeasurement& entry(std::uint32_t index) const {
    return m_entries[index];
  }

  /// The measurements on one module, by index into `entry`
  /// @param module index into `RzLayout::modules`
  /// @return the indices
  std::span<const std::uint32_t> moduleRange(std::uint32_t module) const {
    return {m_order.data() + m_moduleStart[module],
            m_order.data() + m_moduleStart[module + 1]};
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
  const RzLayout* m_layout{};
  std::vector<RzMeasurement> m_entries;
  std::vector<std::uint32_t> m_moduleOf;
  std::vector<std::uint32_t> m_moduleStart;
  std::vector<std::uint32_t> m_order;
};

}  // namespace Acts::Experimental

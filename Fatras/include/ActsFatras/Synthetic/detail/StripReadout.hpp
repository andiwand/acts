// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Reading a crossing of a strip layer out as the stereo pair it really is.
/// See @ref fatras_synthetic_events.

#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/detail/GeneratorScratch.hpp"

#include <array>
#include <optional>
#include <random>
#include <vector>

namespace ActsFatras::Synthetic::detail {

/// What one layer's strip modules look like, resolved once per layout so that
/// the hot loop never looks a subsystem name up. A pixel layer leaves `strip`
/// false and the rest of it alone.
struct StripLayer {
  /// Whether the layer is read out as a stereo pair at all
  bool strip{false};
  /// Cosine of half the stereo angle
  float cosHalfStereo{1.f};
  /// Sine of the same
  float sinHalfStereo{};
  /// Half the distance between the two sensors of a module
  float halfGap{};
  /// Resolution of one strip across itself, the pitch over the square root of
  /// twelve
  float sigma{};
  /// Half the length of a strip
  float halfLength{};
  /// How far outside the strips a crossing may still be recovered
  float gapTolerance{};
  /// Variance of the resolved point along the strip, which is where the
  /// projection error lands: z on a barrel, r on an endcap
  float varianceAlong{};
  /// Variance across it, where the two strips measure nearly one coordinate
  float varianceAcross{};
};

/// What every layer of a layout is read out as, indexed the way
/// `DetectorLayout::layers` is: the trig and the variances its `StripSensor`
/// implies, resolved once so that the hot loop never derives them again.
///
/// @param layout the detector
/// @return one entry per layer, empty if no layer of it carries a sensor
std::vector<StripLayer> stripLayers(const DetectorLayout& layout);

/// Read a crossing of a strip layer out as a stereo pair.
///
/// Two strip lines are laid through the crossing at plus and minus half the
/// stereo angle, one on either sensor of the module, each displaced across
/// itself by its own measurement error. Where the two appear to cross is then
/// resolved the way reconstruction resolves it, from the beam spot rather than
/// from the track -- which is what makes the error a projection error rather
/// than noise, and what leaves the pair able to give the crossing back to
/// anything that knows the direction better.
///
/// @param rng the random engine
/// @param layer the strip parameters of the layer crossed
/// @param cylinder whether the layer is a cylinder, whose strips run along z;
///        an endcap disc's run radially
/// @param position where the track crosses the middle of the module
/// @param direction the direction it crosses in, of unit length
/// @return the space point and the pair it came from, or nothing where the pair
///         does not resolve, which is a strip module measuring nothing
std::optional<StripHit> readStrip(std::mt19937& rng, const StripLayer& layer,
                                  bool cylinder,
                                  const std::array<float, 3>& position,
                                  const std::array<float, 3>& direction);

}  // namespace ActsFatras::Synthetic::detail

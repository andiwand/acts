// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/detail/StripReadout.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/SpacePointFormation/StripSpacePointBuilder.hpp"
#include "ActsFatras/Synthetic/detail/Sampling.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <ranges>

namespace ActsFatras::Synthetic::detail {

namespace {

/// The square root of twelve, a strip measuring its pitch over that.
const float kSqrt12 = std::sqrt(12.f);

/// @return `a` plus `b` scaled by `scale`
constexpr std::array<float, 3> addScaled(const std::array<float, 3>& a,
                                         const std::array<float, 3>& b,
                                         const float scale) {
  return {a[0] + scale * b[0], a[1] + scale * b[1], a[2] + scale * b[2]};
}

/// @return the dot product of two vectors
constexpr float dot(const std::array<float, 3>& a,
                    const std::array<float, 3>& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/// @return the vector as the core library wants it
Acts::Vector3 vector3(const std::array<float, 3>& v) {
  return Acts::Vector3{v[0], v[1], v[2]};
}

}  // namespace

std::vector<StripLayer> stripLayers(const DetectorLayout& layout) {
  const bool any = std::ranges::any_of(
      layout.layers,
      [](const DetectorLayer& l) { return l.sensor.has_value(); });
  if (!any) {
    return {};
  }

  std::vector<StripLayer> out(layout.layers.size());
  for (std::size_t index = 0; index < layout.layers.size(); ++index) {
    const std::optional<StripSensor>& sensor = layout.layers[index].sensor;
    if (!sensor.has_value()) {
      continue;
    }

    StripLayer& strip = out[index];
    strip.strip = true;
    const float halfStereo = 0.5f * sensor->stereoAngle;
    strip.cosHalfStereo = std::cos(halfStereo);
    strip.sinHalfStereo = std::sin(halfStereo);
    strip.halfGap = 0.5f * sensor->moduleGap;
    strip.sigma = sensor->pitch / kSqrt12;
    strip.halfLength = sensor->halfLength;
    // The pair separates a displacement across the strips only to the stereo
    // angle, so a gap seen across them is a gap over its sine along them.
    const float sinStereo = 2.f * strip.sinHalfStereo * strip.cosHalfStereo;
    const float cosStereo = strip.cosHalfStereo * strip.cosHalfStereo -
                            strip.sinHalfStereo * strip.sinHalfStereo;
    strip.gapOverStereo = sensor->moduleGap * cosStereo / sinStereo;

    // Two strips at plus and minus half the stereo angle measure the same
    // point as u1 and u2 across themselves. Solving for the coordinate along
    // the strips and the one across them gives sigma/(sqrt(2) sin) and
    // sigma/(sqrt(2) cos): the pair measures across itself better than one
    // strip does and along itself very much worse.
    const float alongSigma =
        strip.sigma / (std::numbers::sqrt2_v<float> * strip.sinHalfStereo);
    const float acrossSigma =
        strip.sigma / (std::numbers::sqrt2_v<float> * strip.cosHalfStereo);
    strip.varianceAlong = alongSigma * alongSigma;
    strip.varianceAcross = acrossSigma * acrossSigma;
  }
  return out;
}

std::optional<StripHit> readStrip(std::mt19937& rng, const StripLayer& layer,
                                  const float gapParameter, const bool cylinder,
                                  const std::array<float, 3>& position,
                                  const std::array<float, 3>& direction) {
  const float r = std::hypot(position[0], position[1]);
  if (r <= 0.f) {
    return std::nullopt;
  }
  const float cosPhi = position[0] / r;
  const float sinPhi = position[1] / r;

  // The module frame. A barrel module faces outwards with its strips along z,
  // an endcap module faces along the beam with its strips along r; in both the
  // remaining in-plane axis is azimuthal, which is the one the pair measures
  // well.
  const std::array<float, 3> normal =
      cylinder ? std::array<float, 3>{cosPhi, sinPhi, 0.f}
               : std::array<float, 3>{0.f, 0.f, position[2] < 0.f ? -1.f : 1.f};
  const std::array<float, 3> along =
      cylinder ? std::array<float, 3>{0.f, 0.f, 1.f}
               : std::array<float, 3>{cosPhi, sinPhi, 0.f};
  const std::array<float, 3> across{
      normal[1] * along[2] - normal[2] * along[1],
      normal[2] * along[0] - normal[0] * along[2],
      normal[0] * along[1] - normal[1] * along[0]};

  // A track along the module plane resolves nothing, and would send the two
  // sensor crossings off to infinity.
  const float incidence = dot(direction, normal);
  constexpr float minIncidence = 0.05f;
  if (std::abs(incidence) < minIncidence) {
    return std::nullopt;
  }

  // The two sensors sit half the module gap either side of the middle, so the
  // track meets each of them somewhere else.
  Acts::StripSpacePointBuilder::StripEnds ends[2];
  std::array<float, 3> centre[2];
  std::array<float, 3> axis[2];
  for (int side = 0; side < 2; ++side) {
    const float sign = side == 0 ? -1.f : 1.f;
    // rotated about the normal, one sensor each way
    axis[side] = {
        layer.cosHalfStereo * along[0] + sign * layer.sinHalfStereo * across[0],
        layer.cosHalfStereo * along[1] + sign * layer.sinHalfStereo * across[1],
        layer.cosHalfStereo * along[2] + sign * layer.sinHalfStereo * across[2],
    };
    // where the track meets this sensor
    const std::array<float, 3> met =
        addScaled(position, direction, sign * layer.halfGap / incidence);
    // the direction the strip measures in, across itself and in its plane
    const std::array<float, 3> measured{
        layer.cosHalfStereo * across[0] - sign * layer.sinHalfStereo * along[0],
        layer.cosHalfStereo * across[1] - sign * layer.sinHalfStereo * along[1],
        layer.cosHalfStereo * across[2] - sign * layer.sinHalfStereo * along[2],
    };

    // A strip is a strip: it says where the track crossed it to the pitch and
    // says nothing at all about where along itself that was. So the centre is
    // the centre of whichever strip of the grid contains the crossing, not the
    // crossing -- putting the crossing there instead would hand the along-strip
    // coordinate straight back and leave the pair with nothing to resolve.
    const float length = 2.f * layer.halfLength;
    const float here = dot(met, along) / dot(axis[side], along);
    const float middle = (std::floor(here / length) + 0.5f) * length;
    centre[side] = addScaled(met, axis[side], middle - here);
    centre[side] =
        addScaled(centre[side], measured, layer.sigma * sampleNormal(rng));
    ends[side].top =
        vector3(addScaled(centre[side], axis[side], layer.halfLength));
    ends[side].bottom =
        vector3(addScaled(centre[side], axis[side], -layer.halfLength));
  }

  // Resolved from the beam spot, which is all reconstruction knows: the error
  // that leaves is the projection error, not noise.
  Acts::StripSpacePointBuilder::ConstrainedOptions options;
  options.vertex = Acts::Vector3::Zero();
  // How far the beam spot walks the point along the strip at the softest track
  // still wanted. A disc's two sensors are a gap apart along z, so what the
  // pair sees of that gap is what the track's slope makes of it.
  const float slope = cylinder ? 1.f : r / std::abs(position[2]);
  options.stripLengthGapTolerance =
      gapParameter * r * layer.gapOverStereo * slope;
  const Acts::Result<Acts::Vector3> resolved =
      Acts::StripSpacePointBuilder::computeConstrainedSpacePoint(
          ends[0], ends[1], options);
  if (!resolved.ok()) {
    return std::nullopt;
  }
  const Acts::Vector3& point = *resolved;

  StripHit out;
  out.hit.r = static_cast<float>(std::hypot(point.x(), point.y()));
  out.hit.phi = static_cast<float>(std::atan2(point.y(), point.x()));
  out.hit.z = static_cast<float>(point.z());
  // The nominal point sits on the first strip and the calibrated one on the
  // second, which is what `calibrateOuterStripSpacePoint` moves it onto.
  out.strips.outerCenter = centre[1];
  out.strips.innerToOuterSeparation = {centre[1][0] - centre[0][0],
                                       centre[1][1] - centre[0][1],
                                       centre[1][2] - centre[0][2]};
  out.strips.outerHalfVector = {axis[1][0] * layer.halfLength,
                                axis[1][1] * layer.halfLength,
                                axis[1][2] * layer.halfLength};
  out.strips.innerHalfVector = {axis[0][0] * layer.halfLength,
                                axis[0][1] * layer.halfLength,
                                axis[0][2] * layer.halfLength};
  return out;
}

}  // namespace ActsFatras::Synthetic::detail

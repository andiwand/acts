// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// Synthetic ITk-like pixel event generator shared by the seeding benchmarks.
///
/// The layout is a coarse stand-in for the ATLAS ITk pixel detector: five
/// barrel cylinders subdivided into eta modules and nine endcap disks per side
/// subdivided into radial rings.
///
/// The event has two components:
///
///   * *Primaries* -- the charged particles of a pile-up of minimum-bias
///     interactions. They are helices from the luminous region with a flat
///     pseudorapidity distribution and a soft-dominated momentum spectrum.
///   * *Secondaries* -- particles produced when a primary interacts with the
///     detector material. They are generated at the primary's own hits, so
///     their production points follow both the material and the particle flux
///     without needing a separate spatial model, and their yield per crossing
///     is parametrised in (r, z) by the amount of detector already traversed.
///
/// Both components are propagated through the same helix model, so a secondary
/// is simply a soft, displaced track that starts at a radius instead of at the
/// beam line. See `EventConfig` for the numbers the defaults are tuned to.
///
/// The generator is deliberately seeder-agnostic: it produces a plain
/// `Acts::SpacePointContainer` and the layer index each space point belongs to,
/// so it can drive any seeding implementation. Seeder-specific geometry (such
/// as the GBTS layer-connection table) is built from `DetectorLayout` by the
/// individual benchmarks.

#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <numbers>
#include <random>
#include <string>
#include <vector>

namespace ActsTests::SyntheticItk {

/// Whether a layer is a cylinder around the beam axis or a disk normal to it.
enum class LayerType { Barrel, Endcap };

/// One logical layer, i.e. the granularity a seeder reasons about. A barrel
/// cylinder is split into eta modules along z, an endcap disk into rings in r.
struct Layer {
  /// Barrel cylinder or endcap disk
  LayerType type{};
  /// Radius for a barrel layer, signed z position for an endcap disk
  float refCoord{};
  /// Lower bound along the extended coordinate: z for barrel, r for endcap
  float minBound{};
  /// Upper bound along the extended coordinate: z for barrel, r for endcap
  float maxBound{};
  /// +1 for the positive endcap, -1 for the negative endcap, 0 for the barrel
  int side{};
  /// Barrel layer index, or disk index counting outwards from the interaction
  /// point within one endcap
  int layer{};
  /// Eta module index for a barrel layer, ring index for an endcap disk
  int module{};
};

/// A physical detector element together with the logical layers it is split
/// into. Track propagation intersects surfaces; the resulting position then
/// selects one of the surface's modules.
///
/// A surface without modules is passive: it is made of material, so particles
/// interact in it, but it does not measure them. The beam pipe is one.
struct Surface {
  /// Barrel cylinder or endcap disk
  LayerType type{};
  /// Radius for a barrel cylinder, signed z position for an endcap disk
  float refCoord{};
  /// Indices into DetectorLayout::layers of the modules of this surface,
  /// empty for a passive surface
  std::vector<std::uint32_t> modules;
};

/// The full synthetic detector description.
struct DetectorLayout {
  /// All logical layers, in the index space used by the "layerId" column
  std::vector<Layer> layers;
  /// The physical surfaces the layers are grouped into
  std::vector<Surface> surfaces;
};

/// Steering for the synthetic event.
///
/// The defaults describe an ATLAS-like minimum-bias pile-up of 200. The
/// particle-level numbers they reproduce are, per interaction and per unit of
/// pseudorapidity, a charged multiplicity of 6.6 above 100 MeV, of which 45 %
/// are above 500 MeV and 18 % above 1 GeV. Those three numbers fix the
/// momentum spectrum below and give a mean transverse momentum of 0.64 GeV.
struct EventConfig {
  /// Number of minimum-bias interactions overlaid in the event
  std::size_t pileup = 200;
  /// Charged particles per interaction and unit of pseudorapidity above
  /// `minPt`
  float chargedPerUnitEta = 6.6f;

  /// Lowest generated transverse momentum in GeV. Below roughly 90 MeV a
  /// particle curls up before the outermost barrel layer, so little is lost by
  /// cutting here.
  float minPt = 0.1f;
  /// Momentum scale of the spectrum in GeV, see `samplePt`
  float ptScale = 4.7f;
  /// Falloff exponent of the spectrum, see `samplePt`
  float ptExponent = 11.f;

  /// Lower end of the flat pseudorapidity range
  float minEta = -4.f;
  /// Upper end of the flat pseudorapidity range
  float maxEta = 4.f;

  /// Longitudinal spread of the luminous region in mm
  float beamspotSigmaZ = 45.f;
  /// Transverse impact parameter spread of the primaries in mm. Dominated by
  /// the size of the luminous region, widened to stand in for the tail of
  /// tracks from heavy-flavour and strange decays that a single Gaussian
  /// cannot describe.
  float d0Sigma = 0.1f;

  /// Mean number of secondaries a primary produces at a crossing where
  /// `materialDepth` is one, i.e. the overall normalisation of the noise.
  ///
  /// This is an effective yield rather than an interaction probability: it
  /// also absorbs the extra clusters that this deliberately coarse layout
  /// cannot produce by itself, namely module overlaps, the inclined and ring
  /// sections of the real ITk, and particles curling below `minPt`. The
  /// default puts the event at the roughly 200k pixel space points an ITk
  /// ttbar event at a pile-up of 200 produces.
  float secondaryRate = 0.35f;
  /// Radius in mm over which the traversed material doubles in the barrel
  float materialRScale = 300.f;
  /// Longitudinal distance in mm over which the traversed material doubles in
  /// the endcaps
  float materialZScale = 3000.f;
  /// Lowest generated secondary transverse momentum in GeV
  float secondaryMinPt = 0.05f;
  /// Slope of the exponential secondary momentum spectrum in GeV
  float secondaryPtSlope = 0.15f;
  /// RMS opening angle in rad of a secondary against its parent
  float secondaryOpeningAngle = 0.15f;

  /// Gaussian smearing in mm along the measured directions of a sensor, i.e.
  /// the single space point resolution. A barrel cylinder measures r*phi and
  /// z, an endcap disk measures r*phi and r; the remaining direction is normal
  /// to the sensor and is left at the nominal surface position.
  float positionSmearing = 0.015f;
  /// Sensor thickness in mm. The direction normal to a sensor is not measured,
  /// so all the space point knows about it is that it lies within the sensor.
  float sensorThickness = 0.15f;
  /// Solenoid field along z in Tesla
  float bFieldZ = 2.f;
  /// Random seed, so that events are reproducible across runs and seeders
  unsigned int seed = 12345;

  /// Number of primary particles the pile-up settings amount to.
  /// @return the number of primaries to generate
  std::size_t numPrimaries() const {
    return static_cast<std::size_t>(static_cast<float>(pileup) *
                                    chargedPerUnitEta * (maxEta - minEta));
  }
};

/// A generated particle, kept so that the benchmarks can report the truth
/// content of the event without a truth-matching stage.
struct GeneratedParticle {
  /// Transverse momentum in GeV
  float pt{};
  /// Pseudorapidity
  float eta{};
  /// Transverse impact parameter with respect to the beam axis in mm
  float d0{};
  /// Longitudinal position at the perigee in mm
  float z0{};
  /// Radius at which the particle was produced in mm, zero for a primary
  float productionRadius{};
  /// Longitudinal position at which the particle was produced in mm, zero for
  /// a primary
  float productionZ{};
  /// Whether the particle comes from the luminous region or from a material
  /// interaction
  bool primary{};
  /// Number of space points the particle left in the detector
  std::uint32_t numHits{};
};

/// A generated event.
struct Event {
  /// The space points, see `generateEvent` for the columns they carry
  Acts::SpacePointContainer spacePoints;
  /// The generating particles, indexed by the "particleId" column
  std::vector<GeneratedParticle> particles;
};

/// Truth content of an event, so that the benchmarks can report what they are
/// running on.
struct EventSummary {
  /// Number of space points
  std::size_t spacePoints{};
  /// Number of generated primaries
  std::size_t primaries{};
  /// Number of generated secondaries
  std::size_t secondaries{};
  /// Primaries above the momentum threshold that leave enough space points to
  /// be seedable at all
  std::size_t seedablePrimaries{};
  /// Space points left by primaries
  std::size_t primaryHits{};
  /// Space points left by secondaries
  std::size_t secondaryHits{};
};

/// Summarise the truth content of an event.
/// @param event the generated event
/// @param ptThreshold the momentum threshold in GeV a primary must pass
/// @return the summary
inline EventSummary summarize(const Event& event, float ptThreshold) {
  EventSummary summary;
  summary.spacePoints = event.spacePoints.size();
  for (const GeneratedParticle& particle : event.particles) {
    if (particle.primary) {
      ++summary.primaries;
      summary.primaryHits += particle.numHits;
      if (particle.pt >= ptThreshold && particle.numHits >= 3) {
        ++summary.seedablePrimaries;
      }
    } else {
      ++summary.secondaries;
      summary.secondaryHits += particle.numHits;
    }
  }
  return summary;
}

/// What the seeds of an event look like against the generator truth. Cheap
/// enough to run alongside a benchmark, and it catches a seeder configuration
/// that has stopped finding anything long before the timings would.
struct SeedingSummary {
  /// Number of seeds
  std::size_t seeds{};
  /// Seeds whose three space points all come from one primary above the
  /// momentum threshold
  std::size_t trueSeeds{};
  /// Primaries above the threshold that have at least one true seed
  std::size_t matchedPrimaries{};
};

/// Match seeds against the generator truth.
/// @param event the generated event
/// @param seeds the seeds found on it, indexing `event.spacePoints`
/// @param ptThreshold the momentum threshold in GeV a primary must pass
/// @return the summary
inline SeedingSummary evaluateSeeds(const Event& event,
                                    const Acts::SeedContainer& seeds,
                                    float ptThreshold) {
  const auto particleColumn =
      event.spacePoints.column<std::uint32_t>("particleId");

  SeedingSummary summary;
  summary.seeds = seeds.size();
  std::vector<bool> matched(event.particles.size(), false);
  for (const auto seed : seeds) {
    const auto indices = seed.spacePointIndices();
    const std::uint32_t particle =
        event.spacePoints[indices[0]].extra(particleColumn);
    const bool sameParticle =
        std::ranges::all_of(indices, [&](Acts::SpacePointIndex index) {
          return event.spacePoints[index].extra(particleColumn) == particle;
        });
    if (!sameParticle) {
      continue;
    }
    const GeneratedParticle& info = event.particles[particle];
    if (!info.primary || info.pt < ptThreshold) {
      continue;
    }
    ++summary.trueSeeds;
    matched[particle] = true;
  }
  summary.matchedPrimaries =
      static_cast<std::size_t>(std::ranges::count(matched, true));
  return summary;
}

/// Write an event to two CSV files, `<prefix>_spacepoints.csv` and
/// `<prefix>_particles.csv`, for offline inspection of the distributions.
/// @param event the generated event
/// @param layout the detector the event was generated on
/// @param prefix the path prefix of the two files
inline void writeEventCsv(const Event& event, const DetectorLayout& layout,
                          const std::string& prefix) {
  std::ofstream spFile(prefix + "_spacepoints.csv");
  spFile << "x,y,z,r,phi,layer,barrel,particle,primary\n";
  const auto particleColumn =
      event.spacePoints.column<std::uint32_t>("particleId");
  const auto layerColumn = event.spacePoints.column<std::uint32_t>("layerId");
  for (const auto sp : event.spacePoints) {
    const std::uint32_t layer = sp.extra(layerColumn);
    const std::uint32_t particle = sp.extra(particleColumn);
    spFile << sp.x() << ',' << sp.y() << ',' << sp.z() << ',' << sp.r() << ','
           << sp.phi() << ',' << layer << ','
           << (layout.layers[layer].type == LayerType::Barrel ? 1 : 0) << ','
           << particle << ',' << (event.particles[particle].primary ? 1 : 0)
           << '\n';
  }

  std::ofstream particleFile(prefix + "_particles.csv");
  particleFile << "pt,eta,d0,z0,productionRadius,productionZ,primary,"
                  "numHits\n";
  for (const GeneratedParticle& particle : event.particles) {
    particleFile << particle.pt << ',' << particle.eta << ',' << particle.d0
                 << ',' << particle.z0 << ',' << particle.productionRadius
                 << ',' << particle.productionZ << ','
                 << (particle.primary ? 1 : 0) << ',' << particle.numHits
                 << '\n';
  }
}

/// Build the ITk-like pixel layout.
/// @return the detector layout used by the seeding benchmarks
inline DetectorLayout makePixelLayout() {
  // Barrel cylinder radii in mm and their common half-length in z
  constexpr int numBarrel = 5;
  constexpr float barrelRadius[numBarrel] = {34.f, 99.f, 160.f, 228.f, 291.f};
  constexpr float barrelHalfZ = 250.f;
  // The barrel is left unsplit so that the logical layer ids come out as the
  // round 80000, 81000, ... that GBTS special-cases as its innermost layers.
  constexpr int barrelModules = 1;

  // Endcap disk positions in mm and their radial extent
  constexpr int numDisks = 9;
  constexpr float diskZ[numDisks] = {400.f,  600.f,  800.f,  1050.f, 1300.f,
                                     1600.f, 1900.f, 2300.f, 2800.f};
  constexpr float diskRMin = 30.f;
  constexpr float diskRMax = 350.f;
  constexpr int diskRings = 2;

  // The beam pipe carries no readout, but it is the only material in front of
  // the innermost layer and hence the only source of secondaries there.
  constexpr float beamPipeRadius = 25.f;

  DetectorLayout layout;
  layout.surfaces.push_back(
      Surface{LayerType::Barrel, beamPipeRadius, /*modules=*/{}});

  auto addSurface = [&](LayerType type, float refCoord, int side, int layerIdx,
                        float lo, float hi, int nSplit) {
    Surface surface;
    surface.type = type;
    surface.refCoord = refCoord;
    for (int m = 0; m < nSplit; ++m) {
      Layer layer;
      layer.type = type;
      layer.refCoord = refCoord;
      layer.minBound = lo + (hi - lo) * m / nSplit;
      layer.maxBound = lo + (hi - lo) * (m + 1) / nSplit;
      layer.side = side;
      layer.layer = layerIdx;
      layer.module = m;
      surface.modules.push_back(
          static_cast<std::uint32_t>(layout.layers.size()));
      layout.layers.push_back(layer);
    }
    layout.surfaces.push_back(std::move(surface));
  };

  for (int i = 0; i < numBarrel; ++i) {
    addSurface(LayerType::Barrel, barrelRadius[i], 0, i, -barrelHalfZ,
               barrelHalfZ, barrelModules);
  }
  for (int side : {+1, -1}) {
    for (int j = 0; j < numDisks; ++j) {
      addSurface(LayerType::Endcap, side * diskZ[j], side, j, diskRMin,
                 diskRMax, diskRings);
    }
  }

  return layout;
}

/// Draw a transverse momentum in GeV from `dN/dpT ~ (1 + pT/ptScale)^-n`.
///
/// The form is Hagedorn-like: over the first few GeV it is indistinguishable
/// from an exponential of slope `ptScale / n`, about 430 MeV for the defaults,
/// while the power law keeps a tail that a pure exponential would cut off far
/// too early. Its cumulative distribution inverts in closed form, so a single
/// uniform draw suffices.
///
/// @param u a uniform random number in [0, 1)
/// @param minPt the lower end of the spectrum in GeV
/// @param ptScale the momentum scale of the spectrum in GeV
/// @param n the falloff exponent of the spectrum
/// @return the sampled transverse momentum in GeV
inline float samplePt(float u, float minPt, float ptScale, float n) {
  const float exponent = 1.f - n;
  const float lo = std::pow(1.f + minPt / ptScale, exponent);
  // the spectrum reaches to infinity, where the corresponding term vanishes
  return ptScale * (std::pow(lo * (1.f - u), 1.f / exponent) - 1.f);
}

/// How much detector material a particle has traversed to reach (r, z),
/// relative to what it traverses before the beam line.
///
/// Material interactions happen in the detector itself, so the probability of
/// producing a secondary at a crossing grows with the number of surfaces
/// already crossed: with r in the barrel and with |z| in the endcaps. The
/// linear form is the simplest one that keeps the forward region, where a
/// particle crosses many disks, above the barrel.
///
/// @param r the radius in mm
/// @param z the longitudinal position in mm
/// @param cfg the event configuration supplying the two scales
/// @return the traversed material, in units of the material before the beam
///         line
inline float materialDepth(float r, float z, const EventConfig& cfg) {
  return 1.f + r / cfg.materialRScale + std::abs(z) / cfg.materialZScale;
}

/// A helix in a solenoid field along z, parametrised at its perigee with
/// respect to the beam axis.
///
/// The perigee parametrisation covers primaries and secondaries alike: a
/// secondary produced at a radius is just a helix with a large transverse
/// impact parameter that only exists beyond a minimum turning angle.
struct Helix {
  /// Azimuth of the momentum at the perigee
  float phi0{};
  /// Longitudinal position at the perigee in mm
  float z0{};
  /// Ratio of the longitudinal to the transverse momentum
  float cotTheta{};
  /// Radius of curvature in mm
  float radius{};
  /// Charge sign, +1 or -1
  float charge{};
  /// Distance of the centre of the circle from the beam axis in mm
  float rCentre{};
  /// Turning angle in rad below which the particle does not exist yet, i.e.
  /// the angle at which a secondary was produced
  float minGamma{};

  /// Transverse impact parameter with respect to the beam axis in mm.
  /// @return the signed impact parameter
  float d0() const { return charge * (radius - rCentre); }
};

/// Build a helix from perigee parameters.
/// @param d0 the signed transverse impact parameter in mm
/// @param phi0 the azimuth of the momentum at the perigee
/// @param z0 the longitudinal position at the perigee in mm
/// @param cotTheta the ratio of the longitudinal to the transverse momentum
/// @param pt the transverse momentum in GeV
/// @param charge the charge sign, +1 or -1
/// @param bFieldZ the solenoid field in Tesla
/// @return the helix
inline Helix makeHelix(float d0, float phi0, float z0, float cotTheta, float pt,
                       float charge, float bFieldZ) {
  Helix helix;
  helix.phi0 = phi0;
  helix.z0 = z0;
  helix.cotTheta = cotTheta;
  // radius of curvature in mm for a momentum in GeV and a field in Tesla
  helix.radius = pt * 1000.f / (0.3f * bFieldZ);
  helix.charge = charge;
  helix.rCentre = helix.radius - charge * d0;
  return helix;
}

/// The turning angle at which a helix reaches a radius.
/// @param helix the helix to intersect
/// @param r the radius in mm
/// @param gamma set to the turning angle in rad on success
/// @return whether the helix reaches this radius after its production point
inline bool helixGammaAtRadius(const Helix& helix, float r, float& gamma) {
  const float denominator = 2.f * helix.rCentre * helix.radius;
  if (denominator == 0.f) {
    return false;
  }
  const float cosGamma =
      (helix.rCentre * helix.rCentre + helix.radius * helix.radius - r * r) /
      denominator;
  if (cosGamma < -1.f || cosGamma > 1.f) {
    // the particle curls up before this radius, or never comes back to it
    return false;
  }
  gamma = std::acos(cosGamma);
  return gamma > helix.minGamma;
}

/// The azimuth of a helix at a radius.
/// @param helix the helix to intersect
/// @param r the radius in mm
/// @return the azimuth in rad
inline float helixPhiAtRadius(const Helix& helix, float r) {
  const float sinPhi =
      (r * r + helix.rCentre * helix.rCentre - helix.radius * helix.radius) /
      (2.f * r * helix.rCentre);
  return helix.phi0 - helix.charge * std::asin(std::clamp(sinPhi, -1.f, 1.f));
}

/// The radius a helix has reached after a turning angle.
/// @param helix the helix to propagate
/// @param gamma the turning angle in rad
/// @return the radius in mm
inline float helixRadiusAtGamma(const Helix& helix, float gamma) {
  const float squared = helix.rCentre * helix.rCentre +
                        helix.radius * helix.radius -
                        2.f * helix.rCentre * helix.radius * std::cos(gamma);
  return std::sqrt(std::max(0.f, squared));
}

/// Rebuild a helix from a point on it and the momentum direction there.
///
/// Used to launch a secondary from a hit of its parent: the production point
/// and the direction are known there, the perigee parameters are not.
///
/// @param r the radius of the production point in mm
/// @param phi the azimuth of the production point in rad
/// @param z the longitudinal position of the production point in mm
/// @param direction the azimuth of the momentum at the production point in rad
/// @param cotTheta the ratio of the longitudinal to the transverse momentum
/// @param pt the transverse momentum in GeV
/// @param charge the charge sign, +1 or -1
/// @param bFieldZ the solenoid field in Tesla
/// @return the helix, restricted to the part beyond the production point
inline Helix makeHelixFromPoint(float r, float phi, float z, float direction,
                                float cotTheta, float pt, float charge,
                                float bFieldZ) {
  constexpr float pi = std::numbers::pi_v<float>;

  Helix helix;
  helix.cotTheta = cotTheta;
  helix.radius = pt * 1000.f / (0.3f * bFieldZ);
  helix.charge = charge;

  // the centre of the circle sits one radius of curvature to the side of the
  // momentum, on the side the charge bends towards
  const float centreX =
      r * std::cos(phi) + charge * helix.radius * std::sin(direction);
  const float centreY =
      r * std::sin(phi) - charge * helix.radius * std::cos(direction);
  helix.rCentre = Acts::fastHypot(centreX, centreY);
  helix.phi0 = std::atan2(centreY, centreX) + charge * pi / 2.f;

  // walk back from the production point to the perigee, then forbid everything
  // before it
  float gamma = 0.f;
  if (!helixGammaAtRadius(helix, r, gamma)) {
    gamma = 0.f;
  }
  helix.z0 = z - cotTheta * helix.radius * gamma;
  helix.minGamma = gamma;
  return helix;
}

namespace detail {

/// One space point, in the cylindrical coordinates it is measured in, before
/// it is written into the container.
struct Hit {
  float r{};
  float phi{};
  float z{};
  std::uint32_t layer{};
  std::uint32_t particle{};
};

/// Intersection of a track with one surface.
struct Intersection {
  /// Radius in mm
  float r{};
  /// Azimuth in rad
  float phi{};
  /// Longitudinal position in mm
  float z{};
  /// Turning angle in rad, which gives the momentum direction there
  float gamma{};
  /// Index into `DetectorLayout::layers`, only set for a sensitive surface
  std::uint32_t layer{};
  /// Whether the crossed surface measures the particle
  bool sensitive{};
};

/// A track waiting to be propagated.
struct PendingTrack {
  /// The trajectory
  Helix helix;
  /// Transverse momentum in GeV
  float pt{};
  /// Whether it comes from the luminous region
  bool primary{};
  /// Radius at which it was produced in mm
  float productionRadius{};
  /// Longitudinal position at which it was produced in mm
  float productionZ{};
};

/// Find the module of a surface that an intersection falls into.
/// @param layout the detector layout
/// @param surface the intersected surface
/// @param r the radius of the intersection in mm
/// @param z the longitudinal position of the intersection in mm
/// @param layerIndex set to the index into `layout.layers` on success
/// @return whether the intersection is within the bounds of the surface
inline bool selectModule(const DetectorLayout& layout, const Surface& surface,
                         float r, float z, std::uint32_t& layerIndex) {
  const float selector = surface.type == LayerType::Barrel ? z : r;
  for (const std::uint32_t m : surface.modules) {
    const Layer& layer = layout.layers[m];
    if (selector >= layer.minBound && selector < layer.maxBound) {
      layerIndex = m;
      return true;
    }
  }
  return false;
}

/// Propagate a helix through the detector and collect its intersections.
/// @param layout the detector layout
/// @param helix the helix to propagate
/// @param intersections receives the intersections, in surface order
inline void propagate(const DetectorLayout& layout, const Helix& helix,
                      std::vector<Intersection>& intersections) {
  constexpr float pi = std::numbers::pi_v<float>;

  intersections.clear();
  for (const Surface& surface : layout.surfaces) {
    Intersection intersection;

    if (surface.type == LayerType::Barrel) {
      intersection.r = surface.refCoord;
      if (!helixGammaAtRadius(helix, intersection.r, intersection.gamma)) {
        continue;
      }
      intersection.z =
          helix.z0 + helix.cotTheta * helix.radius * intersection.gamma;
    } else {
      intersection.z = surface.refCoord;
      if (helix.cotTheta == 0.f) {
        continue;
      }
      intersection.gamma =
          (intersection.z - helix.z0) / (helix.cotTheta * helix.radius);
      if (intersection.gamma <= helix.minGamma || intersection.gamma >= pi) {
        // the disk is behind the particle, or it curls up before reaching it
        continue;
      }
      intersection.r = helixRadiusAtGamma(helix, intersection.gamma);
      if (intersection.r <= 0.f) {
        continue;
      }
    }

    intersection.sensitive = !surface.modules.empty();
    if (intersection.sensitive &&
        !selectModule(layout, surface, intersection.r, intersection.z,
                      intersection.layer)) {
      continue;
    }
    intersection.phi = helixPhiAtRadius(helix, intersection.r);
    intersections.push_back(intersection);
  }
}

}  // namespace detail

/// Generate a synthetic event.
///
/// The returned container carries the standard X/Y/Z/R/Phi columns, the
/// variances the triplet seeders need, plus three dynamic columns: "layerId"
/// with the index into `layout.layers`, "clusterWidth" and "localPositionY" as
/// expected by the GBTS seeder, and "particleId" with the index into
/// `Event::particles` of the generating particle, so that seeding efficiency
/// can be evaluated without a truth-matching stage.
///
/// @param layout the detector to generate hits on
/// @param cfg steering for the generated event
/// @return the generated event
inline Event generateEvent(const DetectorLayout& layout,
                           const EventConfig& cfg) {
  constexpr float pi = std::numbers::pi_v<float>;

  std::mt19937 rng(cfg.seed);
  std::uniform_real_distribution<float> phiDist(-pi, pi);
  std::uniform_real_distribution<float> etaDist(cfg.minEta, cfg.maxEta);
  std::uniform_real_distribution<float> uniform(0.f, 1.f);
  std::normal_distribution<float> z0Dist(0.f, cfg.beamspotSigmaZ);
  std::normal_distribution<float> d0Dist(0.f, cfg.d0Sigma);
  std::normal_distribution<float> openingDist(0.f, cfg.secondaryOpeningAngle);
  std::normal_distribution<float> smear(0.f, cfg.positionSmearing);

  const std::size_t numPrimaries = cfg.numPrimaries();

  // Tracks are processed in generation order; secondaries are appended while
  // their parent is walked and picked up by the same loop.
  std::vector<detail::PendingTrack> tracks;
  tracks.reserve(numPrimaries * 2);
  for (std::size_t t = 0; t < numPrimaries; ++t) {
    const float pt =
        samplePt(uniform(rng), cfg.minPt, cfg.ptScale, cfg.ptExponent);
    detail::PendingTrack track;
    track.helix = makeHelix(d0Dist(rng), phiDist(rng), z0Dist(rng),
                            std::sinh(etaDist(rng)), pt,
                            uniform(rng) < 0.5f ? -1.f : +1.f, cfg.bFieldZ);
    track.pt = pt;
    track.primary = true;
    tracks.push_back(track);
  }

  Event event;
  event.particles.reserve(tracks.size());

  std::vector<detail::Hit> hits;
  hits.reserve(numPrimaries * 12);
  std::vector<detail::Intersection> intersections;

  for (std::size_t t = 0; t < tracks.size(); ++t) {
    // the vector grows while secondaries are queued, so work on a copy
    const detail::PendingTrack track = tracks[t];
    const std::uint32_t particle = static_cast<std::uint32_t>(t);

    detail::propagate(layout, track.helix, intersections);

    std::uint32_t numHits = 0;
    for (const detail::Intersection& intersection : intersections) {
      if (intersection.sensitive) {
        // smear along the measured directions of the sensor, leaving the
        // normal to it at the nominal surface position
        const bool barrel =
            layout.layers[intersection.layer].type == LayerType::Barrel;
        hits.push_back(
            detail::Hit{intersection.r + (barrel ? 0.f : smear(rng)),
                        intersection.phi + smear(rng) / intersection.r,
                        intersection.z + (barrel ? smear(rng) : 0.f),
                        intersection.layer, particle});
        ++numHits;
      }

      if (!track.primary) {
        // only one generation of secondaries, so that the yield stays the
        // configured one
        continue;
      }

      const float yield = cfg.secondaryRate *
                          materialDepth(intersection.r, intersection.z, cfg);
      std::poisson_distribution<int> secondaryCount(yield);
      const int numSecondaries = secondaryCount(rng);
      for (int s = 0; s < numSecondaries; ++s) {
        // secondaries are boosted along their parent, so they start out close
        // to its direction
        const float parentDirection =
            track.helix.phi0 - track.helix.charge * intersection.gamma;
        const float pt = cfg.secondaryMinPt -
                         cfg.secondaryPtSlope *
                             std::log(std::max(1e-6f, 1.f - uniform(rng)));
        detail::PendingTrack secondary;
        secondary.helix = makeHelixFromPoint(
            intersection.r, intersection.phi, intersection.z,
            parentDirection + openingDist(rng),
            std::sinh(std::asinh(track.helix.cotTheta) + openingDist(rng)), pt,
            uniform(rng) < 0.5f ? -1.f : +1.f, cfg.bFieldZ);
        secondary.pt = pt;
        secondary.primary = false;
        secondary.productionRadius = intersection.r;
        secondary.productionZ = intersection.z;
        tracks.push_back(secondary);
      }
    }

    GeneratedParticle info;
    info.pt = track.pt;
    info.eta = std::asinh(track.helix.cotTheta);
    info.d0 = track.helix.d0();
    info.z0 = track.helix.z0;
    info.productionRadius = track.productionRadius;
    info.productionZ = track.productionZ;
    info.primary = track.primary;
    info.numHits = numHits;
    event.particles.push_back(info);
  }

  Acts::SpacePointContainer spacePoints(
      Acts::SpacePointColumns::CopiedFromIndex | Acts::SpacePointColumns::X |
      Acts::SpacePointColumns::Y | Acts::SpacePointColumns::Z |
      Acts::SpacePointColumns::R | Acts::SpacePointColumns::Phi |
      Acts::SpacePointColumns::VarianceZ | Acts::SpacePointColumns::VarianceR);
  auto layerColumn = spacePoints.createColumn<std::uint32_t>("layerId");
  auto clusterWidthColumn = spacePoints.createColumn<float>("clusterWidth");
  auto localPositionColumn = spacePoints.createColumn<float>("localPositionY");
  auto particleColumn = spacePoints.createColumn<std::uint32_t>("particleId");
  spacePoints.reserve(hits.size());

  const float measuredVariance = cfg.positionSmearing * cfg.positionSmearing;
  // a uniform distribution across the sensor thickness is all the unmeasured
  // direction is known to
  const float normalVariance = cfg.sensorThickness * cfg.sensorThickness / 12.f;
  for (const detail::Hit& hit : hits) {
    const bool barrel = layout.layers[hit.layer].type == LayerType::Barrel;
    // the helix azimuth is not wrapped while it is propagated
    const float phi = std::remainder(hit.phi, 2.f * pi);
    auto sp = spacePoints.createSpacePoint();
    sp.x() = hit.r * std::cos(phi);
    sp.y() = hit.r * std::sin(phi);
    sp.z() = hit.z;
    sp.r() = hit.r;
    sp.phi() = phi;
    sp.varianceZ() = barrel ? measuredVariance : normalVariance;
    sp.varianceR() = barrel ? normalVariance : measuredVariance;
    sp.copiedFromIndex() = sp.index();
    sp.extra(layerColumn) = hit.layer;
    sp.extra(clusterWidthColumn) = 0.f;
    sp.extra(localPositionColumn) = 0.f;
    sp.extra(particleColumn) = hit.particle;
  }

  event.spacePoints = std::move(spacePoints);
  return event;
}

}  // namespace ActsTests::SyntheticItk

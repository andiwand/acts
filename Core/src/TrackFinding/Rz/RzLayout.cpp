// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/Rz/RzLayout.hpp"

#include "Acts/Geometry/ApproachDescriptor.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <set>

namespace Acts::Experimental {

namespace {

/// The RZ geometry of a cylinder or disc surface, or nothing for any other
/// shape.
std::optional<RzSurface> describe(const Surface& surface,
                                  const GeometryContext& gctx) {
  RzSurface out;
  out.geometryId = surface.geometryId();
  out.surface = surface.getSharedPtr();
  const Vector3 center = surface.center(gctx);
  if (surface.type() == Surface::Cylinder) {
    const auto& bounds = static_cast<const CylinderBounds&>(surface.bounds());
    out.shape = RzShape::Cylinder;
    out.refCoord = bounds.get(CylinderBounds::eR);
    const double hz = bounds.get(CylinderBounds::eHalfLengthZ);
    out.minBound = center.z() - hz;
    out.maxBound = center.z() + hz;
    return out;
  }
  if (surface.type() == Surface::Disc) {
    const auto* bounds = dynamic_cast<const DiscBounds*>(&surface.bounds());
    if (bounds == nullptr) {
      return std::nullopt;
    }
    out.shape = RzShape::Disc;
    out.refCoord = center.z();
    out.minBound = bounds->rMin();
    out.maxBound = bounds->rMax();
    return out;
  }
  return std::nullopt;
}

/// Read the material off a surface into bands along its extended coordinate,
/// averaged over azimuth, neighbouring samples merged where they agree.
void sampleMaterial(const Surface& surface, RzSurface& out,
                    const RzLayoutOptions& options) {
  const ISurfaceMaterial* material = surface.surfaceMaterial();
  if (material == nullptr) {
    return;
  }
  const std::uint32_t n = std::max(1u, options.materialSamples);
  const std::uint32_t nPhi = std::max(1u, options.phiSamples);
  const double width = (out.maxBound - out.minBound) / n;

  // A binned lookup takes its two values in the order of its own axes, not
  // in the surface's local frame, so each axis is served what it is binned in
  const std::vector<AxisDirection> axes = material->localAxisDirections();
  auto valueFor = [&](AxisDirection axis, double r, double phi, double z) {
    switch (axis) {
      case AxisDirection::AxisR:
        return r;
      case AxisDirection::AxisPhi:
        return phi;
      case AxisDirection::AxisRPhi:
        return r * phi;
      case AxisDirection::AxisZ:
        return z;
      default:
        return 0.;
    }
  };
  auto localFor = [&](double r, double phi, double z) {
    if (axes.size() >= 2) {
      return Vector2(valueFor(axes[0], r, phi, z),
                     valueFor(axes[1], r, phi, z));
    }
    if (axes.size() == 1) {
      return Vector2(valueFor(axes[0], r, phi, z), 0.);
    }
    // homogeneous: any point will do
    return out.shape == RzShape::Cylinder ? Vector2(r * phi, z)
                                          : Vector2(r, phi);
  };

  std::vector<MaterialSlab> samples;
  samples.reserve(n);
  std::vector<MaterialSlab> ring;
  for (std::uint32_t i = 0; i < n; ++i) {
    const double along = out.minBound + (i + 0.5) * width;
    ring.clear();
    for (std::uint32_t j = 0; j < nPhi; ++j) {
      const double phi =
          -std::numbers::pi + (j + 0.5) * 2. * std::numbers::pi / nPhi;
      const double r = out.shape == RzShape::Cylinder ? out.refCoord : along;
      const double z = out.shape == RzShape::Cylinder ? along : out.refCoord;
      ring.push_back(material->materialSlab(localFor(r, phi, z)));
    }
    MaterialSlab mean = MaterialSlab::combineLayers(ring);
    mean.scaleThickness(1.f / static_cast<float>(nPhi));
    samples.push_back(mean);
  }

  bool anyMaterial = false;
  for (const MaterialSlab& s : samples) {
    anyMaterial |= !s.isVacuum();
  }
  if (!anyMaterial) {
    return;
  }

  // merge runs of samples whose x/X0 agree with the run's first one
  std::uint32_t start = 0;
  while (start < n) {
    std::uint32_t end = start + 1;
    const double ref = samples[start].thicknessInX0();
    while (end < n) {
      const double x = samples[end].thicknessInX0();
      const double scale = std::max(ref, x);
      const bool same =
          scale <= 0. ||
          std::abs(x - ref) <= options.materialBandTolerance * scale;
      if (!same) {
        break;
      }
      ++end;
    }
    std::vector<MaterialSlab> members(samples.begin() + start,
                                      samples.begin() + end);
    MaterialSlab band = MaterialSlab::combineLayers(members);
    band.scaleThickness(1.f / static_cast<float>(members.size()));
    if (out.materialEdges.empty()) {
      out.materialEdges.push_back(out.minBound + start * width);
    }
    out.materialEdges.push_back(out.minBound + end * width);
    out.materialBands.push_back(band);
    start = end;
  }
  // the last edge is the surface end exactly, whatever rounding did above
  out.materialEdges.back() = out.maxBound;
}

/// A sensitive plane as the finder sees it, or nothing for a shape it cannot
/// project onto.
std::optional<RzModule> describeModule(const Surface& surface,
                                       const GeometryContext& gctx) {
  if (surface.type() != Surface::Plane) {
    return std::nullopt;
  }
  RzModule m;
  const Transform3& transform = surface.localToGlobalTransform(gctx);
  m.center = transform.translation();
  m.u = transform.rotation().col(0);
  m.v = transform.rotation().col(1);
  m.normal = transform.rotation().col(2);
  m.geometryId = surface.geometryId();
  m.surface = surface.getSharedPtr();
  if (const auto* planar = dynamic_cast<const PlanarBounds*>(&surface.bounds());
      planar != nullptr) {
    for (const Vector2& vertex : planar->vertices(1)) {
      m.halfU = std::max(m.halfU, std::abs(vertex.x()));
      m.halfV = std::max(m.halfV, std::abs(vertex.y()));
    }
  }
  return m;
}

/// Material a geometry puts on the sensitive modules themselves rather than
/// on the layer, taken as the mean over the modules and stacked onto every
/// band of the layer's surface. The ODD carries none; the generic detector
/// does.
void foldModuleMaterial(const SurfaceArray& array, RzSurface& out) {
  std::vector<MaterialSlab> slabs;
  for (const Surface* surface : array.surfaces()) {
    if (const ISurfaceMaterial* material = surface->surfaceMaterial();
        material != nullptr) {
      slabs.push_back(material->materialSlab(Vector2(0., 0.)));
    }
  }
  if (slabs.empty()) {
    return;
  }
  MaterialSlab mean = MaterialSlab::combineLayers(slabs);
  mean.scaleThickness(static_cast<float>(1. / array.surfaces().size()));
  if (mean.isVacuum()) {
    return;
  }
  if (out.materialBands.empty()) {
    out.materialEdges = {out.minBound, out.maxBound};
    out.materialBands = {mean};
    return;
  }
  for (MaterialSlab& band : out.materialBands) {
    band = MaterialSlab::combineLayers(band, mean);
  }
}

}  // namespace

/// Fill a band's table for a particle
RzMaterialTable tabulate(const MaterialSlab& slab,
                         const ParticleHypothesis& hyp) {
  RzMaterialTable t;
  const float mass = hyp.mass();
  const float absQ = hyp.absoluteCharge();
  const PdgParticle pdg = hyp.absolutePdg();
  t.logThicknessInX0 = std::log(std::max(slab.thicknessInX0(), 1e-12f));
  for (std::uint32_t i = 0; i < RzMaterialTable::kBins; ++i) {
    const double p =
        std::exp(RzMaterialTable::logMinP() + i * RzMaterialTable::logStep());
    const float qOverP = static_cast<float>(hyp.qOverP(p, absQ));
    const float theta0 =
        computeMultipleScatteringTheta0(slab, pdg, mass, qOverP, absQ);
    t.theta0Sq[i] = theta0 * theta0;
    t.energyLoss[i] = computeEnergyLossMean(slab, pdg, mass, qOverP, absQ);
    const float sigma =
        computeEnergyLossLandauSigmaQOverP(slab, mass, qOverP, absQ);
    t.sigmaQOverPSq[i] = sigma * sigma;
  }
  return t;
}

int RzSurface::materialBandAt(double along) const {
  if (materialBands.empty() || along < materialEdges.front() ||
      along >= materialEdges.back()) {
    return -1;
  }
  const auto edge = std::ranges::upper_bound(materialEdges, along);
  return static_cast<int>(edge - materialEdges.begin()) - 1;
}

const MaterialSlab* RzSurface::materialAt(double along) const {
  if (materialBands.empty() || along < materialEdges.front() ||
      along >= materialEdges.back()) {
    return nullptr;
  }
  const auto edge = std::ranges::upper_bound(materialEdges, along);
  return &materialBands[static_cast<std::size_t>(edge - materialEdges.begin()) -
                        1];
}

double RzLayer::phiBinWidth() const {
  return 2. * std::numbers::pi / phiBins;
}

std::uint32_t RzLayout::bin(std::uint32_t layerIndex, double phi,
                            double along) const {
  const RzLayer& layer = layers[layerIndex];
  const int nPhi = static_cast<int>(layer.phiBins);
  int p = static_cast<int>(std::floor(phi / layer.phiBinWidth()));
  p = ((p % nPhi) + nPhi) % nPhi;
  int a = static_cast<int>(
      std::floor((along - layer.alongMin) / layer.alongBinWidth()));
  a = std::clamp(a, 0, static_cast<int>(layer.alongBins) - 1);
  return layer.binOffset + static_cast<std::uint32_t>(p) * layer.alongBins +
         static_cast<std::uint32_t>(a);
}

RzLayout makeRzLayout(const TrackingGeometry& trackingGeometry,
                      const GeometryContext& gctx,
                      const RzLayoutOptions& options, const Logger& logger) {
  RzLayout layout;
  std::set<const Surface*> seen;

  auto addPassive = [&](const Surface& surface) {
    if (surface.surfaceMaterial() == nullptr || !seen.insert(&surface).second) {
      return;
    }
    std::optional<RzSurface> rz = describe(surface, gctx);
    if (!rz.has_value()) {
      ACTS_DEBUG("Skipping material surface " << surface.geometryId()
                                              << " of unsupported shape");
      return;
    }
    sampleMaterial(surface, *rz, options);
    if (rz->materialBands.empty()) {
      return;
    }
    layout.surfaces.push_back(std::move(*rz));
  };

  trackingGeometry.visitVolumes([&](const TrackingVolume* volume) {
    if (const auto* bounds =
            dynamic_cast<const CylinderVolumeBounds*>(&volume->volumeBounds());
        bounds != nullptr) {
      layout.escapeRadius = std::max(layout.escapeRadius,
                                     bounds->get(CylinderVolumeBounds::eMaxR));
      layout.escapeHalfZ = std::max(
          layout.escapeHalfZ, bounds->get(CylinderVolumeBounds::eHalfLengthZ));
    }
    for (const auto& boundary : volume->boundarySurfaces()) {
      addPassive(boundary->surfaceRepresentation());
    }
    if (volume->confinedLayers() == nullptr) {
      return;
    }
    for (const auto& layer : volume->confinedLayers()->arrayObjects()) {
      const Surface& representing = layer->surfaceRepresentation();
      if (const ApproachDescriptor* approach = layer->approachDescriptor();
          approach != nullptr) {
        for (const Surface* surface : approach->containedSurfaces()) {
          addPassive(*surface);
        }
      }
      const SurfaceArray* array = layer->surfaceArray();
      if (array == nullptr) {
        addPassive(representing);
        continue;
      }
      std::optional<RzSurface> rz = describe(representing, gctx);
      if (!rz.has_value()) {
        ACTS_WARNING("Skipping sensitive layer " << representing.geometryId()
                                                 << " of unsupported shape");
        continue;
      }
      std::vector<RzModule> modules;
      for (const Surface* surface : array->surfaces()) {
        if (options.surfaceSelector && !options.surfaceSelector(*surface)) {
          continue;
        }
        std::optional<RzModule> module = describeModule(*surface, gctx);
        if (!module.has_value()) {
          ACTS_WARNING("Skipping module " << surface->geometryId()
                                          << " of unsupported shape");
          continue;
        }
        modules.push_back(*module);
      }
      if (modules.empty()) {
        addPassive(representing);
        continue;
      }
      seen.insert(&representing);
      sampleMaterial(representing, *rz, options);
      foldModuleMaterial(*array, *rz);

      const std::uint32_t layerIndex =
          static_cast<std::uint32_t>(layout.layers.size());
      rz->layer = layerIndex;
      RzLayer rzLayer;
      rzLayer.surface = static_cast<std::uint32_t>(layout.surfaces.size());
      rzLayer.phiBins = std::max(1u, options.phiBins);
      rzLayer.alongMin = rz->minBound;
      rzLayer.alongMax = rz->maxBound;
      rzLayer.alongBins = std::max(
          1u, static_cast<std::uint32_t>(std::ceil(
                  (rz->maxBound - rz->minBound) / options.alongBinWidth)));
      for (RzModule& module : modules) {
        module.layer = layerIndex;
        rzLayer.maxHalfV = std::max(rzLayer.maxHalfV, module.halfV);
        const double offset =
            rz->shape == RzShape::Cylinder
                ? std::hypot(module.center.x(), module.center.y()) -
                      rz->refCoord
                : module.center.z() - rz->refCoord;
        rzLayer.halfThickness =
            std::max(rzLayer.halfThickness, std::abs(offset));
        rzLayer.maxHalfExtent = std::max(
            rzLayer.maxHalfExtent, std::hypot(module.halfU, module.halfV));
        layout.moduleIndex.emplace(
            module.geometryId,
            static_cast<std::uint32_t>(layout.modules.size()));
        layout.modules.push_back(module);
      }
      layout.layers.push_back(rzLayer);
      layout.surfaces.push_back(std::move(*rz));
    }
  });

  // order the navigation lists; the surface indices themselves stay as built
  // so that the layers keep pointing at them
  for (std::uint32_t i = 0; i < layout.surfaces.size(); ++i) {
    (layout.surfaces[i].shape == RzShape::Cylinder ? layout.cylinders
                                                   : layout.discs)
        .push_back(i);
  }
  auto byRef = [&](std::uint32_t a, std::uint32_t b) {
    return layout.surfaces[a].refCoord < layout.surfaces[b].refCoord;
  };
  std::ranges::sort(layout.cylinders, byRef);
  std::ranges::sort(layout.discs, byRef);

  if (options.materialTables) {
    for (RzSurface& surface : layout.surfaces) {
      for (const MaterialSlab& band : surface.materialBands) {
        surface.materialTables.push_back(
            tabulate(band, options.particleHypothesis));
      }
    }
  }

  std::uint32_t bins = 0;
  for (RzLayer& layer : layout.layers) {
    layer.binOffset = bins;
    bins += layer.phiBins * layer.alongBins;
  }
  layout.totalBins = bins;

  // the modules by centre, in the same bins as the measurements will be
  std::vector<std::uint32_t> binOf(layout.modules.size());
  layout.moduleBinStart.assign(bins + 1, 0);
  for (std::uint32_t i = 0; i < layout.modules.size(); ++i) {
    const RzModule& m = layout.modules[i];
    const RzSurface& surface = layout.surfaces[layout.layers[m.layer].surface];
    const double phi = std::atan2(m.center.y(), m.center.x());
    const double along = surface.shape == RzShape::Cylinder
                             ? m.center.z()
                             : std::hypot(m.center.x(), m.center.y());
    binOf[i] = layout.bin(m.layer, phi, along);
    ++layout.moduleBinStart[binOf[i] + 1];
  }
  for (std::size_t i = 1; i < layout.moduleBinStart.size(); ++i) {
    layout.moduleBinStart[i] += layout.moduleBinStart[i - 1];
  }
  layout.moduleOrder.assign(layout.modules.size(), 0);
  std::vector<std::uint32_t> fill(layout.moduleBinStart.begin(),
                                  layout.moduleBinStart.end() - 1);
  for (std::uint32_t i = 0; i < binOf.size(); ++i) {
    layout.moduleOrder[fill[binOf[i]]++] = i;
  }

  ACTS_INFO("RZ layout: " << layout.cylinders.size() << " cylinders, "
                          << layout.discs.size() << " discs, "
                          << layout.layers.size() << " sensitive layers, "
                          << layout.modules.size() << " modules, escape r "
                          << layout.escapeRadius << " |z| "
                          << layout.escapeHalfZ);
  for (const RzSurface& s : layout.surfaces) {
    ACTS_DEBUG((s.shape == RzShape::Cylinder ? "cylinder r=" : "disc z=")
               << s.refCoord << " [" << s.minBound << ", " << s.maxBound
               << "] bands " << s.materialBands.size()
               << (s.layer != kRzNone ? " sensitive" : "") << " from "
               << s.geometryId);
  }
  return layout;
}

}  // namespace Acts::Experimental

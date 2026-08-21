// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Layer.hpp"

#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <algorithm>
#include <vector>

namespace Acts {

Layer::Layer(std::unique_ptr<SurfaceArray> surfaceArray, double thickness,
             std::unique_ptr<ApproachDescriptor> ades, LayerType laytyp)
    : m_nextLayers(NextLayers(nullptr, nullptr)),
      m_surfaceArray(surfaceArray.release()),
      m_layerThickness(thickness),
      m_approachDescriptor(nullptr),
      m_representingVolume(nullptr),
      m_layerType(laytyp),
      m_ssRepresentingSurface(1) {
  if (ades) {
    ades->registerLayer(*this);
    m_approachDescriptor = std::move(ades);
    m_ssApproachSurfaces = 1;  // indicates existence
  }
  // indicates existence of sensitive surfaces
  if (m_surfaceArray) {
    m_ssSensitiveSurfaces = 1;
  }
}

Layer::~Layer() noexcept = default;

const ApproachDescriptor* Layer::approachDescriptor() const {
  return m_approachDescriptor.get();
}

ApproachDescriptor* Layer::approachDescriptor() {
  return const_cast<ApproachDescriptor*>(m_approachDescriptor.get());
}

void Layer::closeGeometry(const IMaterialDecorator* materialDecorator,
                          const GeometryIdentifier& layerID,
                          const GeometryIdentifierHook& hook,
                          const Logger& logger) {
  // set the volumeID of this
  assignGeometryId(layerID);
  // assign to the representing surface
  Surface* rSurface = const_cast<Surface*>(&surfaceRepresentation());
  if (materialDecorator != nullptr) {
    materialDecorator->decorate(*rSurface);
  }
  ACTS_DEBUG("layerID: " << layerID);

  rSurface->assignGeometryId(layerID);

  // also find out how the sub structure is defined
  if (surfaceRepresentation().hasMaterial()) {
    m_ssRepresentingSurface = 2;
  }
  // loop over the approach surfaces
  if (m_approachDescriptor) {
    // indicates the existence of approach surfaces
    m_ssApproachSurfaces = 1;
    // loop through the approachSurfaces and assign unique GeomeryID
    GeometryIdentifier::Value iasurface = 0;
    for (auto& aSurface : m_approachDescriptor->containedSurfaces()) {
      auto asurfaceID = GeometryIdentifier(layerID).withApproach(++iasurface);
      auto mutableASurface = const_cast<Surface*>(aSurface);
      mutableASurface->assignGeometryId(asurfaceID);
      if (materialDecorator != nullptr) {
        materialDecorator->decorate(*mutableASurface);
      }
      // if any of the approach surfaces has material
      if (aSurface->hasMaterial()) {
        m_ssApproachSurfaces = 2;
      }
    }
  }
  // check if you have sensitive surfaces
  if (m_surfaceArray) {
    // indicates the existence of sensitive surfaces
    m_ssSensitiveSurfaces = 1;
    // loop sensitive surfaces and assign unique GeometryIdentifier
    GeometryIdentifier::Value issurface = 0;
    for (auto& sSurface : m_surfaceArray->surfaces()) {
      auto ssurfaceID = GeometryIdentifier(layerID).withSensitive(++issurface);
      ssurfaceID = hook.decorateIdentifier(ssurfaceID, *sSurface);
      auto mutableSSurface = const_cast<Surface*>(sSurface);
      mutableSSurface->assignGeometryId(ssurfaceID);
      if (materialDecorator != nullptr) {
        materialDecorator->decorate(*mutableSSurface);
      }
      // if any of the sensitive surfaces has material
      if (sSurface->hasMaterial()) {
        m_ssSensitiveSurfaces = 2;
      }
    }
  }
}

boost::container::small_vector<NavigationTarget, 10> Layer::compatibleSurfaces(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const NavigationOptions<Surface>& options) const {
  // the list of valid intersection
  boost::container::small_vector<NavigationTarget, 10> surfaceIntersections;

  // fast exit - there is nothing to
  if (!m_surfaceArray || !m_approachDescriptor) {
    return surfaceIntersections;
  }

  double nearLimit = options.nearLimit;
  double farLimit = options.farLimit;

  auto isUnique = [&](const NavigationTarget& b) {
    return std::ranges::none_of(
        surfaceIntersections, [&b](const NavigationTarget& a) {
          return &a.surface() == &b.surface() &&
                 a.intersectionIndex() == b.intersectionIndex();
        });
  };

  // lemma 0 : accept the surface
  auto acceptSurface = [&options](const Surface& sf,
                                  bool sensitive = false) -> bool {
    // surface is sensitive and you're asked to resolve
    if (sensitive && options.resolveSensitive) {
      return true;
    }
    // next option: it's a material surface, and you want to have it
    if (options.resolveMaterial && sf.hasMaterial()) {
      return true;
    }
    // last option: resolve all
    return options.resolvePassive;
  };

  // lemma 1 : check and fill the surface
  auto processSurface = [&](const Surface& surface, bool sensitive = false) {
    // veto if it's start surface
    if (options.startObject == &surface) {
      return;
    }
    // veto if it doesn't fit the prescription
    if (!acceptSurface(surface, sensitive)) {
      return;
    }
    BoundaryTolerance boundaryTolerance = options.boundaryTolerance;
    if (rangeContainsValue(options.externalSurfaces, surface.geometryId())) {
      boundaryTolerance = BoundaryTolerance::Infinite();
    }
    // the surface intersection
    auto [intersection, intersectionIndex] =
        surface.intersect(gctx, position, direction, boundaryTolerance)
            .closestWithIndex();
    NavigationTarget surfaceIntersection(intersection, intersectionIndex,
                                         surface, boundaryTolerance);
    if (intersection.isValid() &&
        detail::checkPathLength(intersection.pathLength(), nearLimit,
                                farLimit) &&
        isUnique(surfaceIntersection)) {
      surfaceIntersections.push_back(surfaceIntersection);
    }
  };

  // (A) approach descriptor section
  //
  // the approach surfaces are in principle always testSurfaces
  // - the surface on approach is excluded via the veto
  // - the surfaces are only collected if needed
  // An approach surface can only be accepted if it carries material and
  // material is being resolved (m_ssApproachSurfaces > 1 flags that at geometry
  // closure), or if everything is being resolved. Otherwise the loop below
  // intersects nothing and only burns the veto.
  if (m_approachDescriptor &&
      ((options.resolveMaterial && m_ssApproachSurfaces > 1) ||
       options.resolvePassive)) {
    // the approach surfaces
    const std::vector<const Surface*>& approachSurfaces =
        m_approachDescriptor->containedSurfaces();
    // we loop through and veto
    // - if the approach surface is the parameter surface
    // - if the surface is not compatible with the collect
    for (const Surface* aSurface : approachSurfaces) {
      processSurface(*aSurface);
    }
  }

  // (B) sensitive surface section
  //
  // check the sensitive surfaces if you have some
  if (m_surfaceArray && (options.resolveMaterial || options.resolvePassive ||
                         options.resolveSensitive)) {
    // get the candidates. The pack of a bin covers everything reachable from
    // anywhere in that bin, so most of it is nowhere near this trajectory --
    // reject those against their precomputed extent before paying for a 3D
    // intersection.
    const SurfaceArray::NeighborQuery query =
        m_surfaceArray->neighborQuery(gctx, position, direction);
    // loop through and veto
    // - if the approach surface is the parameter surface
    // - if the surface is not compatible with the type(s) that are collected
    for (std::size_t i = 0; i < query.surfaces.size(); ++i) {
      if (!query.mayCross(i)) {
        continue;
      }
      processSurface(*query.surfaces[i], true);
    }
  }

  // (C) representing surface section
  //
  // the layer surface itself is a testSurface
  const Surface* layerSurface = &surfaceRepresentation();
  processSurface(*layerSurface);

  return surfaceIntersections;
}

NavigationTarget Layer::surfaceOnApproach(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const NavigationOptions<Layer>& options) const {
  // resolve directive based by options
  // - options.resolvePassive is on -> always
  // - options.resolveSensitive is on -> always
  // - options.resolveMaterial is on
  //   && either sensitive or approach surfaces have material
  const bool resolvePS = options.resolveSensitive || options.resolvePassive;
  const bool resolveMS =
      options.resolveMaterial &&
      (m_ssSensitiveSurfaces > 1 || m_ssApproachSurfaces > 1 ||
       surfaceRepresentation().hasMaterial());

  // The Limits
  const double nearLimit = options.nearLimit;
  const double farLimit = options.farLimit;

  // Approach descriptor present and resolving is necessary
  if (m_approachDescriptor && (resolvePS || resolveMS)) {
    NavigationTarget aSurface = m_approachDescriptor->approachSurface(
        gctx, position, direction, options.boundaryTolerance, nearLimit,
        farLimit);
    return aSurface;
  }

  // Intersect and check the representing surface
  const Surface& layerSurface = surfaceRepresentation();
  const MultiIntersection3D multiIntersection = layerSurface.intersect(
      gctx, position, direction, options.boundaryTolerance);
  for (auto [intersectionIndex, intersection] :
       Acts::enumerate(multiIntersection)) {
    if (intersection.isValid() &&
        detail::checkPathLength(intersection.pathLength(), nearLimit,
                                farLimit)) {
      return NavigationTarget(intersection, intersectionIndex, *this,
                              layerSurface, options.boundaryTolerance);
    }
  }
  return NavigationTarget::None();
}

}  // namespace Acts

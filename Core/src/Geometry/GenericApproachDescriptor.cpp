// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GenericApproachDescriptor.hpp"

#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <boost/container/small_vector.hpp>

namespace Acts {

void GenericApproachDescriptor::registerLayer(const Layer& lay) {
  // go through the surfaces
  for (const Surface* surface : m_surfaceCache) {
    auto* mutableSurface = const_cast<Surface*>(surface);
    mutableSurface->associateLayer(lay);
  }
}

NavigationTarget GenericApproachDescriptor::approachSurface(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryTolerance& boundaryTolerance,
    double nearLimit, double farLimit) const {
  // Almost always 2 surfaces with 2 solutions each, and only the closest one
  // is ever returned, so keep a running minimum instead of materialising every
  // candidate into a container.
  NavigationTarget best = NavigationTarget::None();
  for (const Surface* surface : m_surfaceCache) {
    auto multiIntersection =
        surface->intersect(gctx, position, direction, boundaryTolerance);
    for (auto [intersectionIndex, intersection] :
         Acts::enumerate(multiIntersection)) {
      if (!intersection.isValid() ||
          !detail::checkPathLength(intersection.pathLength(), nearLimit,
                                   farLimit)) {
        continue;
      }
      NavigationTarget candidate(intersection, intersectionIndex,
                                 *surface->associatedLayer(), *surface,
                                 boundaryTolerance);
      if (best.isNone() || NavigationTarget::pathLengthOrder(candidate, best)) {
        best = candidate;
      }
    }
  }
  return best;
}

const std::vector<const Surface*>&
GenericApproachDescriptor::containedSurfaces() const {
  return m_surfaceCache;
}

std::vector<const Surface*>& GenericApproachDescriptor::containedSurfaces() {
  return m_surfaceCache;
}

}  // namespace Acts

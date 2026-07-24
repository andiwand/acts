// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/OverlapChecks.hpp"

namespace Acts::Experimental {

bool overlaps(const GeometryContext& /*gctx*/, const Surface& /*surfaceA*/,
              const Surface& /*surfaceB*/) {
  // TODO not implemented yet
  return false;
}

bool overlaps(const GeometryContext& /*gctx*/, const Volume& /*volumeA*/,
              const Volume& /*volumeB*/) {
  // TODO not implemented yet
  return false;
}

bool overlaps(const GeometryContext& /*gctx*/, const Volume& /*volume*/,
              const Surface& /*surface*/) {
  // TODO not implemented yet
  return false;
}

bool containsFully(const GeometryContext& /*gctx*/, const Volume& /*volume*/,
                   const Surface& /*surface*/) {
  // TODO not implemented yet
  return false;
}

bool containsFully(const GeometryContext& /*gctx*/, const Volume& /*volume*/,
                   const Volume& /*containedVolume*/) {
  // TODO not implemented yet
  return false;
}

}  // namespace Acts::Experimental

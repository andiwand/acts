// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

class GeometryContext;
class Surface;
class Volume;

namespace Experimental {

/// Check whether two surfaces overlap, i.e. whether the two bounded surface
/// areas have a non-empty intersection.
///
/// @note Not implemented yet, always returns false
///
/// @param gctx The geometry context to resolve the surface placements
/// @param surfaceA The first surface
/// @param surfaceB The second surface
/// @return True if the two surfaces overlap
bool overlaps(const GeometryContext& gctx, const Surface& surfaceA,
              const Surface& surfaceB);

/// Check whether two volumes overlap, i.e. whether the intersection of the
/// two bounded volumes is non-empty.
///
/// @note Not implemented yet, always returns false
///
/// @param gctx The geometry context to resolve the volume placements
/// @param volumeA The first volume
/// @param volumeB The second volume
/// @return True if the two volumes overlap
bool overlaps(const GeometryContext& gctx, const Volume& volumeA,
              const Volume& volumeB);

/// Check whether a volume and a surface overlap, i.e. whether any part of
/// the bounded surface area lies inside the bounded volume.
///
/// @note Not implemented yet, always returns false
///
/// @param gctx The geometry context to resolve the placements
/// @param volume The volume
/// @param surface The surface
/// @return True if the surface overlaps with the volume
bool overlaps(const GeometryContext& gctx, const Volume& volume,
              const Surface& surface);

/// Check whether a surface is fully contained inside a volume, i.e. whether
/// the entire bounded surface area lies inside the bounded volume.
///
/// @note Not implemented yet, always returns false
///
/// @param gctx The geometry context to resolve the placements
/// @param volume The volume acting as the container
/// @param surface The surface to be tested for containment
/// @return True if the surface is fully contained in the volume
bool containsFully(const GeometryContext& gctx, const Volume& volume,
                   const Surface& surface);

/// Check whether a volume is fully contained inside another volume, i.e.
/// whether the entire bounded inner volume lies inside the bounded outer
/// volume.
///
/// @note Not implemented yet, always returns false
///
/// @param gctx The geometry context to resolve the placements
/// @param volume The volume acting as the container
/// @param containedVolume The volume to be tested for containment
/// @return True if @p containedVolume is fully contained in @p volume
bool containsFully(const GeometryContext& gctx, const Volume& volume,
                   const Volume& containedVolume);

}  // namespace Experimental

}  // namespace Acts

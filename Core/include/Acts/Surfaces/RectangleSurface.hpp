// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"

namespace Acts {

class RectangleSurface : public PlaneSurface {
  friend class Surface;

 protected:
  RectangleSurface(const RectangleSurface& other);

  RectangleSurface(const GeometryContext& gctx, const RectangleSurface& other,
                   const Transform3& transform);

  RectangleSurface(std::shared_ptr<const RectangleBounds> pbounds,
                   const DetectorElementBase& detelement);

  explicit RectangleSurface(
      const Transform3& transform,
      std::shared_ptr<const RectangleBounds> pbounds = nullptr);

 public:
  RectangleSurface& operator=(const RectangleSurface& other);

  using PlaneSurface::globalToLocal;
  using PlaneSurface::localToGlobal;
  using PlaneSurface::normal;

  const SurfaceBounds& bounds() const override;
  const std::shared_ptr<const RectangleBounds>& boundsPtr() const;
  void assignSurfaceBounds(std::shared_ptr<const RectangleBounds> newBounds);

 protected:
  RectangleBounds m_bounds;

 private:
};

static_assert(RegularSurfaceConcept<RectangleSurface>,
              "RectangleSurface does not fulfill RegularSurfaceConcept");

}  // namespace Acts

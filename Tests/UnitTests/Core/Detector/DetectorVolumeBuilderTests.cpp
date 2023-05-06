// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/DetectorVolumeBuilder.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/interface/IExternalStructureBuilder.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/DetectorVolumeUpdators.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <vector>

using namespace Acts;
using namespace Acts::Experimental;

GeometryContext tContext;

/// @brief Mockup external structure builder
/// @tparam bounds_type the volume bounds type that is constructed
template <typename bounds_type>
class ExternalsBuilder : public IExternalStructureBuilder {
 public:
  ExternalsBuilder(const Transform3& transform, const bounds_type& bounds)
      : IExternalStructureBuilder(),
        m_transform(transform),
        m_bounds(std::move(bounds)) {}

  ExternalStructure construct(
      [[maybe_unused]] const GeometryContext& gctx) const final {
    return {m_transform, std::make_unique<bounds_type>(m_bounds),
            defaultPortalGenerator()};
  }

 private:
  Transform3 m_transform = Transform3::Identity();
  bounds_type m_bounds;
};

/// @brief  Mockup internal surface builder
/// @tparam surface_type the surface type to be constructed
/// @tparam bounds_type the bounds type that is contructed
template <typename surface_type, typename bounds_type>
class InternalSurfaceBuilder : public IInternalStructureBuilder {
 public:
  InternalSurfaceBuilder(const Transform3& transform, const bounds_type& bounds)
      : IInternalStructureBuilder(),
        m_transform(transform),
        m_bounds(std::move(bounds)) {}

  InternalStructure construct(
      [[maybe_unused]] const GeometryContext& gctx) const final {
    auto surface = Surface::makeShared<surface_type>(
        m_transform, std::make_shared<bounds_type>(m_bounds));
    return {{surface}, {}, tryAllPortalsAndSurfaces(), tryNoVolumes()};
  }

 private:
  Transform3 m_transform = Transform3::Identity();
  bounds_type m_bounds;
};

/// @brief  Mockup internal surface builder
/// @tparam surface_type the surface type to be constructed
/// @tparam bounds_type the bounds type that is contructed
template <typename bounds_type>
class InternalVolumeBuilder : public IInternalStructureBuilder {
 public:
  InternalVolumeBuilder(const Transform3& transform, const bounds_type& bounds)
      : IInternalStructureBuilder(),
        m_transform(transform),
        m_bounds(std::move(bounds)) {}

  InternalStructure construct(
      [[maybe_unused]] const GeometryContext& gctx) const final {
    auto bounds = std::make_unique<bounds_type>(m_bounds);
    auto portalGenerator = defaultPortalGenerator();
    auto volume = DetectorVolumeFactory::construct(
        portalGenerator, tContext, "InternalVolume", m_transform,
        std::move(bounds), tryAllPortals());
    return {{}, {volume}, tryAllPortals(), tryRootVolumes()};
  }

 private:
  Transform3 m_transform = Transform3::Identity();
  bounds_type m_bounds;
};

BOOST_AUTO_TEST_SUITE(Detector)

BOOST_AUTO_TEST_CASE(DetectorVolumeBuilder_Misconfigured) {
  // Internal and external structure builder is empty
  DetectorVolumeBuilder::Config dvCfg;
  dvCfg.auxilliary = "*** Test X * Misconfigued ***";
  dvCfg.name = "EmptyCylinder";
  dvCfg.externalsBuilder = nullptr;
  dvCfg.internalsBuilder = nullptr;

  BOOST_CHECK_THROW(auto a = DetectorVolumeBuilder(dvCfg),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(DetectorVolumeBuilder_EmptyVolume) {
  CylinderVolumeBounds cBounds(10, 100, 200);
  auto cBuilder = std::make_shared<ExternalsBuilder<CylinderVolumeBounds>>(
      Transform3::Identity(), cBounds);

  DetectorVolumeBuilder::Config dvCfg;
  dvCfg.auxilliary = "*** Test 0 - Empty Cylinder ***";
  dvCfg.name = "EmptyCylinder";
  dvCfg.externalsBuilder = cBuilder;
  dvCfg.internalsBuilder = nullptr;

  auto dvBuilder = std::make_shared<DetectorVolumeBuilder>(
      dvCfg, getDefaultLogger("DetectorVolumeBuilder", Logging::VERBOSE));

  RootDetectorVolumes roots;
  auto dvComponents = dvBuilder->construct(roots, tContext);

  BOOST_CHECK(roots.volumes.size() == 1u);
  BOOST_CHECK(dvComponents.portals.size() == 4u);

  // Get the volume and check it is empty
  auto volume = roots.volumes.front();
  BOOST_CHECK(volume->surfaces().empty());
  BOOST_CHECK(volume->volumes().empty());
}

BOOST_AUTO_TEST_CASE(DetectorVolumeBuilder_VolumeWithSurface) {
  CylinderVolumeBounds cvBounds(10, 100, 200);
  auto cBuilder = std::make_shared<ExternalsBuilder<CylinderVolumeBounds>>(
      Transform3::Identity(), cvBounds);

  CylinderBounds csBounds(55., 195.);
  auto sBuilder =
      std::make_shared<InternalSurfaceBuilder<CylinderSurface, CylinderBounds>>(
          Transform3::Identity(), csBounds);

  DetectorVolumeBuilder::Config dvCfg;
  dvCfg.auxilliary = "*** Test 1 - Cylinder with internal Surface ***";
  dvCfg.name = "CylinderWithSurface";
  dvCfg.externalsBuilder = cBuilder;
  dvCfg.internalsBuilder = sBuilder;

  auto dvBuilder = std::make_shared<DetectorVolumeBuilder>(
      dvCfg, getDefaultLogger("DetectorVolumeBuilder", Logging::VERBOSE));

  RootDetectorVolumes roots;
  auto dvComponents = dvBuilder->construct(roots, tContext);

  BOOST_CHECK(roots.volumes.size() == 1u);
  BOOST_CHECK(dvComponents.portals.size() == 4u);

  // Get the volume and check it is empty
  auto volume = roots.volumes.front();
  BOOST_CHECK(volume->surfaces().size() == 1u);
  BOOST_CHECK(volume->volumes().empty());
}

BOOST_AUTO_TEST_CASE(DetectorVolumeBuilder_VolumeWithVolume) {
  CylinderVolumeBounds cvBounds(10, 100, 200);
  auto cBuilder = std::make_shared<ExternalsBuilder<CylinderVolumeBounds>>(
      Transform3::Identity(), cvBounds);

  CylinderVolumeBounds ciBounds(15., 95., 195.);
  auto iBuilder = std::make_shared<InternalVolumeBuilder<CylinderVolumeBounds>>(
      Transform3::Identity(), ciBounds);

  DetectorVolumeBuilder::Config dvCfg;
  dvCfg.auxilliary = "*** Test 2 - Cylinder with internal Volume ***";
  dvCfg.name = "CylinderWithVolume";
  dvCfg.externalsBuilder = cBuilder;
  dvCfg.internalsBuilder = iBuilder;
  dvCfg.addToRoot = false;
  dvCfg.addInternalsToRoot = false;

  auto dvBuilder = std::make_shared<DetectorVolumeBuilder>(
      dvCfg, getDefaultLogger("DetectorVolumeBuilder", Logging::VERBOSE));

  RootDetectorVolumes roots;
  auto dvComponents = dvBuilder->construct(roots, tContext);

  BOOST_CHECK(roots.volumes.empty());
  BOOST_CHECK(dvComponents.portals.size() == 4u);
}

BOOST_AUTO_TEST_CASE(DetectorVolumeBuilder_VolumeWithVolumeToRoot) {
  CylinderVolumeBounds cvBounds(10, 100, 200);
  auto cBuilder = std::make_shared<ExternalsBuilder<CylinderVolumeBounds>>(
      Transform3::Identity(), cvBounds);

  CylinderVolumeBounds ciBounds(15., 95., 195.);
  auto iBuilder = std::make_shared<InternalVolumeBuilder<CylinderVolumeBounds>>(
      Transform3::Identity(), ciBounds);

  DetectorVolumeBuilder::Config dvCfg;
  dvCfg.auxilliary =
      "*** Test 3 - Cylinder with internal Volume, adding to root ***";
  dvCfg.name = "CylinderWithVolume";
  dvCfg.externalsBuilder = cBuilder;
  dvCfg.internalsBuilder = iBuilder;
  dvCfg.addToRoot = true;
  dvCfg.addInternalsToRoot = true;

  auto dvBuilder = std::make_shared<DetectorVolumeBuilder>(
      dvCfg, getDefaultLogger("DetectorVolumeBuilder", Logging::VERBOSE));

  RootDetectorVolumes roots;
  auto dvComponents = dvBuilder->construct(roots, tContext);

  BOOST_CHECK(roots.volumes.size() == 2u);
}

BOOST_AUTO_TEST_SUITE_END()

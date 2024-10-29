// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SegmentSeedingAlgorithm.hpp"

#include "Acts/Utilities/AngleHelpers.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <cstddef>
#include <ostream>
#include <stdexcept>

#include <Eigen/src/Core/Matrix.h>

namespace ActsExamples {

namespace {

struct SpacePoints {
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;

  std::vector<float> r;
  std::vector<float> phi;
  std::vector<float> theta;
  std::vector<float> d;

  std::vector<const SimSpacePoint*> sp;

  std::size_t size() const { return x.size(); }

  void push_back(float x_, float y_, float z_, float r_, float phi_,
                 float theta_, float d_, const SimSpacePoint* sp_) {
    x.push_back(x_);
    y.push_back(y_);
    z.push_back(z_);

    r.push_back(r_);
    phi.push_back(phi_);
    theta.push_back(theta_);
    d.push_back(d_);

    sp.push_back(sp_);
  }

  void push_back(float x_, float y_, float z_, const SimSpacePoint* sp_) {
    Eigen::Vector3f xyz(x_, y_, z_);

    push_back(x_, y_, z_, Acts::VectorHelpers::perp(xyz),
              Acts::VectorHelpers::phi(xyz), Acts::VectorHelpers::theta(xyz),
              xyz.norm(), sp_);
  }

  void insert(std::size_t i, float x_, float y_, float z_, float r_, float phi_,
              float theta_, float d_, const SimSpacePoint* sp_) {
    x.insert(x.begin() + i, x_);
    y.insert(y.begin() + i, y_);
    z.insert(z.begin() + i, z_);

    r.insert(r.begin() + i, r_);
    phi.insert(phi.begin() + i, phi_);
    theta.insert(theta.begin() + i, theta_);
    d.insert(d.begin() + i, d_);

    sp.insert(sp.begin() + i, sp_);
  }
};

}  // namespace

SegmentSeedingAlgorithm::SegmentSeedingAlgorithm(
    SegmentSeedingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("SegmentSeedingAlgorithm", lvl), m_cfg(std::move(cfg)) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point input collections");
  }

  for (const auto& spName : m_cfg.inputSpacePoints) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputSpacePoints.emplace_back(
        std::make_unique<ReadDataHandle<SimSpacePointContainer>>(
            this,
            "InputSpacePoints#" + std::to_string(m_inputSpacePoints.size())));
    handle->initialize(spName);
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds output collection");
  }

  m_outputSeeds.initialize(m_cfg.outputSeeds);
}

ProcessCode SegmentSeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  SpacePoints spacePoints;
  for (const auto& isp : m_inputSpacePoints) {
    for (const auto& spacePoint : (*isp)(ctx)) {
      spacePoints.push_back(spacePoint.x(), spacePoint.y(), spacePoint.z(),
                            &spacePoint);
    }
  }

  const float thetaMin = Acts::AngleHelpers::thetaFromEta(3.0f);
  const float thetaMax = Acts::AngleHelpers::thetaFromEta(-3.0f);
  const std::size_t thetaSegments = 101;
  const float thetaStep = (thetaMax - thetaMin) / (thetaSegments - 1);

  const std::size_t phiSegments = 101;
  const float phiStep = 2 * M_PI / phiSegments;

  const float minD = 3 * Acts::UnitConstants::cm;
  const float maxD = 15 * Acts::UnitConstants::cm;

  auto getThetaSegment = [&](float theta) -> int {
    return std::min(static_cast<std::size_t>((theta - thetaMin) / thetaStep),
                    thetaSegments - 1);
  };

  auto getPhiSegment = [&](float phi) -> int {
    return std::min(static_cast<std::size_t>(phi / phiStep), phiSegments - 1);
  };

  std::vector<std::vector<SpacePoints>> segmented;
  segmented.resize(thetaSegments);
  for (auto& segment : segmented) {
    segment.resize(phiSegments);
  }

  for (std::size_t i = 0; i < spacePoints.size(); ++i) {
    const auto& sp = spacePoints;
    const auto thetaSegment = getThetaSegment(sp.theta[i]);
    const auto phiSegment = getPhiSegment(sp.phi[i]);
    for (int j = thetaSegment - 2; j <= thetaSegment + 2; ++j) {
      if (j < 0 || j >= thetaSegments) {
        continue;
      }
      for (int k = phiSegment - 2; k <= phiSegment + 2; ++k) {
        if (k < 0 || k >= phiSegments) {
          continue;
        }
        auto& segment = segmented[j][k];
        auto it = std::lower_bound(segment.d.begin(), segment.d.end(), sp.d[i]);
        const auto idx = std::distance(segment.d.begin(), it);
        segment.insert(idx, sp.x[i], sp.y[i], sp.z[i], sp.r[i], sp.phi[i],
                       sp.theta[i], sp.d[i], sp.sp[i]);
      }
    }
  }

  SimSeedContainer seeds;

  auto findSeeds = [&](const SpacePoints& sp) {
    for (std::size_t i = 0; i < sp.size(); ++i) {
      for (std::size_t j = i + 1; j < sp.size(); ++j) {
        const float d12 = sp.d[j] - sp.d[i];
        if (d12 < minD) {
          continue;
        }
        if (d12 > maxD) {
          break;
        }

        for (std::size_t k = j + 1; k < sp.size(); ++k) {
          const float d23 = sp.d[k] - sp.d[j];
          if (d23 < minD) {
            continue;
          }
          if (d23 > maxD) {
            break;
          }

          seeds.push_back(SimSeed(*sp.sp[i], *sp.sp[j], *sp.sp[k]));
        }
      }
    }
  };

  for (std::size_t i = 0; i < thetaSegments; ++i) {
    for (std::size_t j = 0; j < phiSegments; ++j) {
      const auto& segment = segmented[i][j];
      findSeeds(segment);
    }
  }

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePoints.size() << " space points");

  m_outputSeeds(ctx, std::move(seeds));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

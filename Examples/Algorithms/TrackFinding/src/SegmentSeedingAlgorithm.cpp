// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SegmentSeedingAlgorithm.hpp"

#include "Acts/Utilities/AngleHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <cstddef>
#include <ostream>
#include <stdexcept>

#include <Eigen/src/Core/Matrix.h>
#include <boost/container/flat_set.hpp>

namespace ActsExamples {

namespace {

template <typename Scalar>
struct Circle {
  Scalar x;
  Scalar y;
  Scalar r;
};

template <typename Scalar>
Circle<Scalar> fitCircle(Scalar x1, Scalar y1, Scalar x2, Scalar y2, Scalar x3,
                         Scalar y3) {
  // Calculate the determinants
  Scalar a = x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2;
  if (std::abs(a) < 1e-10f) {
    throw std::runtime_error(
        "The points are collinear and do not define a unique circle.");
  }

  Scalar b1 = (x1 * x1 + y1 * y1);
  Scalar b2 = (x2 * x2 + y2 * y2);
  Scalar b3 = (x3 * x3 + y3 * y3);

  Scalar x_c = (b1 * (y2 - y3) + b2 * (y3 - y1) + b3 * (y1 - y2)) / (2 * a);
  Scalar y_c = (b1 * (x3 - x2) + b2 * (x1 - x3) + b3 * (x2 - x1)) / (2 * a);
  Scalar radius = std::sqrt((x_c - x1) * (x_c - x1) + (y_c - y1) * (y_c - y1));

  return {x_c, y_c, radius};
}

template <typename Scalar>
double distanceToCircle(const Circle<Scalar>& circle, double xp, double yp) {
  // Calculate the distance from the center of the circle to the point
  double dx = xp - circle.x;
  double dy = yp - circle.y;
  double distanceToCenter = std::sqrt(dx * dx + dy * dy);

  // Calculate the distance to the closest point on the circle
  return std::abs(distanceToCenter - circle.r);
}

struct SpacePoints final {
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;

  std::vector<float> r;
  std::vector<float> phi;
  std::vector<float> theta;
  std::vector<float> d;

  std::vector<const SimSpacePoint*> sp;

  template <bool ReadOnly>
  struct ProxyImpl final {
    using Container =
        std::conditional_t<ReadOnly, const SpacePoints, SpacePoints>;

    Container* container;
    std::size_t index;

    ProxyImpl(Container* sp_, std::size_t i_) : container(sp_), index(i_) {}

    template <bool OtherReadOnly>
    ProxyImpl(const ProxyImpl<OtherReadOnly>& other)
      requires(ReadOnly)
        : container(other.container), index(other.index) {}

    float& x()
      requires(!ReadOnly)
    {
      return container->x[index];
    }
    float& y()
      requires(!ReadOnly)
    {
      return container->y[index];
    }
    float& z()
      requires(!ReadOnly)
    {
      return container->z[index];
    }

    float& r()
      requires(!ReadOnly)
    {
      return container->r[index];
    }
    float& phi()
      requires(!ReadOnly)
    {
      return container->phi[index];
    }
    float& theta()
      requires(!ReadOnly)
    {
      return container->theta[index];
    }
    float& d()
      requires(!ReadOnly)
    {
      return container->d[index];
    }

    const SimSpacePoint*& sp()
      requires(!ReadOnly)
    {
      return container->sp[index];
    }

    float x() const { return container->x[index]; }
    float y() const { return container->y[index]; }
    float z() const { return container->z[index]; }

    Eigen::Vector3f xyz() const { return Eigen::Vector3f(x(), y(), z()); }
    Eigen::Vector2f zr() const { return Eigen::Vector2f(z(), r()); }

    float r() const { return container->r[index]; }
    float phi() const { return container->phi[index]; }
    float theta() const { return container->theta[index]; }
    float d() const { return container->d[index]; }

    const SimSpacePoint* sp() const { return container->sp[index]; }

    void derive()
      requires(!ReadOnly)
    {
      r() = Acts::VectorHelpers::perp(xyz());
      phi() = Acts::VectorHelpers::phi(xyz());
      theta() = Acts::VectorHelpers::theta(xyz());
      d() = xyz().norm();
    }
  };

  using Proxy = ProxyImpl<false>;
  using ConstProxy = ProxyImpl<true>;

  template <bool ReadOnly>
  struct IteratorImpl final {
    using Container =
        std::conditional_t<ReadOnly, const SpacePoints, SpacePoints>;

    Container* container;
    std::size_t index;

    IteratorImpl(Container* sp_, std::size_t i_) : container(sp_), index(i_) {}

    Proxy operator*() { return Proxy(container, index); }
    IteratorImpl& operator++() {
      ++index;
      return *this;
    }
    bool operator==(const IteratorImpl& other) const {
      return index == other.index;
    }
  };

  using Iterator = IteratorImpl<false>;
  using ConstIterator = IteratorImpl<true>;

  Iterator begin() { return Iterator(this, 0); }
  Iterator end() { return Iterator(this, size()); }
  ConstIterator begin() const { return ConstIterator(this, 0); }
  ConstIterator end() const { return ConstIterator(this, size()); }

  Proxy operator[](std::size_t i) { return Proxy(this, i); }
  ConstProxy operator[](std::size_t i) const { return ConstProxy(this, i); }

  std::size_t size() const { return x.size(); }

  Proxy emplace_back() {
    x.push_back(0);
    y.push_back(0);
    z.push_back(0);
    r.push_back(0);
    phi.push_back(0);
    theta.push_back(0);
    d.push_back(0);
    sp.push_back(nullptr);
    return Proxy(this, size() - 1);
  }

  Proxy insert(std::size_t i) {
    x.insert(x.begin() + i, 0);
    y.insert(y.begin() + i, 0);
    z.insert(z.begin() + i, 0);
    r.insert(r.begin() + i, 0);
    phi.insert(phi.begin() + i, 0);
    theta.insert(theta.begin() + i, 0);
    d.insert(d.begin() + i, 0);
    sp.insert(sp.begin() + i, nullptr);
    return Proxy(this, i);
  }

  Proxy insert(std::size_t i, const ConstProxy& p) {
    x.insert(x.begin() + i, p.x());
    y.insert(y.begin() + i, p.y());
    z.insert(z.begin() + i, p.z());
    r.insert(r.begin() + i, p.r());
    phi.insert(phi.begin() + i, p.phi());
    theta.insert(theta.begin() + i, p.theta());
    d.insert(d.begin() + i, p.d());
    sp.insert(sp.begin() + i, p.sp());
    return Proxy(this, i);
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
      auto sp = spacePoints.emplace_back();
      sp.x() = spacePoint.x();
      sp.y() = spacePoint.y();
      sp.z() = spacePoint.z();
      sp.sp() = &spacePoint;
      sp.derive();
    }
  }

  const float thetaMin = Acts::AngleHelpers::thetaFromEta(3.5f);
  const float thetaMax = Acts::AngleHelpers::thetaFromEta(-3.5f);
  const std::size_t thetaSegments = 11;
  const float thetaStep = (thetaMax - thetaMin) / (thetaSegments - 1);

  const std::size_t phiSegments = 101;
  const float phiStep = 2 * M_PI / phiSegments;

  const float minZ = -150 * Acts::UnitConstants::mm;
  const float maxZ = 150 * Acts::UnitConstants::mm;
  const float minD = 30 * Acts::UnitConstants::mm;
  const float maxD = 300 * Acts::UnitConstants::mm;
  const float max_l13_d2 = 1 * Acts::UnitConstants::mm;
  const float maxAbsD0 = 5 * Acts::UnitConstants::mm;

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

  for (auto sp : spacePoints) {
    const auto thetaSegment = getThetaSegment(sp.theta());
    const auto phiSegment = getPhiSegment(sp.phi());
    for (int j = thetaSegment - 1; j <= thetaSegment + 1; ++j) {
      if (j < 0 || j >= thetaSegments) {
        continue;
      }
      for (int k = phiSegment - 1; k <= phiSegment + 1; ++k) {
        if (k < 0 || k >= phiSegments) {
          continue;
        }
        auto& segment = segmented[j][k];
        auto it = std::lower_bound(segment.d.begin(), segment.d.end(), sp.d());
        const auto idx = std::distance(segment.d.begin(), it);
        segment.insert(idx, sp);
      }
    }
  }

  SimSeedContainer seeds;
  boost::container::flat_set<const SimSpacePoint*> usedSpacePoints;

  auto findSeeds = [&](const SpacePoints& spacePoints) {
    for (std::size_t i = 0; i < spacePoints.size(); ++i) {
      auto sp1 = spacePoints[i];

      // skip sp 1 if already used
      if (usedSpacePoints.contains(sp1.sp())) {
        continue;
      }

      for (std::size_t j = i + 1; j < spacePoints.size(); ++j) {
        auto sp2 = spacePoints[j];

        // skip sp 2 if already used
        if (usedSpacePoints.contains(sp2.sp())) {
          continue;
        }

        // approximate distance between sp 1 and 2
        const float approx_d12 = sp2.d() - sp1.d();
        if (approx_d12 < minD) {
          continue;
        }
        if (approx_d12 > maxD) {
          break;
        }

        // z0 with sp 1 and 2
        const float z0_12 =
            sp2.z() - sp2.r() * (sp2.z() - sp1.z()) / (sp2.r() - sp1.r());
        if (z0_12 < minZ || z0_12 > maxZ) {
          continue;
        }

        for (std::size_t k = j + 1; k < spacePoints.size(); ++k) {
          auto sp3 = spacePoints[k];

          // skip sp 3 if already used
          if (usedSpacePoints.contains(sp3.sp())) {
            continue;
          }

          // approximate distance between sp 2 and 3
          const float approx_d23 = sp3.d() - sp2.d();
          if (approx_d23 < minD) {
            continue;
          }
          if (approx_d23 > maxD) {
            break;
          }

          // z0 with sp 1 and 3
          const float z0_13 =
              sp3.z() - sp3.r() * (sp3.z() - sp1.z()) / (sp3.r() - sp1.r());
          if (z0_13 < minZ || z0_13 > maxZ) {
            continue;
          }

          // sp 2 distance to sp 1 and 3 line
          // y = m * x + k
          const float l13_k = -z0_13 * sp3.r() / (sp3.z() - z0_13);
          const float l13_m = (sp3.r() - sp1.r()) / (sp3.z() - sp1.z());
          // https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Another_formula
          // d = |k + m * x0 - y0| / sqrt(1 + m^2)
          const float l13_d2 = std::abs(l13_k + l13_m * sp2.z() - sp2.r()) /
                               std::sqrt(1 + l13_m * l13_m);
          if (l13_d2 > max_l13_d2) {
            continue;
          }

          // fit a circle through sp 1, 2, 3
          auto circle =
              fitCircle(sp1.x(), sp1.y(), sp2.x(), sp2.y(), sp3.x(), sp3.y());
          float d0 = distanceToCircle(circle, 0, 0);
          if (d0 > maxAbsD0) {
            continue;
          }

          seeds.push_back(SimSeed(*sp1.sp(), *sp2.sp(), *sp3.sp()));
          usedSpacePoints.insert(sp1.sp());
          usedSpacePoints.insert(sp2.sp());
          usedSpacePoints.insert(sp3.sp());
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

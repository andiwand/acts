// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/Rz/RzTransport.hpp"

#include <algorithm>
#include <numbers>

namespace Acts::Experimental {

namespace {

/// Below this turning rate the transverse motion is a straight line for every
/// length a tracker has.
constexpr double kStraightKappa = 1e-12;

/// Smallest positive path length of the straight line `p + s d` to a radius
std::optional<double> straightPathToCylinder(const RzVector& v, double radius) {
  const double px = v[eRzPos0];
  const double py = v[eRzPos1];
  const double dx = v[eRzDir0];
  const double dy = v[eRzDir1];
  const double a = dx * dx + dy * dy;
  if (a == 0.) {
    return std::nullopt;
  }
  const double b = px * dx + py * dy;
  const double c = px * px + py * py - radius * radius;
  const double disc = b * b - a * c;
  if (disc < 0.) {
    return std::nullopt;
  }
  const double sq = std::sqrt(disc);
  // the root nearer the state first, in the cancellation-free form
  const double q = -(b + std::copysign(sq, b));
  const double s1 = q / a;
  const double s2 = q != 0. ? c / q : s1;
  const double lo = std::min(s1, s2);
  const double hi = std::max(s1, s2);
  if (lo > 0.) {
    return lo;
  }
  if (hi > 0.) {
    return hi;
  }
  return std::nullopt;
}

}  // namespace

std::optional<double> RzHelix::pathToCylinder(const RzVector& v,
                                              double radius) const {
  const double k = kappa(v);
  if (std::abs(k) < kStraightKappa) {
    return straightPathToCylinder(v, radius);
  }
  // Circle against circle: the crossing points of the track's transverse
  // circle with the cylinder's, then the turning angle from the chord to
  // the nearer one ahead. No transcendental but one asin, and no
  // cancellation: d^2 - rho^2 = |p|^2 - 2 p.r0 is formed directly.
  const double px = v[eRzPos0];
  const double py = v[eRzPos1];
  const double dx = v[eRzDir0];
  const double dy = v[eRzDir1];
  const double dT2 = dx * dx + dy * dy;
  if (dT2 == 0.) {
    return std::nullopt;
  }
  const double dT = std::sqrt(dT2);
  const double rho = dT / std::abs(k);
  // from the centre to the state
  const double r0x = -dy / k;
  const double r0y = dx / k;
  const double cx = px - r0x;
  const double cy = py - r0y;
  const double p2 = px * px + py * py;
  const double d2mr2 = p2 - 2. * (px * r0x + py * r0y);
  const double d = std::sqrt(d2mr2 + rho * rho);
  if (d == 0.) {
    return std::nullopt;
  }
  // distance from the axis to the radical line and the half chord on it
  const double x = (d2mr2 + radius * radius) / (2. * d);
  const double y2 = radius * radius - x * x;
  if (y2 < 0.) {
    return std::nullopt;
  }
  const double y = std::sqrt(y2);
  const double ux = cx / d;
  const double uy = cy / d;
  // the two crossings, the one ahead by the shorter chord
  double best = -1.;
  for (const double sign : {1., -1.}) {
    const double qx = ux * x - sign * uy * y;
    const double qy = uy * x + sign * ux * y;
    const double tx = qx - px;
    const double ty = qy - py;
    if (tx * dx + ty * dy <= 0.) {
      continue;
    }
    const double l2 = tx * tx + ty * ty;
    if (best < 0. || l2 < best) {
      best = l2;
    }
  }
  if (best < 0.) {
    // both crossings more than half a turn ahead
    return pathToCylinderClosedForm(v, radius);
  }
  const double half = std::min(1., std::sqrt(best) / (2. * rho));
  return 2. * std::asin(half) / std::abs(k);
}

std::optional<double> RzHelix::pathToCylinderClosedForm(const RzVector& v,
                                                        double radius) const {
  const double k = kappa(v);
  if (std::abs(k) < kStraightKappa) {
    return straightPathToCylinder(v, radius);
  }

  const double px = v[eRzPos0];
  const double py = v[eRzPos1];
  const double dx = v[eRzDir0];
  const double dy = v[eRzDir1];
  // centre of the transverse circle and the state's position relative to it
  const double cx = px + dy / k;
  const double cy = py - dx / k;
  const double r0x = -dy / k;
  const double r0y = dx / k;
  // |p_T(gamma)|^2 = |c|^2 + rho^2 + 2 (A cos(gamma) + B sin(gamma)), the
  // position having turned by gamma = -kappa s about the centre
  const double a = cx * r0x + cy * r0y;
  const double b = r0x * cy - r0y * cx;
  const double c =
      0.5 * (radius * radius - (cx * cx + cy * cy) - (r0x * r0x + r0y * r0y));
  const double amp = std::hypot(a, b);
  if (amp == 0. || std::abs(c) > amp) {
    return std::nullopt;
  }
  const double phase = std::atan2(b, a);
  const double delta = std::acos(std::clamp(c / amp, -1., 1.));

  constexpr double twoPi = 2. * std::numbers::pi;
  // the two turning angles per turn, each carried into the forward direction:
  // s = -gamma / kappa > 0 needs gamma of the opposite sign to kappa
  std::optional<double> best;
  for (const double gamma : {phase + delta, phase - delta}) {
    double g = std::fmod(gamma, twoPi);
    if (k > 0. && g > 0.) {
      g -= twoPi;
    } else if (k < 0. && g < 0.) {
      g += twoPi;
    }
    const double s = -g / k;
    if (s > 0. && (!best.has_value() || s < *best)) {
      best = s;
    }
  }
  if (!best.has_value()) {
    return std::nullopt;
  }

  // one Newton polish on the exact residual takes out the rounding of the
  // closed form, which loses digits for stiff tracks where |c| ~ rho
  double s = *best;
  for (int i = 0; i < 2; ++i) {
    RzVector w = v;
    step(w, s);
    const double f =
        w[eRzPos0] * w[eRzPos0] + w[eRzPos1] * w[eRzPos1] - radius * radius;
    const double df = 2. * (w[eRzPos0] * w[eRzDir0] + w[eRzPos1] * w[eRzDir1]);
    if (df == 0.) {
      break;
    }
    s -= f / df;
  }
  return s;
}

std::optional<double> RzHelix::pathToPlane(const RzVector& v,
                                           const Vector3& point,
                                           const Vector3& normal) const {
  const Vector3 p = v.segment<3>(eRzPos0);
  const Vector3 d = v.segment<3>(eRzDir0);
  const double along = normal.dot(d);
  if (std::abs(along) < 1e-9) {
    return std::nullopt;
  }
  const double straight = normal.dot(point - p) / along;
  // Start from the parabola: the transverse bend over the straight path
  // moves the crossing by half the acceleration times s^2, which leaves an
  // error of the order of kappa^2 s^3, well inside one Newton step of the
  // tolerance for any module distance a layer has.
  const double k = kappa(v);
  const double accel = k * (normal.x() * d.y() - normal.y() * d.x());
  double s = straight;
  if (accel != 0.) {
    // along s + accel s^2 / 2 = along * straight, the root nearest straight
    const double disc = along * along + 2. * accel * along * straight;
    if (disc > 0.) {
      const double q = -(along + std::copysign(std::sqrt(disc), along));
      const double s1 = q / accel;
      const double s2 = -2. * along * straight / q;
      s = std::abs(s1 - straight) < std::abs(s2 - straight) ? s1 : s2;
    }
  }
  for (int i = 0; i < 5; ++i) {
    RzVector w = v;
    step(w, s);
    const double f = normal.dot(w.segment<3>(eRzPos0) - point);
    const double df = normal.dot(w.segment<3>(eRzDir0));
    if (df == 0.) {
      return std::nullopt;
    }
    const double ds = f / df;
    s -= ds;
    // one step from the parabola is a converged one: its own size says so
    if (std::abs(ds) < 1e-5) {
      break;
    }
  }
  return s;
}

double RzHelix::pathToPerigee(const RzVector& v) const {
  const double k = kappa(v);
  const double px = v[eRzPos0];
  const double py = v[eRzPos1];
  const double dx = v[eRzDir0];
  const double dy = v[eRzDir1];
  const double dt2 = dx * dx + dy * dy;
  if (dt2 == 0.) {
    return 0.;
  }
  if (std::abs(k) < kStraightKappa) {
    return -(px * dx + py * dy) / dt2;
  }
  const double cx = px + dy / k;
  const double cy = py - dx / k;
  const double cn = std::hypot(cx, cy);
  if (cn == 0.) {
    return 0.;
  }
  // the point of the circle nearest the axis is the centre pulled in by one
  // radius; the turning angle from here to there, shortest way round
  const double r0x = -dy / k;
  const double r0y = dx / k;
  const double rho = std::hypot(r0x, r0y);
  const double rqx = -rho * cx / cn;
  const double rqy = -rho * cy / cn;
  const double gamma = std::atan2(r0x * rqy - r0y * rqx, r0x * rqx + r0y * rqy);
  return -gamma / k;
}

}  // namespace Acts::Experimental

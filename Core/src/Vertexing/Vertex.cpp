// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {

Vertex::Vertex(const Vector3& position) {
  m_position = position;
  m_seedPosition = position;
}

Vertex::Vertex(const Vector3& position, const SquareMatrix3& covariance,
               std::vector<TrackAtVertex> tracks)
    : m_tracksAtVertex(std::move(tracks)) {
  m_position = position;
  m_seedPosition = position;
  m_covariance.block<3, 3>(ePos0, ePos0) = covariance;
}

Vector3 Vertex::position() const {
  return m_position;
}

SquareMatrix3 Vertex::covariance() const {
  return m_covariance.block<3, 3>(ePos0, ePos0);
}

const std::vector<TrackAtVertex>& Vertex::tracks() const {
  return m_tracksAtVertex;
}

std::pair<double, double> Vertex::fitQuality() const {
  return {m_chiSquared, m_numberDoF};
}

void Vertex::setPosition(const Vector3& position) {
  m_position = position;
}

void Vertex::setCovariance(const SquareMatrix3& covariance) {
  m_covariance.setZero();
  m_covariance.block<3, 3>(ePos0, ePos0) = covariance;
}

void Vertex::setTracksAtVertex(std::vector<TrackAtVertex> tracks) {
  m_tracksAtVertex = std::move(tracks);
}

void Vertex::setFitQuality(double chiSquared, double numberDoF) {
  m_chiSquared = chiSquared;
  m_numberDoF = numberDoF;
}

void Vertex::setFitQuality(std::pair<double, double> fitQuality) {
  m_chiSquared = fitQuality.first;
  m_numberDoF = fitQuality.second;
}

}  // namespace Acts

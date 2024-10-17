// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/KalmanVertexUpdater.hpp"

#include "Acts/Vertexing/Vertex.hpp"

#include <stdexcept>

namespace Acts::KalmanVertexUpdater {

namespace detail {
void updateVertexWithTrackImpl(Vertex& vtx, TrackAtVertex& trk, int sign);

void updateTrackWithVertexImpl(TrackAtVertex& track, const Vertex& vtx);
}  // namespace detail

// The two functions don't contain any of the actual update code, they
// only dispatch into templated functions, effectively doing a
// runtime-to-compile time conversion.

void updateVertexWithTrack(Vertex& vtx, TrackAtVertex& trk) {
  detail::updateVertexWithTrackImpl(vtx, trk, 1);
}

void updateTrackWithVertex(TrackAtVertex& track, const Vertex& vtx) {
  detail::updateTrackWithVertexImpl(track, vtx);
}

}  // namespace Acts::KalmanVertexUpdater

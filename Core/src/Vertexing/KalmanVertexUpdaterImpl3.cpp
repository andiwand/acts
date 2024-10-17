// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/KalmanVertexUpdater.hpp"
#include "Acts/Vertexing/detail/KalmanVertexUpdaterImpl.hpp"

void Acts::KalmanVertexUpdater::detail::updateVertexWithTrackImpl(
    Vertex& vtx, TrackAtVertex& trk, int sign);

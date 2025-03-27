// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/VertexSelector.hpp"

#include <ostream>
#include <stdexcept>

namespace ActsExamples {

VertexSelector::VertexSelector(const Config& config, Acts::Logging::Level level)
    : IAlgorithm("VertexSelector", level), m_cfg(config) {
  if (m_cfg.inputVertices.empty()) {
    throw std::invalid_argument("Missing input vertices collection");
  }
  if (m_cfg.outputVertices.empty()) {
    throw std::invalid_argument("Missing output vertices collection");
  }

  m_inputVertices.initialize(m_cfg.inputVertices);
  m_outputVertices.initialize(m_cfg.outputVertices);

  ACTS_DEBUG("primary vertex ID [" << m_cfg.minPrimaryVertexId << ","
                                   << m_cfg.maxPrimaryVertexId << ")");
}

ProcessCode VertexSelector::execute(const AlgorithmContext& ctx) const {
  // prepare input/ output types
  const SimVertexContainer& inputVertices = m_inputVertices(ctx);

  // helper functions to select tracks
  auto within = [](auto x, auto min, auto max) {
    return (min <= x) && (x < max);
  };

  auto isValidVertex = [&](const SimVertex& v) {
    const bool validPrimaryVertexId =
        within(v.vertexId().vertexPrimary(), m_cfg.minPrimaryVertexId,
               m_cfg.maxPrimaryVertexId);

    return validPrimaryVertexId;
  };

  SimVertexContainer outputVertices;
  outputVertices.reserve(outputVertices.size());

  // copy selected vertices
  for (const auto& inputVertex : inputVertices) {
    if (!isValidVertex(inputVertex)) {
      continue;
    }

    outputVertices.insert(outputVertices.end(), inputVertex);
  }
  outputVertices.shrink_to_fit();

  ACTS_DEBUG("event " << ctx.eventNumber << " selected "
                      << outputVertices.size() << " from "
                      << inputVertices.size() << " vertices");

  m_outputVertices(ctx, std::move(outputVertices));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

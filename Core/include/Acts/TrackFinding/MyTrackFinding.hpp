// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"

namespace Acts {

template <typename source_link_iterator_t>
using SourceLinkAccessorDelegate =
    Delegate<std::pair<source_link_iterator_t, source_link_iterator_t>(
        const Surface&)>;

template <typename propagator_t, typename traj_t>
class MyTrackFinding {
 public:
  struct Config {};

  MyTrackFinding(const Config& cfg) : m_cfg(cfg) {}

  template <typename source_link_iterator_t, typename start_parameters_t,
            typename track_container_t, template <typename> class holder_t,
            typename parameters_t = BoundTrackParameters>
  auto findTracks(
      const start_parameters_t& initialParameters,
      TrackContainer<track_container_t, traj_t, holder_t>& trackContainer) const
      -> Result<std::vector<
          typename std::decay_t<decltype(trackContainer)>::TrackProxy>> {
    using TrackContainer = typename std::decay_t<decltype(trackContainer)>;
    using SourceLinkAccessor =
        SourceLinkAccessorDelegate<source_link_iterator_t>;
  }

 private:
  Config m_cfg;
};

}  // namespace Acts

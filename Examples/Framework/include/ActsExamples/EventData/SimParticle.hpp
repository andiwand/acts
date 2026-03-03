// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/ParticleContainer.hpp"

namespace ActsExamples {

using SimParticleIndex = ::ActsFatras::ParticleIndex;

using SimParticleColumns = ::ActsFatras::ParticleColumns;

using SimParticleContainer = ::ActsFatras::ParticleContainer;
using SimParticle = ::ActsFatras::ConstParticleProxy;
using MutableSimParticle = ::ActsFatras::MutableParticleProxy;

using SimBarcode = ::ActsFatras::Barcode;

class SelectedSimParticles final
    : public Acts::detail::ContainerSubset<
          SelectedSimParticles, SelectedSimParticles, SimParticleContainer,
          SimParticle, SimParticleIndex, true> {
 public:
  /// Base class type
  using Base =
      Acts::detail::ContainerSubset<SelectedSimParticles, SelectedSimParticles,
                                    SimParticleContainer, SimParticle,
                                    SimParticleIndex, true>;

  using Base::Base;
};

}  // namespace ActsExamples

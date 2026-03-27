// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <span>

namespace ActsFatras {

using ParticleIndex = std::uint32_t;
using ParticleIndexSubset = std::span<const ParticleIndex>;

class ParticleContainer;

template <bool>
class ParticleProxy;
using MutableParticleProxy = ParticleProxy<false>;
using ConstParticleProxy = ParticleProxy<true>;

template <bool>
class ParticleSubset;
using ConstParticleSubset = ParticleSubset<true>;
using MutableParticleSubset = ParticleSubset<false>;

template <typename T, bool>
class ParticleColumnProxy;
template <typename T>
using MutableParticleColumnProxy = ParticleColumnProxy<T, false>;
template <typename T>
using ConstParticleColumnProxy = ParticleColumnProxy<T, true>;

using VertexIndex = std::uint32_t;
using VertexIndexSubset = std::span<const VertexIndex>;

class VertexContainer;

template <bool>
class VertexProxy;
using MutableVertexProxy = VertexProxy<false>;
using ConstVertexProxy = VertexProxy<true>;

using HitIndex = std::uint32_t;
using HitIndexSubset = std::span<const HitIndex>;

class HitContainer;

}  // namespace ActsFatras

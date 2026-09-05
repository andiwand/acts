// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>

namespace ActsExamples {

/// The default log level of everything in the examples framework.
constexpr Acts::Logging::Level defaultLogLevel = Acts::Logging::INFO;

/// Make the logger the examples framework uses by default.
///
/// Single point of definition for the default output policy, so the framework
/// can be given a different one without touching every component. Configs
/// default their `logger` to this, and leave it unnamed so the component that
/// receives the config names it.
///
/// @param name Name of the logger. Empty leaves it to be named by its owner.
/// @param level Log level of the logger
/// @return The logger
std::shared_ptr<const Acts::Logger> makeDefaultLogger(
    const std::string& name = "", Acts::Logging::Level level = defaultLogLevel);

/// The logger a component called @p name should use.
///
/// An unnamed logger is cloned under that name, which also gives each copy of
/// a config its own logger instead of sharing one. A logger that already has a
/// name was named deliberately and is used as is. A null logger falls back to
/// the framework default under that name.
///
/// @param logger The logger from the component's config
/// @param name The name of the component
/// @return The logger to use
std::shared_ptr<const Acts::Logger> makeLogger(
    const std::shared_ptr<const Acts::Logger>& logger, const std::string& name);

}  // namespace ActsExamples

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/Logging.hpp"

std::shared_ptr<const Acts::Logger> ActsExamples::makeDefaultLogger(
    const std::string& name, Acts::Logging::Level level) {
  return Acts::getDefaultLogger(name, level);
}

std::shared_ptr<const Acts::Logger> ActsExamples::makeLogger(
    const std::shared_ptr<const Acts::Logger>& logger,
    const std::string& name) {
  if (logger == nullptr) {
    return makeDefaultLogger(name);
  }
  return logger->name().empty() ? logger->clone(name) : logger;
}

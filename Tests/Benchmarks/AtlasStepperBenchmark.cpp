// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/AtlasStepper.hpp"

#include "StepperBenchmarkCommons.hpp"

using namespace Acts;
using namespace ActsTests;

using Stepper = AtlasStepper;

int main(int argc, char* argv[]) {
  BenchmarkStepper benchmark;
  bool gradient = false;
  bool numericalGradient = false;
  bool fieldOff = false;
  double helixStep = 0.;
  po::options_description extra("AtlasStepper options");
  extra.add_options()("gradient", po::bool_switch(&gradient),
                      "include analytic field gradients")(
      "numerical-gradient", po::bool_switch(&numericalGradient),
      "include finite-difference field gradients")(
      "helix-step", po::value<double>(&helixStep)->default_value(0.),
      "constant-field shortcut below this step length in mm")(
      "field-off", po::bool_switch(&fieldOff),
      "disable magnetic deflection and use the straight-line kernel");
  if (auto ret = benchmark.parseOptions(argc, argv, &extra)) {
    return *ret;
  }
  auto bField = benchmark.makeField();
  Stepper::Config config{std::move(bField), {}};
  if (gradient && !numericalGradient) {
    // The benchmark uses a constant field, so the analytic gradient is zero.
    config.fieldGradient = [](const Vector3&, MagneticFieldProvider::Cache&) {
      return Result<SquareMatrix3>::success(SquareMatrix3::Zero());
    };
  }
  Stepper stepper(config);
  benchmark.run(stepper, "AtlasStepper", [&](auto& options) {
    options.useFieldGradient = gradient || numericalGradient;
    options.helixStep = helixStep * UnitConstants::mm;
    options.useMagneticField = !fieldOff;
  });
  return 0;
}

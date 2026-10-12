// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsTests/CommonHelpers/BenchmarkTools.hpp"

#include <iostream>
#include <limits>
#include <utility>

#include <boost/program_options.hpp>

namespace ActsTests {

namespace po = boost::program_options;
using namespace Acts;
using namespace Acts::UnitLiterals;

struct BenchmarkStepper {
  unsigned int toys{};
  double ptInGeV{};
  double BzInT{};
  double maxPathInM{};
  unsigned int lvl{};
  bool withCov{};
  double maxStepInMM{};
  double stepTolerance{};

  std::optional<int> parseOptions(
      int argc, char* argv[],
      const po::options_description* extraOptions = nullptr) {
    try {
      po::options_description desc("Allowed options");
      // clang-format off
      desc.add_options()
        ("help", "produce help message")
        ("toys",po::value<unsigned int>(&toys)->default_value(20000),"number of tracks to propagate")
        ("pT",po::value<double>(&ptInGeV)->default_value(1),"transverse momentum in GeV")
        ("B",po::value<double>(&BzInT)->default_value(2),"z-component of B-field in T")
        ("path",po::value<double>(&maxPathInM)->default_value(5),"maximum path length in m")
        ("cov",po::value<bool>(&withCov)->default_value(true),"propagation with covariance matrix")
        ("verbose",po::value<unsigned int>(&lvl)->default_value(Logging::INFO),"logging level")
        ("max-step",po::value<double>(&maxStepInMM)->default_value(std::numeric_limits<double>::max()),"maximum step length in mm")
        ("tolerance",po::value<double>(&stepTolerance)->default_value(1e-4),"integration error tolerance");
      // clang-format on
      if (extraOptions != nullptr) {
        desc.add(*extraOptions);
      }
      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, desc), vm);
      po::notify(vm);

      if (vm.contains("help")) {
        std::cout << desc << std::endl;
        return 0;
      }
    } catch (std::exception& e) {
      std::cerr << "error: " << e.what() << std::endl;
      return 1;
    }

    return std::nullopt;
  }

  std::unique_ptr<MagneticFieldProvider> makeField() const {
    return std::make_unique<ConstantBField>(
        Vector3{0, 0, BzInT * UnitConstants::T});
  }

  template <typename Stepper>
  void run(Stepper stepper, const std::string& name) const {
    run(std::move(stepper), name, [](auto&) {});
  }

  template <typename Stepper, typename Configure>
  void run(Stepper stepper, const std::string& name,
           Configure configure) const {
    using Propagator = Propagator<Stepper>;
    using PropagatorOptions = typename Propagator::template Options<>;
    using Covariance = BoundMatrix;

    // Create a test context
    GeometryContext tgContext = GeometryContext::dangerouslyDefaultConstruct();
    MagneticFieldContext mfContext = MagneticFieldContext();

    ACTS_LOCAL_LOGGER(getDefaultLogger(name, Logging::Level(lvl)));

    // print information about profiling setup
    ACTS_INFO("propagating " << toys << " tracks with pT = " << ptInGeV
                             << "GeV in a " << BzInT << "T B-field");

    Propagator propagator(std::move(stepper));

    PropagatorOptions options(tgContext, mfContext);
    options.pathLimit = maxPathInM * UnitConstants::m;
    options.stepping.maxStepSize = maxStepInMM * UnitConstants::mm;
    options.stepping.stepTolerance = stepTolerance;
    configure(options.stepping);

    Vector4 pos4(0, 0, 0, 0);
    Vector3 dir(1, 0, 0);
    Covariance cov;
    // clang-format off
    cov << 10_mm, 0, 0, 0, 0, 0,
            0, 10_mm, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0,
            0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1_e / 10_GeV, 0,
            0, 0, 0, 0, 0, 0;
    // clang-format on

    std::optional<Covariance> covOpt = std::nullopt;
    if (withCov) {
      covOpt = cov;
    }
    BoundTrackParameters pars = BoundTrackParameters::createCurvilinear(
        pos4, dir, +1 / ptInGeV, covOpt, ParticleHypothesis::pion());

    double totalPathLength = 0;
    std::size_t numSteps = 0;
    std::size_t numAcceptedSteps = 0;
    std::size_t numStepTrials = 0;
    std::size_t numIters = 0;
    const auto propagationBenchResult = microBenchmark(
        [&] {
          auto state = propagator.makeState(options);
          auto initRes = propagator.initialize(state, pars, nullptr);
          if (!initRes.ok()) {
            ACTS_ERROR("initialization failed: " << initRes.error());
            return;
          }
          auto tmp = propagator.propagate(state);
          auto r = propagator.makeResult(state, tmp, options, true).value();
          if (totalPathLength == 0.) {
            ACTS_DEBUG("reached position "
                       << r.endParameters->position(tgContext).transpose()
                       << " in " << r.steps << " steps");
          }
          totalPathLength += r.pathLength;
          numSteps += r.steps;
          numAcceptedSteps += r.statistics.stepping.nSuccessfulSteps;
          numStepTrials += r.statistics.stepping.nAttemptedSteps;
          ++numIters;
        },
        1, toys);

    ACTS_INFO("Execution stats: " << propagationBenchResult);
    ACTS_INFO("average path length = " << totalPathLength / numIters / 1_mm
                                       << "mm");
    ACTS_INFO("average number of steps = " << 1.0 * numSteps / numIters);
    ACTS_INFO("average accepted integration steps = "
              << 1.0 * numAcceptedSteps / numIters);
    ACTS_INFO("step efficiency = " << 1.0 * numAcceptedSteps / numStepTrials);
  }
};

}  // namespace ActsTests

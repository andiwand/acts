// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <concepts>
#include <memory>
#include <stdexcept>

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/variadic/to_seq.hpp>
#include <pybind11/pybind11.h>

namespace ActsExamples {
struct AlgorithmContext;
enum class ProcessCode;
}  // namespace ActsExamples

namespace ActsPython::Concepts {
template <typename T>
concept has_write_method =
    requires(T& t, const ActsExamples::AlgorithmContext& ctx) {
      { t.write(ctx) } -> std::same_as<ActsExamples::ProcessCode>;
    };

/// A config that carries the component's logger.
template <typename C>
concept LoggingConfig = requires(C& c) {
  { c.logger } -> std::same_as<std::shared_ptr<const Acts::Logger>&>;
};

/// Concept for types usable with declareAlgorithm: the config carries the
/// logger and the algorithm takes only the config.
template <typename A>
concept DeclarableAlgorithm =
    requires { typename A::Config; } && LoggingConfig<typename A::Config> &&
    std::constructible_from<A, const typename A::Config&> &&
    requires(const A& a) {
      { a.config() } -> std::same_as<const typename A::Config&>;
    };
}  // namespace ActsPython::Concepts

#define ACTS_PYTHON_MEMBER(name) \
  _binding_instance.def_readwrite(#name, &_struct_type::name)

#define ACTS_PYTHON_STRUCT_BEGIN(object)               \
  {                                                    \
    [[maybe_unused]] auto& _binding_instance = object; \
    using _struct_type = decltype(object)::type;       \
    do {                                               \
    } while (0)

#define ACTS_PYTHON_STRUCT_END() \
  }                              \
  do {                           \
  } while (0)

/// This macro is needed to use the BOOST_PP_SEQ_FOR_EACH loop macro
#define ACTS_PYTHON_MEMBER_LOOP(r, data, elem) ACTS_PYTHON_MEMBER(elem);

/// Macro that accepts a variadic set of member names that are to be registered
/// into an object as read-write fields
#define ACTS_PYTHON_STRUCT(object, ...)                          \
  do {                                                           \
    ACTS_PYTHON_STRUCT_BEGIN(object);                            \
    BOOST_PP_SEQ_FOR_EACH(ACTS_PYTHON_MEMBER_LOOP, _,            \
                          BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
    ACTS_PYTHON_STRUCT_END();                                    \
  } while (0)

/// Level of the logger a config carries.
template <typename Config>
Acts::Logging::Level configLevel(const Config& cfg) {
  if (cfg.logger == nullptr) {
    throw std::runtime_error("config carries no logger");
  }
  return cfg.logger->level();
}

/// Change the level of the logger a config carries, keeping its output policy
/// and its name.
template <typename Config>
void setConfigLevel(Config& cfg, Acts::Logging::Level level) {
  if (cfg.logger == nullptr) {
    throw std::runtime_error("config carries no logger");
  }
  cfg.logger = cfg.logger->clone(std::nullopt, level);
}

/// Bind the constructor of a component: the config carries its logger.
template <typename T, typename C>
void declareComponentInit(C& c) {
  c.def(pybind11::init<const typename T::Config&>(), pybind11::arg("config"));
}

/// Register the logging members on a bound config class. Done here instead of
/// at every call site so `logger` stays out of the member lists that the
/// declare macros take.
template <typename C>
void declareLoggingConfig(C& c) {
  using Config = typename C::type;
  c.def_readwrite("logger", &Config::logger)
      .def_property("level", &configLevel<Config>, &setConfigLevel<Config>);
}

template <ActsPython::Concepts::DeclarableAlgorithm A, typename B>
auto declareAlgorithm(pybind11::module_& m, const char* name) {
  using Config = typename A::Config;
  namespace py = pybind11;
  auto alg = py::class_<A, B, std::shared_ptr<A>>(m, name)
                 .def(py::init<const Config&>(), py::arg("config"))
                 .def_property_readonly("config", &A::config);
  auto c = py::class_<Config>(alg, "Config");
  if constexpr (std::is_default_constructible_v<Config>) {
    c.def(py::init<>());
  }
  declareLoggingConfig(c);
  return std::tuple{alg, c};
}

/// A macro that uses Boost.Preprocessor to create the python binding for and
/// algorithm and the additional config struct.
#define ACTS_PYTHON_DECLARE_ALGORITHM(algorithm, mod, name, ...)          \
  do {                                                                    \
    auto [alg, c] =                                                       \
        declareAlgorithm<algorithm, ActsExamples::IAlgorithm>(mod, name); \
    ACTS_PYTHON_STRUCT(c, __VA_ARGS__);                                   \
  } while (0)

/// Similar as above for writers
#define ACTS_PYTHON_DECLARE_WRITER(writer, mod, name, ...)                  \
  do {                                                                      \
    using Writer = writer;                                                  \
    using Config = Writer::Config;                                          \
    auto w =                                                                \
        py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>( \
            mod, name)                                                      \
            .def_property_readonly("config", &Writer::config);              \
    declareComponentInit<Writer>(w);                                        \
                                                                            \
    constexpr bool has_write_method =                                       \
        ActsPython::Concepts::has_write_method<Writer>;                     \
                                                                            \
    if constexpr (has_write_method) {                                       \
      w.def("write", &Writer::write);                                       \
    }                                                                       \
    auto c = py::class_<Config>(w, "Config").def(py::init<>());             \
    declareLoggingConfig(c);                                                \
    ACTS_PYTHON_STRUCT(c, __VA_ARGS__);                                     \
  } while (0)

/// Similar as above for readers
#define ACTS_PYTHON_DECLARE_READER(reader, mod, name, ...)                  \
  do {                                                                      \
    using Reader = reader;                                                  \
    using Config = Reader::Config;                                          \
    auto r =                                                                \
        py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>( \
            mod, name)                                                      \
            .def_property_readonly("config", &Reader::config);              \
    declareComponentInit<Reader>(r);                                        \
                                                                            \
    auto c = py::class_<Config>(r, "Config").def(py::init<>());             \
    declareLoggingConfig(c);                                                \
    ACTS_PYTHON_STRUCT(c, __VA_ARGS__);                                     \
  } while (0)

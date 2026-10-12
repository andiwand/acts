#!/usr/bin/env python3
"""Repeat the Atlas/SymPy CPU comparison, serially and in shuffled order.

Build ActsBenchmarkAtlasStepper and ActsBenchmarkSympyStepper in Release mode.
For a before/after comparison, save the Atlas executable before rebuilding and
pass it with --baseline. Each executable performs its own two-second warmup.
Results include exact commands, executable hashes, individual samples and the
median across repetitions. Constant-field helix runs have the same physics as
full RK; gradient runs deliberately measure extra work absent from SymPy.
"""

import argparse
import hashlib
import json
import platform
import random
import re
import statistics
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--toys", type=int, default=50000)
    parser.add_argument("--include-zero-field", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.repeats < 1 or args.toys < 1:
        parser.error("--repeats and --toys must be positive")

    atlas = args.build.resolve() / "bin/ActsBenchmarkAtlasStepper"
    sympy = args.build.resolve() / "bin/ActsBenchmarkSympyStepper"
    variants = [
        ("atlas", atlas, []),
        ("atlas_helix", atlas, ["--helix-step", "1000"]),
        ("sympy", sympy, []),
    ]
    if args.baseline:
        variants.insert(0, ("atlas_baseline", args.baseline.resolve(), []))
    cases = []
    for covariance in (False, True):
        modes = variants.copy()
        if covariance:
            modes += [
                ("atlas_gradient", atlas, ["--gradient"]),
                ("atlas_numerical_gradient", atlas, ["--numerical-gradient"]),
            ]
        for name, exe, flags in modes:
            cases.append((name, exe, flags, covariance, 2))
        if args.include_zero_field:
            for name, exe, flags in (
                ("atlas", atlas, []),
                ("atlas_straight", atlas, ["--field-off"]),
                ("sympy", sympy, []),
            ):
                cases.append((name, exe, flags, covariance, 0))

    cache = (args.build / "CMakeCache.txt").read_text()
    build_settings = {}
    for key in (
        "CMAKE_BUILD_TYPE",
        "CMAKE_CXX_COMPILER",
        "CMAKE_CXX_FLAGS",
        "CMAKE_CXX_FLAGS_RELEASE",
        "Eigen3_DIR",
    ):
        match = re.search(rf"^{key}:[^=]+=(.*)$", cache, re.MULTILINE)
        if match:
            build_settings[key] = match[1]
    result = {
        "platform": platform.platform(),
        "build_settings": build_settings,
        "machine": platform.machine(),
        "toys": args.toys,
        "warmup_seconds_per_process": 2,
        "repeats": args.repeats,
        "executables": {
            str(exe): hashlib.sha256(exe.read_bytes()).hexdigest()
            for _, exe, _, _, _ in cases
        },
        "samples": [],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rng = random.Random(42)
    for repeat in range(args.repeats):
        rng.shuffle(cases)
        for name, exe, flags, covariance, field in cases:
            command = [
                str(exe),
                "--toys",
                str(args.toys),
                "--cov",
                str(covariance).lower(),
                "--pT",
                "1",
                "--path",
                "5",
                "--B",
                str(field),
                *flags,
            ]
            output = subprocess.check_output(
                command, text=True, stderr=subprocess.STDOUT
            )
            if "ERROR" in output or "FATAL" in output:
                raise RuntimeError(output)
            ns = float(re.search(r"([\d.]+)\+/-[\d.]+ns per iteration", output)[1])
            steps = float(re.search(r"average number of steps = ([\d.]+)", output)[1])
            accepted = re.search(
                r"average accepted integration steps = ([\d.]+)", output
            )
            path = float(re.search(r"average path length = ([\d.]+)mm", output)[1])
            if abs(path - 5000) > 1e-3:
                raise RuntimeError(f"Unexpected path length: {output}")
            result["samples"].append(
                {
                    "variant": name,
                    "covariance": covariance,
                    "field_tesla": field,
                    "repeat": repeat,
                    "ns_per_track": ns,
                    "propagator_iterations": steps,
                    "accepted_steps": float(accepted[1]) if accepted else None,
                    "command": command,
                    "output": output,
                }
            )
            args.output.write_text(json.dumps(result, indent=2) + "\n")
            print(
                f"{repeat + 1}/{args.repeats} {name:25} cov={covariance} B={field}: {ns:.0f} ns, {steps:g} propagator iterations",
                flush=True,
            )

    result["summary"] = []
    for name, _, _, covariance, field in sorted(cases):
        values = [
            s["ns_per_track"]
            for s in result["samples"]
            if (s["variant"], s["covariance"], s["field_tesla"])
            == (name, covariance, field)
        ]
        result["summary"].append(
            {
                "variant": name,
                "covariance": covariance,
                "field_tesla": field,
                "median_ns": statistics.median(values),
                "min_ns": min(values),
                "max_ns": max(values),
            }
        )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()

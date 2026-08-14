import argparse
import contextlib
import sys

import sympy as sym
from sympy import MatrixSymbol

from codegen.sympy_common import (
    NamedExpr,
    explicit,
    name_expr,
    find_by_name,
    cxx_printer,
    my_expression_print,
)

step_path_derivatives = explicit(
    MatrixSymbol("step_path_derivatives", 8, 1)
).as_mutable()
step_path_derivatives[7, 0] = 0  # qop

surface_path_derivatives = explicit(
    MatrixSymbol("surface_path_derivatives", 1, 8)
).as_mutable()
surface_path_derivatives[0, 3] = 0
surface_path_derivatives[0, 7] = 0

# M is the bound-to-free jacobian already transported to the current point,
# i.e. what used to be J_t * J_bf.  Its sparsity is everything a sequence of
# steps can produce from a start surface: loc0 and loc1 stay position only,
# phi and theta pick up a direction part (and a position part already at the
# start surface, for a line surface), q/p picks up position, time and
# direction parts, and the time column stays exactly e_time.
J_bf = explicit(MatrixSymbol("M", 8, 6)).as_mutable()
tmp = sym.zeros(8, 6)
tmp[0:3, 0:2] = J_bf[0:3, 0:2]
tmp[0:3, 2:5] = J_bf[0:3, 2:5]
tmp[3, 4] = J_bf[3, 4]
tmp[4:7, 2:5] = J_bf[4:7, 2:5]
tmp[7, 4] = J_bf[7, 4]
tmp[3, 5] = 1
J_bf = tmp

J_fb = explicit(MatrixSymbol("J_fb", 6, 8)).as_mutable()
tmp = sym.zeros(6, 8)
tmp[0:2, 0:3] = J_fb[0:2, 0:3]
tmp[2:4, 4:7] = J_fb[2:4, 4:7]
tmp[5, 3] = 1
tmp[4, 7] = 1
J_fb = tmp


def full_transport_jacobian_generic() -> list[NamedExpr]:
    J_full = name_expr(
        "J_full",
        J_fb * (sym.eye(8) + step_path_derivatives * surface_path_derivatives) * J_bf,
    )

    return [J_full]


def full_transport_jacobian_curvilinear(direction: MatrixSymbol) -> list[NamedExpr]:
    surface_path_derivatives = explicit(
        MatrixSymbol("surface_path_derivatives", 1, 8)
    ).as_mutable()
    surface_path_derivatives[0, 0:3] = -direction.as_explicit().transpose()
    surface_path_derivatives[0, 3:8] = sym.zeros(1, 5)

    J_full = name_expr(
        "J_full",
        J_fb * (sym.eye(8) + step_path_derivatives * surface_path_derivatives) * J_bf,
    )

    return [J_full]


def my_full_transport_jacobian_generic_function_print(name_exprs, run_cse=True):
    printer = cxx_printer
    outputs = [find_by_name(name_exprs, name)[0] for name in ["J_full"]]

    lines = []

    head = (
        "template <typename T> void boundToBoundTransportJacobianImpl("
        "std::span<const T, 48> J_fb, std::span<const T, 48> M,"
        " std::span<const T, 8> step_path_derivatives,"
        " std::span<const T, 8> surface_path_derivatives,"
        " std::span<T, 36> J_full) {"
    )
    lines.append(head)

    code = my_expression_print(
        printer,
        name_exprs,
        outputs,
        run_cse=run_cse,
    )
    lines.extend([f"  {l}" for l in code.split("\n")])

    lines.append("}")

    return "\n".join(lines)


def my_full_transport_jacobian_curvilinear_function_print(name_exprs, run_cse=True):
    printer = cxx_printer
    outputs = [find_by_name(name_exprs, name)[0] for name in ["J_full"]]

    lines = []

    head = (
        "template <typename T> void boundToCurvilinearTransportJacobianImpl("
        "std::span<const T, 48> J_fb, std::span<const T, 48> M,"
        " std::span<const T, 8> step_path_derivatives,"
        " std::span<const T, 3> dir, std::span<T, 36> J_full) {"
    )
    lines.append(head)

    code = my_expression_print(
        printer,
        name_exprs,
        outputs,
        run_cse=run_cse,
    )
    lines.extend([f"  {l}" for l in code.split("\n")])

    lines.append("}")

    return "\n".join(lines)


def check_curvilinear_is_generic_specialised() -> None:
    """Assert the curvilinear jacobian is the generic one at a curvilinear surface.

    The two are printed as separate functions, so nothing otherwise stops them
    drifting apart. A curvilinear surface is just one whose path derivatives
    are -direction for the position part and zero elsewhere, so substituting
    that into the generic form has to reproduce the specialised one exactly.

    Raises AssertionError if they disagree.
    """
    direction = MatrixSymbol("dir", 3, 1)
    generic = full_transport_jacobian_generic()[0].expr
    curvilinear = full_transport_jacobian_curvilinear(direction)[0].expr

    spd = explicit(MatrixSymbol("surface_path_derivatives", 1, 8))
    at_curvilinear = {spd[0, i]: -explicit(direction)[i, 0] for i in range(3)}
    at_curvilinear.update({spd[0, i]: 0 for i in range(3, 8)})

    diff = sym.expand(generic.subs(at_curvilinear) - curvilinear)
    if any(e != 0 for e in diff):
        bad = [
            (i, j)
            for i in range(diff.rows)
            for j in range(diff.cols)
            if diff[i, j] != 0
        ]
        raise AssertionError(
            f"curvilinear jacobian is not the generic one specialised, at {bad}"
        )


def check_sparsity_covariance_transport_assumes() -> None:
    """Assert the shape `transportCovarianceToBoundImpl` is built around.

    That function masks its jacobian down to a diagonal-only q/p row and an
    identity time column, and only ever multiplies the rest. If this one ever
    produced a jacobian outside that shape the covariance transport would
    silently drop the difference, so the assumption is checked where it is
    produced rather than where it is consumed.

    Raises AssertionError if the jacobian reaches outside that shape.
    """
    J = sym.expand(full_transport_jacobian_generic()[0].expr)
    qop, time = 4, 5
    # q/p depends on nothing but itself -- but its own diagonal is free, since
    # energy loss moves it away from one -- and nothing depends on time.
    leaks = [(qop, j) for j in range(6) if j != qop and J[qop, j] != 0]
    leaks += [(i, time) for i in range(6) if J[i, time] != (1 if i == time else 0)]
    if leaks:
        raise AssertionError(
            "bound-to-bound jacobian is not of the shape the covariance "
            f"transport assumes; live entries outside it: {leaks}"
        )


HEADER = """// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Note: This file is generated by generate_sympy_jac.py
//       Do not modify it manually.

#pragma once

#include <cmath>
#include <span>
"""


def main(argv: list[str]) -> None:
    """Generate the transport jacobians.

    @param argv is the command line, argv[0] being the program name
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output",
        nargs="?",
        help="file to write the generated jacobians to; stdout if omitted",
    )
    parser.add_argument(
        "--no-check",
        action="store_true",
        help="skip the symbolic assertions",
    )
    args = parser.parse_args(argv[1:])

    if not args.no_check:
        # If one of these fires the generated code would be wrong, so they
        # guard the generator rather than a test that might not be run.
        check_curvilinear_is_generic_specialised()
        check_sparsity_covariance_transport_assumes()

    with (
        open(args.output, "w") if args.output else contextlib.nullcontext(sys.stdout)
    ) as out:
        out.write(HEADER)
        out.write(
            my_full_transport_jacobian_generic_function_print(
                full_transport_jacobian_generic(), run_cse=True
            )
        )
        out.write("\n")
        out.write(
            my_full_transport_jacobian_curvilinear_function_print(
                full_transport_jacobian_curvilinear(MatrixSymbol("dir", 3, 1)),
                run_cse=True,
            )
        )
        out.write("\n")


if __name__ == "__main__":
    main(sys.argv)

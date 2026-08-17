import argparse
import contextlib
import functools
import random
import sys
from dataclasses import dataclass

import numpy as np

import sympy as sym
from sympy import Symbol, Matrix, ImmutableMatrix, MatrixSymbol
from sympy.codegen.ast import Assignment

from codegen.sympy_common import (
    Derivation,
    explicit,
    make_vector,
    NamedExpr,
    name_expr,
    find_by_name,
    my_subs,
    cxx_printer,
    my_expression_print,
)

# q/p
l = Symbol("l", real=True)

# step length
h = Symbol("h", real=True)

# path length
s = Symbol("s", real=True)

# position
p = make_vector("p", 3, real=True)

# direction
d = make_vector("d", 3, real=True)

# time
t = Symbol("t", real=True)

# mass
m = Symbol("m", real=True, positive=True)

# absolute momentum
p_abs = Symbol("p_abs", real=True, positive=True)

# energy loss per distance
g = Symbol("g", real=True)

# charge
q = Symbol("q", real=True)

# magnetic field
B = make_vector("B", 3, real=True)

# specific magnetic field values
B1 = make_vector("B1", 3, real=True)
B2 = make_vector("B2", 3, real=True)
B3 = make_vector("B3", 3, real=True)

# specific energy loss per distance values
g1 = Symbol("g1", real=True)
g2 = Symbol("g2", real=True)
g3 = Symbol("g3", real=True)
g4 = Symbol("g4", real=True)


@dataclass
class Rk4Step:
    """Everything one plain RK4 step of a second order ODE produces.

    A record rather than nested tuples: the derivatives are as much a result
    as the state is, and positional unpacking three levels deep hid which was
    which.
    """

    #: the integrated state after the step
    new_y: NamedExpr
    #: its derivative after the step
    new_ydot: NamedExpr
    #: the four stage slopes
    k: tuple[NamedExpr, NamedExpr, NamedExpr, NamedExpr]
    #: d(new_y) / d(y, ydot)
    dy: NamedExpr
    #: d(new_ydot) / d(y, ydot)
    dydot: NamedExpr
    #: d(k_i) / d(y, ydot)
    dk: tuple[NamedExpr, NamedExpr, NamedExpr, NamedExpr]
    #: the stage points, kept because the kernels sample the field at them
    x2: NamedExpr
    y2: NamedExpr
    ydot2: NamedExpr
    ydot3: NamedExpr
    x3: NamedExpr
    y3: NamedExpr
    ydot4: NamedExpr

    def named_exprs(self) -> list[NamedExpr]:
        """Every named quantity the step introduced, for resolving.

        @return the stage slopes, stage points and results, in no order
        """
        # dependency order: `resolve` walks this backwards, so a name has to
        # appear before the ones whose bodies mention it
        return [
            *self.k,
            self.x2,
            self.y2,
            self.ydot2,
            self.ydot3,
            self.x3,
            self.y3,
            self.ydot4,
            self.new_y,
            self.new_ydot,
            *self.dk,
            self.dy,
            self.dydot,
        ]


def rk4_subexpr(f, x, y, ydot, h):
    k1 = name_expr("k1", f(1, x, y, ydot))
    x2 = name_expr("x2", x + h / 2)
    y2 = name_expr("y2", y + h / 2 * ydot + h**2 / 8 * k1.name)
    ydot2 = name_expr("ydot2", ydot + h / 2 * k1.name)

    k2 = name_expr("k2", f(2, x2.expr, y2.expr.as_explicit(), ydot2.expr.as_explicit()))
    ydot3 = name_expr("ydot3", ydot + h / 2 * k2.name)

    k3 = name_expr("k3", f(3, x2.expr, y2.expr.as_explicit(), ydot3.expr.as_explicit()))
    x3 = name_expr("x3", x + h)
    y3 = name_expr("y3", y + h * ydot + h**2 / 2 * k3.name)
    ydot4 = name_expr("ydot4", ydot + h * k3.name)

    k4 = name_expr("k4", f(4, x3.expr, y3.expr.as_explicit(), ydot4.expr.as_explicit()))

    new_y = name_expr("new_y", y + h * ydot + h**2 / 6 * (k1.name + k2.name + k3.name))
    new_ydot = name_expr(
        "new_ydot", ydot + h / 6 * (k1.name + 2 * (k2.name + k3.name) + k4.name)
    )

    dk1dyydot = name_expr("dk1dyydot", k1.expr.jacobian([y, ydot]))
    dk2dyydot = name_expr(
        "dk2dyydot",
        k2.expr.jacobian([y, ydot]) + k2.expr.jacobian(k1.name) * dk1dyydot.name,
    )
    dk3dyydot = name_expr(
        "dk3dyydot",
        k3.expr.jacobian([y, ydot]) + k3.expr.jacobian(k2.name) * dk2dyydot.name,
    )
    dk4dyydot = name_expr(
        "dk4dyydot",
        k4.expr.jacobian([y, ydot]) + k4.expr.jacobian(k3.name) * dk3dyydot.name,
    )

    dydyydot = name_expr(
        "dydyydot",
        new_y.expr.as_explicit().jacobian([y, ydot])
        + new_y.expr.as_explicit().jacobian(k1.name) * dk1dyydot.name
        + new_y.expr.as_explicit().jacobian(k2.name) * dk2dyydot.name
        + new_y.expr.as_explicit().jacobian(k3.name) * dk3dyydot.name,
    )
    dydotdyydot = name_expr(
        "dydotdyydot",
        new_ydot.expr.as_explicit().jacobian([y, ydot])
        + new_ydot.expr.as_explicit().jacobian(k1.name) * dk1dyydot.name
        + new_ydot.expr.as_explicit().jacobian(k2.name) * dk2dyydot.name
        + new_ydot.expr.as_explicit().jacobian(k3.name) * dk3dyydot.name
        + new_ydot.expr.as_explicit().jacobian(k4.name) * dk4dyydot.name,
    )

    return Rk4Step(
        new_y=new_y,
        new_ydot=new_ydot,
        k=(k1, k2, k3, k4),
        dy=dydyydot,
        dydot=dydotdyydot,
        dk=(dk1dyydot, dk2dyydot, dk3dyydot, dk4dyydot),
        x2=x2,
        y2=y2,
        ydot2=ydot2,
        ydot3=ydot3,
        x3=x3,
        y3=y3,
        ydot4=ydot4,
    )


SEED_D = Matrix.hstack(sym.eye(3), sym.zeros(3, 1))


@dataclass
class AtlasStages:
    """The named quantities of one ATLAS-form RK4 step.

    Named explicitly rather than collected into a kwargs bag, now that there
    is also a naive form to tell it apart from.
    """

    #: h/3, the position update weight
    S3: NamedExpr
    #: the half-step scaled field at each of the three sample points
    H0: NamedExpr
    H1: NamedExpr
    H2: NamedExpr
    #: the stage slopes, already scaled by h/2
    A0: NamedExpr
    A1: NamedExpr
    A2: NamedExpr
    A3: NamedExpr
    A4: NamedExpr
    A5: NamedExpr
    A6: NamedExpr
    #: dt/ds
    dtds: NamedExpr


def atlas_rk4_stages(deriv: Derivation, taylor_norm: bool) -> AtlasStages:
    """Build the value path of an ATLAS-form RK4 vacuum step.

    ATLAS never evaluates the plain slopes k_i. It carries the half-step
    scaled field H = (h*l/2) * B, so every stage slope comes out as (h/2)*k_i
    directly and no further scaling by h or l appears anywhere in the
    recursion.
    """
    # (h*l/2) is ATLAS' PS2.
    hl2 = deriv.add("hl2", h * l / 2)
    S3 = deriv.add("S3", h / 3)
    S4 = deriv.add("S4", h / 4)

    H0 = deriv.add("H0", hl2.name * B1)
    A0 = deriv.add("A0", d.cross(explicit(H0.name)))  # h/2 * k1
    A2 = deriv.add("A2", explicit(A0.name) + d)  # d + h/2 * k1
    A1 = deriv.add("A1", explicit(A2.name) + d)  # 2d + h/2 * k1
    deriv.add("p2", p + S4.name * explicit(A1.name))

    H1 = deriv.add("H1", hl2.name * B2)
    A3 = deriv.add("A3", d + explicit(A2.name).cross(explicit(H1.name)))  # d + h/2 * k2
    A4 = deriv.add("A4", d + explicit(A3.name).cross(explicit(H1.name)))  # d + h/2 * k3
    A5 = deriv.add("A5", 2 * explicit(A4.name) - d)  # d + h * k3
    deriv.add("p3", p + h * explicit(A4.name))

    H2 = deriv.add("H2", hl2.name * B3)
    A6 = deriv.add("A6", explicit(A5.name).cross(explicit(H2.name)))  # h/2 * k4

    # (A1+A6)-(A3+A4) is h/2 * (k1-k2-k3+k4), hence the leading 2*|h|.
    err_vec = (explicit(A1.name) + explicit(A6.name)) - (
        explicit(A3.name) + explicit(A4.name)
    )
    deriv.add("err", 2 * sym.Abs(h) * err_vec.norm(1))

    deriv.add(
        "new_p",
        p + S3.name * (explicit(A2.name) + explicit(A3.name) + explicit(A4.name)),
    )

    # An = 3 * (d + h/6*(k1 + 2k2 + 2k3 + k4)), three times the unnormalised
    # new direction.
    An = deriv.add(
        "An",
        2 * explicit(A3.name)
        + (explicit(A0.name) + explicit(A5.name) + explicit(A6.name)),
    )

    if taylor_norm:
        # ATLAS replaces 1/|An| by its second order Taylor expansion around
        # |An| = 3, which needs neither a square root nor a division. RK4
        # keeps |An| within ~1e-7 of 3, so the O(u^3) truncation error is far
        # below double precision.
        Dv = deriv.add(
            "Dv", (An.name[0] ** 2 + An.name[1] ** 2) + (An.name[2] ** 2 - 9)
        )
        Dfac = deriv.add(
            "Dfac",
            sym.Rational(1, 3) - sym.Rational(1, 648) * Dv.name * (12 - Dv.name),
        )
        new_d = deriv.add("new_d", Dfac.name * explicit(An.name))
    else:
        inv_norm = deriv.add("inv_norm", 1 / explicit(An.name).norm())
        new_d = deriv.add("new_d", inv_norm.name * explicit(An.name))

    dtds = deriv.add("dtds", sym.sqrt(1 + m**2 / p_abs**2))
    deriv.add("new_t", t + h * dtds.name)

    # B3 is evaluated at p + h*A4, which differs from the step end point by
    # O(h^3). Handing it back lets the caller reuse it as the next step's
    # first field sample and save one lookup per step, as ATLAS does.
    deriv.add("new_B", B3)

    Sl = deriv.add("Sl", 2 / h)  # undoes the h/2 scaling of A6 -> k4
    deriv.add(
        "path_derivatives",
        Matrix.vstack(
            explicit(new_d.name),
            Matrix([dtds.name]),
            Sl.name * explicit(A6.name),
            Matrix([0]),
        ),
    )

    return AtlasStages(
        S3=S3,
        H0=H0,
        H1=H1,
        H2=H2,
        A0=A0,
        A1=A1,
        A2=A2,
        A3=A3,
        A4=A4,
        A5=A5,
        A6=A6,
        dtds=dtds,
    )


def field_contrib(
    deriv: Derivation,
    what: str,
    stage: NamedExpr,
    H: NamedExpr,
    same_as: Matrix,
    seed: Matrix,
) -> Matrix:
    """The term a tangent picks up from H's own dependence on l.

    H is linear in l, so `l * dH/dl == H` exactly; with the tangent's l part
    scaled by l, this term is the H-linear part of the stage, which is already
    available as a named quantity. `same_as` encodes that identity and is
    verified against the plain chain rule before use -- it is what keeps the
    recursion free of any extra cross product.
    """
    contrib = stage.expr.jacobian(H.name) * seed
    deriv.check_same(what, contrib[:, seed.cols - 1], same_as)
    return Matrix.hstack(contrib[:, 0 : seed.cols - 1], same_as)


B2F_LIVE_COLUMNS = [
    ("new_Mp", [(i, 2) for i in (0, 1, 2, 4, 5, 6)]),
    ("new_Mt", [(i, 3) for i in (0, 1, 2, 4, 5, 6)]),
    ("new_Mq", [(i, 4) for i in (0, 1, 2, 3, 4, 5, 6)]),
]

B2F_LIVE = [e for _, g in B2F_LIVE_COLUMNS for e in g]

# A dense step additionally changes q/p itself, so the q/p row of the q/p
# column moves too.
B2F_DENSE_LIVE = (
    [(i, 2) for i in (0, 1, 2, 4, 5, 6)]
    + [(i, 3) for i in (0, 1, 2, 4, 5, 6)]
    + [(i, 4) for i in (0, 1, 2, 3, 4, 5, 6, 7)]
)

# Rows that are structurally zero on input, for the columns that are live at
# all: neither a phi nor a theta perturbation changes the time or the q/p
# coordinate, in vacuum or in matter.
B2F_ZERO_ROWS = {2: (3, 7), 3: (3, 7)}

# bound parameter count -- M is 8 x B2F_COLS, stored column major
B2F_COLS = 6

# The q/p column of M is not stored as the plain bound-to-free jacobian. Its
# free rows are multiplied by q/p and divided by the column's own q/p row:
#
#     M[i, 4] = l * dFree_i/dl_0 / (dl/dl_0)     for i in 0..6
#     M[7, 4] =     dl/dl_0                      unchanged
#
# Both factors are storage convention only, and together they are what makes
# this column cost exactly what the ATLAS stepper's pVector[40] block costs:
# the l makes the term each RK stage picks up from H's own l dependence equal
# to the stage itself instead of a fresh cross product, and the normalisation
# makes the coefficient of that term one instead of a load of the q/p row.
# In vacuum both factors are constant across a step, so the invariant holds by
# itself; rk4_dense changes both and restores it explicitly. The stepper
# converts to and from the plain jacobian where the covariance engine wants it.
B2F_QOP_COLUMN = 4


def b2f_step_update(
    D: Matrix, live: list[tuple[int, int]], l_in=None, l_out=None
) -> Matrix:
    """Apply a free-to-free step jacobian D to the bound-to-free jacobian M.

    Only the live columns are touched, and the structurally zero rows are
    dropped from the input so the products stay sparse. This is the generic
    fallback; the vacuum kernel instead folds the same operation into the RK
    recursion, which is cheaper still because it never builds D.

    The q/p column is stored scaled (see B2F_QOP_SCALING). A step that
    changes q/p changes both factors of that scaling, so it is undone on the
    way in, using the q/p the column was scaled with, and redone on the way
    out with the new one. Pass l_in and l_out to do that; leaving them None
    treats every column as plain.
    """
    M = MatrixSymbol("M", 8, 6)
    qop_col = B2F_QOP_COLUMN
    out = []
    for c in sorted({col for _, col in live}):
        zero = B2F_ZERO_ROWS.get(c, ())
        v = Matrix([[0 if i in zero else M[i, c]] for i in range(8)])
        scaled = c == qop_col and l_in is not None
        if scaled:
            # stored = l_in * plain / plain[7], so plain = stored * stored[7] / l_in
            f = v[7, 0] / l_in
            v = Matrix([[v[i, 0] * f] if i < 7 else [v[i, 0]] for i in range(8)])
        new_v = D * v
        if scaled:
            g = l_out / new_v[7, 0]
            new_v = Matrix(
                [[new_v[i, 0] * g] if i < 7 else [new_v[i, 0]] for i in range(8)]
            )
        out.extend([new_v[i, 0]] for i, col in live if col == c)
    return Matrix(out)


def rk4_vacuum_b2f_atlasexpr(taylor_norm: bool = False) -> list[NamedExpr]:
    """ATLAS-form RK4 vacuum step transporting the bound-to-free jacobian.

    The free-to-free step jacobian D is never built. Instead the columns of
    the bound-to-free jacobian M are pushed through the very same RK recursion
    as the track state itself -- which is why the code below mirrors the value
    path line for line, exactly as ATLAS' d2A/d3A/d4A block mirrors its A0..A6
    block. That removes both the fourth tangent direction (a bound
    parametrisation has two direction degrees of freedom, a free-to-free D
    needs three) and the 8x8 composition.
    """
    deriv = Derivation()
    st = atlas_rk4_stages(deriv, taylor_norm)
    H0, H1, H2 = st.H0, st.H1, st.H2
    A0, A3, A4, A6 = st.A0, st.A3, st.A4, st.A6
    S3, dtds = st.S3, st.dtds

    M = MatrixSymbol("M", 8, 6)

    def col(c, rows):
        return Matrix([[M[i, c]] for i in rows])

    def tangent(name, stage, contribs, field=None):
        T = sym.zeros(3, 1) if field is None else field
        for var, Tvar in contribs:
            T = T + stage.expr.jacobian(var) * Tvar
        return deriv.add(name, T)

    def propagate(tag, c, with_field):
        """Push column `c` of M through the step.

        `with_field` is False for a column with no q/p component, which then
        picks up nothing from H's own l dependence.

        The q/p column is stored in the scaled form described at
        B2F_QOP_SCALING: its free rows carry a factor l and are divided by the
        column's own q/p row. Both factors are pure storage convention, but
        together they are what makes this column cost exactly what the ATLAS
        stepper's pVector[40] block costs. The l turns the term each stage
        picks up from H's l dependence into the stage itself rather than a
        fresh cross product, and the normalisation makes the coefficient of
        that term a literal one instead of a load of the q/p row. Nothing is
        assumed about the q/p row -- rk4_dense restores the invariant whenever
        it changes it.
        """
        seed = col(c, (4, 5, 6))
        pos_fac, dir_fac = S3.name, third.name
        if not with_field:
            fields = [None] * 4
        else:
            one = sym.Integer(1)
            fields = [
                field_contrib(
                    deriv,
                    f"{tag}0",
                    A0,
                    H0,
                    explicit(A0.name) * one,
                    explicit(H0.name) * one,
                ),
                field_contrib(
                    deriv,
                    f"{tag}3",
                    A3,
                    H1,
                    (explicit(A3.name) - d) * one,
                    explicit(H1.name) * one,
                ),
                field_contrib(
                    deriv,
                    f"{tag}4",
                    A4,
                    H1,
                    (explicit(A4.name) - d) * one,
                    explicit(H1.name) * one,
                ),
                field_contrib(
                    deriv,
                    f"{tag}6",
                    A6,
                    H2,
                    explicit(A6.name) * one,
                    explicit(H2.name) * one,
                ),
            ]

        v0 = tangent(f"{tag}0", A0, [(d, seed)], fields[0])
        v2 = deriv.add(f"{tag}2", explicit(v0.name) + seed)
        v3 = tangent(
            f"{tag}3", A3, [(d, seed), (st.A2.name, explicit(v2.name))], fields[1]
        )
        v4 = tangent(
            f"{tag}4", A4, [(d, seed), (A3.name, explicit(v3.name))], fields[2]
        )
        v5 = deriv.add(f"{tag}5", 2 * explicit(v4.name) - seed)
        v6 = tangent(f"{tag}6", A6, [(st.A5.name, explicit(v5.name))], fields[3])

        new_pos = col(c, (0, 1, 2)) + pos_fac * (
            explicit(v2.name) + explicit(v3.name) + explicit(v4.name)
        )
        new_dir = dir_fac * (
            explicit(v0.name)
            + 2 * explicit(v3.name)
            + explicit(v5.name)
            + explicit(v6.name)
        )
        return new_pos, new_dir

    # NOTE: dir_fac must be a Symbol, never a Number. sympy distributes a
    # Number over an Add, so Rational(1,3) * (v0 + 2*v3 + v5 + v6) emits four
    # multiplications per component where ATLAS' equivalent
    # ((d2A0 + 2*d2A3) + (d2A5 + d2A6)) * (1./3.) uses two -- eighteen wasted
    # multiplies once the three columns are counted. A Symbol does not
    # distribute, so the scaling stays factored and no temporary is needed.

    third = deriv.add("third", sym.Rational(1, 3))

    phi_pos, phi_dir = propagate("Tp", 2, False)
    the_pos, the_dir = propagate("Tt", 3, False)
    qop_pos, qop_dir = propagate("Tq", 4, True)

    # dt/ds depends on q/p only, so the time row moves for the q/p column alone.
    # The extra l is the scaling the stored column carries; the division by the
    # q/p row that goes with it has already cancelled, since the plain entry is
    # dtdl times that row.
    #
    # Written over p_abs rather than l: `l * d(h*dtds)/dl` is
    # `h m^2 / (p_abs^2 dtds)`, and the `h m^2 l^2 / dtds` it used to be is
    # only the same thing when the charge is one, since p_abs = |q| / |l|.
    # Same operation count, and no need for `q` as an input.
    dtdl = deriv.add("dtdl", h * m**2 / (p_abs**2 * dtds.name))
    new_time = M[3, 4] + dtdl.name

    deriv.add("new_Mp", Matrix.vstack(phi_pos, phi_dir))
    deriv.add("new_Mt", Matrix.vstack(the_pos, the_dir))
    deriv.add("new_Mq", Matrix.vstack(qop_pos, Matrix([new_time]), qop_dir))

    return deriv.name_exprs


# --- the naive reference: equations of motion, plain RK4, chain-rule
# derivatives. Both kernels are the same derivation with a different right
# hand side, so there is one place where the integrator and the differentiation
# live and two places where the physics does.

# The integrated state is split the way a second order ODE wants it: `y` is
# what is integrated, `ydot` its derivative. Laid out to match FreeVector, so
# the generated code can take the parameter vector as it stands:
#
#   y    = (pos, time)          -> free indices 0,1,2 and 3
#   ydot = (dir, dt/ds, q/p)    -> free indices 4,5,6 and 7, plus dt/ds
#
# dt/ds rides in `ydot` because the dense right hand side evolves it, but it is
# NOT independent: dtds = sqrt(1 + m^2 l^2 / q^2). Every derivative with
# respect to q/p therefore has to pick up the path through it, which is what
# NAIVE_QOP_SEED carries. Getting that wrong is exactly the bug that made the
# dense kernel disagree with the vacuum one on d(time)/d(q/p).


def eom_vacuum(i, x, y, ydot):
    """d(ydot)/ds in vacuum: the field bends the direction, nothing else moves.

    @param i is the Runge-Kutta stage, 1 to 4
    @param x is the path length (unused, the field is sampled per stage)
    @param y is the integrated state (position, time)
    @param ydot is its derivative (direction, dt/ds, q/p)
    @return d(ydot)/ds
    """
    del x, y
    B = [B1, B2, B2, B3][i - 1]
    direction = ydot[0:3, 0]
    qop = ydot[4, 0]
    return Matrix.vstack(
        direction.cross(qop * B),
        Matrix([0]),  # no energy loss, so dt/ds is constant
        Matrix([0]),  # ... and so is q/p
    )


def eom_dense(i, x, y, ydot):
    """d(ydot)/ds through material: the same, plus energy loss on q/p.

    @param i is the Runge-Kutta stage, 1 to 4
    @param x is the path length (unused, the field is sampled per stage)
    @param y is the integrated state (position, time)
    @param ydot is its derivative (direction, dt/ds, q/p)
    @return d(ydot)/ds
    """
    del x, y
    B = [B1, B2, B2, B3][i - 1]
    gloss = [g1, g2, g3, g4][i - 1]
    direction = ydot[0:3, 0]
    dtds = ydot[3, 0]
    qop = ydot[4, 0]
    return Matrix.vstack(
        direction.cross(qop * B),
        Matrix([gloss * m**2 * qop**3 / q**3]),
        Matrix([dtds * qop**2 * gloss / q]),
    )


def propagate_tangent(
    deriv: Derivation, tag: str, step: Rk4Step, y, ydot, u: Matrix
) -> tuple[Matrix, Matrix]:
    """Push one tangent vector through the Runge-Kutta recursion.

    The same chain rule `rk4_subexpr` applies to build the full derivative
    blocks, but contracted with a vector at every stage instead of carried as
    a matrix. A column of the bound-to-free jacobian is exactly such a vector,
    so this transports it without ever forming the 8x8 step jacobian -- which
    is what the vacuum kernel does by hand through the ATLAS stages, done
    generically.

    @param deriv collects the named stage tangents
    @param tag prefixes their names, one set per column
    @param step is the value-path Runge-Kutta step to differentiate
    @param y is the integrated state the step was built over
    @param ydot is its derivative
    @param u is the seed tangent, in (y, ydot) space
    @return the tangents of the new state and of its new derivative
    """
    k = step.k
    ts = []
    for i in range(4):
        contrib = k[i].expr.jacobian([y, ydot]) * u
        if i > 0:
            contrib = contrib + k[i].expr.jacobian(k[i - 1].name) * ts[i - 1]
        ts.append(explicit(deriv.add(f"{tag}k{i + 1}", contrib).name))

    def combine(out):
        total = out.expr.as_explicit().jacobian([y, ydot]) * u
        for i in range(4):
            block = out.expr.as_explicit().jacobian(k[i].name)
            if any(e != 0 for e in block):
                total = total + block * ts[i]
        return total

    return combine(step.new_y), combine(step.new_ydot)


def rk4_dense_tunedexpr() -> list[NamedExpr]:
    big_l = Symbol("big_l", real=True)
    dtds = name_expr("dtds", sym.sqrt(1 + m**2 / p_abs**2))

    step = rk4_subexpr(
        eom_dense,
        s,
        Matrix.vstack(p, Matrix([t, big_l])),
        Matrix.vstack(d, Matrix([dtds.name, l])),
        h,
    )

    k1, k2, k3, k4 = step.k
    new_y, new_ydot = step.new_y, step.new_ydot

    p2 = name_expr("p2", step.y2.expr[0:3, 0])
    l2 = name_expr("l2", step.ydot2.expr[4, 0])
    p3 = name_expr("p3", step.y3.expr[0:3, 0])
    l3 = name_expr("l3", step.ydot3.expr[4, 0])
    l4 = name_expr("l4", step.ydot4.expr[4, 0])
    new_p = name_expr("new_p", new_y.expr[0:3, 0])
    new_t = name_expr("new_t", new_y.expr[3, 0])
    new_d_tmp = name_expr("new_d_tmp", new_ydot.expr[0:3, 0])
    new_d = name_expr("new_d", new_d_tmp.name / new_d_tmp.name.as_explicit().norm())
    new_l = name_expr("new_l", new_ydot.expr[4, 0])
    err = name_expr(
        "err",
        h**2
        * (k1.name[0:3, 0] - k2.name[0:3, 0] - k3.name[0:3, 0] + k4.name[0:3, 0])
        .as_explicit()
        .norm(1),
    )

    new_B = name_expr("new_B", B3)

    path_derivatives = name_expr("path_derivatives", sym.zeros(8, 1))
    path_derivatives.expr[0:3, 0] = new_d.name.as_explicit()
    path_derivatives.expr[3, 0] = new_ydot.name[3, 0]
    path_derivatives.expr[4:7, 0] = k4.name[0:3, 0].as_explicit()
    path_derivatives.expr[7, 0] = new_ydot.name[4, 0]

    # dt/ds rides along as its own component of the integrated state, but it
    # is not an independent variable: dtds = sqrt(1 + m^2 l^2 / q^2), so a
    # derivative with respect to q/p has to pick up the path through it.
    # Differentiating against [t, d, l] alone drops that term -- which is
    # exactly what made this kernel disagree with the vacuum one on
    # d(time)/d(q/p). Differentiate against [t, d, dtds, l] and contract with
    # a seed that carries d(dtds)/dl into the q/p column instead.
    #
    # The system itself confirms the derivative: f gives d(dtds)/ds directly as
    # g m^2 l^3 / q^3, which is d(dtds)/dl times f's dl/ds.
    ddtds_dl = m**2 * l / (q**2 * dtds.name)
    TL = [t, d, dtds.name, l]
    seed_TL = sym.eye(6)[:, [0, 1, 2, 3, 5]]
    seed_TL[4, 4] = ddtds_dl

    # Transport the bound-to-free jacobian the way the vacuum kernel does:
    # push each live column through the same Runge-Kutta recursion as the state,
    # rather than assembling the 8x8 step jacobian and multiplying. Three
    # 5-vectors instead of an 8x8, and the chain rule comes from one shared
    # place rather than being spelled out per kernel.
    M_sym = MatrixSymbol("M", 8, 6)
    deriv = Derivation()
    for ne in [dtds, *step.named_exprs()]:
        deriv.record(ne)

    y = Matrix.vstack(p, Matrix([t, big_l]))
    ydot = Matrix.vstack(d, Matrix([dtds.name, l]))
    seed = free_param_seed()

    columns = []
    for tag, c in (("Tp", 2), ("Tt", 3), ("Tq", B2F_QOP_COLUMN)):
        zero = B2F_ZERO_ROWS.get(c, ())
        v = Matrix([[0 if i in zero else M_sym[i, c]] for i in range(8)])
        scaled = c == B2F_QOP_COLUMN
        if scaled:
            # stored = l * plain / plain[7]; undo it on the way in
            fac = v[7, 0] / l
            v = Matrix([[v[i, 0] * fac] if i < 7 else [v[i, 0]] for i in range(8)])

        ty, tydot = propagate_tangent(deriv, tag, step, y, ydot, seed * v)
        new_v = Matrix(
            [
                [ty[0, 0]],
                [ty[1, 0]],
                [ty[2, 0]],
                [ty[3, 0]],
                [tydot[0, 0]],
                [tydot[1, 0]],
                [tydot[2, 0]],
                [tydot[4, 0]],
            ]
        )
        if scaled:
            # ... and redo it with the q/p the step ended on
            g = new_l.expr / new_v[7, 0]
            new_v = Matrix(
                [[new_v[i, 0] * g] if i < 7 else [new_v[i, 0]] for i in range(8)]
            )
        columns.extend([new_v[i, 0]] for i, col in B2F_DENSE_LIVE if col == c)

    tangent_exprs = deriv.name_exprs[len([dtds, *step.named_exprs()]) :]
    new_M = name_expr("new_M", Matrix(columns))

    return [
        dtds,
        k1,
        step.y2,
        step.ydot2,
        step.ydot3,
        p2,
        l2,
        k2,
        l3,
        k3,
        step.y3,
        step.ydot4,
        p3,
        l4,
        k4,
        err,
        new_B,
        new_y,
        new_ydot,
        new_p,
        new_t,
        new_d_tmp,
        new_d,
        new_l,
        path_derivatives,
        # the per-column stage tangents, in the order they were derived
        *tangent_exprs,
        new_M,
    ]


# Preconditions the derivations rely on but the signatures cannot express.
# Debug only, so they cost nothing in a release build.
#
# `p_abs > 0` is not hypothetical: ParticleHypothesis::extractMomentum returns
# copysign(|q|, q/p) / (q/p), which is 0 for a neutral, and dt/ds is then
# infinite. Better to trip here than to propagate an infinite time.
INPUT_ASSERTS = [
    "  assert(std::abs(d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - 1) < 1e-8 &&",
    '         "direction must be a unit vector");',
    '  assert(p_abs > 0 && "absolute momentum must be positive");',
    '  assert(h != 0 && "step length must be non-zero");',
]


def print_rk4_vacuum_b2f(
    name_exprs: list[NamedExpr], run_cse: bool = False, mode: str = "combined"
) -> str:
    """Print the vacuum kernel.

    `mode` selects which of three shapes is emitted:

    - `jac`/`nojac`: the kernel specialised on covariance transport, which is
      what `SympyStepper::step` calls.  With `run_cse=False` the `nojac` body is
      exactly the `jac` body up to where the jacobian starts, so specialising
      duplicates code and not work.
    - `combined`: one kernel transporting the jacobian only for a non-empty `M`,
      for the dense translation unit's cold vacuum branch, which does not want
      two more kernels inlined into it.

    Do not give `nojac` a CSE pass of its own: sympy picks worse subexpressions
    when it sees fewer uses of them (it drops `B2[i] * l`, shared by two stages),
    which measured 2.9% slower at an equal multiply count.
    """
    printer = cxx_printer

    jac = mode != "nojac"
    output_names = ["p2", "p3", "err", "new_B", "new_p", "new_t", "new_d"]
    if jac:
        output_names += ["path_derivatives", "new_Mp", "new_Mt", "new_Mq"]
    outputs = [find_by_name(name_exprs, name)[0] for name in output_names]

    lines = []

    # not `name`: the expression hooks below bind that in their own loops
    fn_name = {
        "combined": "rk4_vacuum",
        "jac": "rk4_vacuum_jac",
        "nojac": "rk4_vacuum_nojac",
    }[mode]
    jac_params = ", std::span<T, 8> path_derivatives, std::span<T> M" if jac else ""
    lines.append(
        "template <typename T, typename GetB>\n"
        f"int {fn_name}(std::span<const T, 3> p,"
        " std::span<const T, 3> d, const T t, const T h, const T l, const T m,"
        " const T p_abs, std::span<const T, 3> B1, GetB getB, T& err,"
        " const T errTol, std::error_code& fieldErr,"
        " std::span<T, 3> new_p, T& new_t,"
        " std::span<T, 3> new_d, std::span<T, 3> new_B"
        f"{jac_params}) {{"
    )
    lines.extend(INPUT_ASSERTS)
    if jac:
        empty_ok = "M.empty() || " if mode == "combined" else ""
        lines.append(f"  assert({empty_ok}M.size() == {B2F_COLS * 8});")

    def pre_expr_hook(var):
        if str(var) == "p2":
            return "std::array<T, 3> p2;"
        if str(var) == "p3":
            return "std::array<T, 3> p3;"
        for name, group in B2F_LIVE_COLUMNS:
            if str(var) == name:
                return f"std::array<T, {len(group)}> {name};"
        return None

    def post_expr_hook(var):
        if str(var) == "p2":
            return "const auto B2res = getB(std::span<const T, 3>(p2));\n  if (!B2res.ok()) {\n    fieldErr = B2res.error();\n    return 2;\n  }\n  const auto B2 = *B2res;"
        if str(var) == "p3":
            return "const auto B3res = getB(std::span<const T, 3>(p3));\n  if (!B3res.ok()) {\n    fieldErr = B3res.error();\n    return 2;\n  }\n  const auto B3 = *B3res;"
        if str(var) == "err":
            return "if (err > errTol) {\n  return 0;\n}"
        if str(var) == "new_d" and mode == "combined":
            return "if (M.empty()) {\n  return 1;\n}"
        for name, group in B2F_LIVE_COLUMNS:
            if str(var) == name:
                return "\n".join(
                    f"M[{i + 8 * j}] = {name}[{k}];" for k, (i, j) in enumerate(group)
                )
        return None

    code = my_expression_print(
        printer,
        name_exprs,
        outputs,
        run_cse=run_cse,
        pre_expr_hook=pre_expr_hook,
        post_expr_hook=post_expr_hook,
    )
    lines.extend([f"  {l}" for l in code.split("\n")])

    lines.append("  return 1;")

    lines.append("}")

    return "\n".join(lines)


def print_rk4_dense(name_exprs: list[NamedExpr], run_cse: bool = True) -> str:
    printer = cxx_printer
    outputs = [
        find_by_name(name_exprs, name)[0]
        for name in [
            "p2",
            "l2",
            "l3",
            "p3",
            "l4",
            "err",
            "new_B",
            "new_p",
            "new_t",
            "new_d",
            "new_l",
            "path_derivatives",
            "new_M",
        ]
    ]

    lines = []

    lines.append(
        "template <typename T, typename GetB, typename GetG>\n"
        "Acts::Result<bool> rk4_dense(std::span<const T, 3> p,"
        " std::span<const T, 3> d, const T t, const T h, const T l, const T m,"
        " const T q, const T p_abs, std::span<const T, 3> B1, GetB getB,"
        " GetG getG, T& err, const T errTol, std::span<T, 3> new_p, T& new_t,"
        " std::span<T, 3> new_d, T& new_l, std::span<T, 3> new_B,"
        " std::span<T, 8> path_derivatives, std::span<T> M) {"
    )
    lines.extend(INPUT_ASSERTS)
    lines.append(f"  assert(M.empty() || M.size() == {B2F_COLS * 8});")

    lines.append("  const auto g1 = getG(p, l);")

    def pre_expr_hook(var):
        if str(var) == "p2":
            return "std::array<T, 3> p2;"
        if str(var) == "p3":
            return "std::array<T, 3> p3;"
        if str(var) == "l2":
            return "T l2;"
        if str(var) == "l3":
            return "T l3;"
        if str(var) == "l4":
            return "T l4;"
        if str(var) == "new_M":
            return f"std::array<T, {len(B2F_DENSE_LIVE)}> new_M;"
        return None

    def post_expr_hook(var):
        if str(var) == "p2":
            return "const auto B2res = getB(std::span<const T, 3>(p2));\n  if (!B2res.ok()) {\n    return Acts::Result<bool>::failure(B2res.error());\n  }\n  const auto B2 = *B2res;"
        if str(var) == "p3":
            return "const auto B3res = getB(std::span<const T, 3>(p3));\n  if (!B3res.ok()) {\n    return Acts::Result<bool>::failure(B3res.error());\n  }\n  const auto B3 = *B3res;"
        if str(var) == "l2":
            return "const auto g2 = getG(std::span<const T, 3>(p2), l2);"
        if str(var) == "l3":
            return "const auto g3 = getG(std::span<const T, 3>(p2), l3);"
        if str(var) == "l4":
            return "const auto g4 = getG(std::span<const T, 3>(p3), l4);"
        if str(var) == "err":
            return (
                "if (err > errTol) {\n  return Acts::Result<bool>::success(false);\n}"
            )
        if str(var) == "new_l":
            # after new_l, not after new_d: new_l is the energy loss update
            # of q/p and belongs to the step, not to the jacobian transport
            return "if (M.empty()) {\n  return Acts::Result<bool>::success(true);\n}"
        if str(var) == "new_M":
            return "\n".join(
                f"M[{i + 8 * j}] = new_M[{k}];"
                for k, (i, j) in enumerate(B2F_DENSE_LIVE)
            )
        return None

    code = my_expression_print(
        printer,
        name_exprs,
        outputs,
        run_cse=run_cse,
        pre_expr_hook=pre_expr_hook,
        post_expr_hook=post_expr_hook,
    )
    lines.extend([f"  {l}" for l in code.split("\n")])

    lines.append("  return Acts::Result<bool>::success(true);")

    lines.append("}")

    return "\n".join(lines)


def naive_rk4(eom):
    """One plain RK4 step of the equations of motion, nothing rearranged.

    The reference the optimised forms are checked against: no named stage
    shortcuts, no reuse identities, no scaled field. Both kernels use the same
    integrator, so only the right hand side differs.

    @param eom is the right hand side, `f(stage, x, y, ydot) -> d(ydot)/ds`
    @return the Derivation holding every named intermediate, the new state
            and its derivative, and the raw Rk4Step
    """
    dtds_sym = Symbol("dtds", real=True, positive=True)
    # `rk4_subexpr` integrates y with y' = ydot, so the two have to be the same
    # length. `big_l` is the unused integral of q/p that pads y to match; the
    # dense derivation carries the same dummy.
    big_l = Symbol("big_l", real=True)
    y = Matrix.vstack(p, Matrix([t, big_l]))
    ydot = Matrix.vstack(d, Matrix([dtds_sym, l]))
    step = rk4_subexpr(eom, s, y, ydot, h)

    # rk4_subexpr names its stages, so the results refer to k1..k4 rather than
    # spelling them out. Collect everything into a Derivation so the checker
    # can resolve back down to the inputs.
    deriv = Derivation()
    # dt/ds is a component of the integrated state but not an independent one;
    # binding it here is what lets a derivative with respect to q/p find its
    # way through, and what makes the value comparison meaningful at all.
    deriv.record(NamedExpr(dtds_sym, sym.sqrt(1 + m**2 / p_abs**2)))
    for ne in step.named_exprs():
        deriv.record(ne)
    return deriv, step.new_y, step.new_ydot, step


def check_atlas_matches_naive() -> None:
    """Assert the ATLAS-form vacuum step computes the plain RK4 step.

    The ATLAS arrangement cannot be *derived* from the naive one by
    substitution -- once the naive expression is expanded, compound stage
    quantities like A2 = A0 + d no longer appear syntactically, so there is
    nothing for a match to bind to. It can be *proved* though: resolve every
    ATLAS name back down to the inputs and compare. Same guarantee, and it is
    the direction sympy is reliable in.

    Raises AssertionError if the two disagree, which would mean the generated
    code is wrong.
    """
    naive, naive_y, naive_ydot, naive_step = naive_rk4(eom_vacuum)
    naive_y = explicit(naive.resolve(naive_y.expr))
    naive_ydot = explicit(naive.resolve(naive_ydot.expr))

    atlas = Derivation()
    atlas_rk4_stages(atlas, taylor_norm=False)

    def resolved(name):
        return explicit(atlas.resolve(find_by_name(atlas.name_exprs, name)[1]))

    # the naive stage slopes, for the quantities ATLAS expresses over them
    naive_k = [explicit(naive.resolve(k.expr))[0:3, 0] for k in naive_step.k]

    checks = [
        ("new_p", resolved("new_p"), naive_y[0:3, 0]),
        ("new_t", Matrix([resolved("new_t")]), Matrix([naive_y[3, 0]])),
        # An is three times the unnormalised new direction, by construction
        ("An", resolved("An"), 3 * naive_ydot[0:3, 0]),
        # The error estimate, compared before the norm is taken: ATLAS forms
        # (A1+A6)-(A3+A4) and scales by 2|h|, the plain form takes
        # h^2 * |k1-k2-k3+k4|.  Comparing the vectors rather than the assembled
        # norms tests the same algebra without asking expand to reason about
        # absolute values, which it cannot do.
        (
            "err vector",
            (resolved("A1") + resolved("A6")) - (resolved("A3") + resolved("A4")),
            (h / 2) * (naive_k[0] - naive_k[1] - naive_k[2] + naive_k[3]),
        ),
        # the path derivatives ATLAS reaches through Sl*A6 are just k4
        (
            "path_derivatives (dir rows)",
            resolved("path_derivatives")[4:7, 0],
            naive_k[3],
        ),
    ]
    for what, atlas, naive_expr in checks:
        diff = sym.expand(Matrix(atlas) - Matrix(naive_expr))
        if any(e != 0 for e in diff):
            raise AssertionError(
                f"ATLAS form disagrees with plain RK4 on {what}:\n{diff}"
            )


def free_param_seed(at: dict | None = None) -> Matrix:
    """Map derivatives against (y, ydot) onto derivatives against free params.

    `rk4_subexpr` differentiates against (pos, time, big_l) and
    (dir, dt/ds, q/p). Two things have to happen to reach the eight free
    parameters:

    - `big_l`, the dummy that pads y to ydot's length, is not one, so its
      column goes;
    - `dt/ds` is not independent of q/p, so its column folds into the q/p one
      weighted by d(dt/ds)/d(q/p). Leaving that out is what made the dense
      kernel disagree with the vacuum one on d(time)/d(q/p).

    @param at optionally evaluates the seed at a point
    @return the 10 x 8 map
    """
    dtds_sym = Symbol("dtds", real=True, positive=True)
    seed = sym.zeros(10, 8)
    for i in range(4):  # pos, time
        seed[i, i] = 1
    for i in range(3):  # dir
        seed[5 + i, 4 + i] = 1
    seed[8, 7] = m**2 * l / (q**2 * dtds_sym)  # d(dt/ds)/d(q/p)
    seed[9, 7] = 1
    return seed if at is None else seed.subs(at)


#: the point every check evaluates at
SAMPLE_SEED = 20250814


@functools.cache
def sample_point(seed: int = SAMPLE_SEED) -> dict:
    """A rational point to evaluate a symbolic identity at.

    An identity between rational functions either holds everywhere or fails on
    a set of measure zero, so evaluating at one arbitrary rational point
    settles it -- and unlike the closed forms, which run to megabytes for the
    derivative blocks, it stays small. Exact rationals, never floats, so a
    difference of zero means zero.

    Cached, so every check evaluates at the same point and the derived step
    jacobians can be shared between them.

    @param seed fixes the point, so a failure is reproducible
    @return a substitution for every input the kernels take
    """
    rng = random.Random(seed)

    def r():
        return sym.Rational(rng.randint(2, 40), rng.randint(2, 40))

    point = {h: r(), t: r(), l: r(), m: r(), p_abs: r(), s: r()}
    for vec in (p, d, B1, B2, B3):
        for i in range(3):
            point[vec[i]] = r()
    for gi in (g1, g2, g3, g4):
        point[gi] = r()
    point[Symbol("dtds", real=True, positive=True)] = sym.sqrt(
        1 + point[m] ** 2 / point[p_abs] ** 2
    )
    # p = |q| / |q/p| ties the three momentum inputs together
    point[q] = point[p_abs] * point[l]
    return point


@functools.cache
def naive_free_step_jacobian(eom, at_key: int | None = None) -> Matrix:
    """The free-to-free step jacobian D, straight from the chain rule.

    `rk4_subexpr` differentiates against the integrated state and its
    derivative, which is (pos, time, big_l) and (dir, dt/ds, q/p). Two things
    have to happen to turn that into the 8x8 over free parameters:

    - `big_l`, the dummy that pads y to ydot's length, is not a free parameter
      and its column and row are dropped;
    - `dt/ds` is not independent, so its column is folded into the q/p one
      weighted by d(dt/ds)/d(q/p). Leaving that out is what made the dense
      kernel disagree with the vacuum one on d(time)/d(q/p).

    The direction rows are the *unnormalised* direction's, matching what every
    other stepper here transports -- the normalisation derivative is dropped
    by all of them alike.

    Cached on the point, since several checks want the same matrix and it is
    the expensive part of all of them.

    @param eom is the right hand side, as `eom_vacuum` or `eom_dense`
    @param at_key is a `sample_point` seed, or None to stay symbolic
    @return the 8x8 step jacobian over (pos, time, dir, q/p)
    """
    at = None if at_key is None else sample_point(at_key)
    deriv, _, _, step = naive_rk4(eom)
    dy = explicit(deriv.resolve(step.dy.expr, at))
    dydot = explicit(deriv.resolve(step.dydot.expr, at))

    seed = free_param_seed(at)

    D = sym.zeros(8, 8)
    D[0:4, :] = dy[0:4, :] * seed  # position and time
    D[4:7, :] = dydot[0:3, :] * seed  # direction, unnormalised
    D[7, :] = dydot[4, :] * seed  # q/p
    return D


def check_jacobian_matches_naive() -> None:
    """Assert the vacuum kernel's transported columns equal D times M.

    The vacuum kernel never builds D -- it pushes the columns of the
    bound-to-free jacobian through the same Runge-Kutta recursion as the state.
    That is the whole reason it is fast, and also the reason it needs checking:
    nothing about the tangent recursion is obviously the same operation as the
    matrix product.

    Raises AssertionError if they disagree.
    """
    # The closed forms of the derivative blocks run to megabytes, so this one
    # is settled at a rational point rather than symbolically. The point also
    # supplies p = |q| / |q/p|, which ties together the three momentum inputs
    # the two sides spell the time row over.
    at = sample_point()
    D = naive_free_step_jacobian(eom_vacuum, SAMPLE_SEED)
    # q/p does not change in vacuum, so the column's scaling is the same on
    # both sides of the step
    expected = b2f_step_update(D, B2F_LIVE, at[l], at[l]).subs(at)

    deriv = Derivation()
    name_exprs = rk4_vacuum_b2f_atlasexpr()
    for ne in name_exprs:
        deriv.record(ne)
    got = Matrix.vstack(
        *[
            explicit(deriv.resolve(find_by_name(name_exprs, n)[1], at))
            for n in ("new_Mp", "new_Mt", "new_Mq")
        ]
    )

    diff = Matrix(got) - Matrix(expected)
    bad = [i for i, e in enumerate(diff) if sym.simplify(e) != 0]
    if bad:
        # one entry is enough to debug from, and these expressions run to
        # megabytes
        raise AssertionError(
            "transported bound-to-free columns disagree with D * M at "
            f"live entries {[B2F_LIVE[i] for i in bad]}; first difference:\n"
            f"{str(diff[bad[0]])[:2000]}"
        )


def check_dropped_rows_stay_zero() -> None:
    """Assert the rows the bound-to-free update drops can never be populated.

    `B2F_ZERO_ROWS` says the time and q/p rows of the phi and theta columns are
    zero, and the update drops them from the input. That much is a property of
    the bound-to-free jacobian, not something to derive here. What *is*
    derivable is that they stay zero across a step: with those rows zero on the
    way in, (D M)[3, c] is the sum over j not in {3, 7} of D[3, j] M[j, c], so
    it vanishes for arbitrary M only if the time and q/p rows of D have support
    nowhere but the time and q/p columns. If they ever did, dropping the rows
    would silently lose a term on every step.

    Raises AssertionError if either row reaches outside those columns.
    """
    at = sample_point()
    seed = free_param_seed(at)
    for what, eom in (("vacuum", eom_vacuum), ("dense", eom_dense)):
        deriv, _, _, step = naive_rk4(eom)
        # only the time and q/p rows are needed; the direction rows are the
        # expensive part of the jacobian and have nothing to do with this
        rows = {
            3: explicit(deriv.resolve(step.dy.expr[3, :], at)) * seed,
            7: explicit(deriv.resolve(step.dydot.expr[4, :], at)) * seed,
        }
        for row, values in rows.items():
            leaks = [
                col
                for col in range(8)
                if col not in (3, 7) and sym.simplify(values[0, col]) != 0
            ]
            if leaks:
                raise AssertionError(
                    f"{what}: row {row} of the step jacobian reaches columns "
                    f"{leaks}, so B2F_ZERO_ROWS would drop a live term"
                )


def check_dense_reduces_to_vacuum() -> None:
    """Assert the dense equations of motion become the vacuum ones without material.

    The two kernels are alternated freely along a trajectory, so they have to
    agree wherever they overlap. Setting the energy loss to zero is the case
    where they must agree exactly, and it is the check that would have caught
    both the transport sparsity disagreement and the d(time)/d(q/p) one.

    Raises AssertionError if they disagree.
    """
    vac, vac_y, vac_ydot, _ = naive_rk4(eom_vacuum)
    den, den_y, den_ydot, _ = naive_rk4(eom_dense)

    no_loss = [(gi, 0) for gi in (g1, g2, g3, g4)]
    for what, a, b in [
        ("state", vac.resolve(vac_y.expr), den.resolve(den_y.expr)),
        ("derivative", vac.resolve(vac_ydot.expr), den.resolve(den_ydot.expr)),
    ]:
        diff = sym.expand(explicit(a) - explicit(b).subs(no_loss))
        if any(e != 0 for e in diff):
            raise AssertionError(
                f"dense with zero energy loss differs from vacuum on {what}:\n{diff}"
            )


HEADER = """\
// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Note: This file is generated by generate_sympy_stepper.py
//       Do not modify it manually.

#pragma once

#include "Acts/Utilities/Result.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <span>
"""


# The bound-to-free form was measured against the tuned free-to-free one it
# replaced (see git history), clang 17, -O2, arm64, timing SympyStepper::step
# directly with chained steps:
#
#   with covariance transport      ~3-5% faster
#   without covariance transport   ~3-5% slower
#
# The isolated kernel has ~7% fewer instructions, but most of the source level
# saving is something clang already found on its own, and the vacuum path is
# latency bound on the serial chain between the three field lookups rather
# than throughput bound -- pre-scaling the field adds one multiply to exactly
# that chain, which is what costs the no-covariance case.
#
# Two knobs that were measured and left off:
#   taylor_norm=True   ATLAS' sqrt-free direction normalisation. It trades a
#                      square root and a division for a few more
#                      multiplications, and the direction chain is not the
#                      binding one, so it is neutral without covariance
#                      transport and slower with it.
#   run_cse=True       sympy-level CSE on top. Within noise, and it obscures
#                      the correspondence to the ATLAS code.


def main(argv: list[str]) -> None:
    """Generate the Runge-Kutta kernels.

    @param argv is the command line, argv[0] being the program name
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output",
        nargs="?",
        help="file to write the generated kernels to; stdout if omitted",
    )
    parser.add_argument(
        "--no-check",
        action="store_true",
        help="skip the symbolic equivalence assertions (they dominate the run "
        "time; only for iterating on the printers, never for a real build)",
    )
    args = parser.parse_args(argv[1:])

    if not args.no_check:
        # If one of these fires the generated code would be wrong, so they
        # guard the generator rather than a test that might not be run.
        check_atlas_matches_naive()
        check_dense_reduces_to_vacuum()
        check_jacobian_matches_naive()
        check_dropped_rows_stay_zero()

    with (
        open(args.output, "w") if args.output else contextlib.nullcontext(sys.stdout)
    ) as out:
        out.write(HEADER)
        out.write("\n")
        vacuum_exprs = rk4_vacuum_b2f_atlasexpr(taylor_norm=False)
        for vacuum_mode in ("combined", "jac", "nojac"):
            out.write(print_rk4_vacuum_b2f(vacuum_exprs, mode=vacuum_mode))
            out.write("\n\n")
        out.write(print_rk4_dense(rk4_dense_tunedexpr(), run_cse=True))
        out.write("\n")


if __name__ == "__main__":
    main(sys.argv)

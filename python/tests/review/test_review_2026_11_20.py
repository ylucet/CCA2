"""Failing tests for findings discovered on 2026-11-20 by rerunning the convexdb
cross-check after the review 2026-11-19 #3/#4/#12 fixes (task #95).

These are registry entries that the cross-check now ACCEPTS (they are established PLQ
functions with stored, established f* and co_f) but on which the engines DISAGREE with
the registry. The disagreement is an engine defect, not a parser artifact: it is pinned
here so a fix turns the test green. Each test asserts the CORRECT behaviour and is marked
xfail(strict=True) -- the suite stays green while the defect is there.

Finding #95-A -- the thin conjugate is lost: f = y^2/2 + bump(x) with
bump(x) = 1 - x^2 on (-1, 1) and 0 outside (registry `separable_bump_quadratic_2d`).
f is a valid 3-face PLQ. Its conjugate separates: f*(s) = s2^2/2 + h(s1) where
h(s1) = sup_x (s1 x - bump(x)) = 0 for s1 = 0 and +inf for s1 != 0 (the domain is
unbounded on both sides, where f is just y^2/2). So dom f* is the LINE {s1 = 0} and
f*(0, s2) = s2^2/2 -- exactly what the registry stores. Both the Python and the Rust
engine return a 0-face result (+inf everywhere, incl. at s1 = 0): they declare dom f*
empty instead of the line. Verified by hand at (0, 0.5): the truth is 0.125; both
engines give +inf; at (0.5, 0) the truth is +inf and both engines give +inf.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', '..')))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', '..', 'bench')))

from cca2.qua_pol import QuaPol                      # noqa: E402
import cross_convexdb as cc                          # noqa: E402
import sympy as sp                                   # noqa: E402

X, Y = sp.symbols('x y')
INF = float('inf')


def _bump_fixture():
    # f = y^2/2 + bump(x); the cross-check parser's own 3-face subdivision
    return cc.build_fixture_2d(Y ** 2 / 2 + sp.Piecewise((0, X <= -1), (1 - X ** 2, X < 1), (0, True)))


def _bump_quapol():
    fx = _bump_fixture()
    return QuaPol([list(r) for r in fx.V], [list(r) for r in fx.E],
                  [list(r) for r in fx.f], [list(r) for r in fx.F])


@pytest.mark.xfail(strict=True,
                   reason='finding #95-A: conj loses the thin domain (the line s1=0) of the bump conjugate')
def test_conj_bump_thin_domain_is_finite_on_the_line():
    g = _bump_quapol().conj()
    v = g.eval([[0, 0], [0, 0.5], [0, 2], [0, -1]])[0]
    # f*(0, s2) = s2^2/2  (attained at y = s2 and any |x| >= 1, where f = y^2/2)
    assert v == [0, 0.125, 2, 0.5], v


def test_conj_bump_infinite_off_the_line_control():
    # control (the engine is RIGHT off the line -- the defect is the lost line itself):
    # f*(s1, s2) = +inf for s1 != 0, because f is y^2/2 on the unbounded rays x >= 1 and
    # x <= -1, where <s, x> - f = s1 x grows without bound
    g = _bump_quapol().conj()
    v = g.eval([[0.5, 0], [-0.5, 0], [1, 1]])[0]
    assert v == [INF, INF, INF], v


@pytest.mark.xfail(strict=True,
                   reason='finding #95-A: biconj inherits the lost thin domain')
def test_biconj_bump_fenchel_young_at_the_line():
    f = _bump_quapol()
    g = f.conj()
    b = f.biconj()
    s = [0, 0.5]
    # Fenchel-Young: <s, x> - f*(s) is a minorant; at s = (0, 0.5) it is -0.125, finite
    fs = g.eval([s])[0][0]
    assert fs == 0.125
    assert b.eval([s])[0][0] <= 0 + 0.125 + 1e-12

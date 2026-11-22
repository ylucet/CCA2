"""Failing tests for the findings of the 2026-11-19 code review (reviews/codereview-2026-11-19.md).

Each test asserts the CORRECT behaviour and is marked xfail(strict=True): the suite stays green
while the defect is there, and turns red (XPASS) the day it is fixed -- remove the mark then.
"""
import math
import os
import sys
from fractions import Fraction as Fr

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', '..')))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', '..', 'bench')))
sys.path.insert(0, HERE)

from cca2.qua_pol import QuaPol                      # noqa: E402
from cca2 import check                               # noqa: E402
from cca2 import rat_q                               # noqa: E402
from cca2.cell_has_interior import cellHasInterior   # noqa: E402
import plq_oracle as O                               # noqa: E402
from cca2.errors import CCA2Error                    # noqa: E402

INF = math.inf


def _cone_y2():
    # f(x,y) = y^2/2 on the cone {x >= |y|}: rays (1,-1), (1,1)
    return QuaPol([[0, 0], [1, -1], [1, 1]], [[1, 2, 0], [1, 3, 0]],
                  [[0, 0, 0, 0, 0, 0, 1, 0, 0, 0]], [[1, 0], [0, 1]])


# ---------------------------------------------------------------------------------------------
# Finding 1 (CRITICAL): conj on a cone whose zero-curvature direction is strictly inside
# ---------------------------------------------------------------------------------------------
def test_conj_cone_with_interior_null_direction_is_infinite():
    g = _cone_y2().conj()
    v = g.eval([[1, 0], [0.5, 0.3], [0.1, 0], [-1, 0]])[0]
    # along x = (t, 0): <s,x> - f(x) = s1 t -> +inf for s1 > 0
    assert v[:3] == [INF, INF, INF]
    assert v[3] == 0


def test_conj_two_cone_fan_R6_152():
    # piece 1: y^2 - x - 2y on cone((1,-1),(2,1)); piece 2 continuous across the ray (2,1)
    f = QuaPol([[0, 0], [1, -1], [2, 1], [-1, 2]], [[1, 2, 0], [1, 3, 0], [1, 4, 0]],
               [[0, 0, 0, 0, 0, 0, 2, -1, -2, 0], [0, 0, 0, 0, 4, -6, 10, -3, 2, 0]],
               [[1, 0], [2, 1], [0, 2]])
    s = [-0.96875, -2.953125]
    # (t, 0) lies in piece 1 for t > 0, where <s,x> - f = (s1 + 1) t -> +inf
    assert O.conj_exact({'pieces': [{'kind': 'cone', 'apex': (Fr(0), Fr(0)), 'd1': (Fr(1), Fr(-1)),
                                     'd2': (Fr(2), Fr(1)), 'q': (Fr(0), Fr(0), Fr(2), Fr(-1), Fr(-2), Fr(0))}]},
                        (Fr(s[0]), Fr(s[1]))) is None
    assert f.conj().eval([s])[0][0] == INF


# ---------------------------------------------------------------------------------------------
# Finding 2 (MAJOR): cellHasInterior's 'complete' CAD fallback ignores vertical lines
# ---------------------------------------------------------------------------------------------

def test_cell_has_interior_far_vertical_strip():
    EcQ = [rat_q.conic([0, 0, 0, 1, 0, -10]),      # x = 10
           rat_q.conic([0, 0, 0, 2, 0, -21]),      # x = 10.5
           rat_q.conic([1, 0, 1, -200, 0, 9999])]  # circle, centre (100, 0), radius 1
    assert cellHasInterior(EcQ, [[1, 1], [2, -1], [3, 1]])   # (10.25, 0) is interior


def test_cell_has_interior_far_horizontal_strip_control():
    # the horizontal twin is found: the defect is specific to vertical lines
    EcQ = [rat_q.conic([0, 0, 0, 0, 1, -10]), rat_q.conic([0, 0, 0, 0, 2, -21]),
           rat_q.conic([1, 0, 1, 0, -200, 9999])]
    assert cellHasInterior(EcQ, [[1, 1], [2, -1], [3, 1]])


# ---------------------------------------------------------------------------------------------
# Finding 5 (MINOR, test oracle): check.supAt samples recession directions
# ---------------------------------------------------------------------------------------------
def test_definition_oracle_detects_unbounded_cone():
    # R5-053: q = x^2 - 2x - y + 2 on the cone apex (-1,1), rays (-1,0), (1,-2); the null
    # direction (0,-1) of H is strictly inside, and the slope along it is -s2 - 1 > 0 here
    f = QuaPol([[-1, 1], [-2, 1], [0, -1]], [[1, 2, 0], [1, 3, 0]],
               [[0, 0, 0, 0, 2, 0, 0, -2, -1, 2]], [[1, 0], [0, 1]])
    assert check.supAt(f, [-2.96875, -2.953125]) == INF


# ---------------------------------------------------------------------------------------------
# Findings 3 and 4 (MAJOR / MINOR): the convexdb cross-check (#64)
# ---------------------------------------------------------------------------------------------
cc = pytest.importorskip('cross_convexdb')
sp = pytest.importorskip('sympy')
X, Y = sp.symbols('x y')


@pytest.mark.xfail(strict=True, reason='review 2026-11-19 #3: |x|+|y| (PLQ) is refused as non-PLQ')
def test_crosscheck_accepts_norm1():
    fx = cc.build_fixture_2d(sp.Abs(X) + sp.Abs(Y))
    assert len(fx.f) == 4


@pytest.mark.xfail(strict=True, reason='review 2026-11-19 #3: 1-D chain parser ignores the half-plane side')
def test_crosscheck_accepts_bounded_segment_indicator():
    f = sp.Piecewise((0, sp.Eq(X, Y) & (X >= 0) & (X <= 1)), (sp.oo, True))
    fx = cc.build_fixture_1d(f)
    assert sorted(tuple(v) for v in fx.V) == [(0, 0), (1, 1)]


@pytest.mark.xfail(strict=True, reason='review 2026-11-19 #4: int() truncates rational coefficients')
def test_crosscheck_keeps_rational_coefficients():
    row = cc.piece_value_row(X * Y / 3 + X / 2)
    assert row[5] == Fr(1, 3) and row[7] == Fr(1, 2)


@pytest.mark.xfail(strict=True, reason='review 2026-11-19 #4: int() truncates rational half-plane data')
def test_crosscheck_keeps_rational_halfplane():
    (hp,) = cc.atom_to_hp_tuple(Y <= sp.Rational(1, 2))
    assert hp.holds((Fr(0), Fr(1, 4))) and not hp.holds((Fr(0), Fr(3, 4)))


# ---------------------------------------------------------------------------------------------
# reviewer's oracle sanity (not a finding): the exact conjugate on hand-checked cases
# ---------------------------------------------------------------------------------------------
def test_review_oracle_hand_cases():
    q = (Fr(1), Fr(0), Fr(1), Fr(0), Fr(0), Fr(0))
    assert O.conj_exact({'pieces': [{'kind': 'plane', 'q': q}]}, (Fr(1), Fr(2))) == Fr(5, 2)
    tri = [(Fr(0), Fr(0)), (Fr(1), Fr(0)), (Fr(0), Fr(1))]
    zero = (Fr(0),) * 6
    assert O.conj_exact({'pieces': [{'kind': 'poly', 'P': tri, 'q': zero}]}, (Fr(1), Fr(2))) == 2
    cone = {'kind': 'cone', 'apex': (Fr(0), Fr(0)), 'd1': (Fr(1), Fr(0)), 'd2': (Fr(0), Fr(1)), 'q': q}
    assert O.conj_exact({'pieces': [cone]}, (Fr(1), Fr(-2))) == Fr(1, 2)
    y2 = {'kind': 'cone', 'apex': (Fr(0), Fr(0)), 'd1': (Fr(1), Fr(-1)), 'd2': (Fr(1), Fr(1)),
          'q': (Fr(0), Fr(0), Fr(1), Fr(0), Fr(0), Fr(0))}
    assert O.conj_exact({'pieces': [y2]}, (Fr(1), Fr(0))) is None
    assert O.conj_exact({'pieces': [y2]}, (Fr(-1), Fr(0))) == 0


# Finding 2, second mechanism: rows that are DEGENERATE conics sharing a line. conicMeet reports
# the pair 'degenerate' (infinitely many common points) and returns no points, so the other lines'
# crossings never become critical abscissae of the CAD. Reached inside a real conj call (random
# nonconvex hexagon V = (5,3),(2,1),(-1,6),(-6,2),(3,-2),(4,-2), q = (3,-3,3,-1,-2,2), found by
# tests/review/hunt_cell_interior.py); the strict interior point below is checked exactly.
HUNT_ECQ = [[101, 76, -25, 3378, 3000, 9225], [91, 38, -53, -174, -606, -1341],
            [64, 37, -27, 2386, 2113, 6582], [1221, 154, -1067, -42386, -49250, -138147],
            [0, 0, 0, 3, -5, 65], [0, 0, 0, 3, -5, -31], [0, 0, 0, -18, 8, -119],
            [0, 0, 0, 9, -4, 313], [0, 0, 0, 5, -2, 169]]
HUNT_ROWS = [[1, 1], [2, -1], [3, -1], [4, 1], [5, 1], [6, 1], [7, 1], [8, 1], [9, 1]]


def test_hunt_cell_really_has_interior():
    x, y = Fr(-58.425327423549874), Fr(-109.38222339485247)
    for e, side in HUNT_ROWS:
        c = HUNT_ECQ[e - 1]
        v = c[0] * x * x + c[1] * x * y + c[2] * y * y + c[3] * x + c[4] * y + c[5]
        assert side * v > 0


def test_cell_has_interior_shared_line_factor():
    assert cellHasInterior(HUNT_ECQ, HUNT_ROWS)


# ---------------------------------------------------------------------------------------------
# Finding 8 (CRITICAL): conj too HIGH on a bounded nonconvex polygon (PSD-singular quadratic)
# ---------------------------------------------------------------------------------------------
def _L_y2():
    V = [[0, 0], [2, 0], [2, 1], [1, 1], [1, 2], [0, 2]]
    return QuaPol(V, [[i + 1, (i + 1) % 6 + 1, 1] for i in range(6)], [[0, 0, 0, 0, 0, 0, 1, 0, 0, 0]], [[1, 0]] * 6)


def test_conj_L_shape_y2():
    # f = y^2/2 on the L-shape of golden conj record 36. sup at s = (2.53125, -1.453125) is
    # attained at (2, 0): 2 s1 = 5.0625 (y = 0 since s2 < 0). Engines (all three): 6.5137939453125.
    s = [2.53125, -1.453125]
    assert check.supAt(_L_y2(), s) == 5.0625
    assert abs(_L_y2().conj().eval([s])[0][0] - 5.0625) < 1e-12


# ---------------------------------------------------------------------------------------------
# Finding 9 (CRITICAL): biconj values wrong on small continuous inputs (random differential test)
# ---------------------------------------------------------------------------------------------
SQ = ([[0, 0], [1, 0], [1, 1], [0, 1]], [[1, 2, 1], [2, 3, 1], [3, 4, 1], [4, 1, 1]], [[1, 0]] * 4)
LV = [[0, 0], [2, 0], [2, 1], [1, 1], [1, 2], [0, 2]]
LE = [[i + 1, (i + 1) % 6 + 1, 1] for i in range(6)]


@pytest.mark.xfail(strict=True, reason='review 2026-11-19 #9a: f** finite outside conv(dom f)')
def test_biconj_infinite_outside_the_square():
    # f = (y^2 - x^2)/2 on the unit square (E-square-H-101); x = (-0.49, 0.01) is outside the
    # square, so f**(x) = +inf. Engines: 0.24616241455078125.
    f = QuaPol(SQ[0], SQ[1], [[0, 0, 0, 0, -1, 0, 1, 0, 0, 0]], SQ[2])
    assert f.biconj().eval([[-0.4921875, 0.01171875]])[0][0] == INF


@pytest.mark.xfail(strict=True, reason='review 2026-11-19 #9b: f** below <s,x> - f*(s) (Fenchel-Young vs its own conj)')
def test_biconj_fenchel_young_L_shape():
    # f = -xy + y^2/2 on the L-shape (E-L-H0-11). For ANY s, f**(x) >= <s,x> - f*(s); with the
    # engine's OWN conj at s below, the bound is 1.04606, the engine's f**(x) is 0.40871.
    f = QuaPol(LV, LE, [[0, 0, 0, 0, 0, -1, 1, 0, 0, 0]], [[1, 0]] * 6)
    x = [0.2578125, 1.76171875]
    s = [-1.9484626037100141, 1.6789473678974558]
    g = f.conj().eval([s])[0][0]
    assert abs(g - check.supAt(f, s)) < 1e-9          # the engine's conj is right at s
    h = f.biconj().eval([x])
    h = h[0] if isinstance(h, tuple) else h
    assert h[0] >= s[0] * x[0] + s[1] * x[1] - g - 1e-9


def test_biconj_finite_inside_the_domain_R3_072():
    # two continuous pieces on a quadrilateral (R3-072); x is in piece 1, f(x) = -12.29 is finite,
    # so f**(x) <= f(x) < inf. Python returned +inf until #108: realPoint took a spurious ~1e31
    # root along a parabola's axis, so the tie arc that needs the missing atom pair was never
    # sampled (biconj_atoms.realPoint).
    f = QuaPol([[-3, 2], [-1, 0], [2, 0], [-2, 3]], [[1, 2, 1], [2, 3, 1], [3, 1, 1], [3, 4, 1], [4, 1, 1]],
               [[0, 0, 0, 0, 0, 2, -2, 2, 1, 1], [0, 0, 0, 0, 4, 9, 8, -4, -8, 5]],
               [[1, 0], [1, 0], [1, 2], [2, 0], [2, 0]])
    x = [-2.4921875, 1.51171875]
    assert math.isfinite(check.evalQuaPol(f, [x])[0])
    h = f.biconj().eval([x])
    h = h[0] if isinstance(h, tuple) else h
    assert math.isfinite(h[0])


# ---------------------------------------------------------------------------------------------
# Finding 10 (MAJOR): discontinuous input is accepted silently and gives a wrong biconjugate
# ---------------------------------------------------------------------------------------------
def _d_family():
    """The 30 discontinuous two-piece D fixtures of the review generator (plq_gen.py, seed
    20261119), frozen in d_family.json: each is a convex quadrilateral cut by a diagonal, the
    two pieces carrying DIFFERENT quadratics, so the shared edge is a genuine jump."""
    import json
    path = os.path.join(HERE, 'd_family.json')
    return json.load(open(path))


# #108 lifted #96's stopgap refusal: a discontinuous f is read LOWER SEMICONTINUOUS (at a shared
# point, the min over the closed pieces holding it), f* = sup over the closed pieces, f** = cl co f.
# d_family_ref.json (make_d_family_ref.py) holds the Python values the other two engines must match;
# they were graded 0 wrong by compare_engines.py (conj exact, biconj inside the certified bracket).
_D_REF = None
_D_PIECES = None
# Minutes-long in Python (the fast lane keeps the rest; Rust runs them under -- --ignored).
_D_SLOW = {'D-069', 'D-111', 'D-153', 'D-188', 'D-209'}


def _d_ref():
    global _D_REF
    if _D_REF is None:
        import json
        _D_REF = {r['name']: r for r in json.load(open(os.path.join(HERE, 'd_family_ref.json')))}
    return _D_REF


def _d_pieces(name):
    global _D_PIECES
    if _D_PIECES is None:
        import plq_gen
        _D_PIECES = {x['name']: plq_gen.load_pieces(x['pieces'])
                     for x in plq_gen.build(seed=20261119, nrand=210) if x['name'].startswith('D-')}
    return _D_PIECES[name]


def _rel(a, b):
    return abs(a - b) / max(1.0, abs(b))


def _ev(g, P):
    v = g.eval(P)
    return v[0] if isinstance(v, tuple) else v


@pytest.mark.parametrize('name', [pytest.param('D-%03d' % i, marks=pytest.mark.slow)
                                  if 'D-%03d' % i in _D_SLOW else 'D-%03d' % i
                                  for i in range(6, 210, 7)])
def test_discontinuous_input_lsc(name):
    import numpy as np
    from scipy.optimize import minimize
    rec = _d_ref()[name]
    fx = _d_pieces(name)
    f = QuaPol(rec['V'], rec['E'], rec['f'], rec['F'])
    assert not f.isContinuous()[0], name                      # a genuine jump: the family's point
    # ---- conj: the sup over the closed pieces (exact oracle), and the cross-engine reference
    g = f.conj()
    gv = _ev(g, rec['S'])
    for s, v, r in zip(rec['S'], gv, rec['conj']):
        o = O.conj_exact(fx, (Fr(s[0]), Fr(s[1])))
        o = INF if o is None else float(o)
        assert (v == o) if math.isinf(o) else _rel(v, o) < 1e-9, (name, s, v, o)
        r = float(r) if isinstance(r, str) else r
        assert (v == r) if math.isinf(r) else _rel(v, r) < 1e-9
    # ---- biconj
    if 'biconj_err' in rec:
        # a named refusal shared with continuous inputs (curvedCell: a convex piece with a singular
        # Hessian; noWitness), pinned so that a new refusal or a new answer shows here
        with pytest.raises(CCA2Error) as e:
            f.biconj()
        assert e.value.identifier == rec['biconj_err']
        return
    h = f.biconj()
    X = rec['X']
    hv = _ev(h, X)
    for x, v, r in zip(X, hv, rec['biconj']):
        r = float(r) if isinstance(r, str) else r
        assert (v == r) if math.isinf(r) else _rel(v, r) < 1e-7, (name, x, v, r)
    # f** <= a certified upper bound (convex combinations of domain samples)
    samples = O.domain_samples(fx, n=16)
    for x, v in list(zip(X, hv))[::4]:
        hi = O.hull_upper(fx, samples, (Fr(x[0]), Fr(x[1])))
        assert v <= hi + 1e-7 * max(1.0, abs(hi)), (name, x, v, hi)
    # Fenchel-Young against the engine's OWN conjugate, and f** = (f*)*: the sup over s of
    # <s,x> - f*(s) (the five best points of a grid, then Nelder-Mead from each, restarted with a
    # unit simplex -- one short run stalls on the kinks of f*) meets f**(x)
    gs = np.linspace(-16, 16, 65)
    grid = np.array([[a, b] for a in gs for b in gs])
    gg = np.array(_ev(g, grid.tolist()))
    for j, (x, v) in enumerate(zip(X, hv)):
        ph = grid @ np.array(x) - gg
        assert v >= np.max(ph) - 1e-9 * max(1.0, abs(v)), (name, x, v, np.max(ph))
        if not math.isfinite(v) or j % 2:
            continue

        def neg(s, x=x):
            val = _ev(g, [[float(s[0]), float(s[1])]])[0]
            return 1e300 if not math.isfinite(val) else -(s[0] * x[0] + s[1] * x[1] - val)
        cc = -math.inf
        for i in np.argsort(-ph)[:5]:
            s0, c0 = grid[i], ph[i]
            for _ in range(4):
                r = minimize(neg, s0, method='Nelder-Mead',
                             options={'xatol': 1e-12, 'fatol': 1e-14, 'maxiter': 2000,
                                      'initial_simplex': [s0, s0 + [1, 0], s0 + [0, 1]]})
                if -r.fun <= c0:
                    break
                s0, c0 = r.x, -r.fun
            cc = max(cc, c0)
        assert v >= cc - 1e-9 * max(1.0, abs(v)), (name, x, v, cc)        # Fenchel-Young
        assert v - cc <= 1e-6 * max(1.0, abs(v)), (name, x, v, cc)       # f** = conj(conj)


def test_continuous_input_not_refused():
    # R3 (continuous 2-piece split) and R6 (continuous two-cone fan) carry shared edges whose
    # pieces agree along the edge; the check must pass them, so conj still answers.
    import plq_gen
    for fam in ('R3', 'R6'):
        cases = [x for x in plq_gen.build(seed=20261119, nrand=210) if x['name'].startswith(fam + '-')]
        assert cases, fam
        for rec in cases[:5]:
            f = QuaPol(rec['V'], rec['E'], rec['f'], rec['F'])
            try:
                f.conj()
            except CCA2Error as e:
                raise AssertionError('%s: continuous input refused as %s' % (rec['name'], e.identifier))

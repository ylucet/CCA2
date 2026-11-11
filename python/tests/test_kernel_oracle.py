"""The exact kernel against an independent oracle (sympy -- test-only, never in the compute
path): root counting and isolation, signs at algebraic numbers, gcds, square-free parts,
resultant-free conic intersection, linear feasibility."""
import random
from fractions import Fraction

import pytest
import sympy as sp

from cca2 import biv, rat_q
from cca2.alg_row_fails import _oddPart, _relSign, _squareFreeSign
from cca2.isolated_root_sign import isolatedRootSign

x, y = sp.symbols('x y')


def _poly(c):
    return sp.Poly(c, x)


def _rand_poly(rng, deg, lo=-9, hi=9):
    c = [rng.randint(lo, hi) for _ in range(deg + 1)]
    if c[0] == 0:
        c[0] = 1
    return c


def test_sturm_sign_defect_fixed():
    # MATLAB's ratQ.countRootsIn returns 0 here (measured 2026-11-10); x^4 + x - 1 has two real roots
    assert rat_q.countRootsIn([1, 0, 0, 1, -1], -2, 2) == 2


@pytest.mark.parametrize('seed', range(40))
def test_countRootsIn_matches_sympy(seed):
    rng = random.Random(seed)
    p = _rand_poly(rng, rng.randint(1, 6))
    a = Fraction(rng.randint(-40, 0), 8)
    b = a + Fraction(rng.randint(1, 60), 8)
    roots = set(sp.Poly(p, x).real_roots())            # distinct real roots (set drops multiplicity)
    want = sum(1 for r in roots if sp.Rational(a.numerator, a.denominator) < r
               <= sp.Rational(b.numerator, b.denominator))
    assert rat_q.countRootsIn(p, a, b) == want


@pytest.mark.parametrize('seed', range(30))
def test_isolateRoots_and_signAt(seed):
    rng = random.Random(1000 + seed)
    p = _rand_poly(rng, rng.randint(2, 5))
    q = _rand_poly(rng, rng.randint(1, 4))
    roots = sorted(set(sp.Poly(p, x).real_roots()))
    I = rat_q.isolateRoots(p)
    assert len(I) == len(roots)
    for (lo, hi), r in zip(I, roots):
        assert sp.Rational(lo.numerator, lo.denominator) < r <= sp.Rational(hi.numerator, hi.denominator)
        v = sp.Poly(q, x).eval(r)
        want = 0 if sp.simplify(v) == 0 else (1 if sp.N(v, 60) > 0 else -1)
        assert rat_q.signAt(q, p, (lo, hi)) == want
        assert isolatedRootSign(q, rat_q.polySquarefree(p), (lo, hi)) == want


@pytest.mark.parametrize('seed', range(30))
def test_polyGcd_and_squarefree(seed):
    rng = random.Random(2000 + seed)
    g = _rand_poly(rng, rng.randint(0, 2))
    a = rat_q.polyMul1(_rand_poly(rng, rng.randint(1, 3)), g)
    b = rat_q.polyMul1(_rand_poly(rng, rng.randint(1, 3)), g)
    got = sp.Poly(rat_q.polyGcd(a, b), x)
    want = sp.gcd(_poly(a), _poly(b))
    assert sp.degree(got, x) == sp.degree(want, x)
    assert sp.rem(_poly(a), got, domain='QQ') == 0 and sp.rem(_poly(b), got, domain='QQ') == 0
    sq = rat_q.polyMul1(a, g)
    sf = sp.Poly(rat_q.polySquarefree(sq), x)
    assert sf.gcd(sf.diff(x)).degree() == 0                # squarefree
    assert set(sf.real_roots()) == set(_poly(sq).real_roots())


@pytest.mark.parametrize('seed', range(40))
def test_feasible2_matches_vertex_enumeration(seed):
    rng = random.Random(3000 + seed)
    P = [[rng.randint(-5, 5), rng.randint(-5, 5), rng.randint(-6, 6)] for _ in range(rng.randint(1, 6))]
    P = [r for r in P if r[0] or r[1]] or [[1, 0, 0]]
    # strict interior <=> the LP  max t : a.s + c >= t  has a positive optimum (bounded by a box)
    A = P + [[1, 0, 50], [-1, 0, 50], [0, 1, 50], [0, -1, 50]]
    s1, s2, t = sp.symbols('s1 s2 t')
    best = None
    import itertools
    for i, j, k in itertools.combinations(range(len(A)), 3):
        M = sp.Matrix([[A[i][0], A[i][1], -1], [A[j][0], A[j][1], -1], [A[k][0], A[k][1], -1]])
        if M.det() == 0:
            continue
        z = M.solve(sp.Matrix([-A[i][2], -A[j][2], -A[k][2]]))
        if all(r[0] * z[0] + r[1] * z[1] + r[2] - z[2] >= 0 for r in A):
            best = z[2] if best is None else max(best, z[2])
    # the box is far larger than any bounded cell these rows cut out, so it never changes the sign
    assert rat_q.feasible2(P, True) == (best is not None and best > 0)


def _mat_to_expr(M):
    return sum(c * x ** i * y ** j for i, r in enumerate(M) for j, c in enumerate(r))


def _rand_biv(rng, dx, dy):
    M = [[rng.randint(-4, 4) for _ in range(dy + 1)] for _ in range(dx + 1)]
    if not any(any(r) for r in M):
        M[0][0] = 1
    return M


@pytest.mark.parametrize('seed', range(25))
def test_polyGcdBiv_matches_sympy(seed):
    rng = random.Random(4000 + seed)
    g = _rand_biv(rng, rng.randint(0, 1), rng.randint(0, 1))
    A = biv.conv2(_rand_biv(rng, 1, 1), g)
    B = biv.conv2(_rand_biv(rng, 1, 1), g)
    got = sp.expand(_mat_to_expr(biv.polyGcdBiv(A, B)))
    want = sp.gcd(sp.expand(_mat_to_expr(A)), sp.expand(_mat_to_expr(B)))
    # equal up to a rational constant
    assert sp.simplify(sp.cancel(got / want)).free_symbols == set() if want != 0 else got == 0


@pytest.mark.parametrize('seed', range(20))
def test_squareFreeSign_keeps_the_sign(seed):
    rng = random.Random(5000 + seed)
    f1 = _rand_biv(rng, 1, 1)
    f2 = _rand_biv(rng, 1, 1)
    R = biv.conv2(biv.conv2(f1, f1), f2)                   # f1^2 * f2
    S = _squareFreeSign(R)
    eR = _mat_to_expr(R)
    eS = _mat_to_expr(S)
    for _ in range(12):
        px = sp.Rational(rng.randint(-30, 30), 7)
        py = sp.Rational(rng.randint(-30, 30), 11)
        vr = eR.subs({x: px, y: py})
        vs = eS.subs({x: px, y: py})
        if vr != 0:
            assert sp.sign(vr) == sp.sign(vs)
    # and no squared factor is left
    fac = sp.factor_list(sp.expand(eS))[1]
    assert all(m == 1 for _, m in fac)

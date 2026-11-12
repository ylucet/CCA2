"""_poly -- the ONE place the exact engine may call python-flint (FLINT, LGPL), as an optional
accelerator (python/PORT_PLAN.md: "behind a thin facade, added only when a kernel is measured
hot"). Everything here has a pure-Python twin that stays the reference: without python-flint, or
with CCA2_FLINT=0, the callers use it and the answers are the same.

Measured hot (2026-11-11): the odd-multiplicity part of a bivariate integer polynomial (Yun's
square-free decomposition) for algRowFails' surd rows. On erratum Example 2's quartic bridge faces
the pure-Python primitive-PRS gcds ran for over half an hour; FLINT's factor_squarefree is exact
and takes milliseconds.
"""
import os

_ctx = None
if os.environ.get('CCA2_FLINT', '1') != '0':
    try:
        import flint as _flint
        _ctx = _flint.fmpz_mpoly_ctx.get(('x', 'y'), 'lex')
    except ImportError:                             # python-flint absent: the pure-Python path
        _ctx = None


def available():
    return _ctx is not None


def oddPart(S):
    """The product of the square-free factors of ODD multiplicity of the integer polynomial S
    (dense matrix, S[i][j] multiplies x^i y^j), primitive; its sign is left to the caller (as
    Yun's). None when FLINT is not in use."""
    if _ctx is None:
        return None
    d = {(i, j): int(c) for i, r in enumerate(S) for j, c in enumerate(r) if c != 0}
    if not d:
        return None
    _, fac = _ctx.from_dict(d).factor_squarefree()
    O = _ctx.from_dict({(0, 0): 1})
    for f, m in fac:
        if m % 2 == 1:
            O = O * f
    out = {tuple(int(e) for e in k): int(v) for k, v in O.to_dict().items()}
    n1 = max(k[0] for k in out) + 1
    n2 = max(k[1] for k in out) + 1
    M = [[0] * n2 for _ in range(n1)]
    for (i, j), v in out.items():
        M[i][j] = v
    return M


def gcdUni(a, b):
    """gcd of two univariate integer polynomials (lists, highest degree first), up to sign and
    content -- for callers that only count its roots. None when FLINT is not in use."""
    if _ctx is None or not a or not b:
        return None
    g = _flint.fmpz_poly([int(x) for x in a[::-1]]).gcd(_flint.fmpz_poly([int(x) for x in b[::-1]]))
    return [int(x) for x in g.coeffs()[::-1]]

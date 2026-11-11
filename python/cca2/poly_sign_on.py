"""poly_sign_on -- the sign of a bivariate polynomial on the part of a box that meets a region,
by Bernstein enclosures and subdivision (port of polySignOn.m).

P is a coefficient matrix (list of rows): P[i][j] multiplies x^i y^j. box = [x0, x1, y0, y1];
inBox(b) says whether a box meets the region. Returns +1 / -1 when the sign is certain on every
sub-box that meets the region, 0 when undecided at the given depth, and None (MATLAB NaN) when
no sub-box meets the region. Doubles, with a relative tolerance.
"""
import math

import numpy as np

_binomCache = {}
_convCache = {}


def polySignOn(P, box, inBox, depth=6):
    return _signRec(np.array(P, dtype=float), list(box), inBox, depth)


def _signRec(P, b, inBox, depth):
    if not inBox(b):
        return None                                 # misses the region: no constraint from here
    B = _bernstein(P, b)
    m = float(np.max(np.abs(B)))
    tol = 1e-10 * max(m, 2.0 ** -52)
    if np.all(B > tol):
        return 1
    if np.all(B < -tol):
        return -1
    if depth == 0:
        return 0
    xm = (b[0] + b[1]) / 2
    ym = (b[2] + b[3]) / 2
    q = [_signRec(P, [b[0], xm, b[2], ym], inBox, depth - 1),
         _signRec(P, [xm, b[1], b[2], ym], inBox, depth - 1),
         _signRec(P, [b[0], xm, ym, b[3]], inBox, depth - 1),
         _signRec(P, [xm, b[1], ym, b[3]], inBox, depth - 1)]
    q = [v for v in q if v is not None]
    if not q:
        return None
    if all(v == 1 for v in q):
        return 1
    if all(v == -1 for v in q):
        return -1
    return 0


def _bernstein(P, b):
    n1, n2 = P.shape
    A = _shiftMat(n1, b[0], b[1] - b[0]) @ P @ _shiftMat(n2, b[2], b[3] - b[2]).T
    return _convMat(n1) @ A @ _convMat(n2).T


def _binom(n1):
    Cb = _binomCache.get(n1)
    if Cb is None:
        Cb = np.zeros((n1, n1))
        for i in range(n1):
            for k in range(i + 1):
                Cb[i, k] = math.comb(i, k)
        _binomCache[n1] = Cb
    return Cb


_expCache = {}


def _shiftMat(n1, x0, h):
    """S[k, i] = C(i, k) x0^(i-k) h^k: the change of variable x = x0 + h u on coefficients."""
    e = _expCache.get(n1)
    if e is None:
        i = np.arange(n1)
        k = np.arange(n1)[:, None]
        E = i - k
        e = (np.maximum(E, 0), E >= 0, _binom(n1).T.copy(), np.arange(n1)[:, None])
        _expCache[n1] = e
    Epos, mask, CbT, kcol = e
    X = np.where(mask, np.power(float(x0), Epos), 0.0)
    return CbT * X * np.power(float(h), kcol)


def _convMat(n1):
    M = _convCache.get(n1)
    if M is None:
        Cb = _binom(n1)
        n = n1 - 1
        M = np.zeros((n1, n1))
        for i in range(n1):
            for k in range(i + 1):
                M[i, k] = Cb[i, k] / Cb[n, k]       # C(i,k) / C(n,k), lower triangle
        _convCache[n1] = M
    return M

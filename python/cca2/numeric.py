"""numeric -- the few double-precision MATLAB primitives the float screens rely on, spelled so
that they round the way MATLAB's do (roots, linspace, round, polyval).
"""
import math

import numpy as np

EPS = 2.0 ** -52


def roots(p):
    """MATLAB roots(p): complex roots of the polynomial p (highest degree first)."""
    p = [float(x) for x in p]
    while p and p[0] == 0:
        p = p[1:]
    if len(p) < 2:
        return []
    return [complex(z) for z in np.roots(p)]


def polyval(p, x):
    v = 0.0
    for c in p:
        v = v * x + c
    return v


def roundHalfAway(x):
    """MATLAB round(x) for a double: halves away from zero. Returns an int."""
    return int(math.floor(abs(x) + 0.5)) * (1 if x >= 0 else -1)


def linspace(d1, d2, n):
    """MATLAB linspace: d1 + (0:n-1).*(d2-d1)./(n-1), endpoints exact."""
    n1 = n - 1
    if n1 <= 0:
        return [float(d2)]
    y = [d1 + (i * (d2 - d1)) / n1 for i in range(n)]
    if d1 == d2:
        y = [float(d1)] * n
    else:
        y[0] = float(d1)
        y[-1] = float(d2)
    return y


def norm2(v):
    """MATLAB vecnorm / norm of a 2-vector: sqrt(sum(x.^2))."""
    return math.sqrt(v[0] * v[0] + v[1] * v[1])

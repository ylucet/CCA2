"""biv -- exact BIVARIATE integer polynomials as dense coefficient matrices (the 2-D arrays of
ratQ.m and algRowFails.m).

M[i][j] multiplies x^i y^j; a matrix always has at least one row and one column, and its SHAPE
is kept as MATLAB keeps it (padding included), because the Bernstein sign test that consumes
these polynomials depends on the degree it is run at. Entries are Python ints.

polyGcdBiv is the primitive polynomial-remainder-sequence gcd of ratQ.polyGcdBiv, with the
polynomial seen as a polynomial in x whose coefficients are polynomials in y.
"""
from math import gcd

from . import rat_q
from .errors import error


def shape(A):
    return len(A), len(A[0])


def zeros(n1, n2):
    return [[0] * n2 for _ in range(n1)]


def isConst(A):
    """MATLAB numel(A) <= 1: a 1x1 matrix."""
    return len(A) == 1 and len(A[0]) == 1


def isZero(A):
    return not any(x != 0 for r in A for x in r)


def trimN(A):
    """Drop trailing all-zero rows and columns, keeping at least 1x1."""
    A = [list(r) for r in A]
    while len(A) > 1 and not any(A[-1]):
        A.pop()
    n2 = len(A[0])
    while n2 > 1 and not any(r[n2 - 1] != 0 for r in A):
        n2 -= 1
    return [r[:n2] for r in A]


def content(A):
    g = 0
    for r in A:
        for x in r:
            if x != 0:
                g = gcd(g, abs(x))
                if g == 1:
                    return 1
    return g


def contentFree(A):
    g = content(A)
    if g > 1:
        return [[x // g for x in r] for r in A]
    return [list(r) for r in A]


def scal(A, k):
    return [[x * k for x in r] for r in A]


def conv2(A, B):
    a1, a2 = shape(A)
    b1, b2 = shape(B)
    C = zeros(a1 + b1 - 1, a2 + b2 - 1)
    for i in range(a1):
        Ai = A[i]
        for j in range(a2):
            v = Ai[j]
            if v == 0:
                continue
            for k in range(b1):
                Bk = B[k]
                Ck = C[i + k]
                for l in range(b2):
                    if Bk[l]:
                        Ck[j + l] += v * Bk[l]
    return C


def padAdd(A, B):
    a1, a2 = shape(A)
    b1, b2 = shape(B)
    S = zeros(max(a1, b1), max(a2, b2))
    for i in range(a1):
        for j in range(a2):
            S[i][j] = A[i][j]
    for i in range(b1):
        for j in range(b2):
            S[i][j] += B[i][j]
    return S


def pad(A, n1, n2):
    a1, a2 = shape(A)
    if (a1, a2) == (n1, n2):
        return [list(r) for r in A]
    S = zeros(n1, n2)
    for i in range(a1):
        for j in range(a2):
            S[i][j] = A[i][j]
    return S


def transpose(A):
    return [list(c) for c in zip(*A)]


def diff(S, d):
    """Derivative in x (d = 1) or y (d = 2)."""
    if d == 2:
        return transpose(diff(transpose(S), 1))
    if len(S) < 2:
        return [[0] * len(S[0])]
    return [[x * i for x in S[i]] for i in range(1, len(S))]


def equal(A, B):
    return shape(A) == shape(B) and A == B


def toDict(A):
    return {(i, j): x for i, r in enumerate(A) for j, x in enumerate(r) if x != 0}


def fromDict(d):
    if not d:
        return [[0]]
    n1 = max(e[0] for e in d) + 1
    n2 = max(e[1] for e in d) + 1
    M = zeros(n1, n2)
    for e, c in d.items():
        M[e[0]][e[1]] = c
    return M


def polyDivExact(P, D):
    """P / D, exactly; raises ratQ:notDivisible when D does not divide P."""
    from . import poly_n
    return fromDict(poly_n.divExact(toDict(P), toDict(D)))


# ---- the gcd ------------------------------------------------------------------------------------
def _trim1(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def _yGcd(a, b):
    a = _trim1(a)
    b = _trim1(b)
    if not any(a):
        return b
    if not any(b):
        return a
    g = rat_q.polyGcd(a[::-1], b[::-1])[::-1]
    return _trim1(g)


def _yDivExact(a, b):
    a = _trim1(a)
    b = _trim1(b)
    if len(b) == 1:
        if b[0] == 0:
            error('ratQ:divideByZero', 'y-division by zero.')
        if any(x % b[0] != 0 for x in a):
            error('ratQ:notDivisible', 'the y-division is not exact.')
        return [x // b[0] for x in a]
    q = rat_q.polyDivExact1(a[::-1], b[::-1])
    if q is None:
        error('ratQ:notDivisible', 'the y-division is not exact.')
    if not q:
        return [0]
    return _trim1(q[::-1])


def _toCell(A):
    A = trimN(A)
    c = [_trim1(r) for r in A]
    while len(c) > 1 and not any(c[-1]):
        c.pop()
    if len(c) == 1 and not any(c[0]):
        return []
    return c


def _fromCell(c):
    if not c:
        return [[0]]
    w = max(len(r) for r in c)
    return trimN([list(r) + [0] * (w - len(r)) for r in c])


def _primitive(c):
    """(c / cont, cont): the content in y of the x-coefficients."""
    cont = [0]
    for r in c:
        cont = _yGcd(cont, r)
    if len(cont) == 1 and abs(cont[0]) == 1:
        return c, [1]
    return [_yDivExact(r, cont) for r in c], cont


def _conv1(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def _pseudoRem(a, b):
    db = len(b) - 1
    r = [list(x) for x in a]
    guard = 0
    while True:
        while len(r) > 1 and not any(r[-1]):
            r.pop()
        if not r or (len(r) == 1 and not any(r[0])):
            return []
        dr = len(r) - 1
        if dr < db:
            return r
        guard += 1
        if guard > 10000:
            error('ratQ:premNoConverge', 'the bivariate pseudo-remainder did not terminate.')
        lb = b[-1]
        lr = r[-1]
        r = [_conv1(x, lb) for x in r]
        sh = dr - db
        for i in range(db + 1):
            t = _conv1(lr, b[i])
            n = max(len(r[sh + i]), len(t))
            u = list(r[sh + i]) + [0] * (n - len(r[sh + i]))
            for k, v in enumerate(t):
                u[k] -= v
            r[sh + i] = _trim1(u)
        r[-1] = [0]                                 # cancels exactly, by construction
        while len(r) > 1 and not any(r[-1]):
            r.pop()
        if not r or (len(r) == 1 and not any(r[0])):
            return []
        r, _ = _primitive(r)                        # the growth control, INSIDE the loop


def polyGcdBiv(A, B):
    """The gcd of two bivariate integer polynomials, content-free (sign not normalised)."""
    A = contentFree(A)
    B = contentFree(B)
    a = _toCell(A)
    b = _toCell(B)
    if not a:
        return B
    if not b:
        return A
    a, ca = _primitive(a)
    b, cb = _primitive(b)
    cg = _yGcd(ca, cb)
    guard = 0
    while True:
        guard += 1
        if guard > 200:
            error('ratQ:gcdNoConverge', 'the bivariate gcd did not terminate.')
        if len(b) == 1:
            return contentFree(_fromCell([cg]))     # coprime in x
        r = _pseudoRem(a, b)
        if not r:
            b = [_trim1(_conv1(x, cg)) for x in b]
            return contentFree(_fromCell(b))
        a = b
        b, _ = _primitive(r)

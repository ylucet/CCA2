"""corner_wedge_rows -- the two rows that bound a CORNER face's wedge (port of cornerWedgeRows.m).

A corner face sits at a dual point s where a facet of the cell meets the tie; its region is the
wedge p + cone(M) = { x : M^-1 (x - p) >= 0 }. With s irrational (a root sigma of the quadratic
Q), each edge of the wedge is  A(x) + sigma * B(x) = 0  with A, B rational affine forms: a
rational line when B = 0 or B is parallel to A, otherwise ONE LINE of the rational pair
al*A^2 - be*A*B + ga*B^2 = 0, named by a branch. Returns (rows, info); rows is empty when an
edge could not be named (info['refused'] says why).
"""
import math
from fractions import Fraction
from math import gcd

from . import rat_q
from .errors import CCA2Error, error
from .isolated_root_sign import isolatedRootSign
from .rat_q import sign


def cornerWedgeRows(f, n, d, tie, side):
    if not isinstance(f, dict) or not all(k in f for k in ('sigma', 'M', 'p', 's')):
        error('cornerWedgeRows:shape', 'f must be a face of cornerFaceAlg.')
    tie = list(tie)
    n = list(n)
    N = [[n[4], n[5]], [n[5], n[6]]]
    m = [n[7], n[8]]
    sgm = f['sigma']
    al, be, ga = sgm['quad']
    Q = [al, be, ga]
    na = sgm['na']
    qadd = rat_q.qsAdd
    qred = rat_q.qsRed

    def qmul(u, v):
        return rat_q.qsMul(u, v, Q)

    def qneg(u):
        return [-u[0], -u[1], u[2]]

    # the corner s = (s0n + sigma*w)/na, in Q(sigma): elements [n1 n2 D] = (n1 + n2 sigma)/D
    s1 = qred([sgm['s0n'][0], sgm['w'][0], na])
    s2 = qred([sgm['s0n'][1], sgm['w'][1], na])
    g1 = qadd(qadd(qmul([2 * tie[0], 0, 1], s1), qmul([tie[1], 0, 1], s2)), [tie[3], 0, 1])
    g2 = qadd(qadd(qmul([tie[1], 0, 1], s1), qmul([2 * tie[2], 0, 1], s2)), [tie[4], 0, 1])
    a = [f['M'][0][0], f['M'][1][0]]                # the facet normal, RATIONAL and integral
    a = [int(a[0]), int(a[1])]
    M = [[[a[0], 0, 1], qmul([side, 0, 1], g1)],
         [[a[1], 0, 1], qmul([side, 0, 1], g2)]]
    p1 = qred(qadd(qadd(qmul([N[0][0], 0, 1], s1), qmul([N[0][1], 0, 1], s2)), [m[0], 0, 1]))
    p2 = qred(qadd(qadd(qmul([N[1][0], 0, 1], s1), qmul([N[1][1], 0, 1], s2)), [m[1], 0, 1]))
    p1 = qred([p1[0], p1[1], p1[2] * d])
    p2 = qred([p2[0], p2[1], p2[2] * d])

    adjM = [[M[1][1], qneg(M[0][1])], [qneg(M[1][0]), M[0][0]]]

    rows = []
    info = {'witness': None, 'kind': [], 'refused': ''}
    I = rat_q.isolateRoots(Q)                       # OUR root: the larger one when root*al > 0
    ivS = I[(len(I) - 1) if (sgm['root'] * sign(al) > 0) else 0]
    Mf = f['M']
    xw = [f['p'][0] - (Mf[0][0] + Mf[0][1]), f['p'][1] - (Mf[1][0] + Mf[1][1])]   # mu = [1;1]
    info['witness'] = xw
    for i in range(2):
        A, B = _edgeForm(adjM[i], [p1, p2], Q)
        r, kind = _rowFor(A, B, Q, xw, ivS)
        if r is None:
            info['refused'] = kind
            return [], info
        rows.append(r)
        info['kind'].append(kind)
    return rows, info


def _lcm2(a, b):
    return abs(a * b) // gcd(abs(a), abs(b))


def _gcdv(v):
    g = 0
    for x in v:
        if x != 0:
            g = gcd(g, abs(x))
    return g


def _edgeForm(adjRow, p, Q):
    """The edge <adjRow, x - p> = 0 as A(x) + sigma*B(x) = 0, A and B integer affine forms."""
    c1, c2 = adjRow
    k = rat_q.qsRed(rat_q.qsAdd(rat_q.qsMul(c1, p[0], Q), rat_q.qsMul(c2, p[1], Q)))
    D = _lcm2(_lcm2(c1[2], c2[2]), k[2])
    A = [-c1[0] * (D // c1[2]), -c2[0] * (D // c2[2]), k[0] * (D // k[2])]
    B = [-c1[1] * (D // c1[2]), -c2[1] * (D // c2[2]), k[1] * (D // k[2])]
    g = _gcdv(A + B)
    if g > 0:
        A = [x // g for x in A]
        B = [x // g for x in B]
    return A, B


def _linMul(A, B):
    return [A[0] * B[0], A[0] * B[1] + A[1] * B[0], A[1] * B[1],
            A[0] * B[2] + A[2] * B[0], A[1] * B[2] + A[2] * B[1], A[2] * B[2]]


def _rowFor(A, B, Q, xw, ivS):
    if all(x == 0 for x in B):
        g = _gcdv(A)
        if g == 0:
            return None, 'degenerate'
        L = [x // g for x in A]
        v = L[0] * xw[0] + L[1] * xw[1] + L[2]
        if abs(v) <= 1e-9 * max(1, max(abs(x) for x in L)):
            return None, 'witnessOnEdge'
        return [0, 0, 0] + L + [sign(v), 0], 'rational'
    if (A[0] * B[1] - A[1] * B[0] == 0 and A[0] * B[2] - A[2] * B[0] == 0
            and A[1] * B[2] - A[2] * B[1] == 0):
        g = _gcdv(B)
        L = [x // g for x in B]
        v = L[0] * xw[0] + L[1] * xw[1] + L[2]
        if abs(v) <= 1e-9 * max(1, max(abs(x) for x in L)):
            return None, 'witnessOnEdge'
        return [0, 0, 0] + L + [sign(v), 0], 'rationalScaled'
    al, be, ga = Q
    AA = _linMul(A, A)
    AB = _linMul(A, B)
    BB = _linMul(B, B)
    C = [al * x - be * y + ga * z for x, y, z in zip(AA, AB, BB)]
    if all(x == 0 for x in C):
        return None, 'degenerate'
    C = rat_q.conic(C)
    try:
        pd = rat_q.pairLineData(C)
    except CCA2Error as err:
        return None, 'pair:' + err.identifier
    b, ok = _branchOf(pd, A, B, Q, ivS)
    if not ok:
        return None, 'branch'
    v = rat_q.pairLineEval(C, b, [xw])[0]
    if abs(v) <= 1e-9 * max(1, max(abs(x) for x in C)):
        return None, 'witnessOnEdge'
    return C + [sign(v), b], pd['kind']


def _branchOf(pd, A, B, Q, ivS):
    """Which line of the pair is OURS: compare our slope (or level), an element of Q(sigma),
    with a rational strictly between the pair's two roots."""
    if pd['kind'] == 'parallel':
        k = next(i for i, x in enumerate(pd['form']) if x != 0)
        fk = pd['form'][k]

        def num(c):
            return [-A[2] * fk * c[1] - c[0] * A[k], -B[2] * fk * c[1] - c[0] * B[k], 1]
        den = [A[k], B[k], 1]
    else:
        sv = _qsSignBig([A[1], B[1], 1], Q, ivS)
        if pd['vertical']:
            return 1 + (1 if sv == 0 else 0), True  # branch 2 is the vertical line
        if sv == 0:
            return 0, False                         # our line vertical, the pair's not: not ours

        def num(c):
            return [-A[0] * c[1] - c[0] * A[1], -B[0] * c[1] - c[0] * B[1], 1]
        den = [A[1], B[1], 1]
    sd = _qsSignBig(den, Q, ivS)
    if sd == 0:
        return 0, False
    mm = _shortBetween(pd['iv'][0][1], pd['iv'][1][0])
    mm = Fraction(mm)
    st = sd * _qsSignBig(num([mm.numerator, mm.denominator]), Q, ivS)     # sign(t - m)
    return 1 + (1 if st > 0 else 0), True


def _qsSignBig(u, Q, iv):
    if u[1] == 0:
        return sign(u[0])
    return isolatedRootSign([u[1], u[0]], Q, iv) * sign(u[2])


def _shortBetween(a, b):
    """A dyadic rational with few bits strictly between a and b (or a itself when a == b)."""
    a = Fraction(a)
    b = Fraction(b)
    if a == b:
        return a
    k0 = math.ceil(math.log2(float(max(abs(a), abs(b)) + 1)))
    for k in range(k0, -1075, -1):
        step = Fraction(2) ** k
        mm = Fraction(math.ceil(a / step)) * step
        if mm == a:
            mm = mm + step
        if mm < b:
            return mm
    return (a + b) / 2

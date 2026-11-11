"""bridge_face_rat -- the bridge face between an AFFINE dual face and one with a RANK-1 Hessian
(port of bridgeFaceRat.m).

There the tie is a parabola with a rational parametrisation, the sup along it is attained at
t*(x) = -B(x)/(2C(x)), and the face is RATIONAL: a quadratic numerator over a LINEAR denominator
in x. Returns (face, info) with face['kind'] == 'rational'.
"""
import math

from . import rat_q
from .errors import error
from .numeric import roots, roundHalfAway
from .rat_q import sign


def bridgeFaceRat(n1, d1, n2, d2, P=None):
    P = [list(r) for r in (P or [])]
    n1 = _checkFace(n1, 'n1')
    n2 = _checkFace(n2, 'n2')
    d1 = _checkDen(d1, 'd1')
    d2 = _checkDen(d2, 'd2')
    if any(x != 0 for x in n1[4:7]):
        error('bridgeFaceRat:notAffine',
              'n1 must be the AFFINE face; its quadratic part is %s.' % (n1[4:7],))
    N2 = [[n2[4], n2[5]], [n2[5], n2[6]]]
    if all(x == 0 for r in N2 for x in r) or rat_q.detExact(N2) != 0:
        error('bridgeFaceRat:notRankOne',
              'n2 must have a RANK-1 Hessian: rank 0 is two affine faces (no curve to sweep) '
              'and rank 2 is the general case, which bridgeFaceAlg answers.')

    d = rat_q.lcm(d1, d2)
    b1 = [n1[7] * (d // d1), n1[8] * (d // d1)]
    k1 = n1[9] * (d // d1)
    M2 = [[x * (d // d2) for x in r] for r in N2]
    b2 = [n2[7] * (d // d2), n2[8] * (d // d2)]
    k2 = n2[9] * (d // d2)
    dl = [b1[0] - b2[0], b1[1] - b2[1]]
    kp = k1 - k2

    if M2[0][0] != 0:
        u = [M2[0][0], M2[0][1]]
        a = M2[0][0]
    else:
        u = [M2[0][1], M2[1][1]]
        a = M2[1][1]
    if a == 0:
        error('bridgeFaceRat:internal',
              'a rank-1 symmetric matrix has a nonzero diagonal entry; this one appears not to.')
    w = [-u[1], u[0]]                               # along the parabola's axis: u'w = 0
    nu = u[0] * u[0] + u[1] * u[1]
    dw = dl[0] * w[0] + dl[1] * w[1]
    if dw == 0:
        error('bridgeFaceRat:noSeparateFace',
              'the tie is a pair of PARALLEL LINES and this face is AFFINE, so the pair\'s sup '
              'is a linear program on each side of the tie and is attained at a VERTEX -- the '
              'tie acts as an ordinary FACET, and there is no separate bridge face to build.')
    rq = nu
    rl = -2 * a * (dl[0] * u[0] + dl[1] * u[1])
    rc = -2 * a * kp * nu
    rd = 2 * a * dw

    face, info = _buildFace(u, w, nu, rq, rl, rc, rd, b1, k1, d, P)
    info['parabola'] = {'u': u, 'w': w, 'rq': rq, 'rl': rl, 'rc': rc, 'rd': rd}
    return face, info


def _buildFace(u, w, nu, rq, rl, rc, rd, b1, k1, d, P):
    uy = [d * u[0], d * u[1], -(u[0] * b1[0] + u[1] * b1[1])]
    wy = [d * w[0], d * w[1], -(w[0] * b1[0] + w[1] * b1[1])]
    C = [rq * x for x in wy]                        # from r(t)'s t^2 term
    B = [rd * x + rl * y for x, y in zip(uy, wy)]
    A = [rc * x for x in wy]
    A[2] -= nu * rd * k1
    scale = nu * d * rd

    num = [4 * x - y for x, y in zip(polyMulAff(A, C), polyMulAff(B, B))]   # quadratic in x
    den = [4 * scale * x for x in C]                                         # LINEAR in x
    face = {'kind': 'rational'}
    face['fN'] = _quadToWeighted(num)
    face['fD'] = 1
    face['den'] = den
    face['pair'] = {'kind': 'rational', 'u': u, 'w': w, 'nu': nu, 'rq': rq, 'rl': rl, 'rc': rc,
                    'rd': rd, 'B': B, 'C': C, 'scale': scale, 'P': P}
    face['wit'] = _witnessRat(face['pair'])
    face['cells'] = [_regionRows(face['pair'])]
    face['FP'] = _implicitValue(face['fN'], face['den'])
    face['FB'] = -1
    info = {'maxCoeff': max(abs(x) for x in face['fN'] + face['den']), 'scale': scale}
    return face, info


def _witnessRat(pr):
    P = pr['P']
    for t in _candidateParams(pr):
        s = _paraPoint(pr, t)
        if not (math.isfinite(s[0]) and math.isfinite(s[1])):
            continue
        if P and any(r[0] * s[0] + r[1] * s[1] + r[2] <= 0 for r in P):
            continue
        L = [pr['B'][i] + 2 * t * pr['C'][i] for i in range(3)]   # the line t*(x) = t
        g = L[0:2]
        if g[0] == 0 and g[1] == 0:
            continue
        gg = g[0] * g[0] + g[1] * g[1]
        x0 = [-L[2] * g[0] / gg, -L[2] * g[1] / gg]
        dr = [-g[1], g[0]]
        c0 = pr['C'][0] * x0[0] + pr['C'][1] * x0[1] + pr['C'][2]
        cd = pr['C'][0] * dr[0] + pr['C'][1] * dr[1]
        want = -sign(pr['scale']) * max(1.0, abs(c0))   # strictly on the concave side
        if cd != 0:
            k = (want - c0) / cd
            x = [x0[0] + k * dr[0], x0[1] + k * dr[1]]
        elif c0 * pr['scale'] < 0:
            x = x0
        else:
            continue
        for e in range(31):
            xd = 2 ** e
            xn = [roundHalfAway(x[0] * xd), roundHalfAway(x[1] * xd)]
            if abs(xn[0]) >= 2 ** 52 or abs(xn[1]) >= 2 ** 52:
                break
            if not _memberRatAt(pr, xn, xd):
                continue
            return [xn[0], xn[1], xd]
    error('bridgeFaceRat:noWitness',
          'no rational point of this face was found on the search grid. Either the pair admits '
          'no part of the parabola, or the face is thinner than the grid -- which the caller '
          'should be told rather than have guessed at.')


def _candidateParams(pr):
    rts = []
    for r in pr['P']:
        a = r[0:2]
        aw = a[0] * pr['w'][0] + a[1] * pr['w'][1]
        au = a[0] * pr['u'][0] + a[1] * pr['u'][1]
        c2 = aw * pr['rq'] / (pr['rd'] * pr['nu'])
        c1 = (au + aw * pr['rl'] / pr['rd']) / pr['nu']
        c0 = aw * pr['rc'] / (pr['rd'] * pr['nu']) + r[2]
        rts += [z.real for z in roots([c2, c1, c0]) if abs(z.imag) < 1e-12 * max(1.0, abs(z))]
    rts.sort()
    if not rts:
        T = [0.0, 1.0, -1.0]
    else:
        span = max(1.0, rts[-1] - rts[0])
        T = [(a + b) / 2 for a, b in zip(rts[:-1], rts[1:])] + [rts[0] - span, rts[-1] + span]
    return T + [0.0, 1.0, -1.0, 2.0, -2.0, 0.5, -0.5, 4.0, -4.0, 10.0, -10.0]


def _memberRatAt(pr, xn, xd):
    Bv = (pr['B'][0] * xn[0] + pr['B'][1] * xn[1]) / xd + pr['B'][2]
    Cv = (pr['C'][0] * xn[0] + pr['C'][1] * xn[1]) / xd + pr['C'][2]
    if Cv * pr['scale'] >= 0:
        return False                                # not concave along the parabola: no maximum
    t = -Bv / (2 * Cv)
    s = _paraPoint(pr, t)
    return not pr['P'] or all(r[0] * s[0] + r[1] * s[1] + r[2] >= 0 for r in pr['P'])


def _paraPoint(pr, t):
    r = (pr['rq'] * t * t + pr['rl'] * t + pr['rc']) / pr['rd']
    return [(t * pr['u'][0] + r * pr['w'][0]) / pr['nu'],
            (t * pr['u'][1] + r * pr['w'][1]) / pr['nu']]


def _regionRows(pr):
    """The face's region as conic rows: C*scale <= 0, and each row of P pulled back to x."""
    R = []
    s = -sign(pr['scale'])
    _addRow(R, [0, 0, 0] + [s * x for x in pr['C']])
    rs = sign(pr['rd'])                             # nu > 0, so this is rd*nu's sign
    BB = polyMulAff(pr['B'], pr['B'])
    BC = polyMulAff(pr['B'], pr['C'])
    CC = polyMulAff(pr['C'], pr['C'])
    for r in pr['P']:
        a = r[0:2]
        c = r[2]
        aw = a[0] * pr['w'][0] + a[1] * pr['w'][1]
        au = a[0] * pr['u'][0] + a[1] * pr['u'][1]
        al = aw * pr['rq']
        be = au * pr['rd'] + aw * pr['rl']
        ga = aw * pr['rc'] + c * pr['rd'] * pr['nu']
        raw = [al * x - 2 * be * y + 4 * ga * z for x, y, z in zip(BB, BC, CC)]
        _addRow(R, [rs * x for x in raw])
    return R


def _addRow(R, raw):
    if all(x == 0 for x in raw[0:5]):
        return
    nz = next(x for x in raw if x != 0)
    R.append(rat_q.conic(raw) + [sign(nz), 0])


def _implicitValue(fN, den):
    """The face z = num(x)/den(x) as the implicit polynomial den*z - num in (x, y, z)."""
    from fractions import Fraction
    q = [Fraction(fN[4], 2), Fraction(fN[5]), Fraction(fN[6], 2), Fraction(fN[7]),
         Fraction(fN[8]), Fraction(fN[9])]
    den = list(den)
    if any(x.denominator != 1 for x in q):
        q = [2 * x for x in q]
        den = [2 * x for x in den]
    q = [int(x) for x in q]
    FP = {}

    def put(e, c):
        if c != 0:
            FP[e] = c
    put((2, 0, 0), -q[0])
    put((1, 1, 0), -q[1])
    put((0, 2, 0), -q[2])
    put((1, 0, 0), -q[3])
    put((0, 1, 0), -q[4])
    put((0, 0, 0), -q[5])
    put((1, 0, 1), den[0])
    put((0, 1, 1), den[1])
    put((0, 0, 1), den[2])
    return FP


def polyMulAff(f, g):
    return [f[0] * g[0], f[0] * g[1] + f[1] * g[0], f[1] * g[1],
            f[0] * g[2] + f[2] * g[0], f[1] * g[2] + f[2] * g[1], f[2] * g[2]]


def _quadToWeighted(q):
    return [0, 0, 0, 0, 2 * q[0], q[1], 2 * q[2], q[3], q[4], q[5]]


def _checkFace(n, name):
    n = list(n)
    if len(n) != 10:
        error('bridgeFaceRat:shape', '%s must be a 1x10 face numerator.' % name)
    rat_q.chk(n, name)
    if any(x != 0 for x in n[0:4]):
        error('bridgeFaceRat:cubicFace',
              '%s has a CUBIC part; a conjugate face is quadratic.' % name)
    return n


def _checkDen(d, name):
    if not isinstance(d, int) or isinstance(d, bool) or d <= 0:
        error('bridgeFaceRat:denominator', '%s must be a positive integer; got %r.' % (name, d))
    return d

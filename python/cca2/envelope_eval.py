"""envelope_eval -- evaluate an assembled envelope (port of envelopeEval.m).

env = {'quad': [...], 'bridge': [...]} is what envelopeAssemble returns for ONE atom pair: the
value at x is the MAX over the pair's cell pairs of each cell pair's candidates whose region
holds x. A point no candidate of some cell pair covers gets +inf (after nudging it by 1e-7 in
eight directions, to step off a shared boundary). Doubles: evaluation is where a tolerance is
allowed. Returns (val, who): who[i] = {'kind': 'quad' | 'bridge' | 'none', 'idx': 1-based}.
"""
import math

from . import rat_q
from .alg_branch import algBranch
from .bridge_member import bridgeMember
from .conic_meet import _solve2
from .errors import error

SQEPS = math.sqrt(2.0 ** -52)
_R2 = math.sqrt(2.0)
_DIRS = [[1, 0], [-1, 0], [0, 1], [0, -1],
         [1 / _R2, 1 / _R2], [1 / _R2, -1 / _R2], [-1 / _R2, 1 / _R2], [-1 / _R2, -1 / _R2]]


def _pairsOf(q):
    p = q['pair']
    return p if isinstance(p, list) else [p]


def envelopeEval(env, X):
    if not isinstance(env, dict) or 'quad' not in env or 'bridge' not in env:
        error('envelopeEval:shape',
              'env must be what envelopeAssemble returns, with .quad and .bridge.')
    n = len(X)
    val = [math.inf] * n
    who = [{'kind': 'none', 'idx': 0} for _ in range(n)]
    if not env['quad'] and not env['bridge']:
        return val, who
    pairs = set()
    for q in env['quad']:
        pairs.update(_pairsOf(q))
    for f in env['bridge']:
        pairs.add(f['pairIndex'])
    pairs = sorted(pairs)

    for i in range(n):
        x = [float(X[i][0]), float(X[i][1])]
        best = -math.inf
        bestKind = 'none'
        bestIdx = 0
        diverges = False
        for p in pairs:
            pv, pk, pi, cover = _pairValue(env, p, x)
            if not cover:
                h = 1e-7 * max(1.0, abs(x[0]), abs(x[1]))
                for d in _DIRS:
                    _, nk, ni, nc = _pairValue(env, p, [x[0] + h * d[0], x[1] + h * d[1]])
                    if nc:
                        v = _faceValueAt(env, nk, ni, x)
                        if math.isfinite(v):
                            pv, pk, pi, cover = v, nk, ni, True
                        break                       # MATLAB drops the point from `miss` here
            if not cover:
                diverges = True
                break                               # another pair cannot bring it back
            if pv > best:
                best, bestKind, bestIdx = pv, pk, pi
        if not diverges and math.isfinite(best):
            val[i] = best
            who[i] = {'kind': bestKind, 'idx': bestIdx}
    return val, who


def _pairValue(env, p, x):
    cover = False
    pv = -math.inf
    pk = 'none'
    pi = 0
    for k, f in enumerate(env['quad']):
        if p not in _pairsOf(f):
            continue
        if not _inRegion(f['con'], x):
            continue
        v = rat_q.evalFace(f['num'], f['den'], [x])[0]
        if v > pv:
            pv, pk, pi = v, 'quad', k + 1
        cover = True
    for k, face in enumerate(env['bridge']):
        if face['pairIndex'] != p:
            continue
        if not _faceOwns(face, x):
            continue
        v = _bridgeValue(face, x)
        if not math.isfinite(v):
            continue
        cover = True
        if v > pv:
            pv, pk, pi = v, 'bridge', k + 1
    return pv, pk, pi, cover


def _faceValueAt(env, kind, k, x):
    if kind == 'quad':
        f = env['quad'][k - 1]
        return rat_q.evalFace(f['num'], f['den'], [x])[0]
    return _bridgeValue(env['bridge'][k - 1], x)


def _faceOwns(face, x):
    if face['kind'] == 'corner':
        mu = _solve2(face['M'], [face['p'][0] - x[0], face['p'][1] - x[1]])
        if not (math.isfinite(mu[0]) and math.isfinite(mu[1])):
            return False
        lim = -1e-9 * max(1.0, abs(mu[0]), abs(mu[1]))
        return mu[0] >= lim and mu[1] >= lim
    return bridgeMember(face['pair'], [x])[0][0]


def bridgeValue(face, x):
    return _bridgeValue(face, x)


def _bridgeValue(face, x):
    k = face['kind']
    if k == 'corner':
        n = face['src']['n']
        s = face['s']
        g = ((s[0] * (n[4] * s[0] + n[5] * s[1]) + s[1] * (n[5] * s[0] + n[6] * s[1])) / 2
             + n[7] * s[0] + n[8] * s[1] + n[9]) / face['src']['d']
        return x[0] * s[0] + x[1] * s[1] - g
    if k == 'algebraic':
        return algBranch(face['FP'], [x], face['FB'])[0]
    if k == 'rational':
        q = face['fN']
        num = (q[4] / 2 * x[0] ** 2 + q[5] * x[0] * x[1] + q[6] / 2 * x[1] ** 2
               + q[7] * x[0] + q[8] * x[1] + q[9])
        den = face['den'][0] * x[0] + face['den'][1] * x[1] + face['den'][2]
        if den == 0:
            return math.nan
        return num / den / face['fD']
    error('envelopeEval:faceKind',
          'a curved-tie face must be \'algebraic\', \'rational\' or \'corner\'; got \'%s\'.' % k)


def _inRegion(con, x):
    xmax = max(1.0, abs(x[0]), abs(x[1]))
    for r in con:
        c = r[0:6]
        if len(r) >= 8 and r[7] != 0:
            q = rat_q.pairLineEval(c, r[7], [x])[0]
            scale = xmax
        else:
            q = (c[0] * x[0] ** 2 + c[1] * x[0] * x[1] + c[2] * x[1] ** 2
                 + c[3] * x[0] + c[4] * x[1] + c[5])
            scale = max(1.0, max(abs(v) for v in c) * xmax * xmax)
        if not (r[6] * q >= -SQEPS * scale):
            return False
    return True

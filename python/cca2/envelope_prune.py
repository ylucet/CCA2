"""envelope_prune -- drop the quadratic candidates of an assembled envelope that can never be the
max: exact duplicates, and candidates dominated on their whole (polyhedral) region by another
(port of envelopePrune.m). Returns (env, info) -- callers that want only env index [0]... see
`envelopePrune` below, which returns env alone, as every MATLAB caller uses it.
"""
from . import rat_q
from .errors import error


def envelopePruneInfo(env):
    if not isinstance(env, dict) or 'quad' not in env or 'bridge' not in env:
        error('envelopePrune:shape', 'env must be what envelopeAssemble returns.')
    env = {'quad': [dict(q) for q in env['quad']], 'bridge': env['bridge']}
    quad = env['quad']
    n0 = len(quad)
    keep = [True] * n0
    dup = 0
    for i in range(n0):
        if not keep[i]:
            continue
        for j in range(i + 1, n0):
            if not keep[j]:
                continue
            if quad[i]['num'] == quad[j]['num'] and quad[i]['den'] == quad[j]['den'] and \
                    quad[i]['con'] == quad[j]['con']:
                quad[i]['pair'] = _union(quad[i]['pair'], quad[j]['pair'])
                keep[j] = False
                dup += 1
    dom = 0
    idx = [i for i in range(n0) if keep[i]]
    for i in idx:
        if not keep[i]:
            continue
        Li = _linearRows(quad[i]['con'])
        if not Li:
            continue                               # no linear rows: nothing can be decided here
        for j in idx:
            if j == i or not keep[j] or not keep[i]:
                continue
            if not _containsRegion(Li, quad[j]['con']):
                continue
            if not _atLeastAsLarge(quad[j], quad[i], Li):
                continue
            quad[j]['pair'] = _union(quad[j]['pair'], quad[i]['pair'])
            keep[i] = False
            dom += 1
    env['quad'] = [q for q, k in zip(quad, keep) if k]
    return env, {'before': n0, 'after': len(env['quad']), 'duplicates': dup, 'dominated': dom}


def envelopePrune(env):
    return envelopePruneInfo(env)[0]


def _union(a, b):
    a = a if isinstance(a, list) else [a]
    b = b if isinstance(b, list) else [b]
    return sorted(set(a) | set(b))


def _linearRows(con):
    L = []
    for r in con:
        if any(x != 0 for x in r[0:3]):
            continue
        if len(r) >= 8 and r[7] != 0:
            continue                               # one line of a pair: not a half-plane
        L.append([r[6] * x for x in r[3:6]])
    return L


def _containsRegion(Li, conj):
    for r in conj:
        if any(x != 0 for x in r[0:3]):
            return False                           # a curved row: not decided here
        if len(r) >= 8 and r[7] != 0:
            return False
        row = [r[6] * x for x in r[3:6]]
        if rat_q.feasible2(Li + [[-x for x in row]], True):
            return False                           # j's row can fail inside R_i
    return True


def _atLeastAsLarge(fj, fi, Li):
    num = [a * fi['den'] - b * fj['den'] for a, b in zip(fj['num'], fi['num'])]
    if any(x != 0 for x in num[0:4]):
        return False                               # a cubic difference: outside this test's reach
    c = [num[4], 2 * num[5], num[6], 2 * num[7], 2 * num[8], 2 * num[9]]
    if all(x == 0 for x in c):
        return True                                # the two faces are the same function
    nz = next(x for x in c if x != 0)
    s, _ = rat_q.quadSignOnPoly(rat_q.conic(c), Li)
    return s * rat_q.sign(nz) >= 1                 # +1 means v_j - v_i >= 0 on the whole region

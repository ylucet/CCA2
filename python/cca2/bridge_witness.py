"""bridge_witness -- VERIFIED rational points in the interior of a bridge face's region
(port of bridgeWitness.m).

Proposals come from a double-precision sample of the tie curve (the midpoint of the two
gradients at a tie point lies inside the sweep); each is rounded to a coarse rational and then
VERIFIED exactly: Phi has a root mu in [-1, 0] whose base point s(mu) lies in the cell pair P
(signs at a real algebraic number, by Sturm). Returns (W, info): W rows [xn1, xn2, xd], coarsest
first.
"""
import math

import numpy as np

from . import rat_q
from .bridge_phi import bridgePhi, vadd
from .errors import error
from .numeric import EPS, linspace, roundHalfAway


def bridgeWitness(d, N1, m1, N2, m2, Dl, dl, kp, C1, C2, P):
    P = [list(r) for r in (P or [])]
    S = _tieSample(Dl, dl, kp, P)
    tries = 0
    W = []
    nRoots = []
    for s in S:
        x0 = [((N1[0][0] * s[0] + N1[0][1] * s[1] + m1[0])
               + (N2[0][0] * s[0] + N2[0][1] * s[1] + m2[0])) / (2 * d),
              ((N1[1][0] * s[0] + N1[1][1] * s[1] + m1[1])
               + (N2[1][0] * s[0] + N2[1][1] * s[1] + m2[1])) / (2 * d)]
        if not (math.isfinite(x0[0]) and math.isfinite(x0[1])):
            continue
        for xdTry in (1, 2, 4, 8, 16, 32, 64):
            tries += 1
            xnTry = [roundHalfAway(x0[0] * xdTry), roundHalfAway(x0[1] * xdTry)]
            if (C1 and signAtRational(C1, xnTry, xdTry) == 0) or \
               (C2 and signAtRational(C2, xnTry, xdTry) == 0):
                continue
            if not P:
                inR, nR = memberExact(N1, m1, N2, m2, Dl, dl, kp, d, xnTry, xdTry)
            else:
                inR, nR = memberExactArc(N1, m1, N2, m2, Dl, dl, kp, d, xnTry, xdTry, P)
            if inR:
                row = [xnTry[0], xnTry[1], xdTry]
                if row not in W:
                    W.append(row)
                    nRoots.append(nR)
                break                               # this tie point has answered; take the next
        if len(W) >= 12:
            break

    if W:
        # coarsest, then smallest (MATLAB sortrows: stable)
        order = sorted(range(len(W)), key=lambda i: (W[i][2], max(abs(W[i][0]), abs(W[i][1]))))
        W = [W[i] for i in order]
        info = {'nRoots': nRoots[order[0]], 'tries': tries, 'overflows': 0,
                'nWitnesses': len(W)}
        return W, info

    if P and S and len(S) <= 2:
        error('bridgeFaceAlg:tangentialArc',
              'the tie meets this pair in at most a POINT (%d sampled tie points lie in it), so '
              'the pair carries no two-dimensional bridge face. That is the answer for a '
              'TANGENTIAL tie, not a refusal to compute one.' % len(S))
    error('bridgeFaceAlg:noWitness',
          'no interior point of the bridge region could be VERIFIED, after %d attempts over %d '
          'tie points. Either the tie curve has no real points -- the two quadratics never tie, '
          'so the pair carries no bridge and should not have reached here -- or the region is '
          'thinner than the coarsest rational grid this routine will use, which the caller '
          'should be told rather than have guessed at.' % (tries, len(S)))


def _inP(s, P):
    return all(s[0] * r[0] + s[1] * r[1] + r[2] >= 0 for r in P)


def _tieSample(Dl, dl, kp, P):
    """Double-precision points of the tie curve inside P, the least extreme first."""
    if all(x == 0 for r in Dl for x in r):
        D, E = dl[0], dl[1]
        if D == 0 and E == 0:
            return []
        w = [float(E), float(-D)]
        nn = float(D * D + E * E)
        s0 = [-kp * D / nn, -kp * E / nn]
        lo, hi = _arcExtent(s0, w, P)
        if math.isnan(lo):
            return []
        taus = linspace(lo, hi, 43)[1:-1]           # the ENDS are the cut boundaries; stay inside
        S = [[s0[0] + w[0] * t, s0[1] + w[1] * t] for t in taus]
        return [s for s in S if math.isfinite(s[0]) and math.isfinite(s[1])]

    Df = np.array(Dl, dtype=float)
    df = np.array(dl, dtype=float)
    if rat_q.detExact(Dl) != 0:
        c0 = -np.linalg.solve(Df, df)
        v0 = float(c0 @ Df @ c0) / 2 + float(df @ c0) + kp
        ev = float(np.max(np.abs(np.linalg.eigvalsh(Df))))
        L = 4 * math.sqrt(abs(2 * v0) / max(ev, EPS)) + 4
        ctr = float(c0[0])
    else:
        L = 20.0
        ctr = 0.0

    S = []
    for L2 in (L, 4 * L, 20 * L, 100 * L, 500 * L):
        S = _conicPoints(Dl, dl, kp, ctr, L2, 2001)
        if not P or not S:
            break
        keep = [s for s in S if _inP(s, P)]
        if not keep:
            continue
        S = keep
        if len(S) >= 5:
            lo = min(s[0] for s in S)
            hi = max(s[0] for s in S)
            pad = max(1e-9, (hi - lo) / 20)
            S2 = _conicPoints(Dl, dl, kp, (lo + hi) / 2, (hi - lo) / 2 + pad, 2001)
            if S2:
                k2 = [s for s in S2 if _inP(s, P)]
                if k2:
                    S = k2
        break
    if S:                                           # middle first: the least extreme proposals
        mx = sum(s[0] for s in S) / len(S)
        my = sum(s[1] for s in S) / len(S)
        S = sorted(S, key=lambda s: math.sqrt((s[0] - mx) ** 2 + (s[1] - my) ** 2))
    return S


def _arcExtent(s0, w, P):
    lo = -math.inf
    hi = math.inf
    for r in P:
        A = r[0] * w[0] + r[1] * w[1]
        B = r[0] * s0[0] + r[1] * s0[1] + r[2]
        if abs(A) <= 0:
            if B < 0:
                return math.nan, math.nan
            continue
        if A > 0:
            lo = max(lo, -B / A)
        else:
            hi = min(hi, -B / A)
    if lo >= hi:
        return math.nan, math.nan
    span = 24
    if math.isinf(lo):
        lo = -12 if math.isinf(hi) else hi - span
    if math.isinf(hi):
        hi = 12 if math.isinf(lo) else lo + span
    return lo, hi


def _conicPoints(Dl, dl, kp, ctr, L, n):
    S = []
    for t in linspace(-L, L, n):
        s1 = ctr + t
        a = Dl[1][1]
        b = 2 * Dl[0][1] * s1 + 2 * dl[1]
        c = Dl[0][0] * s1 * s1 + 2 * dl[0] * s1 + 2 * kp
        if abs(a) < EPS * max(1.0, abs(b)):
            if abs(b) > 0:
                S.append([s1, -c / b])
            continue
        disc = b * b - 4 * a * c
        if disc < 0:
            continue
        r = math.sqrt(disc)
        S.append([s1, (-b + r) / (2 * a)])
        if r > 0:
            S.append([s1, (-b - r) / (2 * a)])
    return [s for s in S if math.isfinite(s[0]) and math.isfinite(s[1])
            and abs(s[0]) < 1e6 and abs(s[1]) < 1e6]


def memberExactArc(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd, P):
    """Is x = xn/xd in the sweep of the arc of the tie inside P? Exact. Returns (inR, n)."""
    Phi, Dt, U1, U2 = bridgePhi(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd)
    Ph = rat_q.polyTrim(Phi[::-1])
    if len(Ph) <= 1:
        return False, 0
    I = rat_q.isolateRoots(Ph)
    n = 0
    for iv in I:
        lo = max(iv[0], -1)
        hi = min(iv[1], 0)
        if lo > hi:
            continue
        if rat_q.countRootsIn(Ph, lo, hi) != 1 and not (rat_q.polyval(Ph, lo) == 0):
            continue                                # that root is not in [-1,0]
        n += 1
        sD = rat_q.signAt(Dt[::-1], Ph, iv)
        if sD == 0:
            continue                                # det M vanishes there: s(mu) is not a point
        ok = True
        for r in P:
            q = vadd(vadd([r[0] * x for x in U1], [r[1] * x for x in U2]),
                     [r[2] * xd * x for x in Dt])
            sq = rat_q.signAt(q[::-1], Ph, iv)
            if sq * sD <= 0:
                ok = False
                break
        if ok:
            return True, n
    return False, n


def memberExact(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd):
    Phi, _, _, _ = bridgePhi(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd)
    Ph = rat_q.polyTrim(Phi[::-1])                  # highest first
    if len(Ph) <= 1:
        return False, 0
    n = rat_q.countRootsIn(Ph, -1, 0)               # half-open (-1, 0]
    return (n > 0 or rat_q.polyval(Ph, -1) == 0), n


def signAtRational(C, xn, xd):
    v = (C[0] * xn[0] ** 2 + C[1] * xn[0] * xn[1] + C[2] * xn[1] ** 2
         + C[3] * xn[0] * xd + C[4] * xn[1] * xd + C[5] * xd ** 2)
    return rat_q.sign(v)

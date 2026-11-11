"""bridge_member -- is a (double) point in a bridge face's region? (port of bridgeMember.m)

The eval-time PREDICATE of a bridge face: x is in the region when the tie residual Phi(mu) has a
root mu in [-1, 0] whose base point lies in the cell pair. Doubles, with tolerances -- evaluation
is the one place a tolerance is allowed. Returns (tf, info) for a list of points X.
"""
import math

from .bridge_phi import bridgePhi
from .errors import error
from .numeric import EPS, polyval, roots


def bridgeMember(pr, X):
    n = len(X)
    if isinstance(pr, dict) and pr.get('kind') == 'rational':
        tf = []
        for x in X:
            Bv = x[0] * pr['B'][0] + x[1] * pr['B'][1] + pr['B'][2]
            Cv = x[0] * pr['C'][0] + x[1] * pr['C'][1] + pr['C'][2]
            ok = Cv * pr['scale'] < 0
            if ok and pr['P']:
                t = -Bv / (2 * Cv)
                r = (pr['rq'] * t * t + pr['rl'] * t + pr['rc']) / pr['rd']
                S = [(t * pr['u'][0] + r * pr['w'][0]) / pr['nu'],
                     (t * pr['u'][1] + r * pr['w'][1]) / pr['nu']]
                tol = math.sqrt(EPS) * max(1.0, abs(S[0]), abs(S[1]))
                ok = all(S[0] * p[0] + S[1] * p[1] + p[2] >= -tol for p in pr['P'])
            tf.append(bool(ok))
        return tf, {'mu': [math.nan] * n, 'base': [[math.nan, math.nan]] * n}
    if not isinstance(pr, dict) or not all(k in pr for k in
                                           ('d', 'N1', 'm1', 'N2', 'm2', 'Dl', 'dl', 'kp', 'P')):
        error('bridgeMember:pair',
              'pr must carry the cell pair\'s own data (d, N1, m1, N2, m2, Dl, dl, kp, P) -- it '
              'is what `bridgeFaceAlg` puts in `face.pair`.')
    tf = [False] * n
    mus = [math.nan] * n
    base = [[math.nan, math.nan] for _ in range(n)]
    P = pr['P'] or []
    for i, x in enumerate(X):
        xf = [float(x[0]), float(x[1])]
        Phi, Dt, U1, U2 = bridgePhi(pr['N1'], pr['m1'], pr['N2'], pr['m2'], pr['Dl'], pr['dl'],
                                    pr['kp'], pr['d'], xf, 1, False)
        p = [float(v) for v in Phi[::-1]]           # `roots` wants highest power first
        mx = max(abs(v) for v in p)
        nz = next((k for k, v in enumerate(p) if abs(v) > EPS * mx), None)
        if nz is None:
            continue
        p = p[nz:]
        if len(p) < 2:
            continue
        r = [z.real for z in roots(p) if abs(z.imag) <= 1e-9 * max(1.0, abs(z))]
        r = [v for v in r if v >= -1 - 1e-12 and v <= 1e-12]
        dScale = max(abs(v) for v in Dt)
        DtR = [float(v) for v in Dt[::-1]]
        U1R = [float(v) for v in U1[::-1]]
        U2R = [float(v) for v in U2[::-1]]
        for rj in r:
            mu = min(max(rj, -1.0), 0.0)
            D = polyval(DtR, mu)
            if abs(D) <= 1e-9 * max(1.0, dScale):
                continue
            s = [polyval(U1R, mu) / D, polyval(U2R, mu) / D]
            if P:
                lim = -1e-9 * max(1.0, abs(s[0]), abs(s[1]))
                if any(q[0] * s[0] + q[1] * s[1] + q[2] < lim for q in P):
                    continue
            tf[i] = True
            mus[i] = mu
            base[i] = s
            break
        if tf[i] or not any(abs(v + 1) <= 1e-6 for v in r) or polyval(DtR, -1.0) != 0:
            continue
        ok, mu, s = _nearArtefact(pr, xf)
        if ok:
            tf[i] = True
            mus[i] = mu
            base[i] = s
    return tf, {'mu': mus, 'base': base}


def _nearArtefact(pr, x):
    """Near mu = -1 with det M(-1) = 0 the family is re-based at g2 (nu = mu + 1)."""
    N1 = pr['N2']
    N2 = [[2 * a - b for a, b in zip(ra, rb)] for ra, rb in zip(pr['N2'], pr['N1'])]
    m1 = pr['m2']
    m2 = [2 * a - b for a, b in zip(pr['m2'], pr['m1'])]
    Phi, Dt, U1, U2 = bridgePhi(N1, m1, N2, m2, pr['Dl'], pr['dl'], pr['kp'], pr['d'], x, 1, False)
    p = [float(v) for v in Phi[::-1]]
    nz = next((k for k, v in enumerate(p) if v != 0), None)
    if nz is None:
        return False, math.nan, [math.nan, math.nan]
    p = p[nz:]
    if len(p) < 2:
        return False, math.nan, [math.nan, math.nan]
    r = [z.real for z in roots(p) if abs(z.imag) <= 1e-9 * max(1e-12, abs(z))]
    r = [v for v in r if v >= -1e-12 and v <= 1 + 1e-12 and abs(v) <= 1e-6]
    DtR = [float(v) for v in Dt[::-1]]
    U1R = [float(v) for v in U1[::-1]]
    U2R = [float(v) for v in U2[::-1]]
    P = pr['P'] or []
    for rj in r:
        D = polyval(DtR, rj)
        if D == 0:
            continue
        sj = [polyval(U1R, rj) / D, polyval(U2R, rj) / D]
        if not (math.isfinite(sj[0]) and math.isfinite(sj[1])):
            continue
        if P:
            lim = -1e-9 * max(1.0, abs(sj[0]), abs(sj[1]))
            if any(q[0] * sj[0] + q[1] * sj[1] + q[2] < lim for q in P):
                continue
        return True, rj - 1, sj
    return False, math.nan, [math.nan, math.nan]

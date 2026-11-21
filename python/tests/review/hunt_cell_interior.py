"""Does cellHasInterior's 'empty' verdict ever hit a cell that HAS interior inside a real conj call?

Wraps cca2.cell_has_interior.cellHasInterior (in the two modules that import it) with a second,
denser witness search: the engine's own CAD with the missing vertical lines added as critical
abscissae, plus a fine grid per sector. Any cell the engine calls empty but the denser search
finds a STRICT interior point for (every curved row > 0 with a margin, every line row > 0) is
logged with the conj input, then the conj value is checked against plq_oracle at dual points.

    python3 hunt_cell_interior.py N_RANDOM [seed [first last]]   (first <= it < last: one shard of the same N inputs)
"""
import json
import random
import sys
import os
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', '..')))
sys.path.insert(0, HERE)
import numpy as np                                   # noqa: E402
from cca2 import cell_has_interior as CHI           # noqa: E402
from cca2 import assemble_qua_con_cells as AQ       # noqa: E402
from cca2 import cell_interior as CI                # noqa: E402
from cca2.qua_pol import QuaPol                     # noqa: E402
import plq_gen as G                                  # noqa: E402
import plq_oracle as O                               # noqa: E402

ORIG = CHI.cellHasInterior
HITS = []


def strict_point(EcQ, rows):
    """A point where every row holds STRICTLY (exact rational check), searched on a CAD that
    includes vertical lines and every row's x-critical values, with several samples per sector."""
    if any(r[1] == 0 for r in rows):
        return None
    idx = sorted(set(r[0] for r in rows))
    cs = [[float(t) for t in EcQ[i - 1][0:6]] for i in idx]
    xs = []
    for a0, a1, a2, a3, a4, a5 in cs:
        if a2 != 0:
            q = [a1 * a1 - 4 * a2 * a0, 2 * a1 * a4 - 4 * a2 * a3, a4 * a4 - 4 * a2 * a5]
            r = np.roots(q) if any(q) else []
            xs += [z.real for z in r if abs(z.imag) < 1e-9]
        elif a1 != 0:
            xs.append(-a4 / a1)
        elif a0 != 0 or a3 != 0:                     # a vertical line (or two): x-roots
            r = np.roots([a0, a3, a5]) if a0 != 0 else [-a5 / a3]
            xs += [z.real for z in np.atleast_1d(r) if abs(np.imag(z)) < 1e-9]
    for a in range(len(cs)):
        for b in range(a + 1, len(cs)):
            try:
                from cca2.conic_meet import conicMeet
                P, _ = conicMeet(EcQ[idx[a] - 1], EcQ[idx[b] - 1])
                xs += [float(p[0]) for p in P]
            except Exception:                        # noqa: BLE001
                pass
    xs = sorted(set(x for x in xs if np.isfinite(x)))
    if not xs:
        xs = [0.0]
    samp = [xs[0] - 1.0, xs[-1] + 1.0]
    for i in range(len(xs) - 1):
        samp += [xs[i] + (xs[i + 1] - xs[i]) * w for w in (0.25, 0.5, 0.75)]
    for x in samp:
        ys = []
        for a0, a1, a2, a3, a4, a5 in cs:
            B = a1 * x + a4
            C = a0 * x * x + a3 * x + a5
            if a2 != 0:
                d = B * B - 4 * a2 * C
                if d >= 0:
                    ys += [(-B - d ** 0.5) / (2 * a2), (-B + d ** 0.5) / (2 * a2)]
            elif B != 0:
                ys.append(-C / B)
        ys = sorted(set(y for y in ys if np.isfinite(y)))
        if not ys:
            ys = [0.0]
        cand = [ys[0] - 1.0, ys[-1] + 1.0] + [(ys[i] + ys[i + 1]) / 2 for i in range(len(ys) - 1)]
        for y in cand:
            X = (Fr(x), Fr(y))
            ok = True
            for r in rows:
                c = [Fr(t) for t in EcQ[r[0] - 1][0:6]]
                v = c[0] * X[0] ** 2 + c[1] * X[0] * X[1] + c[2] * X[1] ** 2 + c[3] * X[0] + c[4] * X[1] + c[5]
                if r[1] * v <= 0:
                    ok = False
                    break
            if ok:
                return (x, y)
    return None


def wrapped(EcQ, rows):
    v = ORIG(EcQ, rows)
    if not v:
        p = strict_point(EcQ, rows)
        if p is not None:
            HITS.append({'EcQ': [[int(t) for t in c] for c in EcQ], 'rows': [list(r) for r in rows], 'pt': p})
    return v


AQ.cellHasInterior = wrapped
CI.cellHasInterior = wrapped


def main():
    n = int(sys.argv[1])
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    first = int(sys.argv[3]) if len(sys.argv) > 4 else 0
    last = int(sys.argv[4]) if len(sys.argv) > 4 else n
    rng = random.Random(seed)
    found = []
    for it in range(n):
        fam = rng.choice(['R1', 'R2', 'R2', 'far'])
        if fam == 'far':
            k = rng.choice([4, 5, 6])
            base = G.rstar(rng, k) if rng.random() < 0.6 else G.rconvex(rng, rng.choice([3, 4]))
            sh = (rng.randint(-12, 12), rng.randint(-12, 12))
            Pp = [(x + sh[0], y + sh[1]) for (x, y) in base]
        elif fam == 'R1':
            Pp = G.rconvex(rng, rng.choice([3, 4]))
        else:
            Pp = G.rstar(rng, rng.choice([5, 6, 7]))
        q = G.rq(rng, -3, 3)
        pcs = [{'kind': 'poly', 'P': G.P(Pp), 'q': q}]
        if not first <= it < last:
            continue                                 # another shard's input (the draws above keep the stream)
        enc = G.encode(pcs)
        HITS.clear()
        try:
            f = QuaPol(enc['V'], enc['E'], [[int(t) for t in r] for r in enc['f']], enc['F'])
            g = f.conj()
        except Exception as e:                       # noqa: BLE001
            continue
        if not HITS:
            continue
        fx = {'pieces': pcs}
        S = [(Fr(a, 4) + Fr(1, 64), Fr(b, 4) + Fr(3, 128)) for a in range(-60, 61, 3) for b in range(-60, 61, 3)]
        vals = g.eval([[float(s[0]), float(s[1])] for s in S])[0]
        wrong = 0
        ex = None
        for s, v in zip(S, vals):
            o = O.conj_exact(fx, s)
            of = float('inf') if o is None else float(o)
            if (of == float('inf')) != (v == float('inf')) or (of != float('inf') and abs(of - v) > 1e-9 * max(1, abs(of))):
                wrong += 1
                ex = ex or (float(s[0]), float(s[1]), v, of)
        found.append({'it': it, 'P': Pp, 'q': [str(t) for t in q], 'hits': len(HITS), 'wrong': wrong, 'ex': ex,
                      'hit0': HITS[0]})
        print(json.dumps(found[-1]), flush=True)
    print('done', n, 'inputs;', len(found), 'with a missed interior')


if __name__ == '__main__':
    main()

"""Second, independent look at each WRONG verdict of compare_engines.py: re-derive the certificate
with the PROJECT'S OWN definition oracle (cca2.check) and the engine's OWN conjugate.

below the lower bound : find s maximising <s,x> - f*(s) (reviewer's oracle), then report
                        <s,x> - check.supAt(f,s) and <s,x> - (engine conj)(s); if both exceed the
                        engine's f**(x), Fenchel-Young is violated by the engine against itself.
above the upper bound : the LP's convex combination x = sum l_i x_i; report sum l_i f(x_i) with f
                        from check.evalQuaPol.
engine +inf           : x inside the hull of domain samples; report check.evalQuaPol at x and the
                        LP combination value.

    python3 verify_wrong.py FIX.json OUTDIR ID[,ID...]
"""
import json
import math
import os
import sys
from fractions import Fraction as Fr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', '..')))
sys.path.insert(0, HERE)
from cca2.qua_pol import QuaPol   # noqa: E402
from cca2 import check            # noqa: E402
import plq_gen as G               # noqa: E402
import plq_oracle as O            # noqa: E402
import compare_engines as C       # noqa: E402


def main():
    fixes = {x['id']: x for x in json.load(open(sys.argv[1]))}
    outdir = sys.argv[2]
    recs = {e: C.load(outdir, e) for e in ('py', 'rs', 'ml')}
    for cid in [int(t) for t in sys.argv[3].split(',')]:
        fix = fixes[cid]
        fx = G.load_pieces(fix['pieces'])
        f = QuaPol(fix['f'][0]) if not fix['V'] else QuaPol(fix['V'], fix['E'], fix['f'], fix['F'])
        g = f.conj()
        rb = recs['py'].get((cid, 'biconj')) or recs['rs'].get((cid, 'biconj'))
        if 'vals' not in rb:
            rb = recs['rs'].get((cid, 'biconj'))
        vals = [C.fv(v) for v in rb['vals']]
        gs = np.linspace(-16, 16, 65)
        grid = np.array([[a, b] for a in gs for b in gs])
        gv = np.array(C.oracle_conj(fx, grid.tolist()))
        X = fix['X']
        box = ((Fr(-12), Fr(-12)), (Fr(12), Fr(12)))
        samples = O.domain_samples(fx, n=40, box=box)
        shown = 0
        for x, v in zip(X, vals):
            hi = O.hull_upper(fx, samples, (Fr(x[0]), Fr(x[1])))
            lo = C.lower_bound(fx, x, gv, grid)
            tol = 1e-7 * max(1.0, abs(v) if math.isfinite(v) else 1.0)
            if math.isinf(v) and math.isfinite(hi):
                kind = 'engine +inf'
            elif v < lo - tol:
                kind = 'below'
            elif v > hi + tol:
                kind = 'above'
            else:
                continue
            shown += 1
            if shown > 2:
                continue
            print('--- %d %s x=%s engine f**=%r bracket [%r, %r]: %s' % (cid, fix['name'], x, v, lo, hi, kind))
            print('    f(x) by check.evalQuaPol = %r' % check.evalQuaPol(f, [x])[0])
            if kind == 'below':
                from scipy.optimize import minimize
                i = int(np.argmax(grid @ np.array(x) - gv))
                r = minimize(lambda s: -C._phi(fx, x, s) if np.isfinite(C._phi(fx, x, s)) else 1e300,
                             grid[i], method='Nelder-Mead', options={'xatol': 1e-11, 'fatol': 1e-13, 'maxiter': 600})
                s = [float(r.x[0]), float(r.x[1])]
                sx = s[0] * x[0] + s[1] * x[1]
                if max(abs(s[0]), abs(s[1])) > 1e6:
                    print('    s = %r is huge: f** = +inf at x (x outside the closed hull of dom f)' % (s,))
                    continue
                print('    s = %r: <s,x> - f*(s): reviewer %r | check.supAt %r | engine conj %r'
                      % (s, sx - float(O.conj_exact(fx, (Fr(s[0]), Fr(s[1])))), sx - check.supAt(f, s),
                         sx - g.eval([s])[0][0]))
            else:
                from scipy.optimize import linprog
                P = np.array([[float(a), float(b)] for (a, b) in samples])
                fvv = np.array(check.evalQuaPol(f, P.tolist()))
                ok = np.isfinite(fvv)
                A = np.vstack([P[ok].T, np.ones(ok.sum())])
                res = linprog(fvv[ok], A_eq=A, b_eq=[x[0], x[1], 1.0], bounds=(0, None), method='highs')
                lam = res.x
                idx = np.argsort(-lam)[:3]
                print('    LP with check.evalQuaPol values: %r; support %s'
                      % (res.fun, [(P[ok][k].tolist(), round(float(lam[k]), 4), float(fvv[ok][k])) for k in idx]))
        print('=== %d %s: %d flagged points' % (cid, fix['name'], shown))


if __name__ == '__main__':
    main()

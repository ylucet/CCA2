"""Grade the engines' outputs (run_engines.py, reviewDiffDriver.m) against the reviewer's oracle.

    python3 compare_engines.py FIX.json OUTDIR [-j 16] > report.txt

conj   : engine f*(s) vs plq_oracle.conj_exact (exact) at every dual point; rel. tol 1e-9.
biconj : engine f**(x) vs a CERTIFIED bracket lo <= f**(x) <= hi:
           lo = max_s <s,x> - f*(s)   (f* exact; s from a grid, then Nelder-Mead on the concave
                                       function s -> <s,x> - f*(s); any s gives a valid bound)
           hi = LP over convex combinations of dense domain samples (plq_oracle.hull_upper)
         a value below lo or above hi (beyond 1e-7 relative) is a definite defect; a value inside a
         wide bracket is reported as 'loose', never as agreement.
Also: engine-vs-engine (same refusal identifier, values to 1e-9).
"""
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Fr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plq_gen as G      # noqa: E402
import plq_oracle as O   # noqa: E402

TOL_CONJ = 1e-9
TOL_BI = 1e-7


def load(outdir, eng):
    path = os.path.join(outdir, '%s.jsonl' % eng)
    if eng == 'ml':
        recs = []
        for k in range(1, 9):
            p = os.path.join(outdir, 'ml_%d.jsonl' % k)
            if os.path.isfile(p):
                recs += [json.loads(l) for l in open(p) if l.strip()]
    elif os.path.isfile(path):
        recs = [json.loads(l) for l in open(path) if l.strip()]
    else:
        return {}
    out = {}
    for r in recs:
        if r.get('err') == '':
            del r['err']
        out[(r['id'], r['op'])] = r
    return out


def fv(v):
    if isinstance(v, str):
        return {'inf': math.inf, '-inf': -math.inf, 'nan': math.nan}[v]
    return float(v)


def oracle_conj(fx, S):
    out = []
    for s in S:
        v = O.conj_exact(fx, (Fr(s[0]), Fr(s[1])))
        out.append(math.inf if v is None else float(v))
    return out


def _phi(fx, x, s):
    v = O.conj_exact(fx, (Fr(float(s[0])), Fr(float(s[1]))))
    if v is None:
        return -math.inf
    return float(s[0]) * x[0] + float(s[1]) * x[1] - float(v)


def lower_bound(fx, x, grid_vals, grid):
    """max over s of <s,x> - f*(s): grid start, then Nelder-Mead (concave objective)."""
    from scipy.optimize import minimize
    phis = grid @ np.array(x) - grid_vals
    i = int(np.argmax(phis))
    if not np.isfinite(phis[i]):
        return -math.inf
    best = phis[i]
    s0 = grid[i]
    r = minimize(lambda s: -_phi(fx, x, s) if np.isfinite(_phi(fx, x, s)) else 1e300, s0,
                 method='Nelder-Mead', options={'xatol': 1e-11, 'fatol': 1e-13, 'maxiter': 600})
    if np.isfinite(r.fun) and -r.fun > best:
        best = -r.fun
    return best


def grade_fixture(args):
    fix, recs = args
    fx = G.load_pieces(fix['pieces'])
    res = {'id': fix['id'], 'name': fix['name'], 'conj': {}, 'biconj': {}}
    # ---- conj
    ref = oracle_conj(fx, fix['S'])
    for eng, r in recs.items():
        rc = r.get((fix['id'], 'conj'))
        if rc is None:
            res['conj'][eng] = 'missing'
            continue
        if 'err' in rc:
            res['conj'][eng] = 'refuse:%s' % rc['err']
            continue
        vals = [fv(v) for v in rc['vals']]
        bad = []
        for s, a, b in zip(fix['S'], vals, ref):
            if math.isinf(a) or math.isinf(b):
                if a != b:
                    bad.append((s, a, b))
            elif abs(a - b) > TOL_CONJ * max(1.0, abs(b)):
                bad.append((s, a, b))
        res['conj'][eng] = 'ok' if not bad else 'WRONG %d/%d e.g. s=%s engine=%r oracle=%r' % (
            len(bad), len(vals), bad[0][0], bad[0][1], bad[0][2])
    # ---- biconj
    need = any('vals' in r.get((fix['id'], 'biconj'), {}) for r in recs.values())
    if not need:
        for eng, r in recs.items():
            rb = r.get((fix['id'], 'biconj'))
            res['biconj'][eng] = 'missing' if rb is None else 'refuse:%s' % rb.get('err')
        return res
    xs = [v for p in fx['pieces'] for v in (p.get('P') or [])]
    xs += [p['apex'] for p in fx['pieces'] if 'apex' in p]
    X = fix['X']
    lo_b = (min(min(t[0] for t in X), -4) - 4, min(min(t[1] for t in X), -4) - 4)
    hi_b = (max(max(t[0] for t in X), 4) + 4, max(max(t[1] for t in X), 4) + 4)
    box = ((Fr(lo_b[0]), Fr(lo_b[1])), (Fr(hi_b[0]), Fr(hi_b[1])))
    samples = O.domain_samples(fx, n=40, box=box)
    gs = np.linspace(-16, 16, 65)
    grid = np.array([[a, b] for a in gs for b in gs])
    gvals = np.array(oracle_conj(fx, grid.tolist()))
    brackets = []
    for x in X:
        hi = O.hull_upper(fx, samples, (Fr(x[0]), Fr(x[1])))
        lo = lower_bound(fx, x, gvals, grid) if math.isfinite(hi) or True else -math.inf
        brackets.append((lo, hi))
    for eng, r in recs.items():
        rb = r.get((fix['id'], 'biconj'))
        if rb is None:
            res['biconj'][eng] = 'missing'
            continue
        if 'err' in rb:
            res['biconj'][eng] = 'refuse:%s' % rb['err']
            continue
        vals = [fv(v) for v in rb['vals']]
        bad, loose = [], 0
        for x, v, (lo, hi) in zip(X, vals, brackets):
            tol = TOL_BI * max(1.0, abs(v) if math.isfinite(v) else 1.0)
            if math.isinf(v) and v > 0:
                if math.isfinite(hi):
                    bad.append((x, v, lo, hi, 'engine +inf, f** finite'))
                continue
            if v < lo - tol:
                bad.append((x, v, lo, hi, 'below the lower bound'))
            elif v > hi + tol:
                bad.append((x, v, lo, hi, 'above the upper bound'))
            elif math.isfinite(hi) and math.isfinite(lo) and hi - lo > 1e-4 * max(1.0, abs(v)):
                loose += 1
        if bad:
            b = bad[0]
            res['biconj'][eng] = 'WRONG %d/%d e.g. x=%s engine=%r bracket=[%r, %r] (%s)' % (
                len(bad), len(vals), b[0], b[1], b[2], b[3], b[4])
        else:
            res['biconj'][eng] = 'ok' + (' (loose bracket at %d pts)' % loose if loose else '')
    return res


def cross_engine(fixes, recs):
    lines = []
    engs = sorted(recs)
    for fix in fixes:
        for op in ('conj', 'biconj'):
            rs = {e: recs[e].get((fix['id'], op)) for e in engs}
            if any(v is None for v in rs.values()):
                continue
            errs = {e: v.get('err') for e, v in rs.items()}
            if len(set(errs.values())) > 1:
                lines.append('%d %s %s: refusal differs %s' % (fix['id'], fix['name'], op, errs))
                continue
            if any(errs.values()):
                continue
            base = [fv(v) for v in rs[engs[0]]['vals']]
            for e in engs[1:]:
                other = [fv(v) for v in rs[e]['vals']]
                nb = sum(1 for a, b in zip(base, other)
                         if (math.isinf(a) or math.isinf(b)) and a != b
                         or (math.isfinite(a) and math.isfinite(b) and abs(a - b) > 1e-9 * max(1, abs(a))))
                if nb:
                    lines.append('%d %s %s: %s vs %s differ at %d points' % (fix['id'], fix['name'], op, engs[0], e, nb))
    return lines


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('fix')
    ap.add_argument('outdir')
    ap.add_argument('-j', type=int, default=16)
    ap.add_argument('--engines', default='py,rs,ml')
    ap.add_argument('--ids', default=None)
    a = ap.parse_args()
    fixes = json.load(open(a.fix))
    if a.ids:
        keep = {int(t) for t in a.ids.split(',')}
        fixes = [x for x in fixes if x['id'] in keep]
    recs = {e: load(a.outdir, e) for e in a.engines.split(',')}
    recs = {e: r for e, r in recs.items() if r}
    print('engines with output:', sorted(recs), {e: len(r) for e, r in recs.items()})
    with ProcessPoolExecutor(a.j) as ex:
        results = list(ex.map(grade_fixture, [(f, recs) for f in fixes]))
    json.dump(results, open(os.path.join(a.outdir, 'graded.json'), 'w'), indent=0)
    from collections import Counter
    for op in ('conj', 'biconj'):
        for e in sorted(recs):
            c = Counter(r[op].get(e, 'missing').split(':')[0].split(' (')[0].split(' ')[0] for r in results)
            print(op, e, dict(c))
    print()
    for r in results:
        for op in ('conj', 'biconj'):
            for e, v in sorted(r[op].items()):
                if v.startswith('WRONG') or v.startswith('refuse:CRASH') or 'TIMEOUT' in v or 'NOOUTPUT' in v:
                    print('%3d %-22s %-6s %s: %s' % (r['id'], r['name'], op, e, v[:230]))
    print()
    print('engine-vs-engine:')
    for l in cross_engine(fixes, recs):
        print(' ', l)


if __name__ == '__main__':
    main()

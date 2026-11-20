"""Run the Python and Rust engines on the review fixtures (plq_gen.py output), one subprocess per
fixture with a timeout, in parallel. Writes one JSON line per (fixture, op):
  {"id", "op", "err" | "kind", "vals": [...] (inf as "inf")}

    python3 run_engines.py FIX.json OUTDIR [--engines py,rs] [-j 16] [--timeout 300]
    python3 run_engines.py --one-py FIX.json ID           (internal)
"""
import json
import math
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
PYDIR = os.path.normpath(os.path.join(HERE, '..', '..'))
ROOT = os.path.dirname(PYDIR)


def fj(v):
    if isinstance(v, float) and math.isinf(v):
        return 'inf' if v > 0 else '-inf'
    if isinstance(v, float) and math.isnan(v):
        return 'nan'
    return v


def one_py(fixfile, cid):
    sys.path.insert(0, PYDIR)
    from cca2.qua_pol import QuaPol
    from cca2.errors import CCA2Error
    x = [r for r in json.load(open(fixfile)) if r['id'] == cid][0]
    if not x['V']:
        f = QuaPol(x['f'][0])
    else:
        f = QuaPol(x['V'], x['E'], x['f'], x['F'])
    for op, pts in (('conj', x['S']), ('biconj', x['X'])):
        rec = {'id': cid, 'op': op}
        try:
            g = f.conj() if op == 'conj' else f.biconj()
            rec['kind'] = type(g).__name__
            if hasattr(g, 'nf'):
                rec['nf'] = g.nf
            v = g.eval(pts)
            v = v[0] if isinstance(v, tuple) else v
            rec['vals'] = [fj(float(a)) for a in v]
        except CCA2Error as e:
            rec['err'] = e.identifier
        except Exception as e:      # noqa: BLE001 -- an uncaught crash is itself a finding
            rec['err'] = 'CRASH %s: %s' % (type(e).__name__, str(e)[:200])
        print(json.dumps(rec), flush=True)


def run_py(fixfile, cid, timeout):
    try:
        p = subprocess.run([sys.executable, __file__, '--one-py', fixfile, str(cid)],
                           capture_output=True, text=True, timeout=timeout)
        lines = [l for l in p.stdout.splitlines() if l.startswith('{')]
        if not lines:
            return [{'id': cid, 'op': op, 'err': 'NOOUTPUT rc=%d %s' % (p.returncode, p.stderr[-300:])}
                    for op in ('conj', 'biconj')]
        got = [json.loads(l) for l in lines]
        ops = {g['op'] for g in got}
        for op in ('conj', 'biconj'):
            if op not in ops:
                got.append({'id': cid, 'op': op, 'err': 'NOOUTPUT rc=%d' % p.returncode})
        return got
    except subprocess.TimeoutExpired:
        return [{'id': cid, 'op': op, 'err': 'TIMEOUT'} for op in ('conj', 'biconj')]


def golden_bin():
    for c in (os.path.join(ROOT, 'rust', 'target', 'release', 'golden'),):
        if os.path.isfile(c):
            return c
    raise SystemExit('rust golden binary not built')


def run_rs(fix, timeout, gbin):
    cid = fix['id']
    inp = {'nv': len(fix['V']), 'ne': len(fix['E']), 'nf': len(fix['f']),
           'V': fix['V'], 'E': fix['E'], 'f': fix['f'], 'F': fix['F'], 'exact': True}
    out = []
    for op, pts in (('conj', fix['S']), ('biconj', fix['X'])):
        rec = {'fn': 'conjQ' if op == 'conj' else 'biconjQ', 'in': inp, 'out': {'evalX': pts}}
        with tempfile.NamedTemporaryFile('w', suffix='.jsonl', delete=False) as tf:
            tf.write(json.dumps(rec) + '\n')
            tmp = tf.name
        r = {'id': cid, 'op': op}
        try:
            p = subprocess.run([gbin, tmp], capture_output=True, text=True, timeout=timeout)
            line = [l for l in p.stdout.splitlines() if l.startswith('{')]
            if not line:
                r['err'] = 'CRASH rc=%d %s' % (p.returncode, p.stderr[-300:])
            else:
                j = json.loads(line[0])
                if 'err' in j:
                    r['err'] = j['err']
                else:
                    r['kind'] = j.get('kind')
                    if 'nf' in j:
                        r['nf'] = j['nf']
                    if 'evalV' in j:
                        r['vals'] = j['evalV']
                    else:
                        r['err'] = 'EVALERR %s' % j.get('evalErr')
        except subprocess.TimeoutExpired:
            r['err'] = 'TIMEOUT'
        finally:
            os.unlink(tmp)
        out.append(r)
    return out


def main():
    if sys.argv[1] == '--one-py':
        one_py(sys.argv[2], int(sys.argv[3]))
        return
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('fix')
    ap.add_argument('outdir')
    ap.add_argument('--engines', default='py,rs')
    ap.add_argument('-j', type=int, default=16)
    ap.add_argument('--timeout', type=int, default=300)
    ap.add_argument('--ids', default=None)
    a = ap.parse_args()
    fixes = json.load(open(a.fix))
    if a.ids:
        keep = {int(t) for t in a.ids.split(',')}
        fixes = [x for x in fixes if x['id'] in keep]
    os.makedirs(a.outdir, exist_ok=True)
    for eng in a.engines.split(','):
        path = os.path.join(a.outdir, '%s.jsonl' % eng)
        if eng == 'py':
            fn = lambda x: run_py(os.path.abspath(a.fix), x['id'], a.timeout)  # noqa: E731
        else:
            gb = golden_bin()
            fn = lambda x: run_rs(x, a.timeout, gb)  # noqa: E731
        with ThreadPoolExecutor(a.j) as ex, open(path, 'w') as fh:
            for recs in ex.map(fn, fixes):
                for r in recs:
                    fh.write(json.dumps(r) + '\n')
                fh.flush()
        print('%s -> %s' % (eng, path))


if __name__ == '__main__':
    main()

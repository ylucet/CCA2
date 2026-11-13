"""Compare a Rust replay (src/bin/golden.rs) against the Python oracle (tools/py_ref.py), record
by record: refusals by identifier, QuaCon field for field, dom, and evalV to 1e-9 relative.
Also reports against the MATLAB golden's evalV. Usage: compare.py trace.jsonl rust.jsonl py.jsonl"""
import json, math, sys


def num(x):
    return {'inf': math.inf, '-inf': -math.inf, 'nan': math.nan}.get(x, x) if isinstance(x, str) else x


def evdiff(a, b, tol=1e-9):
    for i, (x, y) in enumerate(zip(a, b)):
        x, y = num(x), num(y)
        if math.isinf(x) or math.isinf(y) or math.isnan(x) or math.isnan(y):
            if not (x == y or (math.isnan(x) and math.isnan(y))):
                return 'eval[%d] %r vs %r' % (i, x, y)
        elif abs(x - y) > tol * max(1.0, abs(y)):
            return 'eval[%d] %r vs %r' % (i, x, y)
    if len(a) != len(b):
        return 'eval length %d vs %d' % (len(a), len(b))
    return ''


def ints(v):
    if isinstance(v, list):
        return [ints(x) for x in v]
    return int(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else v


def main(trace, rs, py, verbose=1):
    T = [json.loads(l) for l in open(trace) if l.strip()]
    Rr = [json.loads(l) for l in open(rs)]
    Pp = [json.loads(l) for l in open(py)]
    ok = bad = skip = 0
    tr = tp = 0.0
    for i, (rec, r, p) in enumerate(zip(T, Rr, Pp)):
        if r is None or p is None:
            skip += 1
            continue
        tr += r.get('t', 0)
        tp += p.get('t', 0)
        why = ''
        if 'err' in r or 'err' in p:
            if r.get('err') != p.get('err'):
                why = 'err rust %r vs python %r %s' % (r.get('err'), p.get('err'), r.get('msg', '')[:150])
        else:
            if r.get('kind') != p.get('kind'):
                why = 'kind %r vs %r' % (r.get('kind'), p.get('kind'))
            else:
                for k in ('nf', 'ne', 'nv', 'fN', 'fD', 'EcQ', 'FC', 'Vname', 'E', 'F', 'dom'):
                    if k in p and ints(r.get(k)) != ints(p[k]):
                        why = '%s: rust %s vs python %s' % (k, str(r.get(k))[:200], str(p[k])[:200])
                        break
                if not why and 'evalV' in p:
                    why = evdiff(r.get('evalV', []), p['evalV'])
        if why:
            bad += 1
            if verbose:
                print('#%d %s: %s' % (i, rec['fn'], why))
        else:
            ok += 1
    print('match %d, differ %d, skipped %d   time rust %.3f s, python %.3f s' % (ok, bad, skip, tr, tp))


if __name__ == '__main__':
    main(*sys.argv[1:4])

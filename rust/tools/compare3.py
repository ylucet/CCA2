"""Three-way benchmark summary: MATLAB and Python (CCA2/python/bench/run5, not re-measured) against
Rust (bench/rust_p*.jsonl). Per fixture the BEST time over all passes; groups as in compare.py."""
import json, os, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
B = os.path.join(HERE, '..', '..', 'python', 'bench')


def load(paths):
    out = {}
    for p in paths:
        for line in open(p):
            r = json.loads(line)
            o = out.get(r['id'])
            if o is None:
                out[r['id']] = dict(r)
                continue
            o['err'] = o['err'] or r['err']
            if r['n'] > 0:
                o['best'] = r['best'] if o['n'] == 0 else min(o['best'], r['best'])
                o['n'] += r['n']
    return out


F = {x['id']: x for x in json.load(open(os.path.join(B, 'fixtures.json')))}
M = load([os.path.join(B, 'run5', 'matlab_p%d.jsonl' % k) for k in (1, 2)])
P = load([os.path.join(B, 'run5', 'python_p%d.jsonl' % k) for k in (1, 2)])
Rr = load([os.path.join(HERE, '..', 'bench', 'rust_p%d.jsonl' % k) for k in (1, 2)])
ids = [i for i in F if all(i in X and not X[i]['err'] and X[i]['n'] > 0 for X in (M, P, Rr))]
kinds = [i for i in ids if Rr[i].get('kind') and Rr[i]['kind'] != F[i]['kind']]
print('fixtures timed by all three: %d; result kind differing from the fixture: %s\n' % (len(ids), kinds))
groups = {}
for i in ids:
    x = F[i]
    groups.setdefault((x['bucket'], x['fn'], x['kind']), []).append(i)
print('| bucket | call | result | n | MATLAB (s) | Python (s) | Rust (s) | MATLAB/Python | Python/Rust | MATLAB/Rust |')
print('|---|---|---|---|---|---|---|---|---|---|')
for key in sorted(groups):
    g = groups[key]
    tm = sum(M[i]['best'] for i in g)
    tp = sum(P[i]['best'] for i in g)
    tr = sum(Rr[i]['best'] for i in g)
    print('| %s | %s | %s | %d | %.3f | %.3f | %.3f | %.1f | %.1f | %.0f |' % (key + (len(g), tm, tp, tr, tm / tp, tp / tr, tm / tr)))
tm = sum(M[i]['best'] for i in ids)
tp = sum(P[i]['best'] for i in ids)
tr = sum(Rr[i]['best'] for i in ids)
print('| **all** | | | %d | **%.1f** | **%.1f** | **%.1f** | **%.1f** | **%.1f** | **%.0f** |' % (len(ids), tm, tp, tr, tm / tp, tp / tr, tm / tr))
print('\nSlowest fixtures (best time, seconds):')
print('| id | source | call | result | MATLAB | Python | Rust | Python/Rust |')
print('|---|---|---|---|---|---|---|---|')
for i in sorted(ids, key=lambda i: -M[i]['best'])[:10]:
    x = F[i]
    print('| %d | %s | %s | %s | %.2f | %.2f | %.3f | %.1f |' % (i, x['src'].replace('.jsonl', ''), x['fn'], x['kind'], M[i]['best'], P[i]['best'], Rr[i]['best'], P[i]['best'] / Rr[i]['best']))
rat = [P[i]['best'] / Rr[i]['best'] for i in ids if Rr[i]['best'] > 0]
print('\nPer-fixture Python/Rust ratio: median %.1f, range %.1f - %.1f.' % (statistics.median(rat), min(rat), max(rat)))

#!/usr/bin/env python3
"""Aggregate tools/par_sweep.sh output: per thread count, per fixture the MINIMUM of `best` over the
sweeps, then sums per group. Usage: par_sweep.py RUNDIR"""
import json, glob, re, sys, collections
d = sys.argv[1] if len(sys.argv) > 1 else 'bench/parallel/runs'
best = collections.defaultdict(dict)   # T -> id -> min best
for f in glob.glob(f'{d}/s*_t*.jsonl'):
    t = int(re.search(r'_t(\d+)\.jsonl', f).group(1))
    for l in open(f):
        r = json.loads(l)
        if r['best'] < 0: continue
        i = r['id']; best[t][i] = min(best[t].get(i, 1e9), r['best'])
ts = sorted(best)
heavy = [196, 191, 200, 177, 174, 165]
print('| threads | all 200 fixtures (s) | speedup | ' + ' | '.join(f'#{i} (s)' for i in heavy) + ' | other 194 (s) |')
print('|---|---|---|' + '---|' * (len(heavy) + 1))
base = sum(best[1].values()) if 1 in best else None
for t in ts:
    tot = sum(best[t].values())
    row = [f'{best[t].get(i, float("nan")):.3f}' for i in heavy]
    oth = sum(v for i, v in best[t].items() if i not in heavy)
    print(f'| {t} | {tot:.2f} | {base / tot:.1f}x | ' + ' | '.join(row) + f' | {oth:.2f} |')
if 1 in best:
    print()
    print('per-fixture speedup vs 1 thread (heavy):')
    for i in heavy:
        print(f'  #{i}: ' + '  '.join(f'T{t}={best[1][i] / best[t][i]:.1f}x' for t in ts if t != 1 and i in best[t]))

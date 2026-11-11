"""Compare bench/matlab.jsonl and bench/python.jsonl (same fixtures, same order) and print a
markdown summary: per group the totals of the BEST time per fixture, the median and the range of
the per-fixture ratio Python/MATLAB."""
import json
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    out = {}
    for line in open(os.path.join(HERE, name)):
        r = json.loads(line)
        out[r['id']] = r
    return out


def main():
    F = {x['id']: x for x in json.load(open(os.path.join(HERE, 'fixtures.json')))}
    M = load('matlab.jsonl')
    P = load('python.jsonl')
    groups = {}
    for i, x in F.items():
        if i not in M or i not in P or M[i]['err'] or P[i]['err'] or M[i]['n'] == 0 or P[i]['n'] == 0:
            continue
        key = (x['bucket'], x['fn'], x['kind'])
        groups.setdefault(key, []).append(i)
    print('| bucket | call | result | fixtures | MATLAB total (s) | Python total (s) | Python/MATLAB total | median ratio | ratio range |')
    print('|---|---|---|---|---|---|---|---|---|')
    for key in sorted(groups):
        ids = groups[key]
        tm = sum(M[i]['best'] for i in ids)
        tp = sum(P[i]['best'] for i in ids)
        ratios = [P[i]['best'] / M[i]['best'] for i in ids if M[i]['best'] > 0]
        print('| %s | %s | %s | %d | %.3f | %.3f | %.2f | %.2f | %.2f - %.2f |' % (
            key[0], key[1], key[2], len(ids), tm, tp, tp / tm, statistics.median(ratios),
            min(ratios), max(ratios)))
    print()
    print('Slowest fixtures (best time, seconds):')
    print('| id | source | call | result | MATLAB | Python | ratio | reps M/P |')
    print('|---|---|---|---|---|---|---|---|')
    big = sorted((i for g in groups.values() for i in g), key=lambda i: -max(M[i]['best'], P[i]['best']))[:12]
    for i in big:
        x = F[i]
        print('| %d | %s | %s | %s | %.2f | %.2f | %.2f | %d/%d |' % (
            i, x['src'].replace('.jsonl', ''), x['fn'], x['kind'], M[i]['best'], P[i]['best'],
            P[i]['best'] / M[i]['best'], M[i]['n'], P[i]['n']))
    print()
    spread = [(M[i]['worst'] / M[i]['best'], P[i]['worst'] / P[i]['best'])
              for g in groups.values() for i in g if M[i]['n'] > 1 and P[i]['n'] > 1 and M[i]['best'] > 1e-3]
    if spread:
        print('Spread worst/best over repetitions (fixtures above 1 ms): MATLAB median %.2f max %.2f; '
              'Python median %.2f max %.2f.' % (
                  statistics.median(s[0] for s in spread), max(s[0] for s in spread),
                  statistics.median(s[1] for s in spread), max(s[1] for s in spread)))


if __name__ == '__main__':
    main()

"""Build d_family_ref.json: the D family's conj and biconj values, for the cross-engine tests of #108
(rust/tests/review_2026_11_19.rs, MATLAB/review/reviewFindings20261119Test.m).

The values are the Python engine's, taken from a run_engines.py output directory that
compare_engines.py graded with 0 wrong values (conj against the exact oracle, biconj inside the
certified bracket); a named refusal is recorded as such. Every second primal point and every third
dual point of the generator are kept.

    python3 make_d_family_ref.py FIX.json OUTDIR > d_family_ref.json
"""
import json
import sys


def main():
    fixes = {x['id']: x for x in json.load(open(sys.argv[1]))}
    recs = {}
    for line in open(sys.argv[2] + '/py.jsonl'):
        r = json.loads(line)
        recs[(r['id'], r['op'])] = r
    out = []
    for cid in sorted(i for i in fixes if fixes[i]['name'].startswith('D-')):
        x = fixes[cid]
        S = x['S'][::3]
        X = x['X'][::2]
        rec = {'id': cid, 'name': x['name'], 'V': x['V'], 'E': x['E'], 'f': x['f'], 'F': x['F'],
               'S': S, 'X': X}
        for op, keep in (('conj', 3), ('biconj', 2)):
            r = recs[(cid, op)]
            if 'err' in r:
                if r['err'].startswith(('TIMEOUT', 'CRASH', 'NOOUTPUT')):
                    raise SystemExit('%s %s: %s is not a reference' % (x['name'], op, r['err']))
                rec[op + '_err'] = r['err']
            else:
                rec[op] = r['vals'][::keep]
        out.append(rec)
    json.dump(out, sys.stdout, separators=(',', ':'))
    sys.stdout.write('\n')


if __name__ == '__main__':
    main()

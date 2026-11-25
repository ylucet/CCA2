"""Witness search for TWO BRIDGE FACES OWNING ONE CELL (DESIGN.md "Two bridge faces never share a
cell"). The argument says it cannot happen except on a set of measure zero (or in one degenerate
exception, where it is harmless); this looks for a counterexample two independent ways, per input:

  * DIRECT: for every envelope (one pair of pieces or atoms) with >= 2 bridge faces, draw N random
    points around dom f and ask `bridgeMember` -- the face's DEFINING predicate, independent of the
    arrangement -- how many bridge regions each lies in. A witness is a point in two ALGEBRAIC
    bridges (their membership is exactly the sweep the argument is about). Points in two bridges
    of which one is RATIONAL are counted apart (`inTwoAny`): that membership is an outer bound of
    the sweep (no multiplier-sign condition), so such overlaps are expected and say nothing.
  * ARRANGEMENT: `envelopeSubdivisionAlg`'s own count `info['bridgeOverlaps']` (cells whose majority
    vote puts them in two algebraic bridges).

    python3 hunt_two_bridges.py FIX.json OUT.jsonl [-j 24] [--timeout 600] [--samples 4000]
    python3 hunt_two_bridges.py --one FIX.json ID SAMPLES      (internal)
    python3 hunt_two_bridges.py --summary OUT.jsonl [...]

FIX.json is plq_gen.py output. Every input runs in its own process under a timeout.
"""
import json
import os
import random
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
PYDIR = os.path.normpath(os.path.join(HERE, '..', '..'))


def one(fixfile, cid, nsamp):
    sys.path.insert(0, PYDIR)
    import cca2.envelope_subdivision_alg as E
    from cca2.bridge_member import bridgeMember
    from cca2.qua_pol import QuaPol
    st = {'envs': 0, 'envs2b': 0, 'maxBridges': 0, 'samples': 0, 'inBridge': 0, 'inTwo': 0,
          'inTwoAny': 0, 'algBridges': 0,
          'bridgeOverlaps': 0, 'subdivided': 0}
    box = [-4.0, 4.0, -4.0, 4.0]
    split = E._splitBridge

    def splitBridge(bl):
        corners, bridges = split(bl)
        st['envs'] += 1
        st['maxBridges'] = max(st['maxBridges'], len(bridges))
        if len(bridges) >= 2:
            st['envs2b'] += 1
            rng = random.Random(st['envs'])
            X = [[rng.uniform(box[0], box[1]), rng.uniform(box[2], box[3])] for _ in range(nsamp)]
            M = [bridgeMember(b['pair'], X)[0] for b in bridges]
            alg = [b['kind'] == 'algebraic' for b in bridges]
            st['algBridges'] += sum(alg)
            for k in range(nsamp):
                c = sum(1 for m in M if m[k])
                ca = sum(1 for m, a in zip(M, alg) if a and m[k])
                st['samples'] += 1
                st['inBridge'] += c >= 1
                st['inTwoAny'] += c >= 2
                st['inTwo'] += ca >= 2
        return corners, bridges
    E._splitBridge = splitBridge
    esa = E.envelopeSubdivisionAlg

    def envelopeSubdivisionAlg(env, *a, **k):
        sub, info = esa(env, *a, **k)
        if not isinstance(env, list) or len(env) == 1:
            st['subdivided'] += 1
            st['bridgeOverlaps'] += info.get('bridgeOverlaps', 0)
        return sub, info
    E.envelopeSubdivisionAlg = envelopeSubdivisionAlg
    import cca2.biconj_q as B
    if hasattr(B, 'envelopeSubdivisionAlg'):
        B.envelopeSubdivisionAlg = envelopeSubdivisionAlg
    for mod in list(sys.modules.values()):
        if getattr(mod, 'envelopeSubdivisionAlg', None) is esa:
            mod.envelopeSubdivisionAlg = envelopeSubdivisionAlg
    x = [r for r in json.load(open(fixfile)) if str(r['id']) == str(cid)][0]
    if x['V']:
        xs = [v[0] for v in x['V']]
        ys = [v[1] for v in x['V']]
        box[:] = [min(xs) - .5, max(xs) + .5, min(ys) - .5, max(ys) + .5]
    f = QuaPol(x['f'][0]) if not x['V'] else QuaPol(x['V'], x['E'], x['f'], x['F'])
    try:
        f.biconj()
        err = None
    except Exception as e:      # noqa: BLE001 -- a refusal is recorded, not a finding here
        err = getattr(e, 'identifier', type(e).__name__)
    print(json.dumps({'id': cid, 'name': x.get('name'), 'err': err, **st}), flush=True)


def run(fixfile, out, jobs, timeout, nsamp):
    ids = [r['id'] for r in json.load(open(fixfile))]

    def go(cid):
        try:
            p = subprocess.run([sys.executable, __file__, '--one', fixfile, str(cid), str(nsamp)],
                               capture_output=True, text=True, timeout=timeout)
            lines = [l for l in p.stdout.splitlines() if l.startswith('{')]
            return lines[-1] if lines else json.dumps({'id': cid, 'err': 'CRASH'})
        except subprocess.TimeoutExpired:
            return json.dumps({'id': cid, 'err': 'TIMEOUT'})
    with open(out, 'w') as fo, ThreadPoolExecutor(jobs) as ex:
        for line in ex.map(go, ids):
            fo.write(line + '\n')
            fo.flush()


def summary(files):
    tot = {}
    for fn in files:
        L = [json.loads(l) for l in open(fn) if l.startswith('{')]
        done = [l for l in L if 'envs' in l]
        s = {'inputs': len(L), 'finished': len(done)}
        for k in ('envs', 'envs2b', 'algBridges', 'samples', 'inBridge', 'inTwo', 'inTwoAny',
                  'subdivided', 'bridgeOverlaps'):
            s[k] = sum(l[k] for l in done)
        s['inputs2b'] = sum(1 for l in done if l['envs2b'])
        s['maxBridges'] = max([l['maxBridges'] for l in done] or [0])
        s['witnesses'] = [l['id'] for l in done if l['inTwo'] or l['bridgeOverlaps']]
        print(fn, json.dumps(s))
        for k, v in s.items():
            if isinstance(v, int):
                tot[k] = max(tot.get(k, 0), v) if k == 'maxBridges' else tot.get(k, 0) + v
    print('TOTAL', json.dumps(tot))


if __name__ == '__main__':
    if sys.argv[1] == '--one':
        one(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    elif sys.argv[1] == '--summary':
        summary(sys.argv[2:])
    else:
        import argparse
        ap = argparse.ArgumentParser()
        ap.add_argument('fix')
        ap.add_argument('out')
        ap.add_argument('-j', type=int, default=24)
        ap.add_argument('--timeout', type=int, default=600)
        ap.add_argument('--samples', type=int, default=4000)
        a = ap.parse_args()
        run(a.fix, a.out, a.j, a.timeout, a.samples)

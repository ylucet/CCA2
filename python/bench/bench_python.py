"""Time conjQ / biconjQ on the benchmark fixtures with the Python port. Mirrors bench_matlab.m:
the first call, then `reps` more unless the first took longer than slowCap seconds; wall-clock
inside one process (import time excluded). One JSON line per fixture."""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from cca2 import CCA2Error  # noqa: E402
from cca2 import alg_row_fails, conic_meet, rat_q  # noqa: E402
from cca2.biconj_q import biconjQ  # noqa: E402
from cca2.conj_q import conjQ  # noqa: E402
from cca2.qua_pol import QuaPol  # noqa: E402


def coldPythonOnlyMemos():
    """Empty the memos that MATLAB does NOT have before every timed call, so a repetition is not
    faster than MATLAB's for that reason alone: conicMeet's pairs, the surd corners, the Sturm
    chains and square-free parts, the boundary factors and their screen points, and the float
    arrays of the row matrices. algRowFails' other memos mirror MATLAB's `persistent` ones and
    stay warm on both sides."""
    conic_meet.clearMemo()
    alg_row_fails._memoSurdOf.clear()
    alg_row_fails._memoFactorInts.clear()
    alg_row_fails._memoCurvePts.clear()
    alg_row_fails._arrMemo.clear()
    rat_q._chainMemo.clear()
    rat_q._sqfMemo.clear()


def build(x):
    if not x['V']:
        return QuaPol(x['f'][0] if len(x['f']) == 1 else x['f'])
    return QuaPol(x['V'], x['E'] or [], x['f'], x['F'] or [])


def main(fixFile, outFile, reps=3, slowCap=20.0, only=None):
    F = json.load(open(fixFile))
    conjQ(QuaPol([1, 0, 1, 0, 0, 0]))
    biconjQ(QuaPol([1, 0, 1, 0, 0, 0]))
    with open(outFile, 'w') as fh:
        for x in F:
            if only and x['bucket'] != only:
                continue
            f = build(x)
            call = (lambda: conjQ(f)) if x['fn'] == 'conjQ' else (lambda: biconjQ(f))
            t = []
            err = ''
            try:
                coldPythonOnlyMemos()
                t0 = time.perf_counter()
                call()
                t.append(time.perf_counter() - t0)
                if t[0] <= slowCap:
                    for _ in range(reps):
                        coldPythonOnlyMemos()
                        t0 = time.perf_counter()
                        call()
                        t.append(time.perf_counter() - t0)
            except CCA2Error as e:
                err = e.identifier
            rec = {'id': x['id'], 'first': t[0] if t else -1, 'best': min(t) if t else -1,
                   'worst': max(t) if t else -1, 'n': len(t), 'err': err}
            fh.write(json.dumps(rec) + '\n')
            fh.flush()
            print('%4d %-8s %-11s first %9.4f  best %9.4f  %s' % (
                x['id'], x['fn'], x['kind'], rec['first'], rec['best'], err), flush=True)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], *(int(a) for a in sys.argv[3:4]))

"""Replay a golden trace through the PYTHON port (the oracle) and print one JSON line per record,
in the shape src/bin/golden.rs prints: {"err": id} or the QuaCon fields / kind, "dom" for conjQ,
"evalV" at the record's evalX, and "t" (seconds).

    PYTHONPATH=../python/.pydeps python3 tools/py_ref.py ../python/golden/conj.jsonl > ref/conj.py.jsonl
"""
import json, math, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
PY = os.path.join(HERE, '..', '..', 'python')
sys.path.insert(0, PY)
from cca2 import io_matlab, CCA2Error                      # noqa: E402
from cca2.conj_q import conjQ                              # noqa: E402
from cca2.biconj_q import biconjQ                          # noqa: E402


def fj(x):
    if isinstance(x, float):
        if math.isnan(x):
            return 'nan'
        if math.isinf(x):
            return 'inf' if x > 0 else '-inf'
    return x


def main(path, only=None):
    for line in open(path):
        if not line.strip():
            continue
        rec = json.loads(line)
        if only and rec['fn'] != only:
            print('null')
            continue
        t0 = time.perf_counter()
        try:
            P = io_matlab.quaPolFrom(rec['in'])
            if rec['fn'] == 'conjQ':
                g, dom = conjQ(P)
                out = dict(io_matlab.quaConFields(g), kind='QuaCon', dom=dom)
            else:
                g = biconjQ(P)
                out = {'kind': g.kind()}
                if g.kind() == 'QuaCon':
                    out.update(io_matlab.quaConFields(g))
                elif hasattr(g, 'nf'):
                    out['nf'] = g.nf
            X = io_matlab._mat(rec['out']['evalX']) if rec.get('out') and 'evalX' in rec['out'] else []
            if X:
                v = g.eval(X)
                v = v[0] if isinstance(v, tuple) else v
                out['evalV'] = [fj(float(a)) for a in v]
        except CCA2Error as e:
            out = {'err': e.identifier}
        except Exception as e:  # noqa
            out = {'err': 'python:' + type(e).__name__, 'msg': str(e)[:200]}
        out['t'] = time.perf_counter() - t0
        print(json.dumps(out), flush=True)


if __name__ == '__main__':
    main(*sys.argv[1:])

"""Unit parity against MATLAB traces: the emptiness certificate and the algebraic subdivision."""
import os

import pytest

from cca2 import io_matlab
from cca2.io_matlab import dec, sc, vec
from cca2.alg_row_fails import algRowFails
from cca2.envelope_subdivision_alg import envelopeSubdivisionAlg
from cca2.poly_sign_on import polySignOn

from conftest import GOLDEN

TRACE = os.path.join(GOLDEN, 'trace_biconj.jsonl')


def _load(fn):
    return io_matlab.loadTrace(TRACE, fn) if os.path.exists(TRACE) else []


ARF = _load('algRowFails')
ESA = _load('envelopeSubdivisionAlg')


def _row(a):
    d = {'FP': io_matlab.polyFromTrace(a['FP'], 3), 'FB': sc(a['FB']),
         'num': vec(a['num']) if a.get('num') else None,
         'den': sc(a['den']) if a.get('den') else None,
         'side': sc(a['side']) if a.get('side') else None}
    W = a.get('W')
    d['W'] = {'FP': io_matlab.polyFromTrace(W['FP'], 3), 'FB': sc(W['FB'])} if W else None
    return d


@pytest.mark.parametrize('i', range(len(ARF)))
def test_algRowFails(i):
    """inBox is a MATLAB closure and is not in the trace; it is rebuilt from the region's row
    polynomials (`factors`), which is what it tests when the region has no branch line."""
    rec = ARF[i]
    a = dec(rec['in'])
    rows = [_row(q) for q in (a[0] if isinstance(a[0], list) else [a[0]])]
    box = vec(a[1])
    depth = sc(a[3]) if len(a) > 3 else 6
    factors = a[4] if len(a) > 4 and a[4] else []
    corners = a[5] if len(a) > 5 and a[5] else []

    def inBox(b):
        return not any(polySignOn(P, b, lambda _b: True, 0) == -1 for P in factors)
    assert algRowFails(rows, box, inBox, depth, factors, corners) == bool(sc(dec(rec['out'])[0]))


def _face(w):
    alg = w['alg'] if isinstance(w['alg'], list) else [w['alg']]
    return {'kind': w['kind'], 'num': vec(w['num']) if w['num'] else [], 'den': sc(w['den']),
            'FP': io_matlab.polyFromTrace(w['FP'], 3) if w['FP'] else None,
            'FB': sc(w['FB']) if w['FB'] else 0, 'con': w['con'] or [],
            'alg': [{'FP': io_matlab.polyFromTrace(q['FP'], 3), 'FB': sc(q['FB']),
                     'num': vec(q['num']), 'den': sc(q['den']), 'side': sc(q['side'])}
                    for q in alg if q]}


@pytest.mark.slow
@pytest.mark.parametrize('i', range(len(ESA)))
def test_envelopeSubdivisionAlg(i):
    rec = ESA[i]
    a = dec(rec['in'])
    env = io_matlab.envFromTrace(a[0])
    dom = a[1] if len(a) > 1 else []
    sub, info = envelopeSubdivisionAlg(env, dom or [])
    out = dec(rec['out'])
    want = [_face(w) for w in (out[0] if isinstance(out[0], list) else [out[0]])]
    assert len(sub) == len(want)
    for k, (s, w) in enumerate(zip(sub, want)):
        assert s == w, 'face %d differs in %s' % (k + 1, [f for f in s if s[f] != w[f]])

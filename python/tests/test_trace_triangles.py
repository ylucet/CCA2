"""Unit parity against MATLAB traces: dual points and contact-triangle faces."""
import math
import os

import pytest

from cca2 import CCA2Error, io_matlab
from cca2.io_matlab import dec, sc, vec
from cca2.dual_point import dualPoint, dualPointSign
from cca2.triangle_face import triangleFace, triangleFaceSign

from conftest import GOLDEN

TRACE = os.path.join(GOLDEN, 'trace_tri.jsonl')


def _load(fn):
    return io_matlab.loadTrace(TRACE, fn) if os.path.exists(TRACE) else []


DP = _load('dualPoint')
DPS = _load('dualPointSign')
TF = _load('triangleFace')
TFS = _load('triangleFaceSign')


def _close(a, b, tol=1e-9):
    return all(abs(float(x) - float(y)) <= tol * max(1.0, abs(float(y))) for x, y in zip(a, b))


def _dp(d):
    iv = vec(d['iv'])
    return {'c1': vec(d['c1']), 'c2': vec(d['c2']), 'k': sc(d['k']), 'R': vec(d['R']),
            'iv': (iv[0], iv[1]), 'S1': vec(d['S1']), 'S0': vec(d['S0']), 's': vec(d['s'])}


@pytest.mark.parametrize('i', range(len(DP)))
def test_dualPoint(i):
    rec = DP[i]
    a = dec(rec['in'])
    args = (vec(a[0]), vec(a[1]), vec(a[2]))
    if rec['err']:
        with pytest.raises(CCA2Error) as e:
            dualPoint(*args)
        assert e.value.identifier == rec['err']
        return
    got = dualPoint(*args)
    w = _dp(dec(rec['out'])[0])
    for k in ('c1', 'c2', 'k', 'R', 'S1', 'S0'):
        assert got[k] == w[k], k
    assert [float(v) for v in got['iv']] == [float(v) for v in w['iv']]
    assert _close(got['s'], w['s'])


@pytest.mark.parametrize('i', range(len(DPS)))
def test_dualPointSign(i):
    rec = DPS[i]
    a = dec(rec['in'])
    assert dualPointSign(_dp(a[0]), a[1]) == sc(dec(rec['out'])[0])


def _tf(d):
    tf = {'kind': d['kind'], 'Gn': d['Gn'], 'Gd': vec(d['Gd']), 's': vec(d['s']) if d['s'] else None}
    if d['kind'] == 'rational':
        tf.update(s0n=vec(d['s0n']), s0d=sc(d['s0d']), Tn=d['Tn'], Td=vec(d['Td']))
    elif d.get('dp'):
        tf['dp'] = _dp(d['dp'])
    for k in ('c', 'orient'):
        if d.get(k):
            tf[k] = sc(d[k])
    if d.get('T'):
        tf['T'] = d['T']
    return tf


@pytest.mark.parametrize('i', range(len(TF)))
def test_triangleFace(i):
    rec = TF[i]
    a = dec(rec['in'])[0]
    if not isinstance(a, dict) or 'Gn' not in a:
        with pytest.raises(CCA2Error) as e:
            triangleFace(a)
        assert e.value.identifier == rec['err']
        return
    tri = {'Gn': a['Gn'] or [], 'Gd': vec(a['Gd']), 's': vec(a['s']),
           'G': a.get('G') or []}
    if rec['err']:
        with pytest.raises(CCA2Error) as e:
            triangleFace(tri)
        assert e.value.identifier == rec['err']
        return
    got = triangleFace(tri)
    w = _tf(dec(rec['out'])[0])
    assert got['kind'] == w['kind'] and got['orient'] == w['orient']
    assert _close([x for r in got['T'] for x in r], [x for r in w['T'] for x in r])
    assert abs(got['c'] - w['c']) <= 1e-9 * max(1.0, abs(w['c']))
    if got['kind'] == 'rational':
        assert (got['s0n'], got['s0d'], got['Tn'], got['Td']) == (w['s0n'], w['s0d'], w['Tn'], w['Td'])
    else:
        for k in ('k', 'R', 'S1', 'S0'):
            assert got['dp'][k] == w['dp'][k]


@pytest.mark.parametrize('i', range(len(TFS)))
def test_triangleFaceSign(i):
    rec = TFS[i]
    a = dec(rec['in'])
    assert triangleFaceSign(_tf(a[0]), a[1]) == sc(dec(rec['out'])[0])


# ---- cutting the triangles out of a subdivision ---------------------------------------------------
from cca2.cut_triangles import cutTriangles                      # noqa: E402

CT = _load('cutTriangles')


def _aslist(x):
    if x is None or x == []:
        return []
    return x if isinstance(x, list) else [x]


def _face(w):
    return {'kind': w['kind'], 'num': vec(w['num']) if w['num'] else [], 'den': sc(w['den']),
            'FP': io_matlab.polyFromTrace(w['FP'], 3) if w['FP'] else None,
            'FB': sc(w['FB']) if w['FB'] else 0, 'con': w['con'] or [],
            'alg': [{'FP': io_matlab.polyFromTrace(q['FP'], 3), 'FB': sc(q['FB']),
                     'num': vec(q['num']), 'den': sc(q['den']), 'side': sc(q['side'])}
                    for q in _aslist(w['alg'])]}


@pytest.mark.parametrize('i', range(len(CT)))
def test_cutTriangles(i):
    rec = CT[i]
    a = dec(rec['in'])
    sub = [_face(w) for w in _aslist(a[0])]
    tfs = [_tf(t) for t in _aslist(a[1])]
    mode = a[2] if len(a) > 2 else None
    out = dec(rec['out'])
    got, cuts = cutTriangles(sub, tfs, mode)
    if mode == 'lines':
        assert cuts == (out[1] or [])
        return
    want = _aslist(out[0])
    assert len(got) == len(want)
    for g, w in zip(got, want):
        assert g['kind'] == w['kind']
        assert g['con'] == (w['con'] or [])
        assert g['tri'] == (w['tri'] or [])
        if g['kind'] == 'quad':
            assert g['num'] == vec(w['num']) and g['den'] == sc(w['den'])

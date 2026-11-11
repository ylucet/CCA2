"""Unit parity against MATLAB traces: algBranch and isolatedRootSign."""
import math
import os

import pytest

from cca2 import io_matlab
from cca2.io_matlab import dec, sc, vec
from cca2.alg_branch import algBranch
from cca2.isolated_root_sign import isolatedRootSign

from conftest import GOLDEN

TRACE = os.path.join(GOLDEN, 'trace_biconj.jsonl')


def _load(fn):
    return io_matlab.loadTrace(TRACE, fn) if os.path.exists(TRACE) else []


BR = _load('algBranch')
IRS = _load('isolatedRootSign')


def _same(a, b, tol=1e-8):
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


@pytest.mark.parametrize('i', range(len(BR)))
def test_algBranch(i):
    rec = BR[i]
    a = dec(rec['in'])
    P = io_matlab.polyFromTrace(a[0], 3)
    X = a[1] or []
    b = sc(a[2])
    out = dec(rec['out'])
    if len(out) > 1:
        z, info = algBranch(P, X, b, wantInfo=True)
        want = out[1]['allReal']
        for got, w in zip(info['allReal'], want):
            w = vec(w) if w else []
            assert len(got) == len(w)
            assert all(_same(x, y) for x, y in zip(got, w))
    else:
        z = algBranch(P, X, b)
    wz = vec(out[0])
    assert len(z) == len(wz)
    assert all(_same(x, y) for x, y in zip(z, wz))


def _row(x):
    """A traced polynomial row: a double row, or a cell of BigIntegers."""
    if x and isinstance(x[0], list):
        return [int(v) for v in x[0]] if len(x) == 1 else [int(r[0]) for r in x]
    return [int(v) for v in x]


@pytest.mark.parametrize('i', range(len(IRS)))
def test_isolatedRootSign(i):
    rec = IRS[i]
    N, R, iv = dec(rec['in'])
    iv = vec(iv)
    assert isolatedRootSign(_row(N), _row(R), iv) == sc(dec(rec['out'])[0])


# ---- corner faces ---------------------------------------------------------------------------
from cca2.corner_face_alg import cornerFaceAlg          # noqa: E402
from cca2.corner_wedge_rows import cornerWedgeRows      # noqa: E402

CF = _load('cornerFaceAlg')
CW = _load('cornerWedgeRows')


def _closeList(a, b, tol=1e-9):
    return len(a) == len(b) and all(_same(float(x), float(y), tol) for x, y in zip(a, b))


def _sigma(w):
    return {'quad': vec(w['quad']), 'root': sc(w['root']), 'facet': sc(w['facet']),
            's0n': vec(w['s0n']), 'w': vec(w['w']), 'na': sc(w['na'])}


@pytest.mark.parametrize('i', range(len(CF)))
def test_cornerFaceAlg(i):
    rec = CF[i]
    n, d, P, tie, side = dec(rec['in'])
    faces = cornerFaceAlg(vec(n), sc(d), P or [], vec(tie), sc(side))
    want = dec(rec['out'])[0]
    want = want if isinstance(want, list) else [want]
    assert len(faces) == len(want)
    for f, w in zip(faces, want):
        assert f['FP'] == io_matlab.polyFromTrace(w['FP'], 3)
        assert f['FB'] == sc(w['FB'])
        assert f['sigma'] == _sigma(w['sigma'])
        assert f['con'] == (w['con'] or [])
        assert _closeList(f['s'], vec(w['s'])) and _closeList(f['p'], vec(w['p']))
        assert _closeList([x for r in f['M'] for x in r], [x for r in w['M'] for x in r])


@pytest.mark.parametrize('i', range(len(CW)))
def test_cornerWedgeRows(i):
    rec = CW[i]
    f, n, d, tie, side = dec(rec['in'])
    face = {'sigma': _sigma(f['sigma']), 'M': f['M'], 'p': vec(f['p']), 's': vec(f['s'])}
    rows, info = cornerWedgeRows(face, vec(n), sc(d), vec(tie), sc(side))
    assert rows == (dec(rec['out'])[0] or [])

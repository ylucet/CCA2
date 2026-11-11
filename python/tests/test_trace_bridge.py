"""Unit parity against MATLAB traces: the bridge-face producers, call by call."""
import math
import os

import pytest

from cca2 import CCA2Error, io_matlab
from cca2.io_matlab import dec, sc, vec
from cca2.bridge_phi import bridgePhi
from cca2.bridge_face_alg import bridgeFaceAlg
from cca2.bridge_face_rat import bridgeFaceRat
from cca2.bridge_witness import bridgeWitness
from cca2.bridge_arc_cut import bridgeArcCut
from cca2.bridge_member import bridgeMember

from conftest import GOLDEN

TRACE = os.path.join(GOLDEN, 'trace_biconj.jsonl')


def _load(fn):
    return io_matlab.loadTrace(TRACE, fn) if os.path.exists(TRACE) else []


PHI = _load('bridgePhi')
ALG = _load('bridgeFaceAlg')
RAT = _load('bridgeFaceRat')
WIT = _load('bridgeWitness')
CUT = _load('bridgeArcCut')
MEM = _load('bridgeMember')


def _close(a, b, tol=1e-9):
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(_close(x, y, tol) for x, y in zip(a, b))
    if a is None or b is None:
        return a is b
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    if math.isinf(a) or math.isinf(b):
        return a == b
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _primitive(v):
    """A polynomial coefficient list up to a POSITIVE integer factor (MATLAB removes the content
    only when its double path overflowed)."""
    g = 0
    for x in v:
        g = math.gcd(g, abs(int(x)))
    return [int(x) // g for x in v] if g else [int(x) for x in v]


@pytest.mark.parametrize('i', range(len(PHI)))
def test_bridgePhi(i):
    rec = PHI[i]
    a = dec(rec['in'])
    N1, m1, N2, m2, Dl, dl, kp, d, xn, xd = a[:10]
    exact = a[10] if len(a) > 10 else True
    isBig = exact == 'big'
    if not isBig:
        exact = bool(sc(exact))
    out = dec(rec['out'])
    Phi, Dt, U1, U2 = bridgePhi(N1, vec(m1), N2, vec(m2), Dl, vec(dl), sc(kp), sc(d), vec(xn),
                                sc(xd), True if isBig else exact)
    want = out[0]
    if isinstance(want, dict):                       # MATLAB's big path: {sz, c}, content-free
        want = [0 if (c is None or c == []) else int(c) for c in want['c']]
    else:
        want = vec(want)
    if isBig or exact:
        assert _primitive(Phi) == _primitive(want)
    else:
        assert _close([float(v) for v in Phi], [float(v) for v in want], 1e-12)
    if len(out) > 1 and (isBig or exact):
        got = [Dt, U1, U2]
        for g, w in zip(got, out[1:]):
            if isinstance(w, dict):
                w = w['c']
            assert [int(v) for v in g] == [int(v) for v in vec(w)]


def _sameFace(face, w):
    assert face['kind'] == w['kind']
    assert face['FB'] == sc(w['FB'])
    assert face['wit'] == vec(w['wit'])
    cells = [c or [] for c in w['cells']]
    assert face['cells'] == cells
    assert face['FP'] == io_matlab.polyFromTrace(w['FP'], 3)


@pytest.mark.parametrize('i', range(len(ALG)))
def test_bridgeFaceAlg(i):
    rec = ALG[i]
    a = dec(rec['in'])
    args = [vec(a[0]), sc(a[1]), vec(a[2]), sc(a[3])] + ([a[4] or []] if len(a) > 4 else [])
    if rec['err']:
        with pytest.raises(CCA2Error) as e:
            bridgeFaceAlg(*args)
        assert e.value.identifier == rec['err']
        return
    face, info = bridgeFaceAlg(*args)
    out = dec(rec['out'])
    _sameFace(face, out[0])
    if face['kind'] == 'algebraic':
        pr = out[0]['pair']
        assert face['pair']['d'] == sc(pr['d'])
        assert face['pair']['N1'] == pr['N1'] and face['pair']['N2'] == pr['N2']
        assert face['pair']['m1'] == vec(pr['m1']) and face['pair']['m2'] == vec(pr['m2'])
        assert face['pair']['kp'] == sc(pr['kp'])
        assert face['pair']['P'] == (pr['P'] or [])


@pytest.mark.parametrize('i', range(len(RAT)))
def test_bridgeFaceRat(i):
    rec = RAT[i]
    a = dec(rec['in'])
    args = [vec(a[0]), sc(a[1]), vec(a[2]), sc(a[3])] + ([a[4] or []] if len(a) > 4 else [])
    if rec['err']:
        with pytest.raises(CCA2Error) as e:
            bridgeFaceRat(*args)
        assert e.value.identifier == rec['err']
        return
    face, info = bridgeFaceRat(*args)
    w = dec(rec['out'])[0]
    _sameFace(face, w)
    assert face['fN'] == vec(w['fN']) and face['den'] == vec(w['den'])
    for k in ('nu', 'rq', 'rl', 'rc', 'rd', 'scale'):
        assert face['pair'][k] == sc(w['pair'][k])
    assert face['pair']['B'] == vec(w['pair']['B']) and face['pair']['C'] == vec(w['pair']['C'])


@pytest.mark.parametrize('i', range(len(WIT)))
def test_bridgeWitness(i):
    rec = WIT[i]
    a = dec(rec['in'])
    d, N1, m1, N2, m2, Dl, dl, kp, C1, C2, P = a
    args = (sc(d), N1, vec(m1), N2, vec(m2), Dl, vec(dl), sc(kp), vec(C1), vec(C2), P or [])
    if rec['err']:
        with pytest.raises(CCA2Error) as e:
            bridgeWitness(*args)
        assert e.value.identifier == rec['err']
        return
    W, info = bridgeWitness(*args)
    out = dec(rec['out'])
    assert W == out[0]


@pytest.mark.parametrize('i', range(len(CUT)))
def test_bridgeArcCut(i):
    rec = CUT[i]
    a = dec(rec['in'])
    tie, P, N1, m1, N2, m2, d, wit, C1, C2, kp = a
    args = (vec(tie), P or [], N1, vec(m1), N2, vec(m2), sc(d), vec(wit), vec(C1), vec(C2), sc(kp))
    if rec['err']:
        with pytest.raises(CCA2Error) as e:
            bridgeArcCut(*args)
        assert e.value.identifier == rec['err']
        return
    cells, info = bridgeArcCut(*args)
    out = dec(rec['out'])
    assert cells == [c or [] for c in out[0]]


def _pair(pr):
    if pr.get('kind') == 'rational':
        d = {k: sc(pr[k]) for k in ('nu', 'rq', 'rl', 'rc', 'rd', 'scale')}
        d.update(kind='rational', u=vec(pr['u']), w=vec(pr['w']), B=vec(pr['B']), C=vec(pr['C']),
                 P=pr['P'] or [])
        return d
    return {'d': sc(pr['d']), 'N1': pr['N1'], 'm1': vec(pr['m1']), 'N2': pr['N2'],
            'm2': vec(pr['m2']), 'Dl': pr['Dl'], 'dl': vec(pr['dl']), 'kp': sc(pr['kp']),
            'P': pr['P'] or []}


@pytest.mark.parametrize('i', range(len(MEM)))
def test_bridgeMember(i):
    rec = MEM[i]
    a = dec(rec['in'])
    tf, info = bridgeMember(_pair(a[0]), a[1] or [])
    out = dec(rec['out'])
    assert tf == [bool(v) for v in vec(out[0])]

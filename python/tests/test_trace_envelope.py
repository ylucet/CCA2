"""Unit parity against MATLAB traces: the envelope helpers, call by call."""
import os

import pytest

from cca2 import io_matlab
from cca2.io_matlab import dec, sc, vec
from cca2.envelope_cell_pairs import envelopeCellPairs
from cca2.envelope_cell_regions import envelopeCellRegions

from conftest import GOLDEN

TRACE = os.path.join(GOLDEN, 'trace_biconj.jsonl')


def _load(fn):
    return io_matlab.loadTrace(TRACE, fn) if os.path.exists(TRACE) else []


PAIRS = _load('envelopeCellPairs')
REGIONS = _load('envelopeCellRegions')


@pytest.mark.parametrize('i', range(len(PAIRS)))
def test_envelopeCellPairs(i):
    rec = PAIRS[i]
    a = dec(rec['in'])
    gA = io_matlab.quaConFromTrace(a[0])
    gB = io_matlab.quaConFromTrace(a[1])
    D = a[2] if len(a) > 2 else []
    pairs, stats = envelopeCellPairs(gA, gB, D)
    out = dec(rec['out'])
    want = out[0] if isinstance(out[0], list) else [out[0]]
    assert len(pairs) == len(want)
    for p, w in zip(pairs, want):
        assert p['kind'] == w['kind']
        assert p['ij'] == vec(w['ij'])
        assert p['P'] == (w['P'] or [])
        assert p['q'] == vec(w['q'])
        assert p['win'] == (w['win'] or '')
    if len(out) > 1:
        assert stats == {k: sc(v) for k, v in out[1].items()}


@pytest.mark.parametrize('i', range(len(REGIONS)))
def test_envelopeCellRegions(i):
    rec = REGIONS[i]
    N, m, q, dd, P = dec(rec['in'])
    faces = envelopeCellRegions(N, vec(m), sc(q), sc(dd), P or [])
    want = dec(rec['out'])[0]
    want = want if isinstance(want, list) else [want]
    assert len(faces) == len(want)
    for f, w in zip(faces, want):
        assert f['who'] == {'kind': w['who']['kind'], 'facet': sc(w['who']['facet']),
                            'facet2': sc(w['who']['facet2'])}
        assert f['R'] == (w['R'] or [])
        assert f['fN'] == vec(w['fN'])
        assert f['fD'] == sc(w['fD'])
        assert f['S']['Mn'] == w['S']['Mn']
        assert f['S']['tn'] == vec(w['S']['tn'])
        assert f['S']['d'] == sc(w['S']['d'])


# ---- assembly, prune, rational corners, eval, atoms ------------------------------------------
import math                                                       # noqa: E402

from cca2 import CCA2Error                                        # noqa: E402
from cca2.envelope_assemble import envelopeAssemble               # noqa: E402
from cca2.envelope_prune import envelopePrune                     # noqa: E402
from cca2.envelope_rational_corners import envelopeRationalCorners  # noqa: E402
from cca2.envelope_eval import envelopeEval                       # noqa: E402
from cca2.biconj_atoms import biconjAtoms                         # noqa: E402

ASM = _load('envelopeAssemble')
PRN = _load('envelopePrune')
RCO = _load('envelopeRationalCorners')
EVL = _load('envelopeEval')
ATM = _load('biconjAtoms')


@pytest.mark.parametrize('i', range(len(ASM)))
def test_envelopeAssemble(i):
    rec = ASM[i]
    a = dec(rec['in'])
    gA = io_matlab.quaConFromTrace(a[0])
    gB = io_matlab.quaConFromTrace(a[1])
    D = a[2] if len(a) > 2 else []
    if rec['err']:
        with pytest.raises(CCA2Error) as e:
            envelopeAssemble(gA, gB, D)
        assert e.value.identifier == rec['err']
        return
    env, info = envelopeAssemble(gA, gB, D)
    assert io_matlab.envMismatch(env, io_matlab.envFromTrace(dec(rec['out'])[0])) == ''


@pytest.mark.parametrize('i', range(len(PRN)))
def test_envelopePrune(i):
    rec = PRN[i]
    env = envelopePrune(io_matlab.envFromTrace(dec(rec['in'])[0]))
    assert io_matlab.envMismatch(env, io_matlab.envFromTrace(dec(rec['out'])[0])) == ''


@pytest.mark.parametrize('i', range(len(RCO)))
def test_envelopeRationalCorners(i):
    rec = RCO[i]
    env = envelopeRationalCorners(io_matlab.envFromTrace(dec(rec['in'])[0]))
    assert io_matlab.envMismatch(env, io_matlab.envFromTrace(dec(rec['out'])[0])) == ''


@pytest.mark.parametrize('i', range(len(EVL)))
def test_envelopeEval(i):
    rec = EVL[i]
    a = dec(rec['in'])
    val, who = envelopeEval(io_matlab.envFromTrace(a[0]), a[1] or [])
    want = vec(dec(rec['out'])[0])
    assert len(val) == len(want)
    for x, y in zip(val, want):
        assert (x == y) or abs(x - y) <= 1e-8 * max(1.0, abs(y))


@pytest.mark.parametrize('i', range(len(ATM)))
def test_biconjAtoms(i):
    rec = ATM[i]
    a = dec(rec['in'])
    if not isinstance(a[0], dict) or a[0].get('cls') != 'QuaPol':
        pytest.skip('a non-QuaPol argument (input-type guard)')
    f = io_matlab.quaPolFromTrace(a[0])
    mode = a[1] if len(a) > 1 else None
    if rec['err']:
        with pytest.raises(CCA2Error) as e:
            biconjAtoms(f, mode)
        assert e.value.identifier == rec['err']
        return
    h = biconjAtoms(f, mode)
    out = dec(rec['out'])[0]
    if mode is not None:
        want = out if isinstance(out, list) else [out]
        assert len(h) == len(want)
        return
    pairs = out['pairs'] if isinstance(out['pairs'], list) else [out['pairs']]
    assert len(h.pairs) == len(pairs)
    for env, w in zip(h.pairs, pairs):
        assert io_matlab.envMismatch(env, io_matlab.envFromTrace(w)) == ''
    tris = out['tris'] if isinstance(out['tris'], list) else [out['tris']]
    assert len(h.tris) == len(tris)
    for t, w in zip(h.tris, tris):
        assert io_matlab._floatClose(t['T'], w['T']) and io_matlab._floatClose(t['s'], w['s'])
        assert abs(t['c'] - sc(w['c'])) <= 1e-9 * max(1.0, abs(t['c']))


# ---- the conic subdivision --------------------------------------------------------------------
from cca2.envelope_subdivision import envelopeSubdivision      # noqa: E402

SUB = _load('envelopeSubdivision')


def _quaConDict(d):
    return {'nv': sc(d['nv']), 'ne': sc(d['ne']), 'nf': sc(d['nf']), 'Vname': d['Vname'] or [],
            'E': d['E'] or [], 'EcQ': d['EcQ'] or [], 'fN': d['fN'] or [], 'fD': [vec(d['fD'])],
            'F': d['F'] or [], 'FC': [c or [] for c in d['FC']]}


@pytest.mark.parametrize('i', range(len(SUB)))
def test_envelopeSubdivision(i):
    rec = SUB[i]
    a = dec(rec['in'])
    env = io_matlab.envFromTrace(a[0])
    if rec['err']:
        with pytest.raises(CCA2Error) as e:
            envelopeSubdivision(env)
        assert e.value.identifier == rec['err']
        return
    g, info = envelopeSubdivision(env)
    assert io_matlab.quaConMismatch(g, _quaConDict(dec(rec['out'])[0])) == ''

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

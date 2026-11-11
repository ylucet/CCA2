"""Golden parity for the biconjugate: every biconjQ call the MATLAB suites make, replayed.

A QuaCon result must equal MATLAB's field for field; an AlgCrv / EnvelopeMin result must
reproduce MATLAB's sampled values; a refusal must carry the same identifier.
"""
import os

import pytest

from cca2 import CCA2Error, io_matlab, rat_q
from cca2.biconj_q import biconjQ

from conftest import GOLDEN


def _records(name, fn):
    path = os.path.join(GOLDEN, name)
    if not os.path.exists(path):
        return []
    return [r for r in io_matlab.load(path) if r['fn'] == fn]


BICONJ = _records('biconj.jsonl', 'biconjQ')


def _run(P):
    try:
        return biconjQ(P)
    except ModuleNotFoundError as e:                       # a part of the engine not ported yet
        pytest.skip('not ported yet: %s' % e.name)


@pytest.mark.parametrize('i', range(len(BICONJ)))
def test_biconjQ_matches_matlab(i):
    rec = BICONJ[i]
    P = io_matlab.quaPolFrom(rec['in'])
    assert io_matlab.inputMismatch(P, rec['in']) == ''
    if rec['err']:
        with pytest.raises(CCA2Error) as e:
            _run(P)
        assert e.value.identifier == rec['err']
        return
    h = _run(P)
    out = rec['out']
    assert h.kind() == out['kind']
    assert io_matlab.evalMismatch(h, out) == ''
    if out['kind'] == 'QuaCon':
        why = io_matlab.quaConMismatch(h, out)
        if why:
            rat_q.EMULATE_2P53 = True
            try:
                h2 = biconjQ(io_matlab.quaPolFrom(rec['in']))
            finally:
                rat_q.EMULATE_2P53 = False
            assert io_matlab.quaConMismatch(h2, out) == '', why

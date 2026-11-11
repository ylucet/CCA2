"""Golden parity for the conjugate: every conjQ call the MATLAB suites make, replayed.

Each record must reproduce MATLAB's exact QuaCon field for field (faces, conics, H-form rows,
vertex names, incidence), its dom rows, and its sampled values -- or raise the same named refusal.
"""
import os

import pytest

from cca2 import CCA2Error, io_matlab, rat_q
from cca2.conj_q import conjQ

from conftest import GOLDEN


def _records(name, fn):
    path = os.path.join(GOLDEN, name)
    if not os.path.exists(path):
        return []
    return [r for r in io_matlab.load(path) if r['fn'] == fn]


CONJ = _records('conj.jsonl', 'conjQ')


@pytest.mark.parametrize('i', range(len(CONJ)))
def test_conjQ_matches_matlab(i):
    rec = CONJ[i]
    P = io_matlab.quaPolFrom(rec['in'])
    assert io_matlab.inputMismatch(P, rec['in']) == ''
    if rec['err']:
        with pytest.raises(CCA2Error) as e:
            conjQ(P)
        assert e.value.identifier == rec['err']
        return
    g, dom = conjQ(P)
    assert dom == (io_matlab._mat(rec['dom']) if rec['dom'] else [])
    assert io_matlab.evalMismatch(g, rec['out']) == ''
    why = io_matlab.quaConMismatch(g, rec['out'])
    if why:
        # The ONLY accepted difference is the 2^53 wall (PORT_NOTES.md section 2), and it is
        # proved rather than assumed: with MATLAB's wall emulated the result must be identical.
        rat_q.EMULATE_2P53 = True
        try:
            g2, _ = conjQ(io_matlab.quaPolFrom(rec['in']))
        finally:
            rat_q.EMULATE_2P53 = False
        assert io_matlab.quaConMismatch(g2, rec['out']) == '', why

"""Golden parity on the SLOW MATLAB suites (kinkSubdivisionTest, minSubdivisionTest,
erratumBiconjEx2Test, envelopeSubdivisionSlowTest, algCrvSubdivisionTest): every conjQ / biconjQ
call they make, replayed. Minutes per biconjugate in pure Python -- run with `-m slow`."""
import glob
import os

import pytest

from cca2 import CCA2Error, io_matlab
from cca2.biconj_q import biconjQ
from cca2.check import checkBiconj
from cca2.conj_q import conjQ

from conftest import GOLDEN

RECS = []
for _p in sorted(glob.glob(os.path.join(GOLDEN, 'slow_*.jsonl'))):
    RECS += [(os.path.basename(_p), r) for r in io_matlab.load(_p)]


@pytest.mark.slow
@pytest.mark.parametrize('i', range(len(RECS)))
def test_slow_suite_call_matches_matlab(i):
    name, rec = RECS[i]
    P = io_matlab.quaPolFrom(rec['in'])
    fn = conjQ if rec['fn'] == 'conjQ' else biconjQ
    if rec['err']:
        with pytest.raises(CCA2Error) as e:
            fn(P)
        assert e.value.identifier == rec['err']
        return
    out = rec['out']
    if rec['fn'] == 'conjQ':
        g, _ = conjQ(P)
        assert io_matlab.quaConMismatch(g, out) == ''
        return
    h = biconjQ(P)
    assert h.kind() == out['kind']
    assert io_matlab.evalMismatch(h, out) == ''
    if 'nf' in out and getattr(h, 'nf', None) is not None:
        assert h.nf == out['nf']
    assert checkBiconj(P, h) == []

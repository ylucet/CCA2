"""The definition oracle (MATLAB-free): on every fixture the MATLAB suites exercise, the Python
conjugate must equal sup_x <s,x> - f(x) computed pointwise, and the Python biconjugate must be a
convex minorant of f that equals <s,x> - f*(s) at its own gradient."""
import os

import pytest

from cca2 import CCA2Error, io_matlab
from cca2.biconj_q import biconjQ
from cca2.check import checkBiconj, checkConj
from cca2.conj_q import conjQ

from conftest import GOLDEN


def _records(name, fn):
    path = os.path.join(GOLDEN, name)
    if not os.path.exists(path):
        return []
    return [r for r in io_matlab.load(path) if r['fn'] == fn and not r['err']]


CONJ = _records('conj.jsonl', 'conjQ') + _records('biconj.jsonl', 'conjQ')
# A fixture whose caller ASSERTS fIsConvex is excluded: the engine trusts that flag by design
# (biconjQTest exercises exactly that with a deliberately false flag).
BICONJ = [r for r in _records('biconj.jsonl', 'biconjQ') if not r['in'].get('fIsConvex')]


@pytest.mark.parametrize('i', range(len(CONJ)))
def test_conjugate_is_the_sup(i):
    f = io_matlab.quaPolFrom(CONJ[i]['in'])
    g, _ = conjQ(f)
    assert checkConj(f, g) == []


@pytest.mark.parametrize('i', range(len(BICONJ)))
def test_biconjugate_is_the_convex_envelope(i):
    f = io_matlab.quaPolFrom(BICONJ[i]['in'])
    h = biconjQ(f)
    bad = checkBiconj(f, h)
    assert bad == []

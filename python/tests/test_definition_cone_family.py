"""The definition oracle (checkConj / checkBiconj) on the REVIEW GENERATOR'S CONE FAMILY
(review 2026-11-19, finding 5). test_definition.py runs only on the golden fixtures; the
E-quadrant, R5 and R6 fixtures are exactly where the recession blind spot lives, so the same
oracle checks are run on them. A named refusal of the biconjugate (co f = -infinity) is the
right answer there and is not a defect."""
import os
import sys

import pytest

from cca2 import CCA2Error
from cca2.biconj_q import biconjQ
from cca2.check import checkBiconj, checkConj
from cca2.conj_q import conjQ
from cca2.qua_pol import QuaPol

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'review'))
import plq_gen as G  # noqa: E402

CONEFAM = [x for x in G.build() if x['name'].startswith(('E-quadrant', 'R5', 'R6'))]


def _pol(x):
    if x['V']:
        return QuaPol(x['V'], x['E'], x['f'], x['F'])
    return QuaPol(x['f'][0])


@pytest.mark.parametrize('i', range(len(CONEFAM)))
def test_conj_on_the_cone_family(i):
    x = CONEFAM[i]
    f = _pol(x)
    g, _ = conjQ(f)
    assert checkConj(f, g) == []


@pytest.mark.parametrize('i', range(len(CONEFAM)))
def test_biconj_on_the_cone_family(i):
    x = CONEFAM[i]
    f = _pol(x)
    try:
        h = biconjQ(f)
    except CCA2Error:
        return                                    # a named refusal (co f = -infinity) is the right answer
    bad = checkBiconj(f, h)
    assert bad == []

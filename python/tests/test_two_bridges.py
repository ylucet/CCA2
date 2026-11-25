"""Two BRIDGE faces never own one cell (DESIGN.md "Two bridge faces never share a cell").

The argument: a point x in a bridge region has its base point s in the subdifferential of the pair's
envelope at x, which is a single point for almost every x; two distinct bridges of one envelope have
tie arcs meeting in finitely many points, so two bridge regions meet in a set of measure zero.
Checked here two independent ways on the diagonal kink (c = -2, two bridges): the faces' defining
predicate `bridgeMember` on random points, and `envelopeSubdivisionAlg`'s own count of cells it puts
in two bridges (`info['bridgeOverlaps']`). The witness search over the generator's families is
`tests/review/hunt_two_bridges.py` (slow test below; numbers in DESIGN.md).
Mirrors: MATLAB envelopeSubdivisionAlgTest/twoBRIDGEFacesNeverOwnOneCell, Rust tests/two_bridges.rs.
"""
import os
import random
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..')))

from cca2.bridge_member import bridgeMember                       # noqa: E402
from cca2.conj_q import conjQ                                     # noqa: E402
from cca2.envelope_assemble import envelopeAssemble               # noqa: E402
from cca2.envelope_prune import envelopePrune                     # noqa: E402
import cca2.envelope_subdivision_alg as E                         # noqa: E402
from cca2.qua_pol import QuaPol                                   # noqa: E402


def _tri(V, c):
    q = [0, 0, 0, 0, 2, 0, 2, 0, 0, 0]
    q[7] -= c
    q[8] += c
    return QuaPol(V, [[1, 2, 1], [2, 3, 1], [3, 1, 1]], q, [[1, 0], [1, 0], [1, 0]])


def _kink(c):
    gA = conjQ(_tri([[0, 0], [1, 0], [1, 1]], 0))[0]
    gB = conjQ(_tri([[0, 0], [1, 1], [0, 1]], c))[0]
    return envelopePrune(envelopeAssemble(gA, gB)[0])


def test_two_bridge_regions_meet_in_no_sampled_point():
    env = _kink(-2)
    _, bridges = E._splitBridge(env['bridge'])
    assert len(bridges) == 2, 'the fixture must have two bridges, or the test is vacuous'
    assert all(b['kind'] == 'algebraic' for b in bridges), 'membership exactly the sweep'
    rng = random.Random(20261125)
    X = [[rng.uniform(-2, 3), rng.uniform(-2, 3)] for _ in range(20000)]
    M = [bridgeMember(b['pair'], X)[0] for b in bridges]
    assert all(sum(m) > 300 for m in M), 'each bridge must own a sizeable sample'
    assert sum(1 for k in range(len(X)) if M[0][k] and M[1][k]) == 0


def test_the_subdivision_puts_no_cell_in_two_bridges():
    sub, info = E.envelopeSubdivisionAlg(_kink(-2))
    assert any(f['kind'] == 'bridge' for f in sub)
    assert info['bridgeOverlaps'] == 0


@pytest.mark.slow
def test_witness_search_on_the_review_generator(tmp_path):
    """The R3/R6 continuous two-piece families and E-tri: every envelope with >= 2 bridges is
    sampled; none may put a point (or a cell) in two."""
    sys.path.insert(0, os.path.join(HERE, 'review'))
    import plq_gen
    fx = [r for r in plq_gen.build(20261119, 210) if r['name'].startswith(('E-tri', 'R3', 'R6'))]
    fix = tmp_path / 'fx.json'
    out = tmp_path / 'out.jsonl'
    import json
    fix.write_text(json.dumps(fx))
    subprocess.run([sys.executable, os.path.join(HERE, 'review', 'hunt_two_bridges.py'),
                    str(fix), str(out), '-j', '8', '--timeout', '300', '--samples', '1000'],
                   check=True)
    L = [json.loads(l) for l in out.read_text().splitlines()]
    done = [l for l in L if 'envs' in l]
    assert sum(l['envs2b'] for l in done) >= 10, 'the search must reach envelopes with two bridges'
    assert [l['id'] for l in done if l['inTwo'] or l['bridgeOverlaps']] == []

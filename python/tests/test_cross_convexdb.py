"""Cross-check of the engine against the convexdb registry (outline gap G3) -- fast subset.

Drives the bench script's own machinery (python/bench/cross_convexdb.py) over a fast,
representative subset of the 2-D PLQ registry entries, using the Python engine, and asserts
the cross-check verdicts:

  * no entry may DISAGREE (a disagreement is a real engine-vs-registry defect);
  * the expected non-PLQ entry (norm1_2d = |x|+|y|) is refused by the parser;
  * the plain accepted entries (quadratic / bilinear / rotation / multi-face) agree.

This is the fast regression guard. The full 34-entry set (both engines, the slow/heavy
entries, and the double-protocol / open entries) is the bench run itself; this test keeps
the agreement invariant in the fast bucket without the ~1 min heavy entries or the Rust
binary.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'bench'))

import cross_convexdb as cc  # noqa: E402

cc.bootstrap_pythonpath()

# A fast, representative slice of cc.FAST_SLUGS: every shape family (quadratic triangle,
# bilinear, rotation, multi-face, unbounded) plus the one expected refusal. The heavier
# FAST_SLUGS members (three_convex_edges_split, the erratum hyperbola/circle) are covered
# by the bench run; keeping them out keeps this test in the seconds.
TEST_SLUGS = [
    'quad_psd_triangle',
    'quad_nsd_triangle',
    'bilinear_zero_convex_edges',
    'bilinear_one_convex_edge_rational',
    'general_indefinite_rotation',
    'negative_bilinear_rotation',
    'multi_face_bilinear_square',
    'quad_norm2_half_2d',
    'karmarkar_erratum_unbounded',
    'norm1_2d',
]


@pytest.fixture(scope='module')
def registry():
    try:
        path = cc.find_registry_json()
    except SystemExit:
        pytest.skip('convexdb registry not found (sibling folder missing)')
    _ver, entries = cc.load_entries(path)
    return entries


def _entry(registry, slug):
    return [e for e in registry if e.slug == slug and e.dim == 2][0]


def test_registry_has_fast_slugs(registry):
    slugs = {e.slug for e in registry if e.dim == 2}
    for s in TEST_SLUGS:
        assert s in slugs, 'test slug %s is not in the 2-D registry' % s


@pytest.mark.parametrize('slug', TEST_SLUGS)
def test_no_disagreement(registry, slug):
    rec = cc.process_entry(_entry(registry, slug), ['py'], None)
    for op in ('conj', 'biconj'):
        st = rec['py'][op]['status']
        assert st != 'disagree', (
            'engine disagrees with the registry on %s/%s: %s'
            % (slug, op, rec['py'][op]['detail']))
        assert st in ('agree', 'refuse', 'open'), (slug, op, st)


def test_non_plq_refused(registry):
    rec = cc.process_entry(_entry(registry, 'norm1_2d'), ['py'], None)
    assert not rec['parse'].startswith('ok'), rec['parse']
    assert rec['py']['conj']['status'] == 'refuse'
    assert rec['py']['biconj']['status'] == 'refuse'


def test_accepted_entries_agree(registry):
    for slug in ('quad_psd_triangle', 'quad_nsd_triangle', 'multi_face_bilinear_square',
                 'quad_norm2_half_2d'):
        rec = cc.process_entry(_entry(registry, slug), ['py'], None)
        assert rec['parse'].startswith('ok'), (slug, rec['parse'])
        assert rec['py']['conj']['status'] == 'agree', (slug, rec['py']['conj'])
        assert rec['py']['biconj']['status'] == 'agree', (slug, rec['py']['biconj'])

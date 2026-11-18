"""Star polygons (n = 4..7 here; the Rust test rust/tests/star.rs goes further) (2n vertices, f = (x^2+y^2)/2): conjQ must equal the definition
f*(s) = sup_star <s,x> - f(x), computed independently (interior point x = s, else best boundary point).
Regression for the incomplete cell-interior witness search (DECISIONS.md 2026-11-17)."""
import sys, os
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from star_defect_diagnosis import star, truth      # noqa: E402
from cca2.conj_q import conjQ                      # noqa: E402


@pytest.mark.parametrize('n', [4, 5, 6, 7])
def test_star_conjugate_matches_definition(n):
    g = conjQ(star(n))
    g = g[0]
    rng = np.random.RandomState(n)
    pts = [[float(a), float(b)] for a, b in rng.uniform(-150, 150, (100, 2))]
    pts += [[-36.39, -5.97]]
    vals = list(g.eval(pts)[0])
    for s, v in zip(pts, vals):
        t = truth(n, s)
        assert abs(v - t) <= 1e-6 * max(1.0, abs(t)), (n, s, v, t)

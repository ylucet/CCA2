# cca2 (Python) — exact conjugate and biconjugate of bivariate PLQ functions

Python port of the **exact engine** of the MATLAB toolbox CCA2: `QuaPol.conj` (`conjQ`, returns a
`QuaCon`) and `QuaPol.biconj` (`biconjQ`, returns an `AlgCon`: `QuaCon`, `AlgCrv` or
`EnvelopeMin`). Why Python: `../DECISION/DECISION.md`. Plan: `PORT_PLAN.md`. Every deliberate
difference from MATLAB, and what is not ported: `PORT_NOTES.md`. Measured scope: `PORT_SCOPE.md`.

```python
from cca2.qua_pol import QuaPol

# x*y on the unit triangle: rows of f are (c5 x^2 + 2 c6 xy + c7 y^2)/2 + c8 x + c9 y + c10
f = QuaPol([[0, 0], [1, 0], [0, 1]],            # V
           [[1, 2, 1], [2, 3, 1], [3, 1, 1]],   # E: [v1 v2 flag], 1-based, flag 1 = segment
           [[0, 1, 0, 0, 0, 0]],                # f
           [[1, 0], [1, 0], [1, 0]])            # F: [left face, right face], 0 = outside
g = f.conj()        # QuaCon: exact faces fN/fD, integer edge conics EcQ, named vertices
h = f.biconj()      # the closed convex envelope
print(g.eval([[0.3, 0.4]]), h.eval([[0.25, 0.25]]))
```

Exact data are Python `int`s (no 2^53 wall); the compute path imports no CAS; numpy is used only
for the documented double-precision screens.

## Tests

```
PYTHONPATH=.pydeps python3 -m pytest -q -n 32            # fast bucket (~1.5 min, ~2700 tests)
PYTHONPATH=.pydeps python3 -m pytest -q -n 8 -m slow     # slow suites' fixtures (minutes each)
```

`pip3 install --target .pydeps pytest-xdist` once (the folder is git-ignored). The tests replay
what MATLAB recorded in `golden/` — regenerate with `tools/export_golden.sh` / `tools/export_trace.sh`
(needs MATLAB) — and judge every result against the definition, MATLAB-free (`cca2/check.py`).

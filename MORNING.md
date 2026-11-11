# Morning report — 2026-11-10 overnight run

Branch: overnight/2026-11-10 (from main f62c7dc)

Task: execute `python/PORT_PLAN.md` -- port the exact conjugate / biconjugate engine (conjQ, biconjQ)
from MATLAB to Python. TODO.md "Next up — 2026-11-10 overnight run", items P0..P8.

## What changed
- Checkpoint: `DECISION/` and `PORT/` (the Python-vs-Julia decision and the exactQ probe) were
  untracked on main; committed on this branch. `test.txt` (empty, created in the session) left untracked.

- P0 measured scope: `python/PORT_SCOPE.md` (line coverage of the exact engine under the end-to-end
  suites). The goldens are recorded by `python/tools/export_golden.sh`: a shadow wrapper logs every
  conjQ / biconjQ call the MATLAB suites make (input, exact output or refusal id) as JSON lines.
- P3 DONE: the exact conjugate runs in Python. `pytest tests/test_golden_conj.py`: 120/120 conjQ
  calls reproduce MATLAB -- 119 identical field for field (faces, conics, H-form rows, vertex names,
  incidence, dom rows, sampled values), 1 differs only where MATLAB's 2^53 overflow fallback keeps
  an EMPTY cell; proved by emulating the wall (`rat_q.EMULATE_2P53`), after which it is identical.
- P4 DONE: biconjugate direct cases (convex -> f, edge-concave lower hull incl. recession rays, thin
  domains, chains) in Python; all refusals carry MATLAB's identifiers.
- P5 + P7 DONE: atoms, envelope assembly/prune/eval, bridge faces (algebraic and rational), corner
  faces. Validated function by function against MATLAB: `python/tools/export_trace.sh` records the
  arguments and results of any top-level MATLAB function while suites run (1583 calls of 25
  functions); the Python twin must reproduce each call. ~1500 such calls are green (implicit
  polynomials FP, branch indices, witnesses, cells, wedge rows identical; doubles to 1e-9).
- P6 (most of it): envelopeSubdivision, envelopeSubdivisionAlg (single-pair path), assembleAlgCrv,
  AlgCrv, algBranch, algFaceBelow. `biconj` is END-TO-END in Python: all 61 recorded biconjQ calls
  reproduce MATLAB (43 QuaCon identical field for field, 5 AlgCrv by kind + sampled values, 13
  refusals by identifier).
- `algRowFails` (the 1100-line emptiness certificate) ported, with its bivariate gcd / Yun
  square-free kernel (`biv.py`): 135/135 recorded MATLAB calls agree. With it, 4 of the 6 recorded
  algebraic subdivisions are identical to MATLAB face for face; the 5th fixture (recorded twice)
  differs by 3 faces that MATLAB keeps in a ZERO-WIDTH cell (see What is broken).
- P1b DONE: `cca2/check.py`, a MATLAB-free definition oracle. `supAt` computes sup <s,x> - f(x)
  pointwise from KKT candidates (shares nothing with conjQ). On all 184 fixtures: the Python
  conjugate equals that sup at 49 dual points each, and the Python biconjugate is a convex minorant
  of f equal to <s,x> - f*(s) at its own gradient.
- Whole Python suite: `cd python && PYTHONPATH=.pydeps python3 -m pytest -q -n 32` -> 1525 passed.

## What is broken
- MATLAB (not Python): `ratQ.sturmChain` loses the Sturm property when a chain element has a
  NEGATIVE leading coefficient and the degree drops by 2 (polyPseudoRemDouble scales by
  lc^(gap+1) with its sign; bigPrem uses |lc|). Measured on Frances 2026-11-10:
  `ratQ.countRootsIn([1 0 0 1 -1], -2, 2)` returns 0; x^4 + x - 1 has 2 real roots there. The
  Python port restores the sign and returns 2. See Needs a decision.
- MATLAB (not Python): `envelopeSubdivisionAlg/liveCell` accepts a cell of ZERO WIDTH. On the
  fixture with tie conic (2x-y)^2 = 4, `cellSamplePoint` returns "ok" for the cell
  {(2x-y)^2 >= 4, 2x-y >= 0, branch line 2x-y-2 <= 0}: its sample (1.5026.., 1.0052..) lies exactly
  ON the line 2x-y = 2 and rounding gives a positive slack on both sides (measured on Frances,
  2026-11-10). MATLAB then emits 3 bridge faces for that cell (29 faces; Python 26). The faces own
  no point, so no value changes -- this is very likely the "zero-area faces" item in TODO's
  open list. Python finds no interior point and drops the cell.

## Needs a decision
- **Fix the Sturm sign bug in MATLAB `ratQ.m`?** Options: (a) fix it now in `polyPseudoRemDouble`
  and `polyPseudoRemBig` (scale by |lc|^(gap+1), as `bigPrem` already does) and rerun the suites
  -- recommended, it is a wrong root COUNT, silent; (b) leave MATLAB as is since the port
  supersedes it. Not done tonight: it changes the MATLAB engine, which was not the task.
- **`conjQ/linealityConj` builds the cross-section from `fN` and drops `fD`** (python/PORT_NOTES.md
  3.2). Looks wrong whenever fD ~= 1; the port uses fN/fD. Not yet confirmed by a MATLAB run.

## Where I stopped
P0-P5, P7, P1b done; P6 done except the multi-pair overlay (minOverlay), triangle faces (triangleFace / cutTriangles) and rows against faces of degree >= 3 in z -- all three decline to an EnvelopeMin / "no proof" today, so values are right and only the subdivision is missing.

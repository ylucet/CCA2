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

## What is broken
- MATLAB (not Python): `ratQ.sturmChain` loses the Sturm property when a chain element has a
  NEGATIVE leading coefficient and the degree drops by 2 (polyPseudoRemDouble scales by
  lc^(gap+1) with its sign; bigPrem uses |lc|). Measured on Frances 2026-11-10:
  `ratQ.countRootsIn([1 0 0 1 -1], -2, 2)` returns 0; x^4 + x - 1 has 2 real roots there. The
  Python port restores the sign and returns 2. See Needs a decision.

## Needs a decision
- **Fix the Sturm sign bug in MATLAB `ratQ.m`?** Options: (a) fix it now in `polyPseudoRemDouble`
  and `polyPseudoRemBig` (scale by |lc|^(gap+1), as `bigPrem` already does) and rerun the suites
  -- recommended, it is a wrong root COUNT, silent; (b) leave MATLAB as is since the port
  supersedes it. Not done tonight: it changes the MATLAB engine, which was not the task.
- **`conjQ/linealityConj` builds the cross-section from `fN` and drops `fD`** (python/PORT_NOTES.md
  3.2). Looks wrong whenever fD ~= 1; the port uses fN/fD. Not yet confirmed by a MATLAB run.

## Where I stopped
P0, P3, P4 done; starting P5 (biconjAtoms, envelope*).

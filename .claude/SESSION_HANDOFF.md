# Session Handoff

_2026-09-25_

## Blocked
Nothing.

## State
- Branch `main` @ `6aa5ed9` — "fix: every explicit sum of products under ratQ.chk ..."
- Pushed: pending (12 commits ahead of origin/main)
- Tests: fast 723/0 · normal 24/0 (09-25) · slow 100/0 (09-23)
- Known reds: none; 2 kink tests in biconjCPLQTest quarantined (runtime)

## Next
1. ~118 `chk` calls wrapping a matrix/dot product — same silent-sign risk.
   TODO "exact-sum audit"; `ratQ.chkB` is the tool.
2. biconjQ on poly2 kink stops at `envelopeMesh:bridgeFaces` (corner face,
   no rows). ~90 s repro; TODO, DECISIONS 2026-09-25.
3. Zero-area faces of the c=-2 AlgCrv: drop where the envelope assigns
   quadratic faces. TODO.

## Files
- `ratQ.m` — chkB, mulChk, BigInteger fallbacks
- `algCrvIncidence.m` — arcs traced between vertices
- `exactSumLintTest.m` — guards the fixed pattern
- `algBranch.m` — double-root fibre fix

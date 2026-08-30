# Session Handoff

_2026-08-30_

## Blocked
Nothing external.

## State
- Branch `main` @ `9037a35` — "docs: the vertex-cone fix moves the
  3-piece witness past its old failure into a new one"
- Pushed: yes — `origin/main`. Tags `v0.2`, `v0.2.1` (memoizations +
  scip-09's real-QPLIB validation).
- Tests (2026-08-30): fast 312/0 · normal 12/0 · slow (unboundedFaceTest
  only, 19/0) · earlier full slow 96/0 (2026-08-29, unaffected files).
- Known reds: none silent. G11 quarantined via `assumeFail`.
- scip-09 validated the box envelope on real QPLIB data (0 error,
  837 cases); their remaining wall-clock loss is their own cut mgmt.

## Next
1. **NEW `MATLAB:badsubscript`** in the 3-piece elliptical-edge witness
   (`doc/QuaConExample.md`), reached only after today's vertex-cone fix
   moved it past the old `cplqFailed`. No stack trace captured yet —
   get one first. `SUPPORT_MATRIX.md` 1.2 / G16 still open.
2. Fold-strategy redesign (TODO G4): hub-vertex root cause confirmed
   twice (extra `mergeL` passes: 58->58->58). Needs an N-ARY fan merge
   (new soundness proof) or fixing fragmentation upstream. Not attempted.
3. `splitTwoArcLens`: reproducer FOUND, not wired into a test — target
   cell is an INTERMEDIATE `clipByFace` product, trace forward.

## Files
- `conjConvexOverPiece.m` — REAL BUG FIXED: vertex cone could collapse
  to a point (`edgeDirsAt` probe step vs `tolA`); see DECISIONS 08-30
- `region.m` — THREE memoizations (`unionIsExact`, `getVertices`,
  `simplifyUnboundedRegion`), ~38% faster on the reference fixture
- `unboundedFaceTest.m` — new regression test for the fix above (19/0)
- `DECISIONS.md`, `TODO.md` — extensive, all evidenced

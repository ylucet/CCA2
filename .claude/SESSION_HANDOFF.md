# Session Handoff

_2026-08-30_

## Blocked
Nothing external.

## State
- Branch `main` @ `d49e90c` — "test: fuzz conjConvexPolygon (the
  SCIP-relevant path) across scales -- clean"
- Pushed: yes — `origin/main`. Tags `v0.2`, `v0.2.1` (memoizations +
  scip-09's real-QPLIB validation).
- Tests (2026-08-30): fast 314/0 · normal 12/0 · slow 98/0 (1003s).
  `conjConvexPolygon` (the ACTUAL SCIP box-envelope path, different
  code from today's 2 fixes) independently fuzzed clean, 0/4500.
- Known reds: none silent. G11 quarantined via `assumeFail`.
- scip-09 validated the box envelope on real QPLIB data (0 error,
  837 cases); their remaining wall-clock loss is their own cut mgmt.

## Next
1. Fold-strategy redesign (TODO G4): hub-vertex root cause confirmed
   twice (extra `mergeL` passes: 58->58->58). Needs an N-ARY fan merge
   (new soundness proof) or fixing fragmentation upstream. Not attempted.
2. `splitTwoArcLens`: reproducer FOUND, not wired into a test — target
   cell is an INTERMEDIATE `clipByFace` product, trace forward.
3. `QuaCon`/G16: precondition now genuinely MET (a real elliptical edge
   reaches Step 3 and is handled correctly) — now an EFFICIENCY project
   (H-form), not a correctness gap. Still not started; nothing blocks it.

## Files
- `conjConvexOverPiece.m`, `symbolicFunction.m` — TWO real bugs fixed
  tracing `SUPPORT_MATRIX.md` 1.2 (vertex-cone collapse; tangent on a
  curve missing one ambient var). 3-piece elliptical witness now
  verifies end to end — see DECISIONS 2026-08-30 (three entries).
- `region.m` — THREE memoizations (`unionIsExact`, `getVertices`,
  `simplifyUnboundedRegion`), ~38% faster on the reference fixture
- `DECISIONS.md`, `TODO.md` — extensive, all evidenced

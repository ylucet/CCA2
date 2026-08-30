# Session Handoff

_2026-08-30_

## Blocked
Nothing external.

## State
- Branch `main` @ `92088e9` — "docs: G16's trigger is SUPPORT_MATRIX
  1.2, not G1 -- checked directly, G1 didn't change it"
- Pushed: yes — `origin/main`. Tags `v0.2`, `v0.2.1` (memoizations +
  scip-09's real-QPLIB validation).
- Tests (2026-08-29, no code changes since): fast 312/0 · normal 12/0 ·
  slow 96/0 (493s) · regionTest 27/0 · testcPLQ (verylong) 8/0/1 (G11).
- Known reds: none silent. G11 quarantined via `assumeFail`.
- scip-09 validated the box envelope on real QPLIB data (0 error,
  837 cases); their remaining wall-clock loss is their own cut
  management, not CCA2.

## Next
1. Fold-strategy redesign (TODO G4): hub-vertex root cause confirmed
   TWICE — extra `mergeL` passes on the final result find nothing
   (58->58->58), ruling out a merge-order fix. Real fix needs an N-ARY
   fan merge (new soundness proof) or fixing the fragmentation upstream
   in Steps 1/2. Genuine redesign, not attempted.
2. `splitTwoArcLens`: reproducer FOUND. Target cell is an INTERMEDIATE
   `clipByFace` product — trace forward, not backward.
3. `QuaCon`/G16 blocked on `SUPPORT_MATRIX.md` 1.2 (checked directly,
   not on G1 as its own text claimed) — re-check only after 1.2 moves.

## Files
- `region.m` — THREE memoizations this session (`unionIsExact`,
  `getVertices`, `simplifyUnboundedRegion`); 5 `Inf`/needle bugs fixed
  earlier in `maxAffineOverRegion`
- `maxQuaPar.m` — G1/G4/G10 landed
- `testcPLQ.m` — G11 quarantined via `assumeFail`
- `SCIP_READINESS.md` — gate rewritten, items 1/2/4 met
- `DECISIONS.md`, `TODO.md` — extensive, all evidenced

# Session Handoff

_2026-08-29_

## Blocked
Nothing external.

## State
- Branch `main` @ `82caedf` — "docs: record the unionIsExact memoization
  win in TODO.md's G4 entry"
- Pushed: yes — `origin/main`. Tag `v0.2` cut and pushed.
- Tests (2026-08-29): fast 312/0 · normal 12/0 · slow 96/0 · regionTest
  27/0 · testcPLQ (verylong) 8/0/1 (G11 quarantined).
- Known reds: none silent. G11 quarantined via `assumeFail`.

## Next
1. Fold-strategy redesign (TODO G4): root cause found (high-degree hub
   vertex, 8 cells meet pre-fold). `unionIsExact` memoized as a safe,
   modest mitigation (2186s vs 2226-2546s) — the real fix (resolve each
   hub's fan once, not pairwise) is still open.
2. `splitTwoArcLens`: reproducer FOUND. Target cell is an INTERMEDIATE
   `clipByFace` product — trace forward from operands, not backward.
3. Conic-conic closed form: `conicMeet.m`/`ratQ.m` exist, tested, unused —
   built for `QuaCon`, uncoded.
4. Measure non-box timing against a real QPLIB constraint (spike/SCIP's
   side now — notified via cross-session message, session scip-09).

## Files
- `region.m` — `unionIsExact` memoized (this session); 5 real `Inf`/needle
  bugs fixed in `maxAffineOverRegion` (earlier this session)
- `maxQuaPar.m` — G1/G4/G10 landed
- `testcPLQ.m` — G11 quarantined via `assumeFail`
- `SCIP_READINESS.md` — gate rewritten, items 1/2/4 met
- `DECISIONS.md`, `TODO.md` — extensive, all evidenced

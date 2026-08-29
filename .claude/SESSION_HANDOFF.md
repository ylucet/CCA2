# Session Handoff

_2026-08-28_

## Blocked
Nothing external.

## State
- Branch `main` @ `44d0ee9` — "docs: mark G1/G4/G10 landed in TODO.md"
- Pushed: yes — `origin/main`
- Tests (2026-08-28): fast 312/0 · normal 12/0 · maxQuaParTest 31/0 ·
  conjCPLQTest 33/0 · conjEdgeLowerBoundTest 5/0. Slow/verylong not run in
  full.
- Known reds: `testcPLQ/rectBiconjugateIsAConvexUnderestimator` (verylong,
  SCIP-relevant, see SCIP_READINESS.md Phase C)

## Next
1. G1/G4/G10 LANDED (assemblePieces diff, trade-off in DECISIONS.md
   2026-08-28). Case 21 now hits a known Step-3 legacy gap, not new — don't
   re-derive.
2. Scaling defect (G4): sliver witness has NO valid merge partner, not a
   bug (mergeL confirms different function than its neighbour). Sample 2-3
   more fold-5 cells the same way to see if this generalizes.
3. `splitTwoArcLens`: reproducer FOUND (random search, large curvature
   ratio). Needs a real `maxQuaPar(g1,g2)` regression, not built yet.
4. Conic-conic closed form: `conicMeet.m`/`ratQ.m` exist, tested, unused —
   built for `QuaCon`, uncoded.

## Files
- `maxQuaPar.m` — G1/G4/G10 fix landed (collapseTinyEdges + sagitta test)
- `region.m` — 5 real `Inf`/needle bugs fixed in `maxAffineOverRegion`
- `conjCPLQ.m`, `conjCPLQTest.m`, `conjEdgeLowerBoundTest.m` — fixes + tests
  updated for the case-21 trade-off (same fixture, both files)
- `DECISIONS.md`, `TODO.md` — extensive, all findings evidenced

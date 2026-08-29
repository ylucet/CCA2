# Session Handoff

_2026-08-29_

## Blocked
Nothing external.

## State
- Branch `main` @ `965cfd9` — "docs: item 1 generalized -- 141 refusals
  concentrate into ~12 cell signatures"
- Pushed: yes — `origin/main`
- Tests (2026-08-28): fast 312/0 · normal 12/0 · maxQuaParTest 31/0 ·
  conjCPLQTest 33/0 · conjEdgeLowerBoundTest 5/0. Slow/verylong not run.
- Known reds: `testcPLQ/rectBiconjugateIsAConvexUnderestimator` (verylong,
  SCIP-relevant, see SCIP_READINESS.md Phase C)

## Next
1. G1/G4/G10 LANDED (assemblePieces diff, trade-off in DECISIONS.md
   2026-08-28). Case 21 hits a known Step-3 legacy gap, not new — don't
   re-derive.
2. Scaling defect: the 58-vs-8 surplus traces to ~12 cell signatures
   (89% in top 5), not a diffuse problem. Trace THEIR histories next
   (DECISIONS.md 2026-08-29) — tractable, not a redesign.
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

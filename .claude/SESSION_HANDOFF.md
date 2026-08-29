# Session Handoff

_2026-08-29_

## Blocked
Nothing external.

## State
- Branch `main` @ `aa7b5cd` — "docs: gate item 4 MET -- slow bucket 96/0,
  verylong red now quarantined by name"
- Pushed: yes — `origin/main`
- Tests (2026-08-29): fast 312/0 · normal 12/0 · slow 96/0 · testcPLQ
  (verylong) 8/0/1 (G11 quarantined, not hanging).
- Known reds: none silently red. G11 (`testcPLQ/rectBiconjugateIsAConvexUnderestimator`)
  quarantined via `assumeFail`; delete that line once fold-strategy work lands.

## Next
1. Fold-strategy redesign (TODO G4): root cause found (high-degree hub
   vertex, 8 cells meet pre-fold). Fix resolves each hub's fan once, not
   pairwise — a real algorithm change.
2. `splitTwoArcLens`: reproducer FOUND. Target cell is an INTERMEDIATE
   `clipByFace` product — trace forward from operands, not backward.
3. Conic-conic closed form: `conicMeet.m`/`ratQ.m` exist, tested, unused —
   built for `QuaCon`, uncoded.
4. Measure non-box timing against a real QPLIB constraint, not just the
   one reference fixture used all session.

## Files
- `maxQuaPar.m` — G1/G4/G10 landed this session
- `region.m` — 5 real `Inf`/needle bugs fixed in `maxAffineOverRegion`
- `conjCPLQ.m`, `conjCPLQTest.m`, `conjEdgeLowerBoundTest.m`, `testcPLQ.m`
  — case-21 trade-off + G11 quarantine
- `SCIP_READINESS.md` — gate section fully rewritten, items 1/2/4 met
- `DECISIONS.md`, `TODO.md` — extensive, all evidenced

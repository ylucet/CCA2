# Session Handoff

_2026-08-28_

## Blocked
Nothing external.

## State
- Branch `main` @ `6a8e6dd` — "docs: correct my own overconfident claim --
  the edge-adjacency pre-filter is unsafe"
- Pushed: pending
- Tests (2026-08-28): fast 312/0 · normal 12/0 · regionTest 27/0 ·
  conjCPLQTest 33/0. Slow/verylong not re-run in full this session.
- Known reds: `testcPLQ/rectBiconjugateIsAConvexUnderestimator` (verylong,
  times out — now known SCIP-relevant, see SCIP_READINESS.md Phase C)

## Next
1. Scaling defect (TODO.md G4): real cause is upstream of `unionIsExact` —
   trace `region.merge`'s candidate generation for a sliver's TRUE neighbour
   (shared arc segment, not just curve equation). DECISIONS.md 2026-08-28
   has the witness + a proven-unsafe fix not to retry.
2. `assemblePieces` (G1/G4/G10): the documented globally-consistent
   redesign, not a patch. Acceptance test landed; parked diff's crash is a
   known Step-3 legacy gap, not new — don't re-derive.
3. `splitTwoArcLens`: 4 hand-built reproducer attempts failed. Needs
   `clipPolyByConic`'s real output on a failing fixture, not more guesses.

## Files
- `region.m` — 5 real `Inf`/needle bugs fixed in `maxAffineOverRegion`
- `conjCPLQ.m` — fixed double-conjugating +-inf built a needle, not -inf
- `DECISIONS.md`, `TODO.md` — extensive: fixes, 2 refuted hypotheses, 1
  caught-before-shipping unsafe idea, all with evidence

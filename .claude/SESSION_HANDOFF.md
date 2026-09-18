# Session Handoff

_2026-09-17_

## Blocked
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — delete / track /
  note? Your call. MORNING.md.
- `maxQuaPar:notDegenerate` reachability — an in-contract pair needs the symbolic
  pipeline (hours). DECISIONS 2026-09-06.

## State
- Branch `main`, six commits ahead of `origin/main`, tree clean
- Tests: fast bucket green (bridgeFaceAlgTest 19, ratQTest 80, guardExerciseTest 18,
  QuaConTest/conjQTest/biconjQTest/maxQTest/regionTest 143 together) ·
  normal 13/0 (09-11) · slow 98/0 · verylong 37/0 (09-06/07, not re-run since)
- Known reds: none

## Next
1. **A question about `eval`, before any more geometry.** The bridge face's cut
   boundary is a SEGMENT — a bounded piece of a line — and `QuaCon.FC` can only say
   "this side of this whole curve". Two constructions over sign conditions have now
   been built and measured wrong for exactly that reason. `AlgCon`'s classic mesh
   already has bounded edges. So: can a face's membership be decided from the
   bounded-edge mesh rather than from `FC`? TODO.md 2026-09-17.
2. Then the cut for parabola and hyperbola ties — 23 of the 29 bridge pairs the
   standard fixtures produce.
3. Then the `rank1` / `affine` cell producer, and only then the fold into one
   `AlgCon`.

## Files
- `bridgeFaceAlg.m` — the producer: FP, FB (top branch), region, optional cut
- `bridgeArcCut.m` — the cut; LINE tie exact, conic tie refused with BOTH refutations
- `bridgeWitness.m`, `bridgePhi.m` — shared by the producer and the cut
- `ratQ.pairLine*`, `QuaCon.FC`'s third column — an edge may name one line of a pair

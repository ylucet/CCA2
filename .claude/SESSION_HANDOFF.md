# Session Handoff

_2026-09-17_

## Blocked
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — delete / track /
  note? Your call. MORNING.md.
- `maxQuaPar:notDegenerate` reachability — an in-contract pair needs the symbolic
  pipeline (hours). DECISIONS 2026-09-06.
- **The bridge face's REGION for a conic tie** — needs a type decision only you can
  make: should an `Ec` edge be able to name a BRANCH (which of a degenerate conic's
  two conjugate lines), the way `FB` names a face's branch? Everything left on the
  envelope waits on it. TODO.md 2026-09-17, DECISIONS.md same date.

## State
- Branch `main` @ `3a65c99` — "the CONIC-tie cut was wrong -- measured, and refuted"
- Pushed: NO — five commits ahead of `origin/main`
- Tests: fast 645/0 before the last two commits; bridgeFaceAlgTest 19/0,
  guardExerciseTest 17/0, envelopeCellRegionsTest 8/0, algBranchTest 8/0,
  AlgConTest 8/0 after · normal 13/0 (09-11) · slow 98/0 · verylong 37/0
  (09-06/07, not re-run since)
- Known reds: none

## Next
1. THE TYPE DECISION above. Nothing else on the envelope moves until it is made.
2. Then the cut for parabola and hyperbola ties — 23 of the 29 bridge pairs the
   fixtures produce, against 1 LINE and 2 ellipse.
3. Then the producer for `rank1` and `affine` cells (value side already works;
   region is a fan bounded by the polar conic), and only then the fold into one
   `AlgCon`.

## Files
- `bridgeFaceAlg.m` — the producer: FP (a cubic's discriminant in the multiplier),
  FB (the TOP branch, negative index), region, and the optional cut
- `bridgeArcCut.m` — the cut; LINE tie exact, conic tie refused with the refutation
- `envelopeCellRegions.m` — each candidate carries `.S`, its maximiser map
- `algBranch.m`, `AlgCon.m` — FB < 0 counts from the top; positive unchanged

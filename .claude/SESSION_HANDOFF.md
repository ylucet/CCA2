# Session Handoff

_2026-09-17_

## Blocked
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — delete / track /
  note? Your call. MORNING.md.
- `maxQuaPar:notDegenerate` reachability — an in-contract pair needs the symbolic
  pipeline (hours). DECISIONS 2026-09-06.

## State
- Branch `main` @ `46d0eea` — "envelopeCellRegions returns the MAXIMISER, and two
  measurements that scope the assembly"
- Pushed: NO — two commits ahead of `origin/main`
- Tests: fast 641/1 before this session's work, and the 1 was the guard-inventory
  pin catching the new identifiers; now 654/0 across the suites touched
  (bridgeFaceAlgTest 13, algBranchTest 8, AlgConTest 8, envelopeCellRegionsTest 8,
  envelopeCellPairsTest 10, guardExerciseTest 17) · normal 13/0 (09-11) ·
  slow 98/0 · verylong 37/0 (09-06/07, not re-run since)
- Known reds: none

## Next
1. THE CUT — restrict the bridge face's sweep to the tie arc inside its own pair.
   Measured to be the ONLY thing between the producer and an assembled `biconj`,
   and measured NOT to need a max-fold (so every boundary stays a conic). One open
   question: disambiguating a line-pair conic's sign. TODO.md 2026-09-17 item 1.
2. The DEGENERATE cell pairs — vertex/edge/affine cells are 28 of 29 bridge pairs
   on real fixtures. The value side already works; the region needs the polar
   construction. TODO.md 2026-09-17 item 2.
3. Then the fold into one `AlgCon` and `biconjQ`'s `convexPieceCoupling` route.

## Files
- `bridgeFaceAlg.m` — the producer: FP (a cubic's discriminant in the multiplier),
  FB (the TOP branch), and the region as two conic sign conditions
- `algBranch.m`, `AlgCon.m` — FB < 0 counts from the top; positive is unchanged
- `envelopeCellRegions.m` — each candidate now carries `.S`, its maximiser map
- `envelopeCellPairs.m` — the pair walk and classification, unchanged

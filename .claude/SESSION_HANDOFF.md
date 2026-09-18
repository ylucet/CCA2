# Session Handoff

_2026-09-17_

## Blocked
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — delete / track /
  note? Your call. MORNING.md.
- `maxQuaPar:notDegenerate` reachability — an in-contract pair needs the symbolic
  pipeline (hours). DECISIONS 2026-09-06.

## State
- Branch `main`, ten commits ahead of `origin/main`, tree clean
- Tests: fast bucket green; bridgeFaceAlgTest 22/0 (~50 s, the heaviest in the
  bucket and the one that pins the envelope) · normal 13/0 (09-11) · slow 98/0 ·
  verylong 37/0 (09-06/07, not re-run since)
- Known reds: none

## Next
1. **The SINGULAR PENCIL faces**, 12 of the 29 bridge pairs on the standard
   fixtures: affine x rank-1, where the tie is a PARABOLA and the objective is
   LINEAR, so the face is a RATIONAL function with a linear denominator — the
   `Rat` level, closed form, no elimination. TODO.md / DECISIONS.md 2026-09-17.
2. The fold into one `AlgCon` and `biconjQ`'s `convexPieceCoupling` route. The
   assembly must SKIP tangential pairs (10 of 29 — they carry no 2-D face).
3. The BOUNDARY (a planar arrangement), for figures and the paper's tables.

## Files
- `bridgeFaceAlg.m` — the producer: FP, FB, the outer bound, the pair data
- `bridgeMember.m` — membership from the geometry; `bridgeWitness.m` is its exact
  counterpart, `bridgePhi.m` the one polynomial both rest on
- `bridgeArcCut.m` — exact cut ROWS for a line tie; a conic tie is left to the
  predicate, with both refuted constructions written out
- `ratQ.pairLine*`, `QuaCon.FC`'s third column — an edge may name one line of a pair

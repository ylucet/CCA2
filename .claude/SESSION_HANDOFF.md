# Session Handoff

_2026-09-17_

## Blocked
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — delete / track /
  note? Your call. MORNING.md.
- `maxQuaPar:notDegenerate` reachability — an in-contract pair needs the symbolic
  pipeline (hours). DECISIONS 2026-09-06.

## State
- Branch `main`, eight commits ahead of `origin/main`, tree clean
- Tests: fast bucket green; bridgeFaceAlgTest 20/0 (48 s — the heaviest suite in the
  bucket, and it is the one that pins the envelope) · normal 13/0 (09-11) ·
  slow 98/0 · verylong 37/0 (09-06/07, not re-run since)
- Known reds: none

## Next
1. The producer for `rank1` and `affine` cells — 28 of 29 bridge pairs on the
   standard fixtures still refuse. The VALUE side already works; the region needs
   only an outer bound and a witness now that membership is a predicate.
2. The fold into one `AlgCon`, and `biconjQ`'s `convexPieceCoupling` route.
3. The BOUNDARY (a real planar arrangement) — wanted for figures and the paper's
   topology tables, and deliberately OFF the path to a correct `biconj`.

## Files
- `bridgeFaceAlg.m` — the producer: FP, FB, the region's outer bound, the pair data
- `bridgeMember.m` — membership, decided from the geometry; `bridgeWitness.m` is its
  exact counterpart and `bridgePhi.m` the one polynomial both rest on
- `bridgeArcCut.m` — the exact cut ROWS for a line tie; a conic tie is left to the
  predicate, with both refuted constructions written out
- `ratQ.pairLine*`, `QuaCon.FC`'s third column — an edge may name one line of a pair

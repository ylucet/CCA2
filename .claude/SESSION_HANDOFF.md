# Session Handoff

_2026-09-17_

## Blocked
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — delete / track /
  note? Your call. MORNING.md.
- `maxQuaPar:notDegenerate` reachability — an in-contract pair needs the symbolic
  pipeline (hours). DECISIONS 2026-09-06.

## State
- Branch `main`, eleven commits ahead of `origin/main`, tree clean
- Tests: fast bucket 661/0 over 44 suites · normal 13/0 (09-11) · slow 98/0 ·
  verylong 37/0 (09-06/07, not re-run since)
- Known reds: none

## Next — item 3, and it is now fully specified
1. Fold the per-pair faces into one object and route `biconjQ`'s
   `convexPieceCoupling` through it. The rule per pair: SKIP a tangential one;
   a STRAIGHT tie is a FACET of each side (no bridge face); a CURVED tie takes
   `bridgeFaceAlg`'s face, whose membership is `bridgeMember`. Fold the quadratic
   faces with `maxQ`; a bridge face needs no fold. TODO.md 2026-09-17.
2. The BOUNDARY (a planar arrangement) for figures and the paper's tables — off
   the critical path, because a predicate decides membership.

## Files
- `bridgeFaceAlg.m` — the producer; delegates to `bridgeFaceRat.m` when the pencil
  is singular, so one entry point returns either face kind (`face.kind`)
- `bridgeMember.m` — membership from the geometry; `bridgeWitness.m` exact,
  `bridgePhi.m` the one polynomial both rest on
- `bridgeArcCut.m` — exact cut ROWS for a line tie; both refuted constructions for
  a conic tie are written out in its header
- `ratQ.pairLine*`, `QuaCon.FC`'s third column — an edge may name one line of a pair

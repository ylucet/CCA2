# Session Handoff

_2026-09-18_

## Blocked
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — delete / track /
  note? Your call. MORNING.md.
- `maxQuaPar:notDegenerate` reachability — an in-contract pair needs the symbolic
  pipeline (hours). DECISIONS 2026-09-06.

## State
- Branch `main`, twelve commits ahead of `origin/main`, tree clean
- Tests: fast bucket green over 44 suites · normal gains `envelopeAssembleTest`
  4/0 in ~130 s · slow 98/0 · verylong 37/0 (09-06/07, not re-run since)
- Known reds: none

## What just landed
The convex envelope is ASSEMBLED and verified end to end: exact (0.0e0) against f
on the convex fixture, 8.6e-13 and 3.2e-12 against a numeric sup on the two
non-convex ones, with the bridge faces answering. `envelopeAssemble` returns the
FACES and `envelopeEval` takes the max.

## Next
1. The CORNER candidates — a curved-tie pair's candidates are incomplete where the
   tie meets a facet (measured: 320 and 280 of 600 interior points covered). The
   face is affine with quadratic-irrational coefficients and `FP`/`FB` can hold it.
2. The DOMINATION PRUNE — the pairwise `maxQ` fold explodes (5 → 610 faces over
   seven pairs); the prune is an exact sign test on a polygon.
3. Then a container, and `biconjQ`'s `convexPieceCoupling` route.

## Files
- `envelopeAssemble.m`, `envelopeEval.m` — the assembly and its evaluator
- `bridgeFaceAlg.m` / `bridgeFaceRat.m` — the bridge face, algebraic or rational
- `bridgeMember.m`, `bridgeWitness.m`, `bridgePhi.m` — membership, exact and cached
- `bridgeArcCut.m` — exact cut rows for a line tie; both refuted constructions kept

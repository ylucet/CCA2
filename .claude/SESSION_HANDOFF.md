# Session Handoff

_2026-09-18_

## Blocked
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — delete / track /
  note? Your call. MORNING.md.
- `maxQuaPar:notDegenerate` reachability — an in-contract pair needs the symbolic
  pipeline (hours). DECISIONS 2026-09-06.

## State
- Branch `main`, thirteen commits ahead of `origin/main`, tree clean
- Tests: fast bucket green over 44 suites · `envelopeAssembleTest` 5/0 (~155 s,
  normal) · slow 98/0 · verylong 37/0 (09-06/07, not re-run since)
- Known reds: none

## What just landed
The convex envelope is assembled, and now COMPLETE: the corner candidates (where
the tie meets a facet) close every pair's coverage, so the domain is right as well
as the values — 0.0e0 against f on the convex fixture, 8.6e-13 and 3.2e-12 against
a numeric sup on the kinked ones.

## Next
1. The DOMINATION PRUNE — the pairwise `maxQ` fold explodes (5 → 610 faces over
   seven pairs) against an envelope with a handful of faces. An exact quadratic
   sign test on a polygon decides whether a candidate is beaten on its whole
   region. This is what stands between the face list and a compact mesh.
2. A container for a predicate-decided face — or the observation, once the prune
   lands, that a compacted mesh no longer needs one.
3. Then `biconjQ`'s `convexPieceCoupling` route.

## Files
- `envelopeAssemble.m`, `envelopeEval.m` — the assembly and its evaluator
- `cornerFaceAlg.m` — the corner candidates that complete a curved pair
- `bridgeFaceAlg.m` / `bridgeFaceRat.m` — the bridge face, algebraic or rational
- `bridgeMember.m`, `bridgeWitness.m`, `bridgePhi.m` — membership, exact and cached

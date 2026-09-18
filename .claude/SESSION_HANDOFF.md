# Session Handoff

_2026-09-18_

## Blocked
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — delete / track /
  note? Your call. MORNING.md.
- `maxQuaPar:notDegenerate` reachability — an in-contract pair needs the symbolic
  pipeline (hours). DECISIONS 2026-09-06.

## State
- Branch `main`, fourteen commits ahead of `origin/main`, tree clean
- Tests: fast bucket green over 44 suites · `envelopeAssembleTest` 6/0 (~207 s,
  normal) · slow 98/0 · verylong 37/0 (09-06/07, not re-run since)
- Known reds: none

## Where the envelope stands
It is COMPUTED and VERIFIED -- 0.0e0 against f on the convex fixture, 8.6e-13 and
3.2e-12 against a numeric sup on the kinked ones, domain included -- and now
PRUNED: 120 -> 52, 142 -> 94, 110 -> 38 faces with nothing changed.

## Next — ONE item, where there were two
**THE ARRANGEMENT.** Overlay the pruned faces' regions and take the max per cell.
That is both the mesh `biconj` would return and the boundary the paper's figures
want. `maxQ` cannot do it: a QuaCon is +infinity off its faces while a candidate
that does not apply is -infinity, so folding collapses the list to an empty mesh
(measured). `envelopeEval` already answers correctly from the face list, so
nothing is blocked except the return TYPE.

## Files
- `envelopeAssemble.m`, `envelopePrune.m`, `envelopeEval.m` — assemble, prune, evaluate
- `cornerFaceAlg.m` — the corner candidates that complete a curved pair
- `bridgeFaceAlg.m` / `bridgeFaceRat.m` — the bridge face, algebraic or rational
- `bridgeMember.m`, `bridgeWitness.m`, `bridgePhi.m` — membership, exact and cached

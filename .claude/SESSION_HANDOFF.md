# Session Handoff

_2026-09-18_

## Blocked
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — delete / track /
  note? Your call. MORNING.md.
- `maxQuaPar:notDegenerate` reachability — an in-contract pair needs the symbolic
  pipeline (hours). DECISIONS 2026-09-06.

## State
- Branch `main`, fifteen commits ahead of `origin/main`, tree clean
- Tests: fast bucket 662/0 over 44 suites · `envelopeAssembleTest` 6/0 (normal) ·
  slow 98/0 · verylong 37/0 (09-06/07, not re-run since)
- Known reds: none

## Where the envelope stands
COMPUTED, VERIFIED and PRUNED. Values exact against f on the convex fixture and
8.6e-13 / 3.2e-12 against a numeric sup on the kinked ones, domain included;
120 -> 52, 142 -> 94, 110 -> 38 faces after the prune with nothing changed.
What it is not is a single mesh object.

## Next — TWO routes; both fold orders are now ruled out by measurement
1. **BUILD THE ARRANGEMENT ONCE**: overlay the pruned regions, winner per cell. No
   intermediate function is formed, so nothing blows up on the way — and it is the
   same construction the paper's figures and topology tables need.
2. **Return the FACE LIST** with a type whose eval is the max over faces — a TYPE
   DECISION (a new concrete class, not a QuaCon), and the cheap answer if the
   arrangement proves long. `envelopeEval` already does this correctly.

Folding is ruled out either way: left-to-right reaches 3726 faces by its twelfth
pair, the balanced TREE reaches 570 by level 3 and dies at level 4. The cause is
not the order — a max over any proper subset of the pairs is a partial max, which
carries pieces only the full max deletes. Do not try a third order.

## Files
- `envelopeAssemble.m`, `envelopePrune.m`, `envelopeEval.m` — assemble, prune, evaluate
- `cornerFaceAlg.m` — the corner candidates that complete a curved pair
- `bridgeFaceAlg.m` / `bridgeFaceRat.m` — the bridge face, algebraic or rational
- `bridgeMember.m`, `bridgeWitness.m`, `bridgePhi.m` — membership, exact and cached

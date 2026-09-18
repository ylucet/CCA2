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

## Next — three routes, in order
1. **BUILD THE ARRANGEMENT ONCE**: overlay the pruned regions, winner per cell. No
   intermediate function is formed, so nothing blows up on the way — and it is the
   same construction the paper's figures and topology tables need.
2. Fold in a TREE rather than left to right: same worst case, often far better in
   practice, and a cheap experiment to run first.
3. Return the FACE LIST with a type whose eval is the max over faces — a type
   decision, since it needs a new concrete class rather than a QuaCon.

Folding left to right is ruled out by measurement: semantically right per PAIR
(not per face), and it doubles per fold — 5 to 3726 over twelve pairs.

## Files
- `envelopeAssemble.m`, `envelopePrune.m`, `envelopeEval.m` — assemble, prune, evaluate
- `cornerFaceAlg.m` — the corner candidates that complete a curved pair
- `bridgeFaceAlg.m` / `bridgeFaceRat.m` — the bridge face, algebraic or rational
- `bridgeMember.m`, `bridgeWitness.m`, `bridgePhi.m` — membership, exact and cached

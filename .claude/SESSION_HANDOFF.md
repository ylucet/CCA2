# Session Handoff

_2026-09-18_

## Blocked
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — delete / track /
  note? Your call. MORNING.md.
- `maxQuaPar:notDegenerate` reachability — an in-contract pair needs the symbolic
  pipeline (hours). DECISIONS 2026-09-06.

## State
- Branch `main`, nineteen commits ahead of `origin/main`, tree clean
- Tests: fast 666/0 over 45 suites (09-18) · `envelopeAssembleTest` 8/0 (~394 s, normal) ·
  slow 98/0 · verylong 37/0 (09-06/07, not re-run since)
- Known reds: none

## Where the envelope stands
Assembled, verified, pruned, and now MESHED. `envelopeMesh` builds the arrangement
once — 52 faces -> 13 curves -> 58 cells -> a QuaCon of 321 faces reproducing f
EXACTLY on the convex fixture. `biconjQ` attempts the route for any two-piece
`convexPieceCoupling` input.

## Next
1. **The CORNER ROWS are done and verified** — `cornerWedgeRows`, 6400/6400 points
   agreeing with the KKT predicate, wired into `cornerFaceAlg.con`. What is left is
   the ARRANGEMENT's side: a branch row names ONE line of a pair, so a curve's
   identity becomes (conic, branch) and the emptiness tests (`ratQ.feasible2`,
   `cellHasInterior`) cannot read a branch row as the conic's sign.
   **That is a fork worth one decision** — exact feasibility over Q(sigma) (Fourier-
   Motzkin with `signAt` for the sign tests; clean, but a cell can mix rows from
   several different quadratic fields), or a RELAXED emptiness test that keeps some
   empty cells and lets the winner refinement delete them.
2. **The BRIDGE face's region** stays out (its cut is refuted twice, DECISIONS
   2026-09-17). Either solve that, or carry a bridge as a predicate face in a type
   that allows one — the type decision that was route B. Both halves are needed:
   the arrangement gives the mesh and the boundary, a container gives `biconj` a
   return value on a curved tie.
3. Compaction: 321 faces where a handful is right. The refinement makes many cells
   carrying the same function and `mergeAdjacentCells` only joins single-curve
   neighbours.

## Files
- `envelopeAssemble.m` / `envelopePrune.m` / `envelopeMesh.m` / `envelopeEval.m`
- `cornerFaceAlg.m`, `bridgeFaceAlg.m`, `bridgeFaceRat.m` — the curved-tie faces
- `bridgeMember.m`, `bridgeWitness.m`, `bridgePhi.m` — membership, exact and cached

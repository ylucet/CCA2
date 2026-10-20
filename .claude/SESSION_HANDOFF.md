# Session Handoff

_2026-10-20_

## Blocked
- Push: GitHub rejected the stored credentials ("Invalid username or token") on 2026-10-20; re-auth
  (`gh auth login`) and push main (8+ commits ahead).
- Erratum text changes: Yves to review erratum-lean.tex / erratum-diff.pdf
  (erratum repo D:\AI\67191168115b2b8bc24afaab).
- Example 2's mesh: its triangle faces have degree-4 irrational corners; needs a mesh type (Yves).

## State
- CCA2 `main`: biconj defaults to 'exact' (2026-10-18). biconjQ returns QuaCon / AlgCrv everywhere
  the suites reach except erratum Example 2 (EnvelopeMin, eval only).
- Done this stretch (DECISIONS.md 2026-10-16..20): active-pair prune; one-pair and several-pair
  meshing (min of maxes); the poly2 kink meshes (~3 min); corners carry value; cellInterior's
  curved+branch gap handled in liveCell; certified algebraic-row emptiness (algEmpty/polySignOn).
- Buckets: fast 766/0; slow 106/0 (2026-10-17) plus the new slow suites kinkMeshTest,
  erratumBiconjEx2Test, minMeshTest, envelopeMeshSlowTest; normal 30/0 in 535 s (within budget).
- Lean: unchanged since 2026-10-13 (267 theorems, axiom-clean).

## Next
1. Push (after re-auth).
2. Remaining empty faces (TODO item 3): pieces whose cell the curve R = 0 crosses.
3. Optional: planar Minkowski-Weyl in Lean (lean/README "Not formalized").

## Files
- `DECISIONS.md` 2026-10-16..20 entries -- every change above, with measurements.
- `envelopeMeshAlg.m` -- the mesher (header + minMaxFaces, liveCell, algEmpty).
- `polySignOn.m`, `cellSamplePoint.m`, `envelopeRationalCorners.m` -- new helpers.

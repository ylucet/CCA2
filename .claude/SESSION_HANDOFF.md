# Session Handoff

_Last updated: 2026-08-15_

## Where things stand

- On **`main`**, pushed. (The `overnight/2026-08-13` branch was merged and is history.)
- **`maxQuaParTest` 29 / 0. Fast bucket 204 / 0.** `MAXQP_ASSERT = 1` is now ON for every test in
  `maxQuaParTest`, via a `TestMethodSetup` that restores the previous value on teardown.
- **`maxQuaPar` has NO open case.** The seeded arc-vs-arc sweep (seed 20260803, N=18) is
  **18 exact / 0 wrong / 0 errored**, from 16 / 0 / 2 two sessions ago.
- Slow bucket: the four documented reds
  (`biconjugateOverATwoFaceSubdivisionIsTheEnvelope`, `step3UnboundedAssemblyMatchesTheTruth`,
  two in `unboundedFaceTest`). `testMaxMultiRegion` 24 / 0, `testcPLQ` 8 / 0, `testRegion` 23 / 0
  and `biconjCPLQTest` 10 / 0 were re-run against the current tree and are clean.

## What happened this session

**Bug 5 — FIXED.** `splitTwoArcPiece` found no cut when a piece's two arcs are ADJACENT: its two
candidate chords join the arcs' facing endpoints, which for arcs sharing a vertex ARE the arcs' own
edges, so both chains came out too short and the piece was returned unsplit with one arc flattened
to its chord. Generalised the `nv == 3` shared-vertex fallback to `nv >= 4` with the ordinary
diagonal to a non-adjacent vertex. Pinned by
`arcVsArcSplitsTwoADJACENTArcsOnAPieceWithADiagonal`.

**Bug 1 — three defects fixed, one attempted fix rejected as unsound.** `f**` of `x·y` over the
two-face square is now exact at 5 of 7 probe points; it was 0 of 7.

1. `getEdgeNosInf` numbers an edge by one of its endpoint VERTICES, so a LENS — two edges joining
   the same pair — gets one number for both, and the last-write-wins scatter destroys one.
2. `getNormalConeVertexQ` indexed its second constraint as `j+1` unwrapped, raising
   `badsubscript` on any BOUNDED region — which is why its only caller sent every bounded region
   to the POLYHEDRAL routine, whose cones come from the chord and are wrong for a curved edge.
3. `biconj` handed its second conjugation the curved MESH `conj` has returned since 2026-08-13;
   `quaPolToPlq` refuses a curved domain. It now asks for the symbolic form on purpose.

**Bug 2 — localised, not fixed.** Step 2 is right; the assembled Step 3 cell has lost its `−s1 ≤ 0`
and the quadratics that replaced it are blind to the sign of `s1`.

## Next steps

1. **Bug 1's remainder is a REFACTOR, not a fix.** Give `conjugateOfPiecePoly` an explicit EDGE
   LIST instead of a count with two conventions (`endNv = nv` or `nv−1`; edge `j` at `ineqs(j)` or
   `ineqs(j+1)`). It cannot be done there alone — `j` indexes `getNormalConeEdgeQ`/`Q3`'s output
   and `getSubdiffVertexT2`/`T2Q`'s `subdE` at the same time, so all four move together.
   **Do NOT free a slot by dropping the constraint holding it** — tried, unsound, see
   `DECISIONS.md`.
2. **Bug 2:** dump the cells carrying `s1²/4 + s2²/2` immediately before `functionNDomain.mergeL`.
   Either `region.merge` is unioning two cells whose union is not convex (and `unionIsExact`
   should refuse), or a mirror cell gets that quadratic wrongly earlier. `redundantSubset` is
   already exonerated.
3. **Bugs 3–4 are a missing ALGORITHM, not a defect**, and are worth taking off the bug list:
   `convEnvUnbounded` computes only the AFFINE envelope over an unbounded face and refuses the
   rest by design. Both fixtures have envelopes that are convex and not affine.
4. **Then SCIP/QPLIB**, in the order that bites: wire `biconj` into `SCIP/src/cca2ConvexEnvelope.m`
   → expose value+subgradient off `QuaParCPLQ` → fix diagonal terms over a box → performance.

## How to work on the symbolic layer, learned the hard way

**Build a unit-level reproducer before touching anything.** The half-lens
`{(s1+s2)² ≤ 4·s2, s2 ≤ s1}` carrying `s1`, constructed by hand as a `region` and conjugated
against a brute-force sup over its own boundary, runs in SECONDS and needs no pipeline. Pipeline
runs of the same defect take 10–40 minutes. That one change is what made bug 1 tractable after two
sessions of failed attempts.

**Sign-blind quadratic constraints are a recurring failure mode.** A region bounded by a conic that
cannot tell its two branches apart claims territory on the wrong side. It caused the `maxQuaPar`
defect fixed on 2026-08-14 and it is what bug 2 looks like now.

## Relevant files

- `DECISIONS.md` — dead ends and refuted reasoning. Read before attacking bug 1 or 2; several
  natural-looking approaches are recorded there as unsound.
- `TODO.md` — live action items, with the measurements for each.
- `SUPPORT_MATRIX.md` §4.1–4.3 (arc-vs-arc, the verification tools, the covering proof), §7, §8.
- `maxQuaPar.m` — `splitTwoArcPiece`, `newRaySign`, `pieceIsCurved` (read its header before
  touching `curveAfter` anywhere).
- `functionNDomain.m` — `conjugateOfPiecePoly` and `spreadCollidingEdges`; `region.m` —
  `getNormalConeVertexQ`.

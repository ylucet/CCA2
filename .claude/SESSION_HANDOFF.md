# Session Handoff

_Last updated: 2026-08-15_

## Where things stand

- On **`main`**, pushed. (The `overnight/2026-08-13` branch was merged and is history.)
- **`maxQuaParTest` 29 / 0. `conjCPLQTest` 25 / 0. Fast bucket 204 / 0.** `MAXQP_ASSERT = 1` is now ON for every test in
  `maxQuaParTest`, via a `TestMethodSetup` that restores the previous value on teardown.
- **`maxQuaPar` has NO open case.** The seeded arc-vs-arc sweep (seed 20260803, N=18) is
  **18 exact / 0 wrong / 0 errored**, from 16 / 0 / 2 two sessions ago.
- Slow bucket: **three** reds, down from four —
  `biconjugateOverATwoFaceSubdivisionIsTheEnvelope` and the two in `unboundedFaceTest`.
  `step3UnboundedAssemblyMatchesTheTruth` is GREEN. `testMaxMultiRegion` 24 / 0, `testcPLQ` 8 / 0,
  `testRegion` 23 / 0 and `biconjCPLQTest` 10 / 0 re-run clean.
- **`SUPPORT_MATRIX.md` §8.2's blocker is closed**: Steps 1, 2 and 3 are all done for unbounded
  multi-face domains. What remains there is new item (f), Step 1's CURVED convex envelope.

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

**Bug 2 — FIXED.** `region.removeTangent` built the TANGENT LINE to a quadratic constraint at a
vertex where that quadratic's GRADIENT VANISHES — the apex of a cone, which is exactly where an
unbounded fan's Step 3 split conics meet — and deleted a constraint matching that meaningless
"tangent". It removed `−s1 ≤ 0`, leaving two constraints blind to the sign of `s1`, so the region
went symmetric and claimed the mirror wedge. **This is the SECOND routine to fall into that trap
on that same input**; §8.2(e) records `simplifyUnboundedRegion` doing it, fixed 2026-08-02 by
`region.witnessAwayFrom`.

## Next steps

1. **Bug 1's remainder is a REFACTOR, not a fix.** Give `conjugateOfPiecePoly` an explicit EDGE
   LIST instead of a count with two conventions (`endNv = nv` or `nv−1`; edge `j` at `ineqs(j)` or
   `ineqs(j+1)`). It cannot be done there alone — `j` indexes `getNormalConeEdgeQ`/`Q3`'s output
   and `getSubdiffVertexT2`/`T2Q`'s `subdE` at the same time, so all four move together.
   **Do NOT free a slot by dropping the constraint holding it** — tried, unsound, see
   `DECISIONS.md`.
2. **Bugs 3–4 are a missing ALGORITHM, not a defect**, and are worth taking off the bug list:
   `convEnvUnbounded` computes only the AFFINE envelope over an unbounded face and refuses the
   rest by design. **Both envelopes are already derived in `TODO.md`** —
   `co(x·y + I_K) = y²` on `K = {0 ≤ y ≤ x}` with its proof, and `co(−x²+y²) = −x + y²` on the
   half-strip — together with the pattern they share. Prove the general rule the way
   `convEnvUnbounded`'s header proves the affine case; do not ship a formula that merely matches
   the two fixtures.
4. **Then SCIP/QPLIB**, in the order that bites: wire `biconj` into `SCIP/src/cca2ConvexEnvelope.m`
   → expose value+subgradient off `QuaParCPLQ` → fix diagonal terms over a box → performance.

## How to work on the symbolic layer, learned the hard way

**Build a unit-level reproducer before touching anything.** The half-lens
`{(s1+s2)² ≤ 4·s2, s2 ≤ s1}` carrying `s1`, constructed by hand as a `region` and conjugated
against a brute-force sup over its own boundary, runs in SECONDS and needs no pipeline. Pipeline
runs of the same defect take 10–40 minutes. That one change is what made bug 1 tractable after two
sessions of failed attempts.

**A degenerate geometric object is not a geometric object, and this codebase keeps assuming it is.**
Three separate defects this month, in three different routines, were all one of these: a quadratic
has no TANGENT LINE at its own apex, where the gradient vanishes (`removeTangent`, and
`simplifyUnboundedRegion` before it); two arcs sharing a vertex have no separating CHORD, because
the candidate chords are the arcs' own edges (`splitTwoArcPiece`); a conic cannot tell its two
BRANCHES apart, so a region bounded by one claims territory on the wrong side (`maxQuaPar`, and
bug 2's symptom). When a routine computes a geometric object from a degenerate configuration, check
that the object exists before using it.

## Relevant files

- `DECISIONS.md` — dead ends and refuted reasoning. Read before attacking bug 1 or 2; several
  natural-looking approaches are recorded there as unsound.
- `TODO.md` — live action items, with the measurements for each.
- `SUPPORT_MATRIX.md` §4.1–4.3 (arc-vs-arc, the verification tools, the covering proof), §7, §8.
- `maxQuaPar.m` — `splitTwoArcPiece`, `newRaySign`, `pieceIsCurved` (read its header before
  touching `curveAfter` anywhere).
- `functionNDomain.m` — `conjugateOfPiecePoly` and `spreadCollidingEdges`; `region.m` —
  `getNormalConeVertexQ`.

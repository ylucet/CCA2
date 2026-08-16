# Session Handoff

_Last updated: 2026-08-15_

## What happened this session

**Bug 1 is fixed, and the repository now has no failing test.** It was the last red, marked BLOCKED
after four attempts. The whole thing comes down to one sentence: *both slot conventions identify an
edge by a VERTEX INDEX, and a lens's two edges join the SAME pair of vertices* — so neither
convention can name them apart and no reassignment of slots to constraints ever could. The four
earlier attempts were all variations on reassigning slots. `functionNDomain.edgeIndexList` derives
the edge list from the geometry instead, and the three edge-indexed readers take it as an argument.

## Where things stand

- Branch: `main`, pushed, tree clean.
- **327 pass / 0 fail across all 26 suites** — fast 204 / 0, normal 8 / 0, slow **115 / 0**. The
  repository has no failing test for the first time in weeks.
- `checkBiconjDomainCoverage` re-measures the two-face box — the row that read WRONG — as **OK,
  error 0**, against a ground truth that owes nothing to the conjugate pipeline.
- Measured, unit level: both half-lenses conjugate to **3 cells** (2 vertex cones plus the arc; the
  chord's cell is a ray and drops out), exact against a brute-force sup over the lens at 12 of 12
  points. The old code was `+inf` at `(0,0)` and `(-1,0.5)` where the truth is `0`, and `0` at
  `(2,-1)` where it is `1/2`.
- Two documented cases still ERROR. Both **refuse loudly rather than answering wrongly**, and both
  are unimplemented paths, not defects: the general convex quadrilateral (FIRST conjugation, the
  `nCE == 3` gap) and the parallelogram (SECOND conjugation, `emptyResult`).

## Next steps

1. **cPLQ's Step 2 on a 2-convex-edge triangle is WRONG — a live silent defect, found 2026-08-15
   and the reason the quadrilateral's `nCE == 3` wiring was written and then reverted.** The wiring
   works at Step 1 (4 envelope faces, all `≤ x·y`, no more crash) but the answer is then wrong, and
   a silent wrong answer is worse than a loud refusal. Of the two triangles the test quadrilateral
   splits into, the new `nE = 3` one gets its envelope and **Step 2 returns zero conjugate cells**,
   while the `nE = 2` one — untouched — **carries the whole error**: `f*(0,0) = 0.28647` for a truth
   of `0`, `f*(0.5,1) = 1.00464` for `1`, and a hole in the third quadrant.
   The measurement that localises it, and the trap it avoids: that same triangle conjugated ON ITS
   OWN via `QuaPol.conj` is exact at 7 of 7 — because a single bounded triangle takes the NUMERIC
   route (`conjBoundedPolygon`), not cPLQ. **Checking a suspect piece "on its own" through the
   public API can silently change which implementation runs**; evaluate `p.pieces(k).maxConjugate`
   inside the pipeline instead. Order of attack and the ready-to-re-land branch are in `TODO.md`.
2. **The parallelogram's `emptyResult` — LOCATED, and the error message names the wrong routine.**
   All 12 pieces conjugate (27 cells); it is `functionNDomain.maxOfList` that returns nothing, and
   folding the groups one at a time shows **group 11** emptying the accumulator. Start there, and
   fix `QuaParCPLQ.conj`'s message while you are in it.
3. **Then SCIP/QPLIB**, in the order that bites: wire `biconj` into `SCIP/src/cca2ConvexEnvelope.m`
   → expose value+subgradient off `QuaParCPLQ` → fix diagonal terms over a box (`x²−y²`,
   `(x²+y²)/2` on `[0,1]²` still error in the second conjugation) → performance (~40–60 s/term).

## Two methods that keep paying off

**Build a unit-level reproducer before touching the symbolic layer.** The half-lens conjugates in
~10 s against a brute-force sup; the same piece reached through `biconjugateTest` takes 10–40
minutes. Bug 1's whole fix was developed on a saved copy of that one piece.

**When a routine's index convention is unclear, read its PROBE POINTS.** `getNormalConeVertexQ`'s
two halves probe at vertex `j-1` and vertex `j+1`, which is what settles that they mean the edge
ARRIVING at `j` and the edge LEAVING it. Reading that off the code rather than the comments caught
an off-by-one in the new edge list that no current test would have shown.

## Relevant files

- `DECISIONS.md` — dead ends and refuted reasoning; the newest entry is bug 1's, with the reason
  none of the four earlier attempts could have worked.
- `TODO.md` — live action items; the quadrilateral wiring is first.
- `SUPPORT_MATRIX.md` §7 — the defect table; bug 1's row is struck through with the measurement.
- `functionNDomain.m` — `edgeIndexList`, `edgesStillCollide`, `conjugateOfPiecePoly`, `getInterior`.
- `region.m` — `getNormalConeVertexQ` (optional edge list), `getNormalConeEdgeQE`, `coneNormalAt`.
- `functionNDomainTest.m` — the two unit tests that pin the lens without the pipeline.

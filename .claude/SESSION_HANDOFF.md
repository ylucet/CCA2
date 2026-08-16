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
- **Fast bucket 204 / 0. Normal bucket 8 / 0. Slow bucket green** (`biconjugateTest` 7 / 0 —
  including the test that had been red — and `biconjCPLQTest` 10 / 0; see `MORNING.md` for the
  full run).
- Measured, unit level: both half-lenses conjugate to **3 cells** (2 vertex cones plus the arc; the
  chord's cell is a ray and drops out), exact against a brute-force sup over the lens at 12 of 12
  points. The old code was `+inf` at `(0,0)` and `(-1,0.5)` where the truth is `0`, and `0` at
  `(2,-1)` where it is `1/2`.
- Two documented cases still ERROR. Both **refuse loudly rather than answering wrongly**, and both
  are unimplemented paths, not defects: the general convex quadrilateral (FIRST conjugation, the
  `nCE == 3` gap) and the parallelogram (SECOND conjugation, `emptyResult`).

## Next steps

1. **Wire `nCE == 3` into Case C's Step 1** — the quadrilateral. CCA2 already has the algorithm
   (`convEnvCPLQ`'s `splitThreeConvex`, [COAP] A.5); the vendored `plq_1p.convexEnvelope1` simply
   falls off the end after `nCE == 2`. `TODO.md`'s first item has the route worked out — including
   why it is cheap (`isCanonicalXY` guarantees the piece is exactly `x*y` by then; `plq_1p.conjugate`
   already loops over multiple envelope pieces; `ratPolToPlq.m` shows the construction to copy) and
   the one question it raises about `nCE == 2` and A.4.
2. **The parallelogram's `emptyResult`.** Undiagnosed. Use the per-piece dump that cracked bug 1 —
   conjugate each piece of `f*` on its own and find the one returning no cells — not `biconj`.
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

# Session Handoff

_Last updated: 2026-08-16_

## What happened this session

**Both remaining ERROR cases are fixed** — the parallelogram's `QuaParCPLQ:conj:emptyResult` and the
general quadrilateral's `MATLAB:badsubscript`. The quadrilateral took four attempts, and the one
that worked differed from the one that hung in exactly one respect: **it does the geometry
SYMBOLICALLY.** [COAP] A.4's cevian foot is irrational, so computing it in double precision and
converting gives `2^53` denominators that grow past `1e25` downstream; carried symbolically the
coordinates stay compact surds and the pipeline finishes.

**What it costs is now the open question**, and it is a performance one: surd coordinates make every
symbolic operation downstream work in a quadratic extension, and `testcPLQ` goes from 1542 s to over
3100 s. `CCA2_NO_A45_SPLIT` opts out; fixing Step 3's cell blow-up is what would remove the need.

## Where things stand

- Branch: `main`, pushed, tree clean, no debug instrumentation.
- **330 pass / 0 fail across all 26 suites** — fast 206 / 0, normal 9 / 0, slow **115 / 0**.
- The parallelogram's `biconj` now computes: exact at all four vertices, `+inf` outside the domain,
  8 of 10 interior probe points right against a brute-force double conjugate. The other two are
  about 4% LOW — the residual defect below.
- The general quadrilateral still raises `MATLAB:badsubscript`, loudly, as before.

## Next steps

1. **`getInterior` on a SINGULAR quadratic — the parallelogram's last 4%.** `getInterior` separates
   an edge cell from its neighbours by eliminating `s` between `x = ∂₁f` and `y = ∂₂f`. When `f` is
   a singular convex quadratic the gradient map is not invertible and that elimination returns the
   map's IMAGE LINE, which separates nothing — so two edge cells come out on the SAME region and
   the first one checked wins. Reproduce in about a minute with
   `functionNDomainTest.parallelogramPiece9`.
   **Do NOT attack the `isQuad` chord rewrite for it.** Both alternatives were measured 2026-08-16
   and are in `DECISIONS.md`: chording the vertices the conic actually touches makes that piece
   WORSE (2 wrong of 10 → 3), and skipping the rewrite changes nothing.
2. **Step 3's cross-piece maximum does not scale, and it is now the binding cost.** Measured on the
   quadrilateral: Steps 1 and 2 take about 25 s for all six pieces, `functionNDomain.maxOfList` then
   takes **73 minutes**, and the cell count runs 5, 14, 29, 45, 70, **86** — roughly ten times what
   the answer needs, since `f*` of `x·y` over a convex quadrilateral has a cone per vertex and a
   cell per edge. The surplus is adjacent cells carrying the SAME function that `region.merge` never
   merges. **Start by merging same-function neighbours after each fold rather than only at the end**;
   many of these cells are POLYHEDRAL, where `unionIsExact` already decides exactly. This is what
   stands between the A.4/A.5 split and being affordable without `CCA2_NO_A45_SPLIT`.
3. **Then SCIP/QPLIB**, in the order that bites: wire `biconj` into `SCIP/src/cca2ConvexEnvelope.m`
   → expose value+subgradient off `QuaParCPLQ` → fix diagonal terms over a box (`x²−y²`,
   `(x²+y²)/2` on `[0,1]²` still error in the second conjugation) → performance (~40–60 s/term).

## Three methods that keep paying off

**Build a unit-level reproducer before touching the symbolic layer.** Piece 9 of the
parallelogram's `f*`, hand-built as a `region`, runs in about a minute; reaching the same piece
through `biconj` takes far longer and buries the evidence.

**When a routine's index convention is unclear, read its PROBE POINTS.** `getNormalConeVertexQ`'s
two halves probe at vertex `j−1` and vertex `j+1`, which is what settles that they mean the edge
ARRIVING at `j` and the edge LEAVING it.

**To find which piece of a MAX is wrong, use the property the max must have.** `f**` of a bounded
domain is finite exactly ON that domain and is a max, so EVERY per-piece conjugate must be finite
there. Evaluating all 12 groups at six interior points turned "the max comes out empty" into "these
three pieces are wrong" in one cheap run.

## Relevant files

- `DECISIONS.md` — dead ends and refuted reasoning; the newest entries are the parallelogram's two
  defects and the three failed quadrilateral attempts.
- `TODO.md` — live action items; the quadrilateral leads, with its one-line remaining task.
- `SUPPORT_MATRIX.md` §7 / §7.1 — the defect table and the domain-coverage measurements.
- `region.m` — `simplifyUnboundedRegion` (the witness), `getNormalConeVertexQ`, `getNormalConeEdgeQE`.
- `functionNDomain.m` — `edgeIndexList`, `conjugateOfPiecePoly`, `getInterior`.
- `functionNDomainTest.m` / `regionTest.m` — the unit tests that pin all of this without the pipeline.

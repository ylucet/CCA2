# Session Handoff

_Last updated: 2026-08-17_

## What happened this session

**Every documented ERROR case now has a fix.** Bug 1 (the lens) closed, taking the repository to no
failing test for the first time in weeks; the parallelogram's `QuaParCPLQ:conj:emptyResult` closed,
along with two general defects behind it; and the general convex quadrilateral's
`MATLAB:badsubscript` closed on the fourth attempt. The quadrilateral was done test-first — red,
then green — and the attempt that worked differed from the one that hung in exactly one respect:
**it does the geometry SYMBOLICALLY.** [COAP] A.4's cevian foot is irrational, so computing it in
double precision and converting gives `2^53` denominators that grow past `1e25` downstream; carried
symbolically the coordinates stay compact surds and the pipeline finishes.

**The quadrilateral fix is OPT-IN (`CCA2_A45_SPLIT`) and OFF by default** — the one judgement call
left open. See "Next steps" 1.

## Where things stand

- Branch: `main` @ `ba3457d` — "Make the A.4/A.5 split OPT-IN: with it on, testcPLQ/testRectBiconj
  ERRORS"
- Pushed: yes
- **332 pass / 0 fail across all 26 suites** — fast 206 / 0, normal 11 / 0, slow **115 / 0**, every
  suite at its historical runtime.
- The parallelogram's `biconj` computes: exact at all four vertices, `+inf` outside the domain,
  8 of 10 interior probe points right against a brute-force double conjugate. The other two are
  about 4% LOW — see "Next steps" 2.
- The general quadrilateral raises `MATLAB:badsubscript` by DEFAULT and is exact with
  `CCA2_A45_SPLIT` set (10 of 10 probe points, and 8 of 8 through the full assembly).

## Next steps

1. **AWAITING A DECISION: should the A.4/A.5 split become the default?** It is off because turning
   it on costs `testcPLQ` 1542 s → 4728 s **and** makes `testcPLQ/testRectBiconj` ERROR (a test with
   no assertions, so that is an exception, undiagnosed). Three options were put to the user:
   (a) leave it opt-in and fix Step 3's cost first — recommended, since that is also what the 4%
   error in 2 needs; (b) diagnose the `testRectBiconj` exception first, then flip the default;
   (c) flip it now and accept both costs. Full reasoning in `DECISIONS.md`'s newest entry.
2. **`getInterior` on a SINGULAR quadratic — the parallelogram's last 4%.** It separates an edge
   cell from its neighbours by eliminating `s` between `x = ∂₁f` and `y = ∂₂f`; when `f` is a
   singular convex quadratic the gradient map is not invertible and that elimination returns the
   map's IMAGE LINE, which separates nothing — so two edge cells come out on the SAME region and the
   first one checked wins. Reproduce in about a minute with
   `functionNDomainTest.parallelogramPiece9`.
   **Do NOT attack the `isQuad` chord rewrite for it.** Both alternatives were measured 2026-08-16
   and are in `DECISIONS.md`: chording the vertices the conic actually touches makes that piece
   WORSE (2 wrong of 10 → 3), and skipping the rewrite changes nothing.
3. **Step 3's cross-piece maximum does not scale, and it is the binding cost.** Measured on the
   quadrilateral: Steps 1 and 2 take about 25 s for all six pieces, `functionNDomain.maxOfList` then
   takes **73 minutes**, and the cell count runs 5, 14, 29, 45, 70, **86** — roughly ten times what
   the answer needs, since `f*` of `x·y` over a convex quadrilateral has a cone per vertex and a
   cell per edge. The surplus is adjacent cells carrying the SAME function that `region.merge` never
   merges. **Start by merging same-function neighbours after each fold rather than only at the end**;
   many of those cells are POLYHEDRAL, where `unionIsExact` already decides exactly.
4. **Then SCIP/QPLIB**, in the order that bites: wire `biconj` into `SCIP/src/cca2ConvexEnvelope.m`
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

- `DECISIONS.md` — dead ends and refuted reasoning; newest entries are the opt-in judgement call,
  the quadrilateral's three failed attempts, and the parallelogram's two defects.
- `TODO.md` — live action items; Step 3's cost leads.
- `splitTightTriangleSym.m` — the A.4/A.5 split, with its derivations and why it is symbolic.
- `plq_1p.m` — `appendTriangle` (the `CCA2_A45_SPLIT` gate and what it costs), `triangulate`.
- `SUPPORT_MATRIX.md` §7 / §7.1 — the defect table and the domain-coverage measurements.
- `region.m` — `simplifyUnboundedRegion` (the witness), `getNormalConeVertexQ`, `getNormalConeEdgeQE`.
- `functionNDomain.m` — `edgeIndexList`, `conjugateOfPiecePoly`, `getInterior`.
- `cplqAdapterTest.m` / `functionNDomainTest.m` / `regionTest.m` — the unit tests that pin all of
  this without the pipeline.

# Session Handoff

_Last updated: 2026-08-17_

## What happened this session

**Option (a) was chosen and worked on: the A.4/A.5 split stays opt-in, Step 3's cost is the job.**
The standing rule the user attached to it outranks the flag and settles the whole class of
question: **every computation has to be CORRECT even if it is slow -- a slow correct path gets its
test moved to a slower bucket, never traded away -- with the one exception of a computation so slow
it does not finish, since a timeout helps nobody.**

**THREE double leaks found and fixed; Step 2 is now exact end to end.** `domain.mE`/`cE` were
DOUBLE arrays (an exact slope arrived as `0.6`, an exact zero y-intercept as `-9.06e-72`);
`region.limitOfFAtVertices` let its `limf(j) = 0` branch decide the array's class, silently
rounding every exact gradient limit written after it; and `conjConvexOverPiece` converted `Q, L, c`
and the piece's vertices to double by design. Worst denominator in the quadrilateral's conjugate
cells: **1.4e145 -> 56**.

**And the blow-up is now COUNTED rather than guessed -- which twice refuted what this session
believed.** Exactness turned out not to be what blocks merging, and "nothing merges" was the wrong
framing: of 38 same-function pairs, 21 meet only at a POINT and must not merge, so most refusals
were always correct. What is left is small and splits evenly in two, both measured. See
"Next steps" 1 and 2, and `DECISIONS.md`'s newest three entries for the two retractions.

**332 pass / 0 fail** -- fast 206, normal 11, slow 115 -- with everything above in.

## Where things stand

- Branch: `main` @ `4e23b0d` — the Step 3 measurement work above
- Pushed: yes
- **332 pass / 0 fail across all 26 suites** — fast 206 / 0, normal 11 / 0, slow **115 / 0**, every
  suite at its historical runtime.
- The parallelogram's `biconj` computes: exact at all four vertices, `+inf` outside the domain,
  8 of 10 interior probe points right against a brute-force double conjugate. The other two are
  about 4% LOW — see "Next steps" 2.
- The general quadrilateral raises `MATLAB:badsubscript` by DEFAULT and is exact with
  `CCA2_A45_SPLIT` set (10 of 10 probe points, and 8 of 8 through the full assembly).

## Next steps

1. **DECIDED 2026-08-17 — option (a): the split stays opt-in and STEP 3's COST is the work.** The
   standing rule the user attached to it outranks the flag and settles future versions of the same
   question: **every computation has to be correct even if it is slow; a slow correct path gets its
   test moved to a slower bucket, never traded away — the one exception being a computation so slow
   it does not finish, since a timeout helps nobody.** So `testcPLQ` at 4728 s is a bucket question;
   what blocks the default is the undiagnosed `testcPLQ/testRectBiconj` exception. `DECISIONS.md`'s
   newest entry has it in full.
   **Started, and BLOCKED on the VPN:** instrumentation for the blow-up is written
   (`region.mergeTally`, the refusal reason out of `unionIsExact`, `CCA2_TRACE_MAXP` in
   `functionNDomain.maximumP`) along with the harness `.claude/step3cost.m` — but MATLAB cannot
   check out a licence off the UBC network (`License Manager Error -96`), so **none of it has been
   run, not even for syntax.** `TODO.md`'s newest item has the command to run first and the three
   ranked hypotheses that one run decides between.
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

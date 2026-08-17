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
- Pushed: NO -- six commits sit on local `main`, unpushed (pushing was never asked for)
- **332 pass / 0 fail across all 26 suites** — fast 206 / 0, normal 11 / 0, slow **115 / 0**, every
  suite at its historical runtime.
- The parallelogram's `biconj` computes: exact at all four vertices, `+inf` outside the domain,
  8 of 10 interior probe points right against a brute-force double conjugate. The other two are
  about 4% LOW — see "Next steps" 2.
- The general quadrilateral raises `MATLAB:badsubscript` by DEFAULT and is exact with
  `CCA2_A45_SPLIT` set (10 of 10 probe points, and 8 of 8 through the full assembly).

## Next steps

1. **The 5 pairs whose shared facet `merge` cannot see.** Measured with EXACT arithmetic, so this
   is not rounding: they carry the same hyperplane with opposite orientation and meet in a
   SEGMENT, and `ineqs(i) == -ineqs(j)` still finds nothing. `symbolicFunction.eq` is
   `if (obj1.f == obj2.f)`, a STRUCTURAL test whose own comment says "change to isAlways", so the
   same constraint at a different positive SCALE does not match. `region.normalize1` divides by
   `abs(coeffs(f,vars))(end)` and is supposed to prevent exactly that -- check first whether it
   picks the same term for both operands. Reproduce in about 4 minutes with
   `.claude/step3adjacency.m`.
2. **The 6 pairs `certifiesNonPositive` declines.** They ARE found, reach `unionIsExact`, and are
   refused there -- the fold-1 tally's `lin_exactCurvedTest = 6` matches the count exactly. The
   certificate refuses by design outside its hypothesis: a rational `h`, a non-convex quadratic,
   or a linear relaxation with no vertex. Instrument WHICH of the three fires before extending it;
   the derivation in the method's header says what each would need.
3. **Then re-measure and settle the A.4/A.5 default.** Rerun `.claude/step3cost.m` on the
   quadrilateral once 1 and 2 land (cells ran 5, 14, 29, 45, 70, 86 before any of this work), then
   diagnose `testcPLQ/testRectBiconj`'s exception -- with the cost question answered that is the
   only correctness objection left to turning `CCA2_A45_SPLIT` on.
4. **`getInterior` on a SINGULAR quadratic — the parallelogram's last 4%.** Unchanged by this
   session. It separates an edge cell from its neighbours by eliminating `s` between `x = d1f` and
   `y = d2f`; for a singular convex quadratic the gradient map is not invertible and that
   elimination returns the map's IMAGE LINE, which separates nothing. Reproduce in about a minute
   with `functionNDomainTest.parallelogramPiece9`.
   **Do NOT attack the `isQuad` chord rewrite for it** -- both alternatives were measured
   2026-08-16 and are in `DECISIONS.md`.
5. **Then SCIP/QPLIB**, in the order that bites: wire `biconj` into `SCIP/src/cca2ConvexEnvelope.m`
   -> expose value+subgradient off `QuaParCPLQ` -> fix diagonal terms over a box (`x^2-y^2`,
   `(x^2+y^2)/2` on `[0,1]^2` still error in the second conjugation) -> performance.

## New tools this session, both cheap to rerun

- `.claude/step3cost.m` -- folds the pieces one at a time and reports cells, DISTINCT FUNCTIONS,
  and `region.mergeTally`'s refusal reasons per fold. `CCA2_STEP3_CASE=tri` switches to the
  all-rational control; `CCA2_STEP3_FOLDS` bounds the work.
- `.claude/step3adjacency.m` -- classifies every same-function pair three ways at once: does
  `merge` see a facet, is there a shared hyperplane, do the cells actually meet in a SEGMENT.
  The third column is what showed that 21 of 38 refusals were correct all along.
- `CCA2_TRACE_BIGNUM` in `region`'s constructor prints the stack whenever a region is built from a
  constraint carrying a 15-digit integer. All three double leaks were found with it.

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

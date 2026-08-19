# Session Handoff

_Last updated: 2026-08-18_

## What happened this session

**Two operators got their short-circuits, and `biconj` over a box is now free.** `biconj` IS the
convex envelope, so a convex `f` is its own answer — that plus separability, McCormick and the
diamond case took every row of `checkBoxEnvelopeForSCIP` to **0 s**, from 42–456 s with one
`MATLAB:badsubscript`. All six now return a meshed `QuaPol` the SCIP bridge can read, where they
used to return a mesh-less `QuaParCPLQ`.

**Step 3 got measurably faster by profiling first.** `ptFeasible` gained a numeric filter and
`getVertices` an affine×conic closed form: `solve` −42%, `subs` −17%, one fold 289 s → 265 s.

**The symbolic-removal programme then hit a wall worth knowing about**, and it is not a
performance wall — see "Read this first".

## READ THIS FIRST — what the record says about my own errors

`DECISIONS.md` is the valuable artifact from today. Newest entries, all measured:

1. **`getNormalConeVertexQ` does NOT compute a normal cone.** Tested against the definition
   (`u` in the cone at `v` iff `v` maximises `⟨u,·⟩`), the COMMITTED implementation disagrees on
   4–30 of 72 directions per vertex. A gradient rewrite fails too, in two opposite ways: too big on
   the concave side of a conic, too small at a cusp. **So there is no statement of what a
   replacement must satisfy.** Settle it by reading the CONSUMER, `getSubdiffVertexT1`.
2. **Every `solve()`-for-a-probe-point site is order-dependent.** 42 differential cases, 16
   disagree, every disagreement has both roots real. `solve`'s ordering is unreproducible, so
   substitution changes which probe is used. This is a **caller rewrite**, not a substitution.
3. **A closed-form rewrite was 30% SLOWER** (2.89 → 3.77 s). Extraction cost amortises over uses:
   convert where one extraction serves many operations (`getVertices`), not where it serves one.
4. **A safety check on a hot path must be priced like anything else on it** — a sampled refit
   guard cost more than the `solve` calls it saved.

## Where things stand

- Branch: `main` @ `d3d7454`
- Pushed: **pending** — nothing has been pushed since `51c003d`
- **fast 217 / 0** (89 s, 17 suites), **normal 11 / 0**, `regionTest` 18 / 0
- **SLOW BUCKET RUN CLEAN 2026-08-19 and the verification gap is CLOSED** (`.claude/slowrun.log`,
  ~2 h): 119 / 1 over seven suites, the single red being a STALE type expectation in
  `conjCPLQTest.biconjCoverageByInputCase` that was already red at this session's starting commit
  (verified two ways: `git log -S convexEnough`, and running that test in a worktree at `801ee1f`).
  Fixed; `conjCPLQTest` re-ran alone at **25 / 0**. Everything the `ptFeasible` filter, the
  `getVertices` closed form, Step 0 and the probe rewrite touch is therefore green.

## Next steps

1. **DONE 2026-08-19. The slow bucket is green** — see "Where things stand". Nothing is waiting on
   it any more. Re-run it after the next `region.m` change; a suite that has not been run since a
   change is how the one stale expectation above went unnoticed for a day.
2. **DONE 2026-08-18.** `getNormalConeVertexQ`'s specification is established and is now a test
   (`regionTest.vertexConesMatchTheDefinition`, green, 17/0 in 44 s): the rows' linear parts, read
   as `≤ 0`, are the NORMAL CONE at the vertex; the `s2` coefficient must be `±1` or `0`; the
   constant is free. Given `eIdx` the committed routine is EXACT on all three curved fixtures,
   cusp and three-active-constraint vertex included. Yesterday's "it is not a normal cone" is
   struck in `DECISIONS.md` — it measured the eIdx-less slot fallback on bounded fixtures.
   The fallback is also SOUND on the unbounded layout it was written for (rebuilt from
   `getEdgeNosInf`'s scatter: 0 of 72 wrong; second test). regionTest 18 / 0 in 45 s.
   **Left open, narrow:** a BOUNDED region for which `edgeIndexList` refuses would still reach the
   fallback with the pair off by one; none has been seen.
3. **DONE 2026-08-18. `SUPPORT_MATRIX.md` §0.0.1 is re-derived** from `checkBoxEnvelopeForSCIP`
   on the current tree (log `.claude/boxenvelope.log`): six rows, no ERROR, every one 0 s and
   returning a MESHED `QuaPol`. Three of its four recorded gaps are closed — the diagonal terms
   work (separable / convex short-circuits), the "no mesh" headline is false for every box case,
   and 40–60 s per term is no longer a blocker on that path. `SCIP_READINESS.md`'s A1, A2 and the
   first gate condition are marked met. The caveat that decides whether a QPLIB run is worth
   doing is unchanged: on box+bilinear CCA2 reimplements McCormick, now in 0 s instead of 40.

**Also done 2026-08-18: the first symbolic-removal site.** `getNormalConeVertexQ`'s eight
`solve()` probe calls became four `region.probeOnConstraint` calls — first FEASIBLE root, not
first root — with closed-form quadratic roots (`region.rootsIn`) and `solve()` kept as a fallback.
Live `solve()` in `region.m` 16 → 10. With the edge list the cones are still exactly right
(0 of 72 at every vertex); the eIdx-less fallback moved closer to the definition (32→29, 43→29,
5→5). Three pinned characterization values flipped orientation and are re-pinned with the reason.
One defect went with it: the `cNext` block's second attempt probed `ineqs(cj)`.
4. **DONE 2026-08-18. Phase B is answered, both items.** The map is measured
   (`.claude/phaseBmap.m`, twelve shapes): every single PIECE is closed form and costs ≤ 0.02 s;
   the only symbolic rows are the general polygon (2579 s, Step 3) and a convex MULTI-FACE input
   that pays 26–28 s because `convexEnough` needs the caller's `fIsConvex` flag. B2: the entry
   point is `q.biconj('cplq')` itself — on every shape QPLIB's box terms present it returns a
   meshed `QuaPol` in 0.01 s and never reaches Step 3, so the "40–60 s per term" blocker is stale
   and Phase C must be re-scoped around the NON-box case. The optimisation notes that were filed
   under Phase B are now filed under Phase C, which is what made that gate read like a
   symbolic-removal programme.
5. **DONE 2026-08-18. Step 0 is built** (`mergeSameQuadFaces.m`, called first by both operators).
   The two-triangle unit square is now the one-face square: `biconj` 0.1 s and CORRECT (the
   known-failing `biconjugateOverATwoFaceSubdivisionIsTheEnvelope` is green and untagged), `conj`
   0.8 s. A merge needs a CONVEX union — a reflex face builds and then evaluates to +inf — and an
   edge left separating a face from itself goes with it. `CCA2_NO_STEP0` turns it off, for the
   lens regression only.

## Relevant files

- `DECISIONS.md` — dead ends; the four newest entries are today's and are the reason not to retry
  the normal-cone rewrite blind
- `ALGORITHM.md` — the ORDER of operations for both operators, and why triangulation is the
  indefinite case's apparatus only
- `SCIP_READINESS.md` — the three-phase gate, plus the profile of one Step 3 fold
- `biconjCPLQ.m` — the four short-circuits (convex, separable, McCormick, diamond)
- `region.m` — `ptFeasible`'s numeric filter, `getVertices`' closed forms, `getNormalConeVertexQ`
  (unresolved)
- `regionTest.m` — `normalConesOnCurvedEdgesAreUnchanged`, 13 s, guards the curved cones
- `.claude/step3cost.m`, `.claude/step3adjacency.m` — the Step 3 measurement tools

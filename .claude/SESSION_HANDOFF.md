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

- Branch: `main` @ `35c550b` — "getNormalConeVertexQ does not compute a normal cone"
- Pushed: **pending** — 7 commits unpushed (`05df79d`..`35c550b`); everything up to `51c003d`
  was pushed earlier today
- **fast 206 / 0**, **normal 11 / 0**, `regionTest` 16 / 0 including the new cone test
- **SLOW BUCKET NOT RUN** since the `ptFeasible` filter and the `getVertices` closed form. This is
  the one verification gap. User is running it tonight:
  `CCA2_TEST_TIMEOUT=7200 bash .claude/suite.sh --slow`
- `region.m` is at its last green state — every experimental rewrite this evening was reverted

## Next steps

1. **Read tonight's slow-bucket result.** It covers `region.m`'s two accepted changes. Nothing
   else should be built on them until it is green.
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
3. **Rewrite `SUPPORT_MATRIX.md` §0.0.1.** Dated 2026-08-02 and now wrong in every row: the six
   box cases are all 0 s, and its headline finding ("returns `QuaParCPLQ`, NO MESH") no longer
   holds. It is what a reader consults before planning SCIP work.
4. **Finish Phase B** of `SCIP_READINESS.md` — the direct-formula/symbolic map and the entry point
   a separator should call. Most evidence is now measured.
5. **Step 0 of `ALGORITHM.md`'s `biconj` ordering** — merge adjacent faces carrying the same
   quadratic; the one case in the plan still unbuilt.

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

# Session Handoff

_Last updated: 2026-07-30_

## What happened this session

**Step 3's assembly is fixed and verified.** Both defects reduced to one primitive — maximize a
linear form over a polyhedron, i.e. an LP — and on the reference case `f = xy` over
`conv{(0,0),(3,3),(1,2)}` the assembled partition went from **57 of 289 grid points wrong to 0**,
exact at every one of the fold's seven rounds, and the whole fold got **faster** (1645 s vs
1768 s) because a correct partition carries fewer regions than a damaged one.

That fix made three latent bugs reachable in code that had been shielded by the old
delete-happy behaviour. **All three are now fixed** — the last one, `testcPLQ/testRectBiconj`,
this session.

**`testcPLQ/testRectBiconj` is green.** The cause was NOT the extra constraints as such but a
SILENT COLLISION in `conjugateOfPiecePoly`'s edge-number scatter. `getEdgeNosInf` returns 0 for a
constraint with no vertex (already handled), but for one with EXACTLY ONE vertex it returns that
vertex's own index — the slot the real edge leaving that vertex already claims. The scatter is
last-write-wins, so the intruder EVICTED a genuine edge and left the evicted constraint's old
slot holding a stale duplicate. Piece 23 arrived as the triangle
`{(9s2)/5-s1+5, -s1-7s2-4, s1+2s2-4}` plus a quadratic touching only `(139/44,-45/44)`;
`edgeNo` came out `[3 1 1 2]` and the scatter returned
`[quad, s1+2s2-4, (9s2)/5-s1+5, s1+2s2-4]` — `-s1-7s2-4` gone, `s1+2s2-4` twice. The `isQuad`
branch then chorded `d0.vx(1)`→`d0.vx(2)`, two vertices the quadratic does not join, `solve`
returned a complex pair, and `gtd`'s bare `if (obj1.f>obj2)` could not take it.

The adapter drops a colliding constraint that bounds no edge, from `conjugateOfPiecePoly`'s LOCAL
copy only. Verified beyond "it stops erroring": the dropped quadratic is
`-(s1+7s2)^2-148s1+196s2+684`, decreasing in `u=s1+7s2` over the triangle's `u`-range, hence
maximal on the edge `u=-4` at `s2=-45/44`, where it is exactly 0 — **provably redundant**. A
brute-force max over the ORIGINAL 4-constraint domain matches the symbolic conjugate at 8 sample
points to grid resolution, one piece per point.

**`maxQuaPar` can now split a cell that already carries an arc** (was next-step 2; ~26% of
sampled splits). It needed no conic-conic solver and no multi-arc representation — both of which
this file and `maxQuaPar.m`'s own TODOs previously claimed. Every curved edge is a parabola
(`QuaPar.assertParabolic`), so restricting the splitting conic to the arc via the new
`parabolaArcFrame.conicCoeffs` gives one univariate QUARTIC in the frame's global monotone `u`.
And the splitting curve never CROSSES the arc in this pipeline: measured over the named fixture
plus a 395-quadrilateral sweep, all 22 curved-cell splits left the arc untouched (19) or tangent
(3) — the C1 tangency structure the file already documents — so the arc survives whole in one
half, and ONE ARC PER FACE is kept by subdividing that half with a straight chord
(`splitTwoArcPiece`). Assembled results went **58 → 76 of 395**, every one of the 77 (sweep +
fixture) exact against the closed-form sup (worst **2.8e-14**) and violation-free under
`arrangementViolations`.

## Where things stand

- Branch: **`main`** — `step3-assembly-lp-certificates` was fast-forwarded into it (11 commits,
  from the LP-certificate assembly fix through the `maxQuaPar` arc work) and is no longer the
  place to work. This repo's practice is to commit on `main`; the branch existed only because
  that work was already in flight.
- Pushed: yes — `origin/main` is up to date.
- Working tree: clean
- **Suite: 274 pass / 1 fail** over 24 suites (full sweep, after the `maxQuaPar` arc work),
  against a 270/2 baseline at the start of the session. `testcPLQ` is 8/0, `maxQuaParTest` 16/0,
  the new `functionNDomainTest` 2/0. The single failure is `testRegion/testCreation`
  (longstanding, unrelated — a toolbox-compatibility issue, not the conjugate pipeline).
- **`main` is now unblocked**: the merge gate ("merge only once `testcPLQ/testRectBiconj` is
  green") is satisfied.

## Next steps

1. **Unbounded multi-face `conj`** — re-scoped 2026-07-30, this time MEASURED rather than read.
   Two independent defects, and the ORDER they are fixed in matters.

   **(a) `quaPolToPlq` throws the ray away.** `faceVertexIndices` takes one vertex per edge,
   `E(j,1)`, and never consults QuaPol's ray flag `E(:,3)==0`, for which `E(j,1)` is the base
   point and `E(j,2)` the DIRECTION point. The 4-cone fan `V=[0 0;-1 0;0 1;1 0;0 -1]`,
   `E=[1 2 0;1 3 0;1 4 0;1 5 0]` gives `V=[(0,0);(0,0)]` for the second-quadrant cone, and
   `domain()` — the BOUNDED constructor, which closes the vertex loop — turns that into
   `NaN <= 0` **twice** (not `x <= 0` twice, as recorded before).
   Route: one inequality per edge via `domain.domainEdge`, segments and rays alike being "the
   line through `E(j,1)` and `E(j,2)`", orientation from `F(j,:)`. **Done this session:** the
   loud rejection (`quaPolToPlq:unboundedFace`), so (a) can never be fixed without (b).

   **(b) `plq_1p` reads region's infinity markers as ORDINARY NUMBERS.** This is the real work,
   and it is NOT closed by (a). `domain.domainEdge` genuinely does build the ray encoding —
   the second quadrant comes out as 2 inequalities with vertices `(0,0)`, `(0,intmax)`,
   `(-intmax,0)`, `(-intmax,intmax)`, i.e. source vertex plus a direction vertex per ray — and
   `functionNDomain`/`region` (Step 2/3) handle that correctly, as `testRectBiconj`'s own
   unbounded dual pieces 25–27 show. `plq_1p` does not. Measured on that exact domain with
   `f=(x^2+y^2)/2`, whose conjugate is `min(s1,0)^2/2 + max(s2,0)^2/2` (4 dual pieces):
   - Case C's order (`triangulate` then `maximum`) **errors** — `plq_1p.conjugateFunction` →
     `region.getEdgeNos` → `symbolicFunction.getLinearCoeffs`, "Index exceeds the number of
     array elements".
   - Skipping `triangulate` is WORSE: it **runs and returns garbage** — 8 pieces carrying
     `2147483647*s_2` and the constant `4611686014132420609`, which is exactly `intmax^2`.
     Max error **1.15e18**.

   So the irreducible new work is in `plq_1p` — see (d), which narrows it further to Step 1's
   ENVELOPE. `conjCPLQ.m`'s `isDomBounded` gate comes out **last**, not first.

   **(c) Step 2 needs nothing new — confirmed numerically 2026-07-30.** Two facts, both of which
   (d) then sharpens:
   - `sup` over a union is the `max` of the `sup`s, so the conjugate over an unbounded face
     decomposes over ANY cover of it — and `plq.maximum` already takes a max across pieces. So
     `triangulate` only has to emit triangles PLUS wedges; there is no "conjugate over a general
     unbounded polygon" to write.
   - The conjugate of a convex quadratic over a wedge is the **same T1/T2 active-set
     decomposition cPLQ already builds for a triangle**, with one vertex instead of three and two
     unbounded edges instead of three bounded ones: apex piece + 2 ray pieces + interior piece.
     Confirmed against brute force on `f=(x^2+y^2)/2` over `{x<=0, y>=0}`, whose conjugate is
     `min(s1,0)^2/2 + max(s2,0)^2/2` — exactly 4 pieces, matching to grid resolution (1.7e-5 on a
     1200x1200 grid) at 9 sample points including all four cells and the origin.
     For an AFFINE `f = a.x + c` on a wedge `v + cone(d1,d2)` it collapses further, to one piece:
     `f*(s) = <s-a,v> - c` on `{s : <s-a,d1> <= 0, <s-a,d2> <= 0}`, `+inf` elsewhere — an affine
     function on a translated normal cone, which is the shape `conjugateFunction` already emits.
   **(d) The general conjugate machinery ALREADY does unbounded, for the case the pipeline
   produces. START at `plq_1p`, not at `functionNDomain`.** Measured 2026-07-30; this supersedes
   an earlier note in this file that said the gap was `conjugateOfPiecePoly`'s
   `if obj(i).d.nv > 1` gate on the T2 edge loop. **That was wrong**, and the way it was wrong is
   worth keeping:
   - Run `conjugateOfPiecePoly` on a wedge with an **AFFINE** `f` — the shape Step 1's `nCE==0`
     envelope always produces — and it is **exact**. `f = 2x-3y+1` on `{x<=0, y>=0}`: 1 piece
     (`-1`), err **0.000e+00** at all 8 probes, correctly `+inf` off the apex normal cone. A
     half-strip `{0<=x<=1, y>=0}` with `f = x-2y`: 2 pieces (`0`, `x-1`), err **0.000e+00** at all
     8. Unboundedness is simply not the problem.
   - The missing edge/interior pieces are a **pre-existing, bounded-domain** limitation, nothing
     to do with rays: the same routine on a BOUNDED triangle with the strictly convex
     `f=(x^2+y^2)/2` also returns only the 3 vertex pieces (`2*x-2`, `2*y-2`, `0`) and no edge or
     interior piece, so every interior dual point is uncovered.
   - And that limitation is by DESIGN, not an oversight: for an affine `f`, `grad f` is constant,
     there are no edge/interior pieces to make, and the vertex normal cones tile exactly the
     domain of `f*`. cPLQ always reduces to an affine (or rank-1-PSD, via the `nCE==1/2` special
     formulas) envelope first, so it never needs them.
   So the blocker is entirely `plq_1p`'s front end, and it is **Step 1, not Step 2**:
   `convexEnvelope1`'s `nCE==0` branch builds the affine envelope as the interpolant through
   `vx(1),vx(2),vx(3)` — a formula with no meaning for an unbounded face — and `triangulate`
   picks its fan apex by comparing coordinates (`max(vy)`/`min(vx)`) and rebuilds each triangle
   through the BOUNDED `domain(t,x,y)`. **The open question to answer first is mathematical: what
   IS the convex envelope of a quadratic over an unbounded polygon?** (It can be `-inf`; e.g.
   `x*y` over a wedge containing a direction where it decreases without bound.) Until that is
   settled, do not start coding.
   Decomposition subtlety for when `triangulate` is written: the pieces must be SUBSETS of the
   face whose union is the face (a superset inflates the sup), and fanning an unbounded convex
   polygon from one vertex yields half-strips (`conv{v1,vn} + cone(d)`), not only triangles and
   cones. Both shapes are already fine for Step 2 — see the half-strip result above.

   **Latent bug found while probing (not fixed, not on any path today):**
   `functionNDomain.getInterior` reads `solve()`'s result struct by the literal field names
   `s12.s_1`/`s12.s_2`, so it throws `Unrecognized field name "s_1"` for any region whose vars are
   not named `s_1,s_2`. In-pipeline they always are (it runs on `maxConjugate`), which is why
   nothing hits it. `fieldnames(s12)` would make it var-name agnostic.

2. A native numeric rational branch in `conjPieceCPLQ` would buy **speed, not coverage**.
3. Merge to `main` (now unblocked) and then 0.1 tagging — **do not tag without being asked**.
4. `partialConj` unimplemented for every engine and type (`SUPPORT_MATRIX.md` §2).
5. Longstanding: `'pqp'`/`'graph'` engines; `RatPol.conj`/`biconj`/`add`; the
   `mergeL`/`removeTangent` exact-tie-point bug; `QuaPar.eval` wrong exactly *at* a result
   vertex (~1.4%); `testRegion/testCreation`.

## Do NOT redo these — all tried and reverted, with reasons

- **`merge` guarded by constraint-set EQUALITY.** Provably sound, and makes things WORSE
  (36 → 125 wrong of 289): refusing merges leaves more regions for the other defect to damage.
  The exact condition (`unionIsExact`) is right because it refuses only what it must.
- **`slopes2` SKIPPING a curved constraint with no vertex.** Costs `maxArray` its probe
  directions → `maxEqDom` falls through to `splitmax3` → every undecided region splits →
  compounds round over round. `maximumP 3` went 153 s → >90 min. The finite-vertex CENTROID
  fallback now in the tree runs it in 192 s. Do not "simplify" it back.
- **Replacing `conjugateOfPiecePoly`'s scatter with a sort.** `d.ineqs(edgeNo) = d.ineqs` with
  `edgeNo(i) = j + add` (`j` a VERTEX index) can exceed `numel(ineqs)`, so the scatter
  legitimately GROWS the array, and repeated `edgeNo` values make it last-write-wins. Sorting
  reproduces neither, and it broke a third, previously-passing test
  (`testMaxMultiRegion/testFractional`).
- **PARKING vertexless constraints above the real edge slots.** Worse than useless — the
  `isQuad` branch then builds a chord for such a constraint out of `d0.vx(1),d0.vx(2)`, vertices
  unrelated to it. Dropping them from `conjugateOfPiecePoly`'s local copy is correct: that
  routine is edge-indexed end to end.
- **Concluding cPLQ handles unbounded PRIMAL domains because `region` does.** Half true, and the
  half that is false is the expensive one. `domain.domainEdge` + `region` really do represent a
  ray (source vertex + an `intmax` direction vertex), and Step 2/3 really do consume it. But
  `plq_1p` reads those markers as the NUMBER 2147483647 — see next-step 1(b) for the measurement,
  including a run that returns `intmax^2` as a coefficient instead of erroring. Do not relax
  `conjCPLQ`'s `isDomBounded` gate, or `quaPolToPlq`'s new `unboundedFace` rejection, until
  `plq_1p` has a wedge case. "It ran without erroring" is not evidence here.
- **Deciding "bounds no edge" by the VERTEX COUNT alone, gated on the region being BOUNDED.**
  Tried; it crashes `poly2orderUnbounded:312` on piece 24 of the same test. Two separate reasons,
  either fatal: (a) a RAY edge of an unbounded region is also active at exactly one finite vertex
  and is load-bearing — `getEdgeNosInf` keeps it precisely because slot 1 is reserved for it, so
  it never collides; (b) there is no cheap boundedness test to gate on — piece 24 is genuinely
  unbounded (recession direction `(2,-1)`) yet carries **no vertex at infinity** for `removeInfV`
  to find, so `nv`-before vs `nv`-after says "bounded". Key the drop on the COLLISION instead: it
  is the actual damage, and it leaves every piece whose `edgeNo` has no repeats untouched.
- **Offering ALL constraints to `deleteIfRedundant`** (not just the heuristic's candidates).
  The candidate set already contains every constraint missing a finite vertex; widening it only
  adds constraints that ARE active at vertices, and deleting a redundant-but-active constraint
  changes the vertex list `getVertices` computes, which the whole pipeline reads.

## Environment / harness notes

- MATLAB R2024b here is **network-licensed only** (`SERVER SLMS-SMATLABP1.ead.ubc.ca`, no local
  `license.dat`). Off the UBC VPN, `matlab -batch` dies with `Licensing error: -96,7 /
  System Error: 11002`. **Connect the VPN before any session that needs to run anything.**
- Harness in `.claude/`: `run.m` (round-by-round fold vs the per-piece max, reads `$CCA2DIR`),
  `suite.m` (per-suite pass/fail table), `smoke.m` (direct unit checks of the LP helpers).
  Make a pristine baseline with `git archive HEAD | tar -x -C baseline`.
- Timings: the harness prints per-round times — **trust those, not wall-clock between polls.**
  A contended machine (16 cores at ~43% from other work) made me wrongly conclude twice that a
  round had stalled when it had not.
- Long waiters must watch a **per-run** log filename. Three `until grep -q "^TOTAL" fix.log`
  waiters spun overnight because `fix.log` was later overwritten by a re-run, making their exit
  condition unsatisfiable.

## Relevant files

- `region.m` — LP-certificate block at the top (after `probeAlong`/`probePerp`); instance
  helpers `linearForm`/`redundantSubset`/`deleteIfRedundant`/`unionIsExact` just before
  `finiteVertices`; `slopes2`; `getEdgeNosInf`; `simplifyOpenRegion1`;
  `simplifyUnboundedRegion`; `merge` (header records the correctness argument and the 36 → 125
  history).
- `functionNDomain.m` — `conjugateOfPiecePoly`: the two empty-domain guards, and the entry
  ADAPTER just before the scatter (both parts — the vertexless drop and the collision drop — with
  the piece-23 trace in the comment). Its `if obj(i).d.nv > 1` gate on the T2 edge loop is where
  next-step 1(d) begins.
- `maxQuaPar.m` — `splitCell` (now handles a curved cell), `splitTwoArcPiece` (the one-arc-per-face
  chord subdivision), `arcHasStrictCrossing` (the tangency-vs-crossing check).
- `parabolaArcFrame.m` — `conicCoeffs`, the quartic restriction of a second conic to the parabola;
  `lineCoeffs`' companion.
- `regionTest.m` — 9 tests covering `maxLinear`'s three answers, `linearForm`'s affine flags,
  four redundancy cases, and three merge cases.
- `functionNDomainTest.m` — NEW, 2 tests. Reproduces `testcPLQ/testRectBiconj`'s piece 23
  directly, so the collision defect is pinned in ~26 s instead of that test's ~22 min. Verified
  to FAIL on pristine `HEAD~1` with the original error. It asserts the collision still exists as
  a precondition — if `getEdgeNosInf` ever stops producing it, the test stops covering anything.
- `conjCPLQ.m` — `assertStep3MatchesPieces` (the gate stays; it is a real invariant);
  `conjCPLQ.m:103`'s `isDomBounded` gate for next-step 1, which comes out LAST.
- `quaPolToPlq.m` — header now records both halves of next-step 1 and the measurements behind
  them; the new `quaPolToPlq:unboundedFace` rejection is at the top of the function body.
- `plq_1p.m` — `triangulate` (finite vertex fan, picks its start by `max(vy)`/`min(vx)` and
  rebuilds each triangle through the BOUNDED `domain(t,x,y)`), `convexEnvelope1` and
  `conjugateFunction` (closed-form triangle formulas indexing `vx(1..3)`). This is where
  next-step 1(b)'s wedge case has to go.
- `maxQuaPar.m` — `splitCell`'s `pieceIsCurved` guard and `insertPassthroughVertices` for
  next-step 2; also the overstated "conic-conic" TODOs at lines 145/150/194.
- `SUPPORT_MATRIX.md` §1.2 (4-face row now OK), §8 (blocker list).

## Still true from before

- **There is no "3-convex-edge case" to implement.** [COAP] Appendix A.5's split reduces such a
  triangle to 2-convex-edge sub-triangles and Step 1 already applies it. Describe these cases by
  the **envelope's face count**, not the input's edge count.
- **The rational-face `0/0` did not need exact arithmetic** — it was cleanly removable, and
  `symbolicFunction.limitDirectional` resolves it by a directional limit.
- **Which half of cPLQ is reusable:** its **Step 2/3**, on **CCA2's** Step 1 output.
- The reverted `maxArray` `intmax`-vs-`inf` "fix" — don't reintroduce it; the comment at the
  site explains why a vertical constraint must fall into the arithmetic mean.

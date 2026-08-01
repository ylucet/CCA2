# Session Handoff

_Last updated: 2026-08-01_

## What happened this session (2026-08-01)

Three named tasks, all done, plus one new blocker found and made loud.

**1. Curved convex envelope over an unbounded face — DONE.** `conjConvexOverPiece.m`. When a
face's `q` is convex, `co q = q`, so Step 1 has nothing to do and Step 2 must conjugate a CURVED
function — which cPLQ has no branch for (its Step 1 always hands Step 2 an affine or rank-1-PSD
envelope, so `conjugateOfPiecePoly` returns only the VERTEX cells and leaves the interior dual
uncovered). The construction is the KKT active set of `max{<s,x> - q(x) : x in P}`: one cell per
vertex (`s - grad q(v)` in the normal cone), one per edge with `d'Qd > 0`, one interior cell when
`Q` is nonsingular. It covers the bounded TRIANGLE and the unbounded WEDGE and HALF-STRIP with
the same code, because a ray contributes its direction exactly as a bounded edge does — which is
the "very specific geometries" restriction Yves asked for: those three shapes are all
`fanUnboundedFace` emits, and no general polyhedron case was written. The handoff's headline case
`(x²+y²)/2` over `{x<=0,y>=0}` returns `min(s1,0)²/2 + max(s2,0)²/2`, exact at 10 probes.

**2. The sidesteps are REPAIRED, not worked around — DONE.** The ray representation was never at
fault, as Yves said. Two ordinary bugs around it were:
- `region.getEdges` raised `MATLAB:unassignedOutputs` instead of returning an empty list when NO
  constraint is active at a point — and the box-clip corner `(intmax,intmax)` of the first
  quadrant is exactly such a point, since it lies on the implicit `±intmax` box rather than on
  either constraint. `edges` was only ever assigned inside the loop.
- `region.poly2orderUnbounded` then indexed `edges(1)`/`edges(2)` unconditionally and THREW on
  the simplest unbounded region there is.
`getNormalConeVertex` walks CONSECUTIVE vertices and so needs boundary cyclic order; with the
ordering routine repaired it returns the true cone at every finite vertex (first-quadrant apex:
`{s1<=0, s2<=0}`, previously `{s1+s2<=0, s1<=0}`). **`conjAffineOverPiece.m`, the parallel
half-plane construction added as a stopgap, is DELETED** — the native vertex path now gives
identical, correct answers, and every other vertex-based consumer benefits from the repair.

**3. Case C for a general quadratic — DONE, by rotating the domain.** Step 1 now classifies by
the SIGNS OF THE EIGENVALUES of `Q`, never by `nCE` (which counts edges of positive finite SLOPE
and therefore only classifies `x*y`):
- **convex / affine** -> `co q = q`, no envelope work; Step 2 uses `conjConvexOverPiece`.
- **concave** -> the envelope is affine; `convEnvUnbounded` builds it from the ACTUAL values
  of `q`, for a triangle, wedge or half-strip alike.
- **indefinite** -> `xyFrame.m` moves the problem into the frame where `q` IS `x*y`, cPLQ's own
  closed forms run there unchanged, and the conjugate comes back via `f*(s) = h*(M's - a) - c0`.
  The congruence is two steps: eigen-rescale to `diag(1,-1)`, then the 45-degree map to `u*v`.
Verified exact against brute force on the unit square for `(x²+y²)/2`, `x*y`, `x²-y²`,
`3xy+7x-2y+5` and `-(x²+y²)/2`. `f*(0.3,0.4)` was `0.4` against a truth of `0.125`; it is now
`0.125`. Pinned by `conjCPLQTest.caseCValuesAreCorrectForAGeneralQuadratic` — Case C previously
had NO value assertion at all, only type assertions, which is why this survived so long.

**THE NEW BLOCKER: Step 3's CROSS-PIECE maximum drops cells.** Newly reachable, because (1) makes
unbounded convex faces produce multi-cell conjugates. On the 4-cone fan with convex faces, each
of the 4 faces produces a CORRECT 4-cell conjugate and the per-piece maximum preserves all 4 —
then `plq.maximumConjugate` keeps only 4 of the 16, losing face 1's `s2²/2` cell on
`{s1<=0, s2>=0}`, so `f*(-0.5,2)` returns `1.125` for a truth of `2`. Localized by printing the
conjugates before and after each stage, so the attribution is measured, not inferred.
It is now LOUD rather than silent: `conjCPLQ`'s `assertStep3MatchesPieces` is applied to Case C
(it previously guarded only `conjEnvelopeViaCPLQ`), raising `PLQ:conjCPLQ:cplqFailed`.
**That gate itself had to be fixed first**, and the fix matters for anyone touching it: it built
its reference as the max of the FINITE per-piece values, skipping any piece that is `+inf` at the
sample point. That is right only when every piece's conjugate is finite everywhere — true when
every PRIMAL face is bounded, false the moment one is not. Unfixed it fired on a CORRECT answer
(`max(0,x,y)` as three wedges, whose conjugate really is the indicator of the simplex).

**A CHAIN of latent bugs in `conjugateOfPiecePoly`, unmasked one at a time.** Case C's
BICONJUGATE has never worked (ZERO pieces on pristine `HEAD`, i.e. `f** = +inf` everywhere, for a
convex `f`). It is not caused by this session's work — but this session's work makes the FIRST
conjugation richer (11 pieces rather than 9), which carries the second one further in and exposes
the next bug down each time one is fixed. All are the same shape: an index or an output used
outside the guard that was meant to protect it.
- `region.getEdges` — output only assigned inside the loop, so "no active constraint" raised
  `MATLAB:unassignedOutputs` instead of returning `[]`. **Fixed.**
- `region.splitmax3` — same shape: `r` unassigned when NO vertex has `f1 >= f2`, i.e. when `f2`
  is the max on the whole region and there is nothing to split. That case has a well-defined
  answer in the routine's own convention (`ineqs(1)` delimits where `f(1)` wins), so it returns
  `[-ineq, ineq]`: first half empty, second half everything. **Fixed.**
- `region.getNormalConeVertexQ` — `py = py(1)` placed BEFORE its own `if isempty(py)` guard,
  making the guard dead code. The identical block a few lines above has the right order.
  **Fixed.**
- `functionNDomain.getInterior` — indexes `c2(2)` under a guard that only tests `size(c1,2)`.
  **NOT fixed**; this is where the chain currently stops. Fixing it will likely expose the next.
`conjCPLQTest.biconjCoverageByInputCase` pins the limitation as "errors" rather than by
identifier, because the identifier moves as the chain is peeled; the invariant that matters is
that it does not silently return a wrong `f**`.

**Performance note.** `xyFrame` is built with EXACT symbolic arithmetic, not doubles. A double
becomes, under `sym`, a binary rational with a 2^52 denominator, and cPLQ's `solve`/`simplify`
then work with 17-digit numerators — measured to take `q = x²+3xy-2y²+x` from seconds to over 20
minutes. Built symbolically it stays a clean surd. A `Q` that is ALREADY `[[0,1],[1,0]]` skips
the eigendecomposition entirely (only the affine part needs stripping), which keeps cPLQ's own
inputs on rational arithmetic. Quadratics whose eigenvalue RATIO is not a rational square still
carry nested surds and are slow; that irrationality is inherent (reducing an indefinite form to
`u*v` needs a square root), not a defect.

---

## Previous session
## Previous session

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

**Merged to `main`, and next-step 1's `-inf` gate is implemented.** The branch was
fast-forwarded in (11 commits) once `testRectBiconj` went green. Then Yves answered the open
envelope question, which unblocked the first real code on next-step 1:
`region.quadUnboundedBelow(Q,L)` decides whether `conv q` over an unbounded face is `-inf` —
see 1(e) below for the closed form and the three things it is deliberate about.

## Where things stand

- Branch: **`main`**, building on `d82fe99` ("Add region.quadUnboundedBelow"). This repo's
  practice is to commit on `main`; `step3-assembly-lp-certificates` is fully merged and is no
  longer the place to work.
- **Not tagged.** 0.1 tagging is still deliberately not done — do not tag without being asked.
- **Suite: 292 pass / 1 fail over 25 suites** (full sweep, final code), against 274/1 over 24 at
  the start of this arc. The only failure is `testRegion/testCreation` (longstanding, toolbox
  compatibility, unrelated to the conjugate pipeline). `conjCPLQTest` 22/0, `unboundedFaceTest`
  12/0, `regionTest` 13/0, `biconjCPLQTest` 10/0, `testcPLQ` 8/0.
- **PERFORMANCE REGRESSED, and it is not noise.** `conjCPLQTest` alone now takes roughly an hour:
  `indefiniteTriangleThreeConvexEdgesUsesStep3` is **2002 s** and
  `caseCValuesAreCorrectForAGeneralQuadratic` **375 s**. Cause: conjConvexOverPiece emits 4-7
  cells per piece where the old path emitted 1-3, and plq.maximumConjugate overlays pieces
  pairwise, so the cell count going into each maximumP grows multiplicatively. This is the
  cost of getting the right answer, not a stall -- but it makes the suite painful and is worth
  attacking before anyone depends on it.
- Working tree: **not committed**. Modified: `region.m`, `plq_1p.m`, `quaPolToPlq.m`,
  `conjCPLQ.m`, `QuaParCPLQ.m`, `regionTest.m`, `conjCPLQTest.m`, `biconjCPLQTest.m`,
  `SUPPORT_MATRIX.md`, `.claude/SESSION_HANDOFF.md`. New: `fanUnboundedFace.m`,
  `convEnvUnbounded.m`, `conjAffineOverPiece.m`, `unboundedFaceTest.m`.
- **Three tests changed their expectation, and each change is a finding, not an accommodation.**
  Re-read these before assuming a test was merely "updated":
  * `conjCPLQTest/multiFaceUnboundedDomainStillNotImplemented` -> `...WithCurvedEnvelope...`: the
    blanket refusal became a specific one, which is what lets the affine case through.
  * `conjCPLQTest/biconjCoverageByInputCase`: was passing on an EMPTY result. Measured on
    pristine `HEAD` via `git archive`: `caseC.conj()` gives 9 pieces, `caseC.biconj()` gives
    ZERO -- `f** = +inf` everywhere for a convex `f`. `.kind()` is `'QuaParCPLQ'` either way,
    which is why nothing noticed. Now pinned as an error.
  * `biconjCPLQTest/unsupportedShapesStillErrorAsBefore`: the refusal moved off boundedness,
    which was never the real precondition, onto the Step 2 gap that actually causes it.

## Next steps

1. **Unbounded multi-face `conj`** — Steps 1 and 2 DONE 2026-07-31; the blocker is now **Step 3**.
   The `conjCPLQ` `isDomBounded` gate and the `quaPolToPlq:unboundedFace` rejection are both
   removed, and an unbounded QuaPol runs the whole pipeline. What each part needed, and what is
   left:

   **(a) DONE — `quaPolToPlq` no longer throws the ray away.** `faceVertexIndices` takes one
   vertex per edge, `E(j,1)`, and never consults QuaPol's ray flag `E(:,3)==0`, for which
   `E(j,1)` is the base point and `E(j,2)` the DIRECTION point; two rays off a shared apex both
   reported the apex and the bounded `domain()` constructor turned that into `NaN <= 0`. Faces
   carrying a ray are now built from HALF-PLANES (`faceDomainFromHalfPlanes`), one inequality per
   edge, segments and rays alike being "the line through `E(j,1)` and `E(j,2)`". Orientation is
   read off `P{k}`'s own sign convention (`orderEdges` stores `-j` when the face is on the LEFT),
   NOT off a centroid probe — an unbounded face has no bounded centroid to probe with.

   **(b) DONE for an AFFINE envelope — `plq_1p` no longer reads infinity markers as numbers.**
   Three separate fixes, and all three were needed:
   - `triangulate` sent unbounded faces down a fan that reads `vx/vy` as coordinates and rebuilds
     each triangle through the BOUNDED `domain(t,x,y)`. It now routes them to
     `fanUnboundedFace`, which works from the half-planes.
   - `convexEnvelope1`'s `nCE==0` branch was a closed form in three vertex coordinates. It now
     calls `convEnvUnbounded`.
   - `conjugateFunction` dispatched on `nCE = obj.d.nE`, which counts edges of positive finite
     SLOPE computed by walking consecutive `vx/vy` pairs — meaningless when some of those pairs
     are `intmax` markers, so an unbounded piece fell into whichever quadratic branch the garbage
     selected. It now dispatches on the ENVELOPE. A bounded piece with `nCE ~= 0` is deliberately
     left untouched and unprobed: the `nCE==1` envelope is RATIONAL and evaluating it at `(0,0)`
     divides by zero (this is live in `conjCPLQTest`, not hypothetical — it cost a suite run).

   **(c) DONE — the dual REGIONS were wrong too, and NOT for the reason first recorded.** The
   first diagnosis this session said Step 3 over-merged because `maximumP` assumes every piece's
   conjugate has full domain. That was WRONG, and the way it was wrong is worth keeping. The real
   defect was one step earlier: `region.getNormalConeVertex` reads the VERTEX LIST, and
   `getVertices` FABRICATES corner vertices at `(±intmax, ±intmax)` for an unbounded region — it
   tests `ptFeasible(vars,[intmax,intmax])` and appends whichever pass. Those corners are not
   vertices, and they displace real facets of the normal cone. On the first quadrant the corner
   `(intmax,intmax)` supplied the direction `(1,1)`, so the cone at the apex came out
   `{s1+s2 <= 0, s1 <= 0}` instead of `{s1 <= 0, s2 <= 0}` — a different set that agrees with the
   truth on the obvious probes and reports `f*(-10,5) = 0` where the answer is `+inf`.
   **This nearly shipped as "verified exact":** the first probe set for the wedge missed the
   wedge of error entirely, and only widening it caught the bug. New file `conjAffineOverPiece.m`
   builds the conjugate of an affine envelope from the HALF-PLANES instead, using
   `N_P(v) = {u : <u,e> <= 0 for each edge direction e leaving v}` — no vertex list, and a ray
   contributes its direction exactly as a bounded edge contributes its own. With it the 4-cone
   fan with `f = |x|+|y|` returns the indicator of `[-1,1]^2`, exact at all 10 probes; before, it
   returned the single piece `0` everywhere. The bounded path deliberately stays on
   `getNormalConeVertex`: there are no fabricated corners there, and it is what the whole
   existing suite exercises.

   **(c') REMAINING, and narrower than it first looked.** `f = |xy|` on the same fan (`x*y` on
   quadrants I/III, `-x*y` on II/IV) still errors: `MATLAB:unassignedOutputs`, `objR3` not
   assigned in `functionNDomain.maximumP`. Note what that answer IS: every face's envelope is 0,
   so every face's conjugate is the indicator of a quadrant and the max is the indicator of the
   single POINT `{0}` — a zero-dimensional region, which `region` reports as `region.empty`
   (`getVertices` returns `nv == 0`). So this reads as a DEGENERATE-INTERSECTION gap in Step 3,
   not the structural "assumes full domain" problem first recorded here. Confirm that before
   building anything on it. Fixture: the 4-cone fan in `unboundedFaceTest`.

   **(d) ALSO REMAINING — a CURVED convex envelope over an unbounded face.** A strictly convex
   `q` has `co q = q`, which is not affine, so the support-function construction does not apply;
   it needs cPLQ's T1/T2 active-set branch generalized off the bounded triangle. Currently
   REFUSED, loudly, as `plq_1p:conjugateFunction:unboundedNonAffine` — do not let this become a
   silent wrong answer, it is how the 1.15e18 errors arose. This is the headline
   `f=(x²+y²)/2` over `{x<=0,y>=0}` case whose answer is `min(s1,0)²/2 + max(s2,0)²/2`.
   NOTE: an earlier version of this file said this was "confirmed against brute force, exactly 4
   pieces". That confirmed the mathematical CLAIM, not any routine — no code in the tree produces
   those 4 pieces. `conjugateOfPiecePoly` on a bounded triangle with a strictly convex `f` returns
   only the 3 vertex pieces, so it would not either.

   **On reading [GARDINER].** Done, and the useful part is small. [GARDINER-10] is univariate.
   [GARDINER-13] §3.4.1 is the transferable idea: a primal face maps to the dual entity
   `F* = co{grad f(x_i) : F = co{x_i}}`, and "the formula extends to the case F is unbounded by
   considering points at infinity" — i.e. a ray contributes a DIRECTIONAL condition where a
   vertex contributes a point condition. That is exactly the 3-conditions-for-3-degrees-of-freedom
   structure `convEnvUnbounded` is built on (triangle = 3 vertex values; wedge = 1 value + 2
   directional derivatives; half-strip = 2 values + 1 directional derivative, its two rays being
   parallel). Its Fact 3 / Formula (3) — `f*(s) = (s-b)'Q^-(s-b)/2 + I_{b+Im Q}(s)` with `Q^-`
   the pseudo-inverse — is the closed form to reach for when (d) gets done. None of this required
   implementing engine `'pqp'`.

   **Latent bug found while probing (not fixed, not on any path today):**
   `functionNDomain.getInterior` reads `solve()`'s result struct by the literal field names
   `s12.s_1`/`s12.s_2`, so it throws `Unrecognized field name "s_1"` for any region whose vars are
   not named `s_1,s_2`. In-pipeline they always are (it runs on `maxConjugate`), which is why
   nothing hits it. `fieldnames(s12)` would make it var-name agnostic.

2. A native numeric rational branch in `conjPieceCPLQ` would buy **speed, not coverage**.
3. 0.1 tagging — **do not tag without being asked**. (The merge half of this item is DONE.)
4. `partialConj` unimplemented for every engine and type (`SUPPORT_MATRIX.md` §2).
5. Longstanding: `'pqp'`/`'graph'` engines; `RatPol.conj`/`biconj`/`add`; the
   `mergeL`/`removeTangent` exact-tie-point bug; `QuaPar.eval` wrong exactly *at* a result
   vertex (~1.4%); `testRegion/testCreation`.

## Do NOT redo these — all tried and reverted, with reasons

- **Returning recession directions as `(cos t, sin t)`.** `region.recessionRays` finds its
  candidates as angles, and rebuilding the direction from the angle is the obvious next line. It
  is wrong: a direction of a rational half-plane is rational, and the round trip puts `6.123e-17`
  where a `0` belongs. That is not cosmetic — the direction goes on to BUILD the sub-face
  half-planes in `fanUnboundedFace`, so `x <= 0` came back as
  `x - 4967757600021511/81129638414606681695789005144064*y <= 0`, a half-plane that is no longer
  pointed and whose fan then missed a whole boundary line (45 grid points). Keep the raw
  `[-a2, a1]` vector alongside its angle and return that. Pinned by
  `regionTest/recessionRaysReturnsExACTdirectionsNotTrigRoundTrips`.
- **Probing every envelope for affineness in `conjugateFunction`.** Deciding the dispatch needs
  the envelope's affine coefficients, and evaluating at `(0,0)`/`(1,0)`/`(0,1)` is the natural way
  to get them. But the `nCE==1` envelope is RATIONAL and its denominator vanishes at `(0,0)`, so
  an unconditional probe throws `symbolic:kernel:DivisionByZero` and took out 5 tests across
  `conjCPLQTest` and `cplqAdapterTest`. A bounded piece with `nCE ~= 0` must not be probed at
  all — its dispatch did not change, so it does not need the answer.

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
  `finiteVertices`; **`quadUnboundedBelow`** (the `-inf` gate) and **`recessionRays`** (the
  extreme rays, read from the INEQUALITIES and kept RATIONAL — read the two together, they share
  the same closed-form "the recession cone is an arc" argument) immediately before
  `redundantSubset`; `slopes2`; `getEdgeNosInf`; `simplifyOpenRegion1`;
  `simplifyUnboundedRegion`; `merge` (header records the correctness argument and the 36 → 125
  history).
- `functionNDomain.m` — `conjugateOfPiecePoly`: the two empty-domain guards, and the entry
  ADAPTER just before the scatter (both parts — the vertexless drop and the collision drop — with
  the piece-23 trace in the comment). Do NOT go looking at its `if obj(i).d.nv > 1` T2 gate for
  next-step 1 — an earlier version of this file sent a session there and it was a dead end. Its
  `getInterior` has the `s12.s_1` field-name bug noted in 1. Its `maximumP` is where next-step
  1(c), the Step 3 blocker, lives.
- `maxQuaPar.m` — `splitCell` (now handles a curved cell), `splitTwoArcPiece` (the one-arc-per-face
  chord subdivision), `arcHasStrictCrossing` (the tangency-vs-crossing check).
- `parabolaArcFrame.m` — `conicCoeffs`, the quartic restriction of a second conic to the parabola;
  `lineCoeffs`' companion.
- `regionTest.m` — 13 tests: `maxLinear`'s three answers, `linearForm`'s affine flags, four
  redundancy cases, three merge cases, two `quadUnboundedBelow` cases (the gate's answers, incl.
  the ray case that shows why eigenvectors alone are not enough; and the curved-facet refusal),
  and two `recessionRays` cases — the four cone shapes, and a REGRESSION pinning that a zero
  component comes back as exactly 0 rather than 6.1e-17.
- `unboundedFaceTest.m` — NEW, 9 tests, ~25 s. The fan's cover-and-no-overclaim property on a
  grid; the wedge and half-strip envelopes; the `f`-dependence regression; two exact end-to-end
  conjugates; and BOTH refusals (`-inf`, and a curved envelope over an unbounded face). The
  refusal tests are the load-bearing ones — they are what stops the 1.15e18 failure mode from
  coming back silently.
- `functionNDomainTest.m` — NEW, 2 tests. Reproduces `testcPLQ/testRectBiconj`'s piece 23
  directly, so the collision defect is pinned in ~26 s instead of that test's ~22 min. Verified
  to FAIL on pristine `HEAD~1` with the original error. It asserts the collision still exists as
  a precondition — if `getEdgeNosInf` ever stops producing it, the test stops covering anything.
- `fanUnboundedFace.m` — NEW. Covers a pointed unbounded face with triangles + one half-strip +
  one wedge; its header carries the proof that the fan covers (follow the ray from the apex
  through a point: it either leaves through a boundary piece or the point is a recession
  direction) and the reason a COVER suffices while a SUPERSET would not.
- `convEnvUnbounded.m` — NEW. The affine envelope for each of those three shapes, general in `f`.
  Its header carries the argument that matters: an affine minorant touching `q` is NOT on its own
  the envelope, since the sup of affine minorants is generally piecewise affine. The proof is
  that the gap-minimizing minorant comes out INDEPENDENT of the evaluation point, so one affine
  function is best everywhere at once. Also owns the `-inf` gate and the convex short-circuit.
- `conjAffineOverPiece.m` — NEW. The conjugate of an AFFINE envelope over one fan piece, built
  from the half-planes. Its header records the fabricated-corner defect it exists to avoid, and
  the one-line normal-cone formula it uses instead. Called only for UNBOUNDED pieces; bounded
  pieces stay on `getNormalConeVertex`.
- `conjCPLQ.m` — `assertStep3MatchesPieces` (the gate stays; it is a real invariant). The
  `isDomBounded` gate is GONE; unbounded domains take the same Case C route.
- `quaPolToPlq.m` — `faceDomainFromHalfPlanes` builds a ray-carrying face from half-planes, with
  orientation from `P{k}`'s sign. The header still records both original defects and the
  measurements behind them.
- `plq_1p.m` — `pieceVars` (the domain is the authority on which plane a piece lives in;
  `f.getVars` reports only the variables that OCCUR, so `q = -x^2` used to break every
  `vars(2)`); `affineParts` (the safe affine probe — guards the rational envelope);
  `triangulate`'s unbounded branch; `convexEnvelope1`'s unbounded branch and its now
  `f`-dependent `nCE==0`; `conjugateFunction`'s envelope-keyed dispatch, whose comment records
  both defects it replaced.
- `SUPPORT_MATRIX.md` §7 (new row: the `nCE==1/2` `x*y`-only envelope), §8 item 2 (rewritten:
  (a) and (b) done, (c) Step 3 is the blocker).

## Still true from before

- **There is no "3-convex-edge case" to implement.** [COAP] Appendix A.5's split reduces such a
  triangle to 2-convex-edge sub-triangles and Step 1 already applies it. Describe these cases by
  the **envelope's face count**, not the input's edge count.
- **The rational-face `0/0` did not need exact arithmetic** — it was cleanly removable, and
  `symbolicFunction.limitDirectional` resolves it by a directional limit.
- **Which half of cPLQ is reusable:** its **Step 2/3**, on **CCA2's** Step 1 output.
- The reverted `maxArray` `intmax`-vs-`inf` "fix" — don't reintroduce it; the comment at the
  site explains why a vertical constraint must fall into the arithmetic mean.

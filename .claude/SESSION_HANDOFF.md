# Session Handoff

_Last updated: 2026-07-29T00:00:00Z (superseded by the 2026-07-29 later session below)_

## STATE OF THE TREE: uncommitted work in progress

`region.m` and `regionTest.m` carry an **unverified but well-evidenced** fix for Step 3's
assembly (previous next-step 1). **Nothing was committed this session, deliberately** — see
"Why not committed" below. `git status` should show `M region.m` and `?? regionTest.m`.

Reference case throughout: `f = xy` over `T = conv{(0,0),(3,3),(1,2)}`, Step 1 envelope = 4 faces.

## Environment note (cost me the first 20 minutes)

MATLAB R2024b here is **network-licensed only** (`SERVER SLMS-SMATLABP1.ead.ubc.ca ... USE_SERVER`,
no local `license.dat`). Off the UBC VPN, `matlab -batch` dies with
`Licensing error: -96,7 / System Error: 11002`. **Connect the VPN before any session that needs to
run anything.**

## What was done: the Step 3 assembly fix

Both defects reduce to ONE primitive — **maximize a linear form over a polyhedron**, i.e. an LP,
which is exact and answers unboundedness/infeasibility as first-class results (these regions are
routinely unbounded). New static helpers `region.maxLinear` / `region.impliedBy`, new instance
helpers `region.linearForm` / `redundantSubset` / `deleteIfRedundant` / `unionIsExact`.

1. **`simplifyUnboundedRegion` dropping non-redundant constraints.** Replaced the proxy "delete any
   constraint not passing through a finite vertex" with the real test
   `max{g_i : g_j <= 0, j != i} <= 0`. Both deletion sites (`simplifyOpenRegion1`'s
   `obj.ineqs(mark)=[]` and `simplifyUnboundedRegion`'s `obj.ineqs(markF0)=[]`) now go through
   `deleteIfRedundant`, so a heuristic upstream can only ever *propose* a deletion. Erring is
   one-directional by construction: a non-affine `g_i` is never called redundant, and non-affine
   `g_j` are dropped from the set tested against (a relaxation ENLARGES the feasible set, so
   redundancy in the relaxation implies redundancy here). Candidates are judged against the
   constraints still standing, so two copies of one constraint cannot certify each other away.

2. **`merge` over-claiming.** With `A = A' n {g<=0}`, `B = B' n {g>=0}`, merge returns `M = A' n B'`.
   `M` never LOSES a point (any `x` in `A' n B'` has `g<=0`, so lies in `A`, or `g>=0`, so lies in
   `B`). It can GAIN points belonging to neither. `M = A u B` **exactly** when `A subset B'` and
   `B subset A'` — equivalently, when `A u B` is convex. `unionIsExact` decides that by LP before
   any facet is deleted. Also refuses when TWO facets are shared, where the "never loses a point"
   argument fails outright (a point with `g1<=0, g2>=0` is in neither operand yet survives `M`).
   The constraints TESTED must be affine; the region tested OVER enters only as its linear
   relaxation (a superset — sound), which lets a pair meeting along a shared PARABOLIC facet still
   be decided. That is the shape this codebase actually produces (arc between two rays).

   Note this is much WEAKER than the constraint-set-equality guard the previous session tried and
   reverted. That one was sound but refused far too much; the exact condition refuses only what it
   must. Its `36 -> 125` history is recorded in `merge`'s header — keep it.

3. **`region.slopes2` latent bug, made reachable by (1).** It takes a curved constraint's tangent at
   a region vertex lying on that curve, and a curved constraint need not have one — `pt` was left
   unassigned (`Unrecognized function or variable 'pt'`), or worse, retained the PREVIOUS
   iteration's point and silently computed a tangent at a point not on the curve. Unreachable while
   such constraints were being deleted. **Two variants were tried; see the open question below.**

## Measured results (17x17 dual grid, first-match evaluation, vs the pointwise max of the
## per-piece Step 2 conjugates — the same f* computed the other way)

Steps 1 and 2 remain **exact**: the per-piece max matches `sup_{x in T} <s,x> - xy` at all 289
points, in every run.

**BASELINE at HEAD (complete, 1768 s):**

| round | nRegions | wrong |
|---|---|---|
| seed | 5 | 0 |
| `* piece 2` | 14 | 0 |
| `maximumP 2` | 12 | **8** |
| `* piece 3` | 29 | **15** |
| `maximumP 3` | 22 | **8** |
| `* piece 4` | 47 | **8** |
| `maximumP 4` | 33 | **57** |

So the final assembled partition is wrong at **57 of 289** points. (The old handoff's "36" was
measured differently; 57 is this harness's number and the two are not comparable.)

**FIX, variant 1 (`slopes2` SKIPS a curved constraint with no vertex on it):**

| round | nRegions | wrong |
|---|---|---|
| seed | 5 | 0 |
| `* piece 2` | 14 | **0** |
| `maximumP 2` | 12 | **0**  (baseline 8) |
| `* piece 3` | 28 | **0**  (baseline 15) |
| `maximumP 3` | — | **stalled >90 min** (baseline 153 s) |

Every round measured is corrected. Test suite with variant 1: 8 suites, **0 failures**
(`PLQVCTest` 47, `QuaParTest` 10, `RatParTest` 12, `RatPolTest` 9, `addQuaParTest` 4,
`addQuaPolTest` 6, `biconjCPLQTest` 10, `clipArcByHalfPlaneTest` 7) before it too stalled inside
`conjCPLQTest` (which contains the 4-face case). New `regionTest`: **9/9**.

**FIX, variant 2 (`slopes2` falls back to the finite-vertex CENTROID) — currently in the tree, and
the one to keep:**

| round | nRegions | wrong | time (baseline) |
|---|---|---|---|
| seed | 5 | **0** | — |
| `* piece 2` | 14 | **0** | 51 s (48 s) |
| `maximumP 2` | 12 | **0** (baseline 8) | 96 s (88 s) |
| `* piece 3` | 28 | **0** (baseline 15) | 144 s (141 s) |
| `maximumP 3` | 25 | **0** (baseline 8) | 192 s (153 s) |
| `* piece 4` | 40 | **0** (baseline 8) | 273 s (302 s) |
| `maximumP 4` | 40 | **0** (baseline **57**) | 427 s (554 s) |
| TOTAL | | | **1645 s (1768 s)** |

**COMPLETE CLEAN FOLD.** Every round exact -- 0 disagree, 0 gap, 0 extra -- so the final assembled
partition goes from **57 of 289 wrong to 0 of 289**. And the whole fold is FASTER than baseline
(1645 s vs 1768 s), because a correct partition carries fewer regions than a damaged one: the last
two rounds are where it wins (273 vs 302, 427 vs 554). Keep variant 2; discard variant 1.

The `slopes2` fallback choice was therefore load-bearing for PERFORMANCE, not just for avoiding the
unassigned-`pt` crash: variant 1 (skip the constraint) stalled `maximumP 3` past 90 minutes, variant
2 (centroid tangent) runs it in 192 s. The mechanism is as hypothesised — fewer probe directions
make `region.maxArray` undecided more often, `maxEqDom` falls through to `splitmax3`, and every
undecided region splits, compounding round over round. Do not "simplify" `slopes2` back to skipping.

CAUTION when reading timings mid-run: the machine was heavily contended all session (16 cores at
~43% from other work, MATLAB getting ~6-12% of one core), so wall-clock between CHECK lines badly
overstates cost. Trust the per-round times the harness prints, not elapsed time between polls — an
earlier reading of this same run wrongly concluded `maximumP 3` had stalled when it had not.

## THE ONE OPEN QUESTION — resolve this first

**Why is `maximumP 3` slow under the fix, when rounds 1-3 are at baseline speed?** Baseline 153 s;
variant 1 >90 min; variant 2 >20 min. The fold's other five steps track baseline within 10% under
identical contention, so this is NOT "extra constraints cost more everywhere" — it is specific to
this one call.

Candidate mechanism, now only partly supported: fewer usable probe directions -> `region.maxArray`
returns undecided more often -> `maxEqDom` falls through to `splitmax3` -> every undecided region
splits -> region count compounds. Variant 1 (skip the constraint) fits this well; variant 2
(centroid tangent) was written to test it and DID help substantially (>90 min -> >20 min) without
fully curing it, so something else contributes too. Secondary suspect: `region()` -> `getVertices`
runs an O(n^2) symbolic `solve` over facet pairs, and the logs are dominated by
`symbolic:sym:isAlways:TruthUnknown` warnings, each a MuPAD round trip printing a stack trace.

A tempting "minimization" idea was considered and **rejected** — do not redo it: offering ALL
constraints (not just the heuristic's candidates) to `deleteIfRedundant` does not help, because the
candidate set already contains every constraint missing a finite vertex; widening it would only add
constraints that ARE active at vertices, and deleting a redundant-but-active constraint changes the
vertex list `getVertices` computes, which the whole pipeline reads.

## Commit status

The fold is COMPLETE and CLEAN (0/289 against baseline's 57/289, and faster).

**Full suite: 269 pass / 3 fail / 2 incomplete over 23 suites.** All 18 CCA2 suites are green,
including `conjCPLQTest` **20/0** (the rewritten 4-face test passes) and the new `regionTest` 9/0.
Against the 262/1 baseline that is **two NEW failures**, both in the vendored cPLQ suites and both
the SAME latent bug — the usual pattern here, one more index assumption that only became reachable
once irredundant constraints stopped being deleted:

    testMaxMultiRegion/testMax     ERRORED  MATLAB:badsubscript
    testcPLQ/testRectBiconj        ERRORED  MATLAB:badsubscript
      region.getEdgeNosInf:2917  <- functionNDomain.conjugateOfPiecePoly:1001 <- plq.biconjugateF:208

`vertexOfEdge(i)` returns `nv == 0` for a constraint with no vertex on the region, so `vx` is empty
and `vx(start)` overruns. `testRegion/testCreation` is the longstanding pre-existing failure and is
unrelated.

**Emitting `edgeNo(i) = 0` was TRIED AND REVERTED — do not redo it.** It is the right meaning (no
vertex => bounds no edge => no edge number) but `conjugateOfPiecePoly:1002` uses `edgeNo` as a
PERMUTATION, `d.ineqs(edgeNo) = d.ineqs`, so a 0 only moves the crash one line down ('Array indices
must be positive integers'). Both suites still failed. The comment at `getEdgeNosInf` records this.

**UPDATE — one of the two is FIXED; `testMaxMultiRegion` is back to 24/0.** Three layers, each
found by reading the actual contract rather than guessing at the data shape:

1. `getEdgeNosInf` reports **0** for a constraint with no vertex on the region (it bounds no edge).
2. `conjugateOfPiecePoly` **drops** those from its own local copy of `d` before the scatter. The
   whole routine is edge-indexed — the `isQuad` chord rewrite, `getNormalConeVertex`,
   `getSubdiffVertexT1/T2/T2Q` all address `d.ineqs` BY EDGE — so a vertexless facet has no
   representation there at all. (Parking them above the real slots was tried first and is worse
   than useless: the `isQuad` branch then builds a chord for such a constraint out of
   `d0.vx(1),d0.vx(2)`, vertices unrelated to it.) Dropping gives that routine exactly the
   information it had before `redundantSubset` started preserving these constraints.
3. Empty-domain guards in **both** of its loops — `removeTangent` can return `region.empty`, and
   `obj.nv` / `d0.ineqs` on a 0x0 region array is a comma-separated list with 0 values. The second
   guard must still set `ia(i+1) = size(pc,2)+1` before `continue`, since `ia` is a running index
   into `pc` with one entry per INPUT piece (`ia(i+1)==ia(i)` = "contributed nothing", the same
   convention `maxEqDom` uses).

**STILL OPEN: `testcPLQ/testRectBiconj`.** Fails in the same `conjugateOfPiecePoly` `isQuad`
branch, at `if isAlways(nineq.subsF(vars,[mx,my])>0)` -> `symbolicFunction.gtd:754`, with
'Conversion to logical from sym is not possible' (note `gtd` does a bare `if (obj1.f>obj2)`, so it
cannot take an undecidable sym at all -- `isAlways` never gets a chance to help).

What that means, and it is the useful part: dropping the vertexless constraints was NOT enough, so
the difference between old and new constraint sets is **not only** those. `simplifyUnboundedRegion`'s
second phase also deleted constraints that DO touch a vertex but fail its probe test, and
`redundantSubset` keeps those too unless they are provably redundant. So regions arriving at this
edge-indexed vendored model are richer in more than one way, and the chord it builds from
`d0.vx(1),d0.vx(2)` stops being meaningful.

The principled fix is an ADAPTER at `conjugateOfPiecePoly`'s entry: reduce `d` to the edge
description this model requires, rather than teaching the model about richer regions. Note it
already does shape normalization right there (`poly2order`/`poly2orderUnbounded`), so that is the
natural place. Do NOT instead weaken `redundantSubset` — those constraints are load-bearing for the
assembly (57 -> 0 wrong points).

**EARLIER, TRIED AND REVERTED — read this before attempting again.** Emitting `edgeNo(i)=0` AND
replacing the call site's scatter with a sort:

    hasEdge = edgeNo > 0;  [~,ord] = sort(edgeNo(hasEdge));  iHas = find(hasEdge);
    d.ineqs = [d.ineqs(iHas(ord)), d.ineqs(~hasEdge)];

on the reasoning that sorting by `edgeNo` equals scattering through it whenever `edgeNo` is a
permutation. **That reasoning is WRONG and it broke a THIRD test** (`testMaxMultiRegion/
testFractional`, which had been passing): `edgeNo(i) = j + add` with `j` a VERTEX index and
`add ∈ {0,1}` can exceed `numel(ineqs)`, so `d.ineqs(edgeNo) = d.ineqs` legitimately **GROWS** the
array — and duplicate `edgeNo` values make it last-write-wins. A sort reproduces neither. Whatever
replaces the scatter has to preserve both behaviours.

So the real question is what `conjugateOfPiecePoly` actually wants that array to BE afterwards
(it is indexed positionally downstream), which needs reading beyond the scatter line — that was
not done. Do that first; do not guess again.

Both attempts are reverted; `region.m` and `functionNDomain.m` are at the committed state, which
has the 2-failure regression and the explanatory comment at `getEdgeNosInf`.

## THIS IS THE STATE AT SESSION END — READ BEFORE TOUCHING ANYTHING

- The Step 3 assembly fix is **verified and correct**: 0/289 on the reference case (baseline
  57/289), all 18 CCA2 suites green.
- It **regresses two vendored cPLQ biconjugate tests**, both the single `getEdgeNosInf` /
  `conjugateOfPiecePoly` issue above. 269/3 against a 262/1 baseline.
- **Nothing was committed.** Do not commit until those two are green (expected then: 271/1).
- Do NOT "fix" the regression by weakening `redundantSubset` back toward the old delete-happy
  behaviour — that would re-break the assembly. The constraints it now keeps are load-bearing.

Already updated for the fix, so the commit is otherwise ready:
- `conjCPLQTest.indefiniteTriangleThreeConvexEdgesUsesStep3` — rewritten: it used to PIN the loud
  `PLQ:conjCPLQ:cplqFailed`, and now checks the answer against `supBilinearOverPoly` at 8 dual
  points, including `s=(1,1)` (baseline gave 1.0 for 1.125) and `s=(-3,-3)` (baseline gave 9 for 0).
- `conjCPLQ.m` — `assertStep3MatchesPieces`'s header and error text: the gate STAYS (a wrong
  partition fails silently by nature), but now points at the LP certificates as the first thing to
  check rather than describing the assembly as broken.
- `SUPPORT_MATRIX.md` — §1.2's 4-face row is now **OK**; the assembly-defect block is rewritten as
  resolved, keeping the `36 -> 125` history and the `slopes2` performance finding; §8's blocker
  list renumbered (unbounded multi-face is now blocker 2).

## Next steps

1. **Finish verifying the fix in the tree.** On an idle machine: one full `run.m` fold (expect all
   rounds 0 wrong, ending 0/289 against baseline's 57/289) and one full suite (expect the 262/1
   baseline plus 9 new `regionTest` cases = 271/1). Decide variant 1 vs variant 2 of `slopes2` on
   the measured timing. Then update `SUPPORT_MATRIX.md` section 1.2 (the two assembly defects are
   the "Update (2026-07-29)" block) and section 8 blocker 2, rewrite
   `conjCPLQTest.indefiniteTriangleThreeConvexEdgesUsesStep3` (it currently PINS the loud
   `PLQVC:conjCPLQ:cplqFailed`; it must become a correctness check), and soften
   `conjCPLQ.assertStep3MatchesPieces`'s message (keep the gate — it is a good invariant).
   The harness lives in the scratchpad: `run.m` (setup + round-by-round fold, reads `$CCA2DIR`),
   `suite.m`, `smoke.m`, and a pristine HEAD copy under `baseline/` made with
   `git archive HEAD | tar -x -C baseline`.

2. **Unbounded multi-face `conj`** — scoped this session, and it breaks EARLIER and more quietly
   than `conjCPLQ.m:103`'s guard suggests. `quaPolToPlq` feeds ray DIRECTION POINTS to `domain()`
   as if they were vertices, so for the 4-cone fan `V=[0 0;-1 0;0 1;1 0;0 -1]`, `E=[1 2 0;1 3 0;
   1 4 0;1 5 0]` the second-quadrant cone comes out as a degenerate 2-vertex "polygon" whose
   inequality list is `x <= 0` **twice** — silently wrong, not an error. Behind that,
   `plq_1p.triangulate` is a pure vertex fan and `plq_1p.convexEnvelope1`/`conjugateFunction`
   index `vx(1),vx(2),vx(3)` directly, so cPLQ has no unbounded-piece case at all in Step 1 or
   Step 2. Route: build face domains from half-planes via `domain.domainEdge` (which takes
   inequalities, so it handles unbounded regions) instead of from a vertex list; then the
   irreducible new work is the conjugate of a quadratic over a **wedge** (two rays from a common
   apex), since an unbounded convex face fans into bounded triangles plus one such wedge.
   *Cheap first increment, worth doing on its own:* make `quaPolToPlq` REJECT an unbounded face
   loudly instead of silently corrupting it.

3. **`maxQuaPar`: split a cell that already carries an arc** (~26%, 30 of 115 sampled splits) —
   also scoped, and the cheap route does NOT work. On the fixture that trips the guard
   (`maxQuaParTest.maxQuaParRejectsSplittingACellThatAlreadyCarriesAnArc`), 15 of the 18 candidate
   splitting curves `{f1row-f2row=0}` are pure straight lines, BUT every curved cell comes from the
   one curved face, and the three splitting curves meeting it are exactly the non-line ones (one
   pair-of-lines, two genuine parabolas). So the common case really does need conic-conic
   intersection. Parametrizing the arc by `parabolaArcFrame`'s global monotone `u` and substituting
   into the splitting conic gives a quartic in `u` — tractable. The harder half is representation:
   each half can need TWO curved edges (a sub-arc of the original plus the new one), which the
   `pieces` struct's single `curveAfter`/`curveEc` slot and `facePoly`'s one-curved-edge assertion
   both forbid. QuaPar's own per-edge `Ec` CAN hold several, so the limit is maxQuaPar-internal.
   Multi-arc pieces is the natural unit, as previously noted.

4. `partialConj` unimplemented for every engine and type (`SUPPORT_MATRIX.md` section 2).
5. A native numeric rational branch in `conjPieceCPLQ` would buy **speed, not coverage**.
6. Then 0.1 tagging — **do not tag without being asked**.
7. Longstanding: `'pqp'`/`'graph'` engines; `RatPol.conj`/`biconj`/`add`; the
   `mergeL`/`removeTangent` exact-tie-point bug; `QuaPar.eval` wrong exactly *at* a result vertex
   (~1.4%); `testRegion/testCreation`.

## Relevant files

- `region.m` — the LP-certificate block at the top (after `probeAlong`/`probePerp`), the instance
  helpers just before `finiteVertices`, `slopes2`, `simplifyOpenRegion1`, `simplifyUnboundedRegion`,
  and `merge` (whose header now records both the correctness argument and the `36 -> 125` history).
- `regionTest.m` — NEW, 9 tests: `maxLinear`'s three answers, `linearForm`'s affine flags, four
  redundancy cases (including "never delete both copies" and the unbounded irredundant constraint),
  and three merge cases (union-is-convex merges, the L-shape that must be refused, two shared
  facets refused).
- `conjCPLQ.m` — `assertStep3MatchesPieces`; `conjCPLQ.m:103` for next-step 2.
- `maxQuaPar.m` — `splitCell`'s `pieceIsCurved` guard for next-step 3.

## Still true from before

- **There is no "3-convex-edge case" to implement.** [COAP] Appendix A.5's split reduces such a
  triangle to 2-convex-edge sub-triangles and Step 1 already applies it. Describe these cases by
  the **envelope's face count**, not the input's edge count.
- **The rational-face `0/0` did not need exact arithmetic** — it was cleanly removable, and
  `symbolicFunction.limitDirectional` resolves it by a directional limit.
- **Which half of cPLQ is reusable:** its **Step 2/3**, on **CCA2's** Step 1 output.
- The reverted `maxArray` `intmax`-vs-`inf` "fix" — don't reintroduce it; the comment at the site
  explains why a vertical constraint must fall into the arithmetic mean.

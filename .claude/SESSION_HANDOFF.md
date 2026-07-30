# Session Handoff

_Last updated: 2026-07-30_

## What happened this session

**Step 3's assembly is fixed and verified.** Both defects reduced to one primitive — maximize a
linear form over a polyhedron, i.e. an LP — and on the reference case `f = xy` over
`conv{(0,0),(3,3),(1,2)}` the assembled partition went from **57 of 289 grid points wrong to 0**,
exact at every one of the fold's seven rounds, and the whole fold got **faster** (1645 s vs
1768 s) because a correct partition carries fewer regions than a damaged one.

That fix made three latent bugs reachable in code that had been shielded by the old
delete-happy behaviour; two are fixed, one is still open (below).

## Where things stand

- Branch: `step3-assembly-lp-certificates` @ `cd985c2` — "Fix testMaxMultiRegion: teach
  conjugateOfPiecePoly about vertexless facets" (3 commits ahead of `main`)
- Pushed: no upstream configured for this branch (`origin` exists:
  `https://github.com/ylucet/CCA2.git`)
- Working tree: clean
- **Suite: 270 pass / 2 fail** over 23 suites, against a 262/1 baseline. All 18 CCA2 suites are
  green (`conjCPLQTest` 20/0, new `regionTest` 9/0). The 2 failures are
  `testcPLQ/testRectBiconj` (NEW — see below) and `testRegion/testCreation` (longstanding,
  unrelated).

## Next steps

1. **`testcPLQ/testRectBiconj` — the one open regression.** Fails in
   `functionNDomain.conjugateOfPiecePoly`'s `isQuad` branch at
   `if isAlways(nineq.subsF(vars,[mx,my])>0)` → `symbolicFunction.gtd:754`, "Conversion to
   logical from sym is not possible". Note `gtd` does a bare `if (obj1.f>obj2)`, so it cannot
   take an undecidable sym at all — `isAlways` never gets a chance to help.

   **Diagnosis (established, not guessed):** dropping the vertexless constraints was NOT enough,
   so the difference between the old and new constraint sets is *not only* those.
   `simplifyUnboundedRegion`'s second phase also deleted constraints that DO touch a vertex but
   fail its probe test, and `redundantSubset` keeps those unless provably redundant. Regions
   arriving at this edge-indexed vendored model are therefore richer in more than one way, and
   the chord it builds from `d0.vx(1),d0.vx(2)` stops being meaningful.

   **Route:** an ADAPTER at `conjugateOfPiecePoly`'s entry that reduces `d` to the edge
   description that model requires — it already does shape normalization right there
   (`poly2order`/`poly2orderUnbounded`), so that is the natural place.
   **Do NOT weaken `redundantSubset` instead** — its extra constraints are what took the
   assembly from 57 wrong points to 0.

2. **Unbounded multi-face `conj`** — scoped this session. Breaks EARLIER and more quietly than
   `conjCPLQ.m:103`'s guard suggests: `quaPolToPlq` feeds ray DIRECTION POINTS to `domain()` as
   if they were vertices, so for the 4-cone fan `V=[0 0;-1 0;0 1;1 0;0 -1]`,
   `E=[1 2 0;1 3 0;1 4 0;1 5 0]` the second-quadrant cone comes out as a degenerate 2-vertex
   "polygon" whose inequality list is `x <= 0` **twice** — silently wrong, not an error. Behind
   that, `plq_1p.triangulate` is a pure vertex fan and `convexEnvelope1`/`conjugateFunction`
   index `vx(1),vx(2),vx(3)` directly, so cPLQ has no unbounded-piece case at all.
   Route: build face domains from half-planes via `domain.domainEdge` (it takes inequalities, so
   it handles unbounded regions); the irreducible new work is then the conjugate of a quadratic
   over a **wedge**. *Cheap first increment, worth doing alone:* make `quaPolToPlq` REJECT an
   unbounded face loudly instead of silently corrupting it.

3. **`maxQuaPar`: split a cell that already carries an arc** (~26%, 30 of 115 sampled splits) —
   also scoped, and the cheap route does NOT work. On the guard-tripping fixture
   (`maxQuaParTest.maxQuaParRejectsSplittingACellThatAlreadyCarriesAnArc`), 15 of 18 candidate
   splitting curves are pure straight lines, BUT every curved cell comes from the one curved
   face, and the three curves meeting it are exactly the non-line ones (one pair-of-lines, two
   genuine parabolas). So the common case really does need conic-conic intersection.
   Parametrizing the arc by `parabolaArcFrame`'s global monotone `u` gives a quartic in `u` —
   tractable. The harder half is representation: each half can need TWO curved edges, which the
   `pieces` struct's single `curveAfter`/`curveEc` slot and `facePoly`'s one-curved-edge
   assertion both forbid. Multi-arc pieces is the natural unit.

4. `partialConj` unimplemented for every engine and type (`SUPPORT_MATRIX.md` §2).
5. A native numeric rational branch in `conjPieceCPLQ` would buy **speed, not coverage**.
6. Merge to `main` and then 0.1 tagging — **do not tag without being asked**.
7. Longstanding: `'pqp'`/`'graph'` engines; `RatPol.conj`/`biconj`/`add`; the
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
- `functionNDomain.m` — `conjugateOfPiecePoly` (the open regression, and the two empty-domain
  guards + the vertexless-facet drop).
- `regionTest.m` — NEW, 9 tests covering `maxLinear`'s three answers, `linearForm`'s affine
  flags, four redundancy cases, and three merge cases.
- `conjCPLQ.m` — `assertStep3MatchesPieces` (the gate stays; it is a real invariant);
  `conjCPLQ.m:103` for next-step 2.
- `maxQuaPar.m` — `splitCell`'s `pieceIsCurved` guard for next-step 3.
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

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

## Where things stand

- Branch: `step3-assembly-lp-certificates` — 4 commits ahead of `main`
- Pushed: the first 3 commits only. The `testRectBiconj` fix is committed locally, **not pushed**.
- Working tree: clean
- **Suite: 271 pass / 1 fail** over 23 suites, against a 270/2 baseline. `testcPLQ` is 8/0. The
  single failure is `testRegion/testCreation` (longstanding, unrelated). Plus the NEW
  `functionNDomainTest` 2/0, created after that sweep's glob — 273/1 all told.
- **`main` is now unblocked**: the merge gate ("merge only once `testcPLQ/testRectBiconj` is
  green") is satisfied.

## Next steps

1. **Unbounded multi-face `conj`** — scoped this session. Breaks EARLIER and more quietly than
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

2. **`maxQuaPar`: split a cell that already carries an arc** (~26%, 30 of 115 sampled splits) —
   also scoped, and the cheap route does NOT work. On the guard-tripping fixture
   (`maxQuaParTest.maxQuaParRejectsSplittingACellThatAlreadyCarriesAnArc`), 15 of 18 candidate
   splitting curves are pure straight lines, BUT every curved cell comes from the one curved
   face, and the three curves meeting it are exactly the non-line ones (one pair-of-lines, two
   genuine parabolas). So the common case really does need conic-conic intersection.
   Parametrizing the arc by `parabolaArcFrame`'s global monotone `u` gives a quartic in `u` —
   tractable. The harder half is representation: each half can need TWO curved edges, which the
   `pieces` struct's single `curveAfter`/`curveEc` slot and `facePoly`'s one-curved-edge
   assertion both forbid. Multi-arc pieces is the natural unit.

3. A native numeric rational branch in `conjPieceCPLQ` would buy **speed, not coverage**.
4. Merge to `main` (now unblocked) and then 0.1 tagging — **do not tag without being asked**.
5. `partialConj` unimplemented for every engine and type (`SUPPORT_MATRIX.md` §2).
6. Longstanding: `'pqp'`/`'graph'` engines; `RatPol.conj`/`biconj`/`add`; the
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
  the piece-23 trace in the comment).
- `regionTest.m` — 9 tests covering `maxLinear`'s three answers, `linearForm`'s affine flags,
  four redundancy cases, and three merge cases.
- `functionNDomainTest.m` — NEW, 2 tests. Reproduces `testcPLQ/testRectBiconj`'s piece 23
  directly, so the collision defect is pinned in ~26 s instead of that test's ~22 min. Verified
  to FAIL on pristine `HEAD~1` with the original error. It asserts the collision still exists as
  a precondition — if `getEdgeNosInf` ever stops producing it, the test stops covering anything.
- `conjCPLQ.m` — `assertStep3MatchesPieces` (the gate stays; it is a real invariant);
  `conjCPLQ.m:103` for next-step 1.
- `maxQuaPar.m` — `splitCell`'s `pieceIsCurved` guard for next-step 2.
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

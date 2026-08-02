# Session Handoff

_Last updated: 2026-08-02_

## What happened this session

**The biconjugate works, including the convex envelope over a box.** `x*y` on `[0,1]^2` now
returns `max(0, x+y-1)` — the McCormick envelope — exact to 1.1e-16. Two defects stood in the
way: the second pass was missing the algorithm's final MAX entirely, and
`conjugateOfPiecePoly` returned wrong regions for pieces whose domain carried a redundant
constraint. Also tagged and pushed **v0.1** so the SCIP spike can run its benchmark, after
finding that the `QuaPoly` -> `QuaPol` rename had silently broken SCIP's bridge.

## Where things stand

- Branch: `main` @ `64d5665` — "Make the biconjugate work: add the missing max, fix
  conjugateOfPiecePoly's regions"
- Pushed: pending (see Step 6 outcome in the session summary; `v0.1` and everything up to
  `2288dc1` are already on `origin/main`)
- **Tag `v0.1` is PUSHED and public** (`github.com/ylucet/CCA2`). Moving it is no longer cheap —
  anything further should be `v0.1.1`.
- **Suite: incomplete.** A full sweep was still running at session end (9 of 26 suites, all green
  except the one known failure below). The last COMPLETE sweep was 293/0 over 25 suites, at the
  tag, before this session's changes. **Re-run `.claude/suite.sh` first thing.**

### THE ALGORITHM (get this right; it was mis-stated repeatedly)

```
f  -> triangulate -> convex envelope of f+I_T per triangle T (splitting as needed)
                  -> conjugate of g+I_T per piece
                  -> MAX of all those conjugates                       = f*

h = f* on pieces P (P polyhedral, or with a single parabolic edge)
   -> conjugate of h+I_P per piece
   -> MAX of all those conjugates                                      = f**
```

Two things that cost hours by being confused:
- **Step 1's output is NOT the convex envelope of f over its domain.** It is a per-triangle
  intermediate. The envelope over a box/polygon is `f**`, i.e. two passes. `convEnvCPLQ`'s own
  header says so; do not read its result as "the envelope" and conclude a shape is unsupported.
- **The max is the same code in both passes**, and must not be duplicated. It is
  `functionNDomain.maxOfList`, called by `plq.maximumConjugate`, `plq.biconjugateF` and
  `QuaParCPLQ.conj`. The conjugate's version is the general one (operands may carry parabolic
  edges); the biconjugate's is a strict special case (polyhedral only).
- `ia` from `conjugateOfPiecePoly` is only the block delimiter (`pc(ia(k):ia(k+1)-1)` is piece
  k's own conjugate). It exists for grouping; nothing deeper.
- `region.plus` is INTERSECTION — it concatenates the inequality lists. "Adding equations" IS
  intersection. `functionNDomain.mtimes` (`*`) does NOT multiply: it intersects each pair of
  regions and stores BOTH functions there so `maximumP` has two values to compare.

## Next steps

1. **Re-run the full suite** (`CCA2_TEST_TIMEOUT=5400 bash .claude/suite.sh`) and record the
   total. Expect 1.5–2 h. Nothing below should be trusted until this is green.

2. **The two-face / PARABOLIC defect — the one failing test**, and the sharpest lead available.
   `biconjugateTest/biconjugateOverATwoFaceSubdivisionIsTheEnvelope`. The unit square with
   `f = x*y` is CORRECT as one face and WRONG as two. Both first conjugates are verified correct
   against brute force (one-face `f*` 6 regions, two-face `f*` 9 regions, identical values at 7
   dual points). The two-face `f*` carries parabolic pieces, and there the second pass returns
   the PER-FACE envelope — at every probe on the shared diagonal `y = x` it returns `x*y` itself,
   which is exactly what a diagonal split fails to relax.
   **The lead:** blocks 1 and 4 of that `f*` differ ONLY in the sign of the conic
   `2*s1*s2 - 4*s1 + s1^2 + s2^2` — opposite sides of one parabola — yet produce IDENTICAL
   conjugate cells. At least one is wrong. Branch: `getNormalConeVertexQ` /
   `getSubdiffVertexT2Q` / the `isQuad` chord rewrite. Already ruled out: replacing the chord's
   hardcoded `vx(1),vx(2)` endpoints with the vertices the conic actually touches — changed
   nothing, and was reverted. Full diagnosis is in the test file at the failure site.

3. **Performance — rank this high.** Cost scales with the number of `f*` regions, which grows
   fast with vertex count. A PENTAGON's `f*` alone takes **885 s and yields 41 regions**, after
   which the second pass must conjugate and max all 41; it was tried as a test and removed.
   `conjCPLQTest` is ~1 h; the full suite 1.5–2 h. This caps the toolbox at small polygons.

4. **Step 3's cross-piece max drops cells on SOME unbounded assemblies.** Data-dependent: the
   4-cone fan with convex faces is exact under `F=[1 2;2 3;3 4;4 1]` and loses 12 of 16 cells
   under `F=[3 2;2 1;1 4;4 3]`. Caught by `conjCPLQ`'s `assertStep3MatchesPieces` and raised, not
   returned. Pinned by `conjCPLQTest/step3DropsCellsOnSomeUnboundedAssemblies`.

5. **`QuaPar.eval` at a vertex — fix is UNVERIFIED.** A relative tolerance replaced the exact
   `all(vals <= 0, 2)`; it matches the recorded mechanism and does not regress, but the ~1.4%
   measurement's generator was a throwaway script with no recorded seed and could not be
   reproduced. Needs a committed, SEEDED sweep per `SUPPORT_MATRIX.md` §0.1.

6. **Reproducibility debt.** `SUPPORT_MATRIX.md` §0.1 now requires any quoted measurement to name
   a committed, seeded generator, and lists three existing claims that fail that bar (the
   1.4%/0.8%, the 115/85 splits, the 58→76). The committed suite itself has no randomness at all.

7. **`mergeL`/`removeTangent` exact-tie-point** — the remaining §7 defect that returns a wrong
   number rather than erroring.

8. Unimplemented: `partialConj` (all engines/types), `'pqp'`/`'graph'` engines,
   `RatPol.conj`/`biconj`/`add`.

9. **`QuaParCPLQ` data-shape wart.** Uniform in TYPE (all `conj` results are `RatPar`, with
   `kind()`), but it stores `fnd` and leaves `V/E/Ec/f/F` empty, and its `add` refuses a
   `QuaPol`/`QuaPar` operand, so cross-case `infConv` errors. `RETURN_TYPE.md` records this.

## SCIP (`AI/spike/SCIP`) — read before renaming anything

SCIP bridges in through the MATLAB Engine and calls **exactly one entry point, `convEnvCPLQ`**,
via its own glue `SCIP/src/cca2ConvexEnvelope.m`, consuming the returned `RatPol`'s
`V/E/F/f(:,5:10)/den` as plain arrays. It never calls `conj`, `biconj`, `partialConj`, or any
`QuaPar` method — so **none of the open items above is on its path**.

Its bridge was broken by the 2026-07-27 `QuaPoly` -> `QuaPol` rename, which left no shim on the
stated grounds that "nothing external can depend on it". `QuaPoly.m` is now an alias (as `PLQVC`
already was) and the bridge is verified end to end against SCIP's own Phase 2 reference instance:
`conv(xy)` over `conv{(1,1),(0,0),(2,0)}` returns `2y^2/(y-x+2)`, the published [COAP] A.3.3
closed form. See `SUPPORT_MATRIX.md` §0.0. **Check it before renaming any public name.**

## Do NOT redo these — tried, with reasons

- **Reading `convEnvCPLQ`'s output as "the convex envelope" and concluding a box is unsupported.**
  It is Step 1's per-triangle intermediate. The box's envelope is `f**`.
- **Replacing `addEq`'s intersection with a union.** Adding equations IS intersection, and that is
  correct: `f**` is finite only where EVERY piece's conjugate is finite. `addEq` is not the bug.
- **Adding a second, private max inside `biconjugateF`.** Use `functionNDomain.maxOfList`.
- **The `isQuad` chord endpoints from the conic's own vertices** (instead of `vx(1),vx(2)`).
  Reasonable, but changed nothing on the two-face case; reverted. Needs new evidence first.
- **Returning recession directions as `(cos t, sin t)`.** A direction of a rational half-plane is
  rational; the round trip puts `6.123e-17` where a `0` belongs, and that direction goes on to
  BUILD sub-face half-planes. Pinned by `regionTest/recessionRaysReturnsExACTdirectionsNotTrigRoundTrips`.
- **Probing every envelope for affineness in `plq_1p.conjugateFunction`.** The `nCE==1` envelope
  is RATIONAL and its denominator vanishes at `(0,0)`; an unconditional probe threw
  `symbolic:kernel:DivisionByZero` across 5 tests. `plq_1p.affineParts` guards it with
  `isPolynomial`.
- **`merge` guarded by constraint-set EQUALITY** — provably sound and makes things WORSE
  (36 → 125 wrong of 289). `unionIsExact` refuses only what it must.
- **`slopes2` SKIPPING a curved constraint with no vertex** — `maximumP 3` went 153 s → >90 min.
- **Replacing `conjugateOfPiecePoly`'s scatter with a sort** — the scatter legitimately GROWS the
  array and is last-write-wins; sorting reproduces neither and broke
  `testMaxMultiRegion/testFractional`.
- **PARKING vertexless constraints above the real edge slots** — the `isQuad` branch then chords
  them from unrelated vertices.

## Environment / harness

- MATLAB R2024b is **network-licensed only**. Off the UBC VPN, `matlab -batch` dies with
  `Licensing error: -96,7`. **Connect the VPN before any session that runs anything.**
- **`.claude/suite.sh`** — one MATLAB per suite under `timeout`; default 3600 s,
  `CCA2_TEST_TIMEOUT` to override, `0` disables, `TIMEOUT` reported separately from failure.
  Use this, not `.claude/suite.m`, which cannot time anything out.
  Also takes suite names: `bash .claude/suite.sh biconjugateTest conjCPLQTest`.
- Pristine baseline for A/B: `git archive HEAD | tar -x -C /tmp/cca2base`.
- Long waiters must watch a **per-run** log filename; a reused name makes the exit condition
  unsatisfiable.

## Relevant files

- `functionNDomain.m` — `maxOfList` (the shared max, in a `methods (Static)` block at the end);
  `conjugateOfPiecePoly`'s preprocessing: the redundancy drop (0), the vertexless drop (1) and the
  collision drop (2) including the tie rule, all on the LOCAL copy; `maximumP` (its `objR3`
  initializer, and the Step 3 cell-dropping of next-step 4); `addEq` (correct as-is);
  `getInterior` (indexes `c2(2)` under a guard testing only `size(c1,2)` — next in the latent
  chain, NOT fixed).
- `plq.m` — `maximumConjugate` and `biconjugateF`, both now calling `maxOfList`.
- `QuaParCPLQ.m` — `conj`, also on `maxOfList`; the `emptyResult` guard.
- `biconjugateTest.m` — NEW, 6 tests, one failing by design. Read the failure-site comment before
  touching the parabolic branch.
- `region.m` — `getEdges`, `splitmax3`, `getNormalConeVertexQ`, `poly2orderUnbounded` (all four
  had an output/index used outside its guard; all fixed); `redundantSubset`/`unionIsExact` (LP
  certificates); `recessionRays` and `quadUnboundedBelow` (read together).
- `plq_1p.m` — `quadKind`/`quadParts` (classify by eigenvalue signs, never by `nCE`);
  `convexEnvelope1`; `conjugateFunction`'s envelope-keyed dispatch; `affineParts`; `pieceVars`.
- `xyFrame.m` / `transformDomain.m` / `substituteFrame.m` — the indefinite-quadratic frame change,
  built with EXACT symbolic arithmetic (as doubles it took one case from seconds to >20 min).
- `fanUnboundedFace.m`, `convEnvUnbounded.m`, `conjConvexOverPiece.m`, `unboundedFaceTest.m` —
  the unbounded-face work.
- `QuaPoly.m` — the alias SCIP depends on.
- `SUPPORT_MATRIX.md` — §0.0 downstream consumers, §0.1 the reproducibility rule, §7 known
  wrong-answer defects, §8 the blocker list.

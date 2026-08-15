# Morning report — 2026-08-15 overnight run

Branch: **`main`** (merged, and pushed as the run went — you authorised pushing for this run)

Task as given: housekeeping (merge the previous branch, turn the piece invariants on, retire the
stale TODO entry, commit and push), then fix bugs 1–5. Later, explicitly: fix bug 5. Standing
instructions: parallel agents permitted, do not stop after three failures, do not wait for input.

## Headline

**FOUR of the five bugs are FIXED — 2, 3, 4 and 5 — and bug 1 went from 0 to 5 of 7 probe
points.** The whole repository is down to **ONE red**.

`maxQuaParTest` 29 / 0, `conjCPLQTest` 25 / 0, `unboundedFaceTest` 18 / 0, fast bucket 204 / 0,
slow bucket **114 / 1**. `maxQuaPar` has no open case (seeded sweep **18 / 0 / 0**, from 16 / 0 / 2
two sessions ago), and **`SUPPORT_MATRIX.md` §8's second release blocker — unbounded multi-face
domains — is closed**: Steps 1, 2 and 3 are all done, the last two pieces on the same day.

Three of the four fixes were one lesson, in three different routines:

> **A degenerate geometric object is not a geometric object, and this codebase keeps assuming it
> is.** Two arcs sharing a vertex have no separating CHORD, because the candidate chords are the
> arcs' own edges (bug 5). A quadratic has no TANGENT LINE at its own apex, where the gradient
> vanishes (bug 2) — and §8.2(e) already recorded a *different* routine falling into that exact
> trap on that exact input. A lens's two edges join the same two vertices, so numbering an edge by
> an ENDPOINT cannot tell them apart (bug 1). In each case a routine computed a geometric object
> from a degenerate configuration and used it without checking that the object exists.

## What changed

- **Housekeeping — DONE and pushed** (`5f22486`). Previous branch merged into `main`.
  `MAXQP_ASSERT = 1` is now on for every `maxQuaParTest` test, via a `TestMethodSetup` that
  restores the previous value on teardown. The stale `twoCurvedWhereTheSplitCurveCrossesAnArc`
  entry is retired.

- **BUG 5 — FIXED** (`db4a188`). `splitTwoArcPiece` found no cut when a piece's two arcs are
  ADJACENT: its two candidate chords join the arcs' facing endpoints, which for arcs sharing a
  vertex ARE the arcs' own edges, so both chains came out too short and the piece was returned
  unsplit with one arc flattened to its chord — surfacing three stages later as an orphan
  boundary edge. Generalised the `nv == 3` shared-vertex fallback to `nv >= 4` with the ordinary
  diagonal to a non-adjacent vertex. **Sweep 17 / 0 / 1 → 18 / 0 / 0.** Pinned by
  `arcVsArcSplitsTwoADJACENTArcsOnAPieceWithADiagonal`, which exists because the sweep test
  asserts `nWrong == 0` and would have counted this input in `nErr`.

- **BUG 1 — three defects fixed; `f**` exact at 5 of 7 probe points, was 0 of 7**
  (`b105814`, `508b336`). In the order they had to be found:
  1. `getEdgeNosInf` numbers an edge by one of its endpoint VERTICES, so a LENS — two edges
     joining the SAME pair — gets one number for both and the last-write-wins scatter destroys
     one. That, not the `size(ineqs,2) == nv` count the record blamed, is the "conjugate of the
     chord".
  2. `getNormalConeVertexQ` indexed its second constraint as `j+1` unwrapped, raising
     `badsubscript` on any BOUNDED region — which is why its only caller sent every bounded region
     to the POLYHEDRAL routine, whose cones come from the CHORD and are wrong for a curved edge.
  3. `biconj` handed its second conjugation the curved MESH `conj` has returned since 2026-08-13;
     `quaPolToPlq` refuses a curved domain, so it died before reaching any of the above. It now
     asks `conjCPLQ` for the symbolic form on purpose.
  Unit level: the half-lens conjugates to 3 cells exact against a brute-force sup at all 10 probe
  points (2 identical wrong cells before).

- **BUG 2 — FIXED** (`1d04c75`). `region.removeTangent` builds the TANGENT LINE to a quadratic
  constraint at a vertex and deletes any constraint matching it. At the APEX OF A CONE the
  gradient vanishes, there is no tangent, and what it computes is meaningless — and that apex is
  exactly where an unbounded fan's Step 3 split conics meet. It deleted `−s1 ≤ 0`, leaving two
  constraints **blind to the sign of `s1`**, so the region went symmetric and claimed the mirror
  wedge. `f*(-3,-2.4)` is now **4.5**, the truth; it was 5.130. **`conjCPLQTest` 25 / 0.**
  Cleared on the way, so nobody re-checks them: `redundantSubset`, `simplifyUnboundedRegion`, and
  Step 2 itself. The old `step3DropsCellsOnSomeUnboundedAssemblies` pinned the GATE firing; it no
  longer does, so it is renamed and rewritten to pin what the gate protects.

- **BUGS 3 and 4 — FIXED** (`1b98e30`). These were a missing ALGORITHM rather than a defect:
  `convEnvUnbounded` computed only the AFFINE envelope and refused the rest by design. Two shapes
  are now derived and implemented, each with its proof in the source:
  a **wedge** with one flat and one convex ray, whose envelope is `q` with its CROSS TERM deleted;
  and a **half-strip** convex along the ray whose base edge is Q-orthogonal to it, where `q`
  separates into the affine interpolant along the concave base plus the convex part along the ray.
  A negative cross term on the wedge means the envelope is `−inf` — now reported, not answered —
  and a non-orthogonal half-strip is refused loudly rather than approximated.
  **`unboundedFaceTest` 18 / 0**, from 16 / 2. The formulas were checked against the fixtures, not
  fitted to them: they reproduce the `y²` and `−x + y²` those tests derive by hand.

## What is broken

**One red in the whole repository**: `biconjugateOverATwoFaceSubdivisionIsTheEnvelope`, now 5
of 7 points right instead of 0. Slow bucket 114 / 1, fast bucket 204 / 0.

## Needs a decision

**Bug 1's remainder is a refactor, and I stopped rather than force it.** The lens needs the
BOUNDED index layout and cannot get it: it has 2 vertices and 2 genuine edges, but arrives with 5
constraints, and the layout is chosen by that COUNT. `TODO.md` now maps the full indexing
contract — two layouts across four routines, with the loop variable indexing `NCE`'s rows,
`subdE`'s rows and `d.ineqs` at once — and records the two ways to give the lens the right layout
together with the trap in each. Freeing the slots by DROPPING the non-edge constraints was tried
and is **unsound**; it would be sound behind a redundancy PROOF, which `redundantSubset` cannot
supply past a conic.

I would take that refactor next, from the hand-built lens probe rather than from the pipeline.

## Where I stopped

All work is committed and pushed on `main`. `DECISIONS.md` carries the unsound approaches so they
are not re-derived; `TODO.md` and `.claude/SESSION_HANDOFF.md` carry the next concrete step for
each open item.

The transferable lesson, now in the handoff: **build a unit-level reproducer before touching the
symbolic layer.** The hand-built half-lens runs in seconds against a brute-force sup; the same
defect through the pipeline takes 10–40 minutes per run. That is what made bug 1 tractable after
two sessions of failed attempts on it.

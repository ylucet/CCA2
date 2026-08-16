# Morning report — 2026-08-15 overnight run

Branch: **`main`** (pushed as the run went — you authorised pushing for this run)

Task as given: fix the remaining bugs. Standing instructions: parallel agents permitted, do not
stop after three failures, do not wait for input.

## Headline

**BUG 1 IS FIXED, and with it the repository has NO failing test: 327 pass / 0 fail** across all
26 suites — fast 204 / 0, normal 8 / 0, slow **115 / 0**. `biconjugateTest` is 7 / 0, including
`biconjugateOverATwoFaceSubdivisionIsTheEnvelope`, which had been the single red for weeks and was
marked BLOCKED after four attempts.

Confirmed independently: `checkBiconjDomainCoverage` re-measures the two-face box — the row that
read **WRONG** — as **OK with error 0**, against a ground truth that owes nothing to the conjugate
pipeline (the lower convex hull of the sampled graph).

The four earlier attempts all failed for one reason, and stating it is the fix:

> **Both slot conventions identify an edge by a VERTEX INDEX** — "edge `j` is at `ineqs(j)`", or
> "at `ineqs(j+1)`" — and a **LENS's two edges join the SAME pair of vertices**. So neither
> convention can name them apart, and *no* reassignment of slots to constraints ever could. The
> question was never which constraint gets which slot; it is that slots are the wrong addressing
> scheme for this shape.

## What changed

- **BUG 1 — FIXED** (`ddb01c6`, `d06453c`). `functionNDomain.edgeIndexList` derives `eIdx(j)` — the
  constraint bounding the edge from vertex `j` to vertex `j+1` — from the geometry: both endpoints
  on the constraint AND its own curve between them inside the region. The second half is not a
  refinement; on the pipeline's own lens BOTH conics pass through both vertices and the redundant
  one meets the region nowhere else. The three edge-indexed readers take the list as an argument.

  **The refactor was much smaller than the mapped-out plan predicted**, and the reason is worth
  keeping: `getNormalConeEdgeQ` and `getNormalConeEdgeQ3` are the SAME routine under two slot
  conventions — both build the perpendicular to the edge's own constraint at each endpoint,
  oriented by the other endpoint — and collapse into one; `getSubdiffVertexT2` and `T2Q` are
  identical on these inputs. **Three of the "four routines that move together" were two routines
  wearing four names.** And both loops are inside `conjugateOfPiecePoly`, so the list needs no
  field and no signature change on the pieces.

  **Scope:** the new path is entered only when two constraints that each bound a genuine edge still
  share an edge number after `spreadCollidingEdges` — the lens signature, reached by nothing else —
  and where it is entered, today's slot is preferred whenever it is still geometrically valid. No
  constraint is dropped, so the two unsound drops recorded in `DECISIONS.md` stay ruled out.

  **Measured:** both half-lenses conjugate to **3 cells** (2 vertex cones plus the arc; the chord's
  cell is a ray and drops out), exact against a brute-force sup over the lens at **12 of 12** points.
  The old code was `+inf` at `(0,0)` and `(-1,0.5)` where the truth is `0`, and `0` at `(2,-1)`
  where it is `1/2`.

- **Two defects in the new code, found by RE-READING it rather than by a test** (`d06453c`), plus a
  third that was already on the open list:
  1. The two halves of `getNormalConeVertexQ` were handed `eIdx(j)` and `eIdx(j+1)` where the
     routine's own probe points say they mean `eIdx(j-1)` and `eIdx(j)`. Invisible on a lens
     (`nv = 2`, predecessor and successor coincide), wrong for any larger region.
     **When a routine's index convention is unclear, read its PROBE POINTS — they are the geometry
     it is actually using.**
  2. `edgeIndexList` built `nv-1` entries for an unbounded region while its consumer walks all
     `nv` vertices: an index used outside the guard that produced it, for the fourth time in this
     codebase. It now refuses an unbounded region outright (a lens is bounded by construction).
  3. `functionNDomain.getInterior` indexed `c2(2)` under a guard testing only `size(c1,2)` —
     `SUPPORT_MATRIX.md` §7 had this recorded as the next link in Case C's biconjugate chain. The
     two sizes are independent, and any `f` with an `s1²` term but no `s1·s2` term raised
     `badsubscript`.

- **New unit tests** (`functionNDomainTest`, ~10 s where reaching the same piece through
  `biconjugateTest` takes 10–40 minutes): `twoEdgesOnOneVertexPairAreBothKept` pins the precondition
  (the numbering does collide, the geometry does settle it, and one of the two edges is the arc) and
  `halfLensConjugateIsFiniteEverywhereAndExact` pins the conjugate being finite everywhere and exact.

## What is broken

**Nothing that has a test.** Every suite is green: `biconjCPLQTest` 10 / 0, `biconjugateTest` 7 / 0,
`conjCPLQTest` 25 / 0, `testMaxMultiRegion` 24 / 0, `testRegion` 23 / 0, `testcPLQ` 8 / 0,
`unboundedFaceTest` 18 / 0.

Two documented cases still ERROR, and both are **unimplemented paths, not wrong answers** — they
refuse loudly rather than answering:

- **General convex quadrilateral, one face.** Fails in the FIRST conjugation, at
  `plq_1p.conjugateFunction` indexing an EMPTY envelope: `plq_1p.convexEnvelope1` branches on
  `nCE == 0, 1, 2` and falls off the end at `nCE == 3`. CCA2's own `convEnvCPLQ` HAS that case
  ([COAP] A.5, `splitThreeConvex`); it is simply not reachable from `conj`/`biconj`. A concrete
  wiring plan is in `TODO.md`.
- **Parallelogram, one face.** Fails in the SECOND conjugation with
  `QuaParCPLQ:conj:emptyResult` — **and the message names the wrong routine**. It says
  "conjugateOfPiecePoly returned no pieces"; it did not, all **12** pieces conjugate, to **27**
  cells. What returns nothing is `functionNDomain.maxOfList`: folding the groups one at a time, the
  accumulator runs 2, 3, 4, 3, 3, 3, 3, 1, 1 cells and then group 11 empties it. A max's domain is
  the INTERSECTION of the domains, so shrinking is legitimate; empty is not, since that asserts
  `f** = +inf` everywhere. Next attempt starts at group 11, not in Step 2.

## Needs a decision

Nothing is blocked on you. The next item is the quadrilateral wiring above, and `TODO.md` records
the route I would take (build the sub-envelopes with `convEnvCPLQ` and install them as two
`plq_1p` envelope pieces — `conjugate` already loops over them, and `ratPolToPlq` shows the exact
construction) together with the one question it raises: `plq_1p`'s `nCE == 2` branch does not apply
[COAP] A.4's tightness split, which `convEnvCPLQ` does, so routing through `convEnvCPLQ` makes the
two Step 1s agree rather than just filling a hole.

## Where I stopped

All work committed and pushed on `main`; tree clean; no debug instrumentation left. `TODO.md`,
`DECISIONS.md` and `SUPPORT_MATRIX.md` §7 are updated — bug 1's row is struck through with the
measurement, and the four failed attempts are kept against it in `DECISIONS.md` with the one-line
reason none of them could have worked.

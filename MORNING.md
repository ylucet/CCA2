# Morning report — 2026-08-16

Branch: **`main`** (pushed as the run went — you authorised pushing for this run)

Tasks as given: (1) the sub-triangle split that fixes the vendored Step 1's minorant defect and
unblocks the general quadrilateral, and (2) the parallelogram's empty max. Then, separately: write
a test for the quadrilateral, let it fail, and fix it — symbolic first, explicit formulas after.

## Headline

**BOTH ARE FIXED.** The general quadrilateral works, on the fourth attempt, and the method you
prescribed is exactly what made the difference: **the third attempt failed because it computed the
same geometry in double precision.** Doing it symbolically is what kept the numbers small enough
for the pipeline to finish.

- **The quadrilateral:** `f*` of `x·y` over `conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}` is exact at 10 of
  10 probe points, and the fully assembled answer at 8 of 8. Was `MATLAB:badsubscript`. **The fix
  is OPT-IN** (`CCA2_A45_SPLIT`) — see "Needs a decision".
- **The parallelogram:** `emptyResult` is gone, and two general defects went with it. Its
  biconjugate is exact at all four vertices, `+inf` outside the domain, and right at 8 of 10
  interior points.

**One thing to decide, and it is a cost, not a bug** — see "Needs a decision".

## What changed

- **The parallelogram — TWO defects, both FIXED, both general** (`e528943`). Its worst piece went
  from **6 of 10** probe points wrong or uncovered to **2**, against a brute-force sup.

  1. **`region.simplifyUnboundedRegion` declared any region with no finite VERTEX empty** — a
     half-plane, a slab, the whole plane. And **a half-plane is exactly what a TANGENT vertex
     produces**: the cone at a vertex is built from the two edges meeting there, and when those are
     tangent — an arc and its chord touching, which is how a curvilinear piece ends — both
     half-planes are the SAME one. On the offending piece the cone at `(1/4,7/8)` is
     `{2x/3 + y ≥ 4/3}` carrying `x/4 + 7y/8 − 1/2`, and it was being deleted. Now refuted by a
     WITNESS: a feasible point certifies non-emptiness, failing to find one proves nothing, so the
     old verdict stands whenever no witness turns up.
  2. **The edge list, in bug 1's other form.** The piece is bounded with 3 vertices and 3 genuine
     edges plus a conic touching one vertex: 4 constraints for 3 vertices, so `size(ineqs,2) == nv`
     called a bounded region unbounded, `endNv` came out `nv−1`, and it was built one edge cell
     short. The edge list now covers any bounded piece the count mislabels, not only a lens.

  Pinned by `regionTest/aHalfPlaneIsNotEmpty` (with its converse, so a genuinely empty region is
  still reported empty) and
  `functionNDomainTest/aBoundedPieceWithATangentVertexConjugatesOntoTheWholePlane`.

- **The general quadrilateral — FIXED** (`eb7d11d`, `309b0f7`), after a red test written first
  (`6505077`). `splitTightTriangleSym` splits a triangle into sub-triangles on each of which
  cPLQ's own closed form for THAT sub-triangle's convex-edge count IS the convex envelope, and
  `triangulate` emits them as pieces — which keeps every one of them on a path Step 2 already has.
  Two tests: one asserting only what must hold of an envelope (it exists, it is `≤ x·y`, and where
  `x·y ≥ 0` it is `≥ 0`), one on `f*` against the vertex-attained sup.

  **Why symbolic mattered, since it is the whole difference from the attempt that hung:** A.4's
  cevian has slope `−sqrt(mh·mw)`, so its foot is irrational. Taking it from the existing
  double-precision routine gives `2^53` denominators that grow past `1e25` downstream, and the
  symbolic engine then cannot decide anything. Carried symbolically the coordinates stay compact
  surds — `5/2 − sqrt(5)/2`, `3/2 − 3·sqrt(5)/10`.

  **Two results fell out of writing it, both now in the source.** The `±` branch of A.4's form is
  vestigial: substituting the edge into either sign gives the same three coefficients, so the
  touching test cannot separate them and `+1` always wins. And the "no split needed" test is a
  perfect square, `(sqrt(mh·mw)·dx + dy)²/(sqrt(mh)+sqrt(mw))²`, so it vanishes exactly when the
  weak edge is PARALLEL to the cevian direction — no `q1` and no branch selection needed.

- **The A.4/A.5 domain split taken from double precision — built, measured, REVERTED**
  (`b0a36de`). The attempt the above replaced. See below.

## What is broken

Nothing that has a test. **330 pass / 0 fail across all 26 suites** — fast 206 / 0, normal 9 / 0,
slow **115 / 0**, every suite green on the committed tree.

Two things are known-imperfect and documented rather than hidden:

- **The parallelogram is 4% LOW at 2 of 10 interior points** — `0.986` for `1.031`, `1.913` for
  `1.950`. Cause found: on that piece `f` is a SINGULAR convex quadratic (constant along one whole
  edge, `∇f = 0` at two of three vertices), so `getInterior` — which separates an edge cell from
  its neighbours by eliminating `s` between `x = ∂₁f` and `y = ∂₂f` — returns the gradient map's
  image LINE rather than a separating curve, two edge cells come out on the same region, and the
  wrong one is checked first. **Do not attack the `isQuad` chord rewrite for it**: chording the
  vertices the conic actually touches makes the piece WORSE (2 wrong of 10 → 3), and skipping the
  rewrite changes nothing at all. Both measured today.
- **The general quadrilateral still raises `MATLAB:badsubscript` BY DEFAULT.** The fix is written
  and tested; it is opt-in. See "Needs a decision".

## Needs a decision

**The quadrilateral fix is correct but OPT-IN (`CCA2_A45_SPLIT`), and I left it off by default.**
That is the one judgement call in this run, and here is what it rests on.

Turning it on costs two things, both measured:

- **`testcPLQ` goes from 1542 s to 4728 s** (its historical time is 1427 s). A.4's cevian foot is
  IRRATIONAL, so a split sub-triangle has SURD coordinates and every symbolic operation downstream
  works in a quadratic extension instead of the rationals. Only two of that suite's six domains even
  gain a piece (2 → 3 and 1 → 2), so this is the algebraic degree of the coordinates, not the piece
  count. Nothing else moves — the fast bucket is unchanged at 206 / 0, because the no-split path
  costs 20 ms.
- **`testcPLQ/testRectBiconj` then ERRORS.** That test has no assertions; it runs the pipeline, so
  the failure is an exception. Undiagnosed.

So switching it on trades a documented, LOUD failure on one domain shape for a new one on another,
and I did not think that was mine to decide silently. **Say the word and I will flip the default**,
or fix the `testRectBiconj` error first — that is the smaller of the two remaining jobs.

**The larger one is Step 3, now the top `TODO.md` item with its numbers.** Assembling the
cross-piece maximum for the quadrilateral takes 73 minutes while Steps 1 and 2 take 25 seconds, and
the cell count runs 5, 14, 29, 45, 70, **86**. Eighty-six is about ten times what the answer needs —
`f*` of `x·y` over a convex quadrilateral has a cone per vertex and a cell per edge — and the
surplus is adjacent cells carrying the same function that `region.merge` never merges. Fix that and
the split becomes affordable everywhere.

## Where I stopped

All work committed and pushed on `main`; tree clean; no debug instrumentation left. `TODO.md`,
`DECISIONS.md` and `SUPPORT_MATRIX.md` §7 / §7.1 are updated — the parallelogram's row is struck
through with the measurement, the quadrilateral's carries all three failed attempts and what not to
re-try, and the residual `getInterior` defect has its own row and its own one-minute reproducer.

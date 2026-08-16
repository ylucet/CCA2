# Morning report — 2026-08-16

Branch: **`main`** (pushed as the run went — you authorised pushing for this run)

Task as given: do (1) the sub-triangle split that fixes the vendored Step 1's minorant defect and
unblocks the general quadrilateral, and (2) the parallelogram's empty max. Do not wait for input.

## Headline

**(2) is FIXED — the parallelogram's `emptyResult` is gone, and two general defects went with it.**
Its biconjugate now computes: exact at all four vertices, `+inf` outside the domain, and 8 of 10
interior probe points correct against a brute-force double conjugate.

**(1) is NOT fixed, and the third attempt is what finally names the blocker: ARITHMETIC, not
structure.** The split was built exactly as the previous write-up prescribed, it worked at Step 1,
and it turned the quadrilateral's crash into a **HANG**. Reverted.

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

- **The A.4/A.5 domain split — built, measured, REVERTED** (`b0a36de`). See below.

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
- **The general quadrilateral still raises `MATLAB:badsubscript`**, as before. Loudly.

## Needs a decision

Nothing is blocked on you. The quadrilateral's remaining work is now one line — implement [COAP]
A.4's cevian and A.5's smooth-fit line SYMBOLICALLY and have `triangulate` emit the sub-triangles
as pieces — and `TODO.md` leads with it.

Why the third attempt failed, since it is the useful part: the split was taken from `convEnvCPLQ`'s
own faces, which is double precision, and `sym` of a double is EXACT — a denominator near `2^53`.
Snapping the new vertices to the simplest rational within `1e-10` bounds the VERTEX denominators
but not the downstream ones, because the conjugate is a rational function of those coordinates: a
few squarings carry `1e5` to `1e25`, and MuPAD's `isAlways` then decides nothing. The conjugate ran
45+ minutes with no output behind 3.8 MB of `TruthUnknown` warnings. A.4's cevian has slope exactly
`−sqrt(mh·mw)`, so its foot is an exact algebraic number and `sqrt` is something the symbolic layer
keeps small — that is the way in.

## Where I stopped

All work committed and pushed on `main`; tree clean; no debug instrumentation left. `TODO.md`,
`DECISIONS.md` and `SUPPORT_MATRIX.md` §7 / §7.1 are updated — the parallelogram's row is struck
through with the measurement, the quadrilateral's carries all three failed attempts and what not to
re-try, and the residual `getInterior` defect has its own row and its own one-minute reproducer.

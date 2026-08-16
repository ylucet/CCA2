# Session Handoff

_Last updated: 2026-08-16_

## What happened this session

**The parallelogram's `QuaParCPLQ:conj:emptyResult` is fixed**, and the two defects behind it are
general, not special cases. **The general quadrilateral is not fixed**, and the third attempt at it
is the one that finally names the blocker: the split has to be carried SYMBOLICALLY, because taking
it from `convEnvCPLQ`'s double-precision faces makes the exact symbolic arithmetic downstream
explode.

## Where things stand

- Branch: `main`, pushed, tree clean, no debug instrumentation.
- **330 pass / 0 fail across all 26 suites** — fast 206 / 0, normal 9 / 0, slow **115 / 0**.
- The parallelogram's `biconj` now computes: exact at all four vertices, `+inf` outside the domain,
  8 of 10 interior probe points right against a brute-force double conjugate. The other two are
  about 4% LOW — the residual defect below.
- The general quadrilateral still raises `MATLAB:badsubscript`, loudly, as before.

## Next steps

1. **`getInterior` on a SINGULAR quadratic — the parallelogram's last 4%.** `getInterior` separates
   an edge cell from its neighbours by eliminating `s` between `x = ∂₁f` and `y = ∂₂f`. When `f` is
   a singular convex quadratic the gradient map is not invertible and that elimination returns the
   map's IMAGE LINE, which separates nothing — so two edge cells come out on the SAME region and
   the first one checked wins. Reproduce in about a minute with
   `functionNDomainTest.parallelogramPiece9`.
   **Do NOT attack the `isQuad` chord rewrite for it.** Both alternatives were measured 2026-08-16
   and are in `DECISIONS.md`: chording the vertices the conic actually touches makes that piece
   WORSE (2 wrong of 10 → 3), and skipping the rewrite changes nothing.
2. **The general quadrilateral — implement [COAP] A.4's cevian and A.5's smooth-fit line
   SYMBOLICALLY**, and have `plq_1p.triangulate` emit the sub-triangles as PIECES. The cevian's
   slope is exactly `−sqrt(mh·mw)`, so its foot is an exact algebraic number and `sqrt` is something
   the symbolic layer keeps small. Bounded work — two line intersections and a curvature test — but
   new code, not wiring. **Three attempts are recorded in `DECISIONS.md` with what not to re-try**;
   the shortest version is that envelope-level splitting cannot work (Step 2 has no
   rational-envelope branch) and domain-level splitting taken from `convEnvCPLQ` hangs (irrational
   cevian foot → `1e25` denominators → `isAlways` decides nothing).
3. **Then SCIP/QPLIB**, in the order that bites: wire `biconj` into `SCIP/src/cca2ConvexEnvelope.m`
   → expose value+subgradient off `QuaParCPLQ` → fix diagonal terms over a box (`x²−y²`,
   `(x²+y²)/2` on `[0,1]²` still error in the second conjugation) → performance (~40–60 s/term).

## Three methods that keep paying off

**Build a unit-level reproducer before touching the symbolic layer.** Piece 9 of the
parallelogram's `f*`, hand-built as a `region`, runs in about a minute; reaching the same piece
through `biconj` takes far longer and buries the evidence.

**When a routine's index convention is unclear, read its PROBE POINTS.** `getNormalConeVertexQ`'s
two halves probe at vertex `j−1` and vertex `j+1`, which is what settles that they mean the edge
ARRIVING at `j` and the edge LEAVING it.

**To find which piece of a MAX is wrong, use the property the max must have.** `f**` of a bounded
domain is finite exactly ON that domain and is a max, so EVERY per-piece conjugate must be finite
there. Evaluating all 12 groups at six interior points turned "the max comes out empty" into "these
three pieces are wrong" in one cheap run.

## Relevant files

- `DECISIONS.md` — dead ends and refuted reasoning; the newest entries are the parallelogram's two
  defects and the three failed quadrilateral attempts.
- `TODO.md` — live action items; the quadrilateral leads, with its one-line remaining task.
- `SUPPORT_MATRIX.md` §7 / §7.1 — the defect table and the domain-coverage measurements.
- `region.m` — `simplifyUnboundedRegion` (the witness), `getNormalConeVertexQ`, `getNormalConeEdgeQE`.
- `functionNDomain.m` — `edgeIndexList`, `conjugateOfPiecePoly`, `getInterior`.
- `functionNDomainTest.m` / `regionTest.m` — the unit tests that pin all of this without the pipeline.

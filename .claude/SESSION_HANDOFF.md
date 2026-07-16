# Session Handoff

_Last updated: 2026-07-15T01:00:00Z_

## What happened this session

**Part 1 (complete, verified): fixed the prior session's #1-priority envelope-tightness gap** in
`convEnvCPLQ.m` (COAP Appendix A.4, 2-convex-edge triangle case), diagnosed but not fixed last
time. Derived and implemented the genuine sub-partition that diagnosis anticipated: the triangle
is split by a cevian from one weak-edge endpoint into the opposite convex edge, giving two
sub-triangles -- the one containing the two convex edges' shared vertex keeps the original
`twoEdgeQuadPlain` quadratic unchanged (still exactly tight there); the other uses a new quadratic
(`buildEdgeAffinePiece`) that touches `f=u1u2` along its remaining convex edge and the AFFINE
CHORD (not `f` itself) along the now-fully-contained weak edge. The cevian's direction is forced
(found via the same "factor the difference of two same-touching quadratics" trick as
`splitThreeConvex`), and exactly one of its two candidate directions is geometrically valid,
verified across ~6500 random triangles. A separate, measure-zero case (e.g. the mirror-symmetric
triangle `(0,0),(2,1),(1,2)`) needs no split at all and is detected directly and cheaply.

Verified via: the paper's own Appendix A.4.3 example (previously wrong at the flagged bad dual
point and at the "dip" point `x0=(0.474343,0)`, both now correct); a new regression test
(`convEnvCPLQTest.bilinearTwoConvexEdgesSplitIsTight`); and a ~60-triangle MATLAB-side randomized
stress test comparing each split's two sub-triangle conjugates (`conjPieceCPLQ`) against exact
closed-form ground truth (`supBilinearOverPoly`) -- max error `5e-14` (machine precision) across
every triangle that produced a genuine split. **This part stands on its own and is not in
question**, regardless of Part 2 below.

While implementing, found and fixed a bug in my own first MATLAB port of the derivation (an
algebra transcription error in the analytic null-space formula inside `buildEdgeAffinePiece`,
caught by a randomized stress test producing wrong `convEnvCPLQ:internal` errors and one
PSD-classification crash before the fix).

**Part 2 (investigated at length, PARTIALLY fixed, not resolved): the downstream `maxQuaPar` gap**
Fixing Step 1 means `conjCPLQ`'s existing Step-3 fallback (`conjMaxOfSubTriangles`, already used
for the 3-convex-edge case) now also gets exercised for the 2-convex-edge case, and `maxQuaPar`
(Step 3) frequently -- ~52% of random split cases in a 138-triangle stress test -- errors
combining the two sub-triangles' conjugates (`maxQuaPar:internal`, an unmatched boundary
half-edge). Spent this session's remaining time digging into WHY:

- **First, confirmed this is not a correctness question at all.** `max(g1,g2)` is *always*
  exactly `f*` for any domain partition (`sup` over a union is the max of the `sup`s over the
  pieces -- a trivial identity, true whether the seam is smooth or kinked). Verified directly:
  bypassing `maxQuaPar` entirely and taking a plain `max` of the two sub-triangles' own
  `conjPieceCPLQ` conjugates matches ground truth to `5e-14` in every case, including the ~52%
  where `maxQuaPar`'s ASSEMBLY throws. So the bug is entirely in `maxQuaPar` building a
  self-consistent `QuaPar` object, never in any value it manages to produce (0 mismatches among
  the cases that *did* succeed).
- **Root cause, as far as diagnosed**: `splitThreeConvex`'s (3-convex-edge case) sub-triangles are
  built so their conjugates paste together SMOOTHLY (C1) along the shared seam -- confirmed
  `maxQuaPar` has only ever been exercised against that case. `splitTwoConvexEdges`'s new cevian
  (this session's Part 1 fix) only guarantees the two PRIMAL pieces agree in VALUE along the seam,
  not gradient -- confirmed by direct computation (gradients of the two pieces match at one seam
  endpoint, diverge increasingly toward the other -- a real, continuously-varying kink). A primal
  kink gives the two pieces' conjugates `g1,g2` a genuine POSITIVE-AREA tie region in dual space
  (`g1=g2` over a whole 2D set, not just a curve) near the dual image of the shared seam vertex --
  expected because both `g1` and `g2` independently produce their own vertex-cone face for that
  shared point, and unlike the always-smooth case, these can overlap or sit adjacent to several of
  the other side's faces (a "fan"). `maxQuaPar`'s per-`(k,l)`-pair loop isn't built to expect this.
- **Implemented `dedupPieces`** (new function in `maxQuaPar.m`, called right before
  `assemblePieces`): detects groups of pieces with IDENTICAL geometry produced by different
  `(k,l)` pairs (one concrete symptom of the tie phenomenon) and collapses each group to one,
  reconciling disagreeing winners by evaluating every candidate row at an interior point and
  keeping the largest (sound, since both `g1`'s and `g2`'s rows are independently exact
  everywhere). Verified correct on the paper's own example (where two `(k,l)` pairs produce an
  identical cell with a wrong vs. right winner) and causes zero regressions (146/146 still pass).
  **However, it empirically never independently fixes a single failing triangle** (checked with
  and without it across 138 random split triangles: identical pass/fail outcome, `dedupPieces`
  literally never changes which triangles succeed) -- every triangle hitting the exact-duplicate
  pattern ALSO hits a second, more complex pattern (several small cone/strip faces near the kink
  vertex whose RAYS fail to pair with any counterpart, not exact duplicates) that isn't understood
  well enough yet to fix. Kept anyway: it's correct, safe, and necessary infrastructure even if not
  sufficient on its own.
- Fails LOUDLY and cleanly throughout (never a silent wrong answer) -- this is a real, fairly
  common gap, but not a blocker for Part 1's correctness fix, which stands independently verified.

## Where things stand

- Branch: `main` @ `719943d` -- "Fix convEnvCPLQ envelope-tightness bug; partial fix for maxQuaPar
  assembly gap".
- Pushed: yes.
- Files changed: `convEnvCPLQ.m` (Part 1 fix), `convEnvCPLQTest.m` (+1 regression test + 2 static
  helpers), `conjPieceCPLQTest.m` (1 test updated for the new 2-face Step-1 output), `maxQuaPar.m`
  (`dedupPieces`, new -- Part 2 partial fix), `DESIGN.md` (documented both parts). No other files
  touched.
- Full test suite: **146/146 PASS**.

## Next steps

- **Top priority now**: finish diagnosing and fix the residual `maxQuaPar` gap (Part 2 above).
  `dedupPieces` handles the exact-duplicate flavor of the "positive-area tie region" phenomenon;
  the remaining ~52% failure rate is a second, adjacent-but-not-identical flavor of the same root
  phenomenon (small cone/strip faces near the kink vertex's dual image, unmatched rays) that needs
  its own diagnosis. Recommended direction: give the shared seam-vertex's dual-side "fan" of faces
  first-class handling BEFORE the general `(k,l)` double loop (detect, for any vertex `g1`/`g2`
  inherit from a shared primal point, whether their own faces there overlap/complement, and
  resolve before the generic clip loop runs) -- likely needs vertex-PROVENANCE tracking (tying
  each dual face back to the specific primal vertex/edge it came from), a tool this file's own
  HISTORY comments already anticipated for a different, now-resolved ambiguity.
  - Repro triangles (regenerate via `rand(3,2)*10-5`, keep ones where `classifyConvexEdges` gives
    exactly 2 rows, and `convEnvCPLQ(...).nf==2`): the paper's own `V=[2 1;0 0;1 0]`, and
    `V=[1.508518 2.818371; 2.687354 4.499057; -1.870671 3.524095]`. Both throw
    `maxQuaPar:internal` ("piece 1" always has the specific unmatched edge in these traces),
    *even with `dedupPieces` active*.
  - Ground-truth check to validate any fix: compare `q.conj('cplq').eval(s)` against
    `convEnvCPLQTest.supBilinearOverPoly(s, V)` (exact closed form, already in the codebase) at
    several random `s` -- the same pattern `bilinearTwoConvexEdgesSplitIsTight` already uses at the
    per-piece (not end-to-end) level.
  - Once fixed: strengthen `bilinearTwoConvexEdgesSplitIsTight` (or add a new test) to call
    `q.conj('cplq')` end-to-end instead of manually taking `max` of the two per-piece conjugates,
    and add a dedicated `maxQuaParTest`/`conjCPLQTest` regression using one of the repro triangles.
  - This session's diagnostic scratch files (`stress2edge.m`, `widersearch.m`, `checkg2.m`,
    `maxQuaPar_backup.m`) are in the session scratchpad, NOT committed -- likely gone next session;
    easy to regenerate from the descriptions above.
- Related, still unconfirmed whether same root cause: `QuaPar.orderEdges`/`createP` rejects the
  face topology for some near-degenerate/thin triangles with a DIFFERENT error ("Face k has ..." /
  "expected 2 but got N").
- Lower priority / untouched: `partialConj` for `'cplq'`/`'pqp'`; `add` for `RatPol` and the
  `RatPar` parent class; conjugate engines `'pqp'`/`'graph'` (unimplemented, error explicitly); the
  standalone `RatPol.conj` gap.

## Relevant files

- `convEnvCPLQ.m` -- `splitTwoConvexEdges` (new), `buildEdgeAffinePiece` (new), `seamPoint` (new):
  Part 1's fix. `envelopeFromClassified`/`twoEdgeQuadPlain`/`buildTwoEdge` themselves are UNCHANGED
  (still used as-is for the 0/1-convex-edge cases and, deliberately unmodified, inside
  `splitThreeConvex`'s 3-convex-edge sub-triangle loop -- that path was already empirically
  correct per its own extensive `maxQuaParTest` coverage, so it was left alone to avoid regressing
  it).
- `maxQuaPar.m` -- `dedupPieces` (new), called from the main `maxQuaPar` function right before
  `assemblePieces`. Part 2's partial fix; see its own header and this file's HISTORY-style
  commentary for the full rationale. The residual gap (adjacent-but-not-identical tie cells) is
  still inside this file, not yet isolated to a specific function.
- `DESIGN.md` -- the `maxQuaPar.m` Implementation-status bullet now has the full writeup for both
  Part 1 (the fix + verification) and Part 2 (root-cause diagnosis, the partial fix, and what's
  still open), appended after last session's diagnosis.
- `convEnvCPLQTest.m` -- new test `bilinearTwoConvexEdgesSplitIsTight` (the paper's own example,
  checked for tightness via each sub-triangle's own conjugate, not just primal value -- NOT via
  `q.conj('cplq')` end-to-end, since that still throws for this triangle per Part 2) plus two new
  static helpers (`extractTriFace`, `supBilinearOverPoly`) reused from `maxQuaParTest.m`'s pattern.
- `conjPieceCPLQTest.m` -- `psdRank1QuadraticEndToEnd` updated to extract one face from Step 1's
  now-2-face output before feeding it to Step 2 (documented at the call site).
- `conjCPLQ.m` -- NOT modified this session; its existing `conjMaxOfSubTriangles`/Step-3 fallback
  already generically handles a multi-face Step-1 envelope (originally written for the 3-edge
  case) and needed no changes to correctly dispatch the 2-edge case too -- the remaining gap is
  entirely inside `maxQuaPar.m`.

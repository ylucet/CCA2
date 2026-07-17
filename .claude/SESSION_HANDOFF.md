# Session Handoff

_Last updated: 2026-07-16T00:00:00Z_

## What happened this session

Continued the prior session's #1-priority item: the `maxQuaPar` assembly gap that makes
`conjCPLQ`/`maxQuaPar` throw `maxQuaPar:internal` on roughly half of random 2-convex-edge triangle
splits (the case fixed structurally, but not yet fully wired end-to-end, last session). **Found and
fixed two more concrete bugs within the "vertex fan" phenomenon diagnosed last session; both are
real, verified, zero-regression fixes, but neither moves the observed aggregate crash rate** --
established via a properly-controlled experiment, described below. The underlying gap remains open.

Diagnostic method: reproduced the paper's own `V=[2 1;0 0;1 0]` example end to end, then
temporarily instrumented `assemblePieces` (via `assignin('base',...)`, removed before every
committed state) to dump the internal `pieces`/half-edge list right before the crash, and directly
tested pairs of pieces for literal point-set overlap via a coverage/sampling check (not just
trusting the crash message) -- this is what surfaced both bugs below as genuine geometric defects,
not just symptoms.

1. **`matchHalfEdges`' ray matching accepted a pairing on `(apex, direction)` equality alone, with
   no check that the two candidate pieces are on OPPOSITE sides of that ray.** When a `g1` (or `g2`)
   face is cut into several `(k,l)` sub-pieces by different opposing faces, every sub-piece
   independently inherits the parent face's own boundary rays -- not because each is adjacent to
   some other sub-piece there, but simply because they all descend from the same ancestor. The old
   code paired up the first two such inheritors it found as if mutual neighbours; confirmed via
   direct point-sampling that two of this triangle's pieces genuinely overlapped over a
   positive-area region, not merely touched a shared boundary. **Fixed** by `oppositeSides` (new, in
   `maxQuaPar.m`): tests each ray candidate's own adjacent geometry via a 2D cross product against
   the shared ray direction, accepting a pairing only when the two candidates fall on opposite
   sides.
2. **A related but distinct redundancy**: two different `(k,l)` pairs sharing the identical winning
   row `f` can produce pieces that are not exact geometric duplicates (`dedupPieces` only collapses
   those) yet still overlap, one wholly containing the other. **Fixed** by `dropSubsumedPieces` (new
   in `maxQuaPar.m`, called right after `dedupPieces`): drops whichever piece's entire geometry
   (vertices, and recession directions if unbounded) is contained in another piece agreeing on `f`.
   Confirmed: piece count for the same triangle drops from 12 to 10 with zero change in any resolved
   value.

**Both fixes are individually correct** (each addresses a genuine, hand-confirmed geometric defect)
**and cause zero regressions** (full suite 147/147: 146 prior + 1 new regression test). **Important,
initially counter-intuitive finding, established carefully after an early false positive**: neither
fix changes the pass/fail outcome for ANY triangle in a large random sample. A same-triangle-set A/B
comparison -- pre-generating every triangle candidate AND every dual sample point up front so the
"with fixes" and "without fixes" runs are byte-for-byte reproducible regardless of which specific
triangles happen to succeed -- across 1494 randomly generated 2-convex-edge triangles found **zero**
outcome flips in either direction (931/1494, 62.3%, succeed with or without the fixes). (An earlier,
naive attempt at this same comparison drew random numbers lazily inside the success/failure
branches, which silently desynchronised the two runs' RNG streams the moment they first behaved
differently -- producing a spurious-looking "52% -> 36%" improvement that evaporated once the
methodology was corrected. Recorded here so it isn't mistaken for a real result if seen in old
scratch output.) So: both fixes are real and worth keeping (the overlap bug especially -- it is a
genuine correctness hazard, not just a crash-avoidance nicety, even though it happened not to flip
any outcome in this sample), but the ~38% aggregate crash rate on generic random triangles is driven
by a THIRD, still-undiagnosed pattern within the same "fan" phenomenon.

**One dead end recorded so it isn't retried**: broadening `insertPassthroughVertices`' candidate
points from just the current cell's own two parent faces to EVERY vertex of both full inputs
`g1`/`g2` (reasoning: a missing "T-junction" split could in principle come from a third, unrelated
face) neither fixed the repro triangle nor helped the aggregate rate, and caused 7 new regressions
in the existing suite -- reverted immediately, not present in the final diff.

## Where things stand

- Branch: `main` @ `8828375` -- "Fix two maxQuaPar assembly bugs (ray same-side pairing, subsumed
  pieces)".
- Pushed: pending (see end-of-session summary for the actual outcome).
- Files changed: `maxQuaPar.m` (`oppositeSides` + `raySideVector`, called from `matchHalfEdges`;
  `dropSubsumedPieces` + `isSubsumed`, called from the main `maxQuaPar` function right after
  `dedupPieces`), `maxQuaParTest.m` (+1 regression test,
  `matchHalfEdgesRejectsSameSideRayPairingAndDropsSubsumedPieces`, pinning the CURRENT documented
  failure mode -- `maxQuaPar:internal` on the paper's own example -- as a concrete target for
  finishing the fix), `DESIGN.md` (full write-up appended after last session's diagnosis). No other
  files touched.
- Full test suite: **147/147 PASS** (146 prior + 1 new).
- The `maxQuaPar:internal` crash rate on generic random 2-convex-edge triangles is **unchanged by
  this session's fixes**: still ~38% (measured on this session's own 1494-triangle sample; the
  ~52% figure from the immediately-prior session's differently-generated 138-triangle sample is not
  directly comparable methodology, so don't read a trend into the two numbers differing).

## Next steps

- **Top priority, unchanged in substance from last session, now with two known-NOT-it patterns ruled
  out**: the dominant remaining cause of `maxQuaPar:internal` on 2-convex-edge splits is still a
  third, undiagnosed pattern -- confirmed NOT the same-side ray mismatch or the subsumed-piece
  redundancy this session fixed (both present and fixed on the repro triangle, which still fails).
  Diagnostic groundwork now in place for whoever continues:
  - Reproduce via `maxQuaParTest.buildG1G2ForTriangle(V)` (already in the test file) + `maxQuaPar`,
    e.g. `V=[2 1;0 0;1 0]` (the paper's own example; verified it still throws after both this
    session's fixes) or `V=[1.508518 2.818371; 2.687354 4.499057; -1.870671 3.524095]` (from the
    prior session, not re-checked this session but likely still failing too -- same underlying
    pattern).
  - To inspect the internal `pieces`/half-edge state right before the crash (this session's own
    diagnostic method, not committed): temporarily add
    `assignin('base','DBG_pieces',pieces); assignin('base','DBG_HE',HE); assignin('base','DBG_opp',opp); assignin('base','DBG_rootOf',rootOf);`
    right after the `[V, rootOf] = buildGlobalVertices(...)` line in `assemblePieces` (inside
    `maxQuaPar.m`), call `maxQuaPar` from a **script** (not a function -- so `assignin('base',...)`
    lands in the same workspace you inspect from) wrapped in try/catch, then dump `DBG_pieces`'
    vertices/dirIn/dirOut/f and `DBG_HE`'s `opp` column to see exactly which half-edges are
    unmatched and why. Remove the `assignin` line before committing anything.
  - A useful independent check for any candidate fix: a direct point-sampling coverage test (sample
    a grid or random points across the region near the unmatched pieces, check via half-plane
    membership tests how many of the CURRENT `pieces` claim each point -- 0 means a genuine gap,
    ≥2 means a genuine overlap) is more reliable than reasoning about the geometry by hand; this is
    what surfaced both of this session's bugs and is worth writing fresh each time rather than
    trusting intuition about which pieces "should" be adjacent.
  - For the specific repro triangle, the remaining unmatched pieces (after this session's fixes)
    are a small cluster near the dual image of primal vertex `A=(0,0)` (5-6 tiny pieces that don't
    close up combinatorially) PLUS one larger, seemingly unrelated orphaned edge on the big
    quadrilateral piece -- worth checking whether these are really the same root cause or two
    separate remaining bugs.
  - Recommended direction, still unchanged from last session: vertex-provenance tracking (tying
    each dual face/edge back to the specific primal vertex/edge that produced it), to give the
    shared seam-vertex's dual-side "fan" first-class handling before the generic `(k,l)` clip loop
    runs, rather than trying to patch the generic pairwise-matching machinery further -- this
    session's two fixes were exactly that kind of generic patch, and while both were genuine,
    correct bugs, neither touched the dominant failure pattern, suggesting further generic patches
    are unlikely to be a good use of time versus the structural fix.
  - Ground-truth check for any fix: `maxQuaParTest.supBilinearOverPoly(s, V)` (exact closed form,
    already in the test file) vs. `maxQuaPar(g1,g2).eval(s)` at several random `s`.
  - Once genuinely fixed: strengthen or replace this session's
    `matchHalfEdgesRejectsSameSideRayPairingAndDropsSubsumedPieces` test (which currently only
    documents the known failure via `verifyError`) with a real end-to-end value check against
    `supBilinearOverPoly`, and consider re-running this session's 1494-triangle A/B methodology
    (script not committed, easy to regenerate from the description above) to get a real
    before/after aggregate number for whatever the actual fix turns out to be.
- Related, still unconfirmed whether same root cause: `QuaPar.orderEdges`/`createP` rejects the face
  topology for some near-degenerate/thin triangles with a DIFFERENT error ("Face k has ..." /
  "expected 2 but got N").
- Lower priority / untouched: `partialConj` for `'cplq'`/`'pqp'`; `add` for `RatPol` and the `RatPar`
  parent class; conjugate engines `'pqp'`/`'graph'` (unimplemented, error explicitly); the
  standalone `RatPol.conj` gap.

## Relevant files

- `maxQuaPar.m` -- `oppositeSides`/`raySideVector` (new, called from `matchHalfEdges`) and
  `dropSubsumedPieces`/`isSubsumed` (new, called from the main `maxQuaPar` function right after
  `dedupPieces`): this session's two fixes. Both have full header comments with the diagnosis. The
  residual gap is still somewhere in this file's `(k,l)`-pair generation or matching machinery, not
  yet isolated to a specific function.
- `maxQuaParTest.m` -- new test
  `matchHalfEdgesRejectsSameSideRayPairingAndDropsSubsumedPieces`, documenting both fixes and
  pinning the current (still-failing) status of the paper's own repro triangle.
- `DESIGN.md` -- the `maxQuaPar.m` Implementation-status bullet now has this session's full write-up
  (root causes, fixes, and the controlled A/B methodology/result) appended after the prior two
  sessions' diagnosis.
- `convEnvCPLQ.m` -- NOT modified this session; last session's Part 1 fix (envelope tightness) is
  unaffected and still independently verified.
- `conjCPLQ.m` -- NOT modified this session.

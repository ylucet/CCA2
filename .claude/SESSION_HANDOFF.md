# Session Handoff

_Last updated: 2026-07-18T02:15:00Z_

## What happened this session

Continued the Part 2 (2-convex-edge tightness) investigation. Found and fixed the actual root
cause (Part 2c): the correct cevian is where `q1` and the neighboring Appendix A.3 piece are
C1-TANGENT (a perfect-square factor when clearing denominators, not just equal in value) along the
line through the anchor vertex with slope `-sqrt(mh*mw)` -- provably exact, independent of `beta`
and the affine shift (both cancel out of the derivation completely). Replaces the session's earlier,
still-imperfect `edgeClipCevian` criterion. Verified against ground truth; stress-test gap rate
dropped from ~45% (session start) to 23.1% (first fix) to 9.3% (this fix).

Diagnosing that remaining 9.3% found the SAME tightness issue also affects the separate `nCE==3`
(three-convex-edge, `splitThreeConvex`) code path, which never checked it. Generalized
`convEnvCPLQ.m` to apply the check recursively wherever a 2-convex-edge sub-triangle appears. This
closes the gap in the one repro checked, but has a much bigger footprint than intended: of 9
"3-convex-edge" triangles sampled (a small, non-random sample -- see DESIGN.md's hedge), all 9 now
produce 4 pieces (2 quadratic + 2 rational) instead of 2, and `conjPieceCPLQ` can't yet conjugate a
rational piece. This breaks 6 existing tests whose actual purpose is regression-testing
`maxQuaPar`'s own assembly logic (not Step 1's tightness), via a shared helper
(`maxQuaParTest.buildG1G2ForTriangle`) that can no longer build test fixtures the way it used to.
**Committed as WIP, undecided** -- see Next steps.

## Where things stand

- Branch: `main` @ `e6e936a` -- "WIP: convEnvCPLQ 2CE split -- exact tangency cevian, recursive 3CE
  generalization (test-breaking, undecided)".
- Pushed: no (not asked this session; see Next steps).
- Full MATLAB suite: **73/79 passing** (6 failures, all attributable to the `nCE==3`
  generalization described above -- `conjCPLQTest.indefiniteTriangleThreeConvexEdgesUsesStep3` and
  5 `maxQuaParTest` cases: `splitCellAcceptsGenuineNonDegenerateParabola`,
  `dedupHitsMergesCrossingsAtACellCorner`, `assemblePiecesResolvesNearDuplicateApexCluster`,
  `checkOrphanHalfEdgesDropsProvablyDegenerateOrphanEdges`,
  `insertPassthroughVerticesDropsNearDuplicateCrossingPoint`).
- Untracked, not part of this repo: `cPLQ/`; not touched.

## Next steps

- **Decide, deliberately, before doing anything else**: keep the `nCE==3` generalization (and
  follow through on fixing the 6 broken tests -- likely by freezing each test's current `g1,g2` as
  hardcoded `QuaPar` fixtures computed from the PRE-session code, decoupling them from the live
  Step1+Step2 pipeline, since they test `maxQuaPar`'s assembly logic, not Step 1's tightness), OR
  scope it back out (keep only `tangentCevian` replacing `edgeClipCevian` for the standalone
  `nCE==2` path; revert `nCE==2`'s own sub-triangle loop and the `nCE==3` loop to calling
  `envelopeFromClassified` directly as before; remove `solveTriangleBF`/`assemblePiecesBF`). Neither
  path has been done. See DESIGN.md's `convEnvCPLQ.m` entry, "Part 2c", for the full derivation,
  what's confirmed vs. hedged, and more detail on both options.
- Whichever way that goes, the full MATLAB suite must be back to a passing count before
  considering this session's work "done" -- do not push with 6 known failures.
- Fill in the exact `[LOCATELLI]` citation in `DESIGN.md`'s reference list -- still not done.
- Still open, unchanged from before: the 2/741 (0.3%) residual `maxQuaPar:internal` crashes (see
  DESIGN.md history for repro triangles) -- not investigated this session.
- Lower priority / untouched: `QuaPar.orderEdges`/`createP`'s "Face k has ... expected 2 but got
  N" error on near-degenerate/thin triangles; `partialConj` for `'cplq'`/`'pqp'`; `add` for
  `RatPol`/`RatPar`; conjugate engines `'pqp'`/`'graph'`; the standalone `RatPol.conj` gap.

## Relevant files

- `convEnvCPLQ.m` -- `tangentCevian` (new, replaces `edgeClipCevian`, itself replaced
  `seamPoint`/`buildEdgeAffinePiece` earlier this session) is the verified-correct Part 2c fix, no
  known issues. `solveTriangleBF`/`assemblePiecesBF` (new) are the `nCE==3` generalization whose
  scope is the open decision above -- reverting them is a self-contained, mechanical undo (see
  DESIGN.md for exactly what to revert to).
- `DESIGN.md` -- full derivation trail; search "Part 2c" under `convEnvCPLQ.m`'s entry for
  everything from this session, including the hedged 9/9-sample finding and both paths forward.
- `conjCPLQTest.m`, `maxQuaParTest.m` -- contain the 6 currently-failing tests; NOT yet touched
  this session (the fix decision above determines what changes, if any, they need).
- `cPLQ/` -- untracked clone of the original reference implementation; not touched this session.

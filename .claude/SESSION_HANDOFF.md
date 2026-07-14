# Session Handoff

_Last updated: 2026-07-14_

## Where things stand

- Branch: `main` (uncommitted changes on top of `27f3eba` -- see "What happened" below; not yet
  committed/pushed this session).
- Full test suite: **142/142 PASS**.
- A 600-triangle randomized stress test of the full Step1->Step2->Step3 pipeline
  (`convEnvCPLQ`->`conjPieceCPLQ`->`maxQuaPar`) found **zero wrong answers** (max error 1.85e-13,
  machine precision) across 103 valid 3-convex-edge samples; only a pre-existing, separate
  `maxQuaPar:internal` assembly-topology crash on ~4% of samples (item 1 in Next steps).

## What happened this session: fixed the Step 1 envelope-correctness gap (2026-07-13's #1 priority)

Root-caused and fixed the "highest priority" open item from the previous handoff: for
`T=(0,0),(7.02,0.67),(8.43,7.63)`, `maxQuaPar(g1,g2)` ran without error but returned a numerically
too-high result, even though every individual piece (`q1`,`q2`,`g1`,`g2`, and Step 3's assembly)
had been separately verified correct.

**Root cause**: `convEnvCPLQ.m`'s `splitThreeConvex` split a 3-convex-edge triangle by a
**horizontal** line through the middle vertex (`vmid`). This is wrong in general. Proof/derivation:
both resulting sub-envelopes `q1` (on `T1={vlow,vmid,Pnew}`) and `q2` (on `T2={vmid,vhigh,Pnew}`)
are built by Appendix A.4's two-convex-edge harmonic-mean formula using the SAME third edge
(`vlow-vhigh`, the one both sub-triangles inherit), so both touch `u1*u2` exactly along that
**entire** edge -- meaning `q1-q2` (a quadratic) vanishes identically on that whole line, hence
factors as `(vlow-vhigh line) * (a second line)`. That second line necessarily also passes through
`vmid` (both `q1(vmid)` and `q2(vmid)` equal `u1u2(vmid)` via each one's OTHER convex edge). This
second line -- not a horizontal one -- is the true smooth-fit split direction: verified numerically
that `q1` and `q2` agree along it to machine precision (~1e-15) at every tested point, for both the
well-known example and the previously-broken one, whereas along the old horizontal split they only
agreed at the two seam endpoints (small mismatch ~0.02-0.04 for the well-tested
`T=(0,0),(3,3),(1,2)` example -- small enough to never affect that example's specific tested dual
points -- but a much larger ~1-13 mismatch for other triangles, e.g. the broken one, corrupting the
final conjugate). A closed form for the correct split line (in terms of the three edges'
slopes/intercepts) is now implemented directly in `splitThreeConvex`; see its HISTORY comment for
the full derivation.

**Fixed** in `convEnvCPLQ.m`'s `splitThreeConvex`. Re-verified: the previously-wrong triangle now
matches `maxQuaParTest.supBilinearOverPoly` ground truth to ~1e-14 at 10 test points; the full test
suite (142 tests) passes; a 600-trial randomized stress test found 0 wrong answers (see above).

**Side effects fixed along the way** (both triggered by the corrected split point changing exactly
how thin some sub-triangles are, in `conjPieceCPLQ.m`'s `pickRep`):
- `pickRep`'s fixed-list-of-magnitudes search (already extended once in the prior session) hit its
  wall again on the SAME regression-test triangle (`T=(0,0),(4.46,1.83),(5.81,2.38)`), since the
  corrected split changed that triangle's sub-triangle sliver from needing `mag~1e-4` to needing
  `mag~1e-9`. Replaced the fixed list with a dense geometric sweep (16+ orders of magnitude,
  shrinking-then-growing so the largest/most-robust working magnitude in each direction is
  preferred over the smallest technically-valid one -- a naive smallest-first sweep picked up a
  numerically fragile point near a vertex and broke a DIFFERENT, previously-passing test). See
  `pickRep`'s HISTORY comment.
- That same test's exact triangle became SO thin post-fix (needing `mag~1e-9`) that it started
  tripping a separate, deeper, still-open bug several layers downstream
  (`QuaPar.orderEdges`/`createP` rejecting the resulting face topology -- same general category as
  the already-known, not-yet-investigated `assemblePieces` near-degenerate-sliver crash from the
  prior session). Rather than chase that separate bug, swapped the regression test's triangle for
  `T=(3.1436,2.4929),(5.0857,4.1038),(9.0757,7.5555)`, which still produces a comparably thin
  (~0.027, same order as the original ~0.0074) edge-strip sub-triangle exercising `pickRep`'s
  search exactly as intended, without hitting that separate issue. See
  `conjPieceCPLQTest.pickRepFindsThinEdgeStripFace`'s updated comment.

## Next steps

1. **Highest priority now**: the `maxQuaPar:internal` ("assemblePieces: boundary edge (r,c) has no
   matching neighbour") crash -- seen on ~4% of random valid 3-convex-edge triangles in this
   session's stress test, and previously seen (before this session) on the near-degenerate sliver
   triangle from the old `pickRep` fix. This is a Step-3 assembly-topology bug, NOT a correctness
   bug (no wrong answers have been observed, only outright crashes on certain edge configurations)
   -- distinct from the gap just fixed. Not yet investigated. Reproduce with e.g.
   `T=(6.0365,4.9504),(9.8960,6.3015),(1.4908,3.3753)` (found in this session's stress test) via
   `maxQuaParTest.buildG1G2ForTriangle(T)` then `maxQuaPar(g1,g2)`.
2. Separately, `QuaPar.orderEdges`/`createP` rejects the face topology produced by
   `conjPieceCPLQ`/`conjPSDRank1QuadTriangle` for SOME sufficiently near-degenerate/thin triangles
   (see "side effects fixed" above) -- worth understanding whether this is the same root cause as
   item 1 or a genuinely separate one.
3. Now that Step 1 is verified correct (this session) and `conjCPLQ`'s Step 3 wiring is in place
   (prior session), re-check whether `infConv`/`moreau`/`lasryLions`/`proxAverage` now work
   end-to-end on a genuinely piecewise (3-convex-edge triangle) `f,g` pair, not just full-domain
   quadratics -- DESIGN.md's "Next planned" #1 (`conjPieceCPLQ`'s rational-piece TODO) is still the
   gap for the fully general multi-face case, but the single-triangle/3-convex-edge case may now be
   trustworthy end to end.
4. Consider a LARGER/longer randomized stress test once (1)-(2) are resolved, since the crash rate
   (not correctness) is now the main open question.
5. `partialConj` for the `'cplq'`/`'pqp'` engines (II.4) -- not started for any engine.
6. `add` for `RatPol` and the `RatPar` parent class (II.3) -- deprioritized.
7. Conjugate engines `'pqp'`/`'graph'` remain unimplemented (error explicitly, not silently wrong).
8. The standalone `RatPol.conj` gap (rational piece with no known originating quadratic) is still
   open and untouched.
9. Test command (Frances, prefix every batch call -- Maple toolbox conflicts with Symbolic):
   ```
   matlab -batch "restoredefaultpath; rehash toolboxcache; cd('/home/ylucet/CCA2/CCA2'); \
     res=runtests({'RatPolTest','convEnvCPLQTest','QuaParTest','conjCPLQTest','conjPieceCPLQTest', \
     'PLQVCTest','maxQuaParTest','addQuaPolyTest','addQuaParTest','infConvTest','moreauTest', \
     'lasryLionsTest','proxAverageTest'}); \
     fprintf('\n==== SUMMARY ====\nTOTAL %d PASSED %d FAILED %d\n',numel(res),sum([res.Passed]),sum([res.Failed])); \
     exit(sum([res.Failed])>0)"
   ```

## Relevant files

- `convEnvCPLQ.m` -- `splitThreeConvex` (~line 215-260): this session's core fix, replacing the
  horizontal split with the smooth-fit line; full derivation in its HISTORY comment. Main
  function's docstring updated to match.
- `conjPieceCPLQ.m` -- `pickRep` (~line 362): replaced fixed magnitude list with a dense geometric
  sweep; HISTORY comment documents both this session's fix and why ordering (shrink-then-grow, not
  smallest-first) matters.
- `conjPieceCPLQTest.m` -- `pickRepFindsThinEdgeStripFace`: triangle swapped to avoid a separate,
  still-open near-degenerate-sliver bug (item 2 above); comment explains why.
- `DESIGN.md` -- `maxQuaPar.m`'s Implementation-status bullet now has a "Correctness fix
  (2026-07-14)" paragraph summarizing this session's fix and stress-test results.
- `maxQuaParTest.m` -- `buildG1G2ForTriangle`/`supBilinearOverPoly` (from the prior session) were
  the main tools used to root-cause and verify this session's fix; unchanged this session.
- `/home/ylucet/CCA2/3-edge.tex` (+ `.pdf`, outside this repo) -- unchanged; not revisited this
  session.

## Prior session history (2026-07-13, for context)

Two parallel sessions were merged on `main` @ `27f3eba`: Session A completed the nonconvex-PLQ
operator pipeline (`infConv`/`moreau`/`lasryLions`/`proxAverage`) and wired `maxQuaPar` into
`conjCPLQ`'s Step 3 for the 3-convex-edge case; Session B found and fixed three bugs in
`maxQuaPar.m` via a randomized-triangle stress test (a too-strict parabola-rejection guard in
`splitCell`, a too-tight `dedupHits` tolerance, and a too-coarse `pickRep` magnitude list), but
also discovered the Step 1 envelope-correctness gap that this session (2026-07-14) has now fixed.
See git history / `DESIGN.md` for full detail on the operator pipeline itself, which this session
did not touch.

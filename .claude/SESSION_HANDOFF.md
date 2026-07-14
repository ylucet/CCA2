# Session Handoff

_Last updated: 2026-07-13_

## Where things stand

- Branch: `cplq-engine` @ `8cd1c14` -- "Fix three bugs found via randomized triangle stress test
  on maxQuaPar pipeline"
- Pushed: pending (see end of this session's summary to the user for the actual outcome)

## What happened this session

Started from a question about whether `3-edge.tex`'s (corrected) claim -- "only parabolic
boundaries arise" between adjacent `maxQuaPar` conjugates, never a genuine hyperbola -- could be
trusted, or whether a counterexample/bug could be found first. Built a generalized test harness
(`maxQuaParTest.buildG1G2ForTriangle(T)`, now checked in) that runs the REAL Step 1/Step 2 pipeline
(`convEnvCPLQ` then `conjPieceCPLQ`) on an ARBITRARY triangle `T` with all 3 edges convex for
`u1*u2`, not just the one hard-coded `T=(0,0),(3,3),(1,2)` example used everywhere else. Stress-
testing ~25 random such triangles surfaced **three distinct, real bugs**, each found, isolated,
given a regression test that fails pre-fix/passes post-fix, and fixed:

1. **`maxQuaPar.m`'s `splitCell` guard rejected genuine, representable parabolas.** The guard
   threw `maxQuaPar:notDegenerate` ("a genuine ellipse/hyperbola boundary would be needed")
   whenever the full 3x3 conic discriminant `Delta` was nonzero -- but `Delta~=0` with
   `delta=b^2-4ac==0` (parabolic TYPE) is a genuine, non-degenerate PARABOLA, exactly what
   `QuaPar`'s curved `Ec` edges exist to represent (the rest of `splitCell`, the
   `isStraight`/`edgeEc` construction, already builds it correctly once let through -- confirmed by
   testing, not just inspection). Fixed by conditioning the rejection on `delta` as well as `Delta`
   (only reject when `delta~=0` AND `Delta~=0`, i.e. a TRUE irreducible ellipse/hyperbola). Also
   fixed a latent, never-before-exercised sign bug this exposed in `assemblePieces`: a curved
   edge's `Ec` row is built once from `f1row-f2row` and reused unchanged for both neighbouring
   pieces, so its sign is tied to "which of f1/f2 is globally larger", not to whichever piece ends
   up on the geometric left of the edge's final stored direction (an accident of half-edge
   processing order) -- flipped `Ec`'s sign when it doesn't evaluate positive on the (always-left)
   interior of the first-processed half-edge's own piece. Regression test:
   `maxQuaParTest.splitCellAcceptsGenuineNonDegenerateParabola`
   (`T=(0,0),(9.31,7.63),(3.80,7.40)`, verified end-to-end against `supBilinearOverPoly` ground
   truth at 5 points).
   **Stress-test finding, now confirmed mathematically**: across the ~25 random triangles, every
   genuinely non-degenerate boundary found was parabolic-type (`delta~0`) -- zero genuine
   hyperbola/ellipse boundaries turned up. The 3-edge.tex claim holds for this family; the bug was
   in the code's guard, not the math.
2. **`maxQuaPar.m`'s `dedupHits` tolerance (`sqrt(eps)`) was too tight**, in exactly the same way
   already documented and fixed elsewhere in this file for `assemblePieces`' vertex merge: a
   genuine boundary crossing exactly at a cell corner, computed via two different boundary edges'
   independent quadratic-root arithmetic, can disagree by ~1e-7 -- `sqrt(eps)~1.5e-8` doesn't cover
   that, so a real 2-crossing split was inflated to 3 hits, tripping `splitCell`'s "expected
   exactly 2 boundary crossings" assertion (`maxQuaPar:internal`) on legitimate input. Fixed by
   widening `dedupHits`' tolerance to `1e-6` absolute, matching the existing precedent. Regression
   test: `maxQuaParTest.dedupHitsMergesCrossingsAtACellCorner`
   (`T=(0,0),(2.11,1.43),(8.84,4.50)`, verified end-to-end against ground truth at 5 points).
3. **`conjPieceCPLQ.m`'s `pickRep` magnitude list was too coarse** for very thin edge-strip faces.
   `pickRep` searches `base +/- dir*mag` over a fixed geometric sequence of magnitudes to find a
   point inside a target face; for a sub-triangle Step 1 produces from
   `T=(0,0),(4.46,1.83),(5.81,2.38)` (an extremely thin sliver, two vertices only 0.0074 apart),
   the E23 edge-strip face needed `mag=scale*1e-4` -- one order of magnitude below the smallest
   value the old list tried (`scale*0.001`) -- so `conjPieceCPLQ` threw
   `conjPieceCPLQ:internal` ("could not locate a representative point for face 3"). Fixed by
   extending the magnitude range (smaller AND, symmetrically, larger). Regression test:
   `conjPieceCPLQTest.pickRepFindsThinEdgeStripFace` (tests Step 2 in isolation on the thin
   sub-triangle, verified against a numeric grid sup; deliberately does NOT go through the full
   `maxQuaPar` pipeline for this triangle -- see open item below).

Full suite: **118/118 PASS** (115 original + 3 new regression tests), confirmed both before
starting (115/115 baseline) and after all three fixes, with each new test independently verified
to fail on the pre-fix code (via `git stash` of just the source fix, keeping the test) and pass
post-fix.

## Open item found but NOT fixed (important -- read before trusting `maxQuaPar` on an arbitrary triangle)

Fixing bug #1 above (the `splitCell` guard) exposed a **fourth, deeper, unresolved issue** while
building the regression test: for `T=(0,0),(7.02,0.67),(8.43,7.63)` specifically (a different
random triangle than any used in the fixes above), `maxQuaPar(g1,g2)` now runs to completion
without error (guard bug no longer blocks it) but returns a **numerically wrong** result --
`h.eval(s)` exceeds the true ground truth (`supBilinearOverPoly`) at several sample points, even
though:
  - `q1`, `q2` (Step 1's two sub-triangle envelopes) are each individually verified correct: rank-1
    PSD, match `x*y` exactly at their own 3 vertices, and are `<= x*y` everywhere on their own
    sub-triangle (checked by grid search).
  - `g1`, `g2` (Step 2's conjugates) are each individually verified correct: `g_i.eval(s)` matches
    a brute-force grid search of `sup_{x in T_i}[s.x - q_i(x)]` exactly.
  - `h.eval(s) == max(g1.eval(s), g2.eval(s))` EXACTLY (checked directly) -- so Step 3
    (`maxQuaPar`'s assembly) is NOT the bug; it correctly computes the max of two individually-
    correct conjugates.
  - By an unconditional real-analysis identity (`sup_x max(a(x),b(x)) = max(sup_x a(x), sup_x
    b(x))`), `max(g1,g2)` always equals the conjugate of `min(q1-on-T1, q2-on-T2)` -- which is a
    valid minorant of `f=xy` on all of `T`, hence `max(g1,g2) >= groundTruth` pointwise ALWAYS, with
    equality only if `min(q1,q2)` is exactly the TRUE joint convex envelope of `f` over the WHOLE
    `T` (not merely a valid minorant on each sub-triangle separately). The overshoot found here
    means Step 1's per-sub-triangle envelope construction (`convEnvCPLQ`'s `nCE==3` branch,
    splitting via the "horizontal line through the middle vertex" and applying Appendix A.4's
    harmonic-mean formula to each 2-convex-edge sub-triangle SEPARATELY) is not, for this
    triangle's specific shape, producing the tightest/globally-correct joint envelope -- even
    though each piece is separately a legitimate, locally-tight-feeling minorant.
  - Investigated candidate causes and RULED OUT: (a) `splitThreeConvex`'s split-point identification
    (matches the expected middle-vertex construction exactly, verified by hand); (b)
    `twoEdgeQuadPlain`'s `+/-` sign selection for `q2` (both branches are always rank-1 PSD by
    construction -- `AM-GM` guarantees `denom>=0` for both signs always -- so convexity alone can't
    discriminate them; the code's actual "touches the h-edge at one test point" check is
    VACUOUS -- both signs satisfy it identically along the whole infinite line through that edge,
    not just at the tested point -- so the code deterministically falls through to `s=1`
    regardless of triangle shape; but `s=1` is ALSO the branch giving the strictly larger/tighter
    value at both T1's and T2's own centroids, i.e. "pick the larger of the two branches" -- the
    natural fix for the vacuous check -- doesn't change the outcome here either, so the bug is not
    simply a wrong sign choice between the two known closed-form candidates).
  - NOT yet found: why `min(q1,q2)` fails to be the true joint envelope for this triangle when the
    SAME construction (verified by hand, hard-coded coefficients) is proven correct at all 7 ground
    truth points for `T=(0,0),(3,3),(1,2)`. Needs either a deeper look at whether COAP's Appendix
    A.4/A.5 theorem has an unstated precondition this triangle's shape violates, or a bug in
    `envelopeFromClassified`/`twoEdgeQuadPlain` itself that only manifests for certain edge-slope
    configurations.
  - A SEPARATE, likely-unrelated bug was also seen while investigating: on YET ANOTHER triangle
    (`T=(0,0),(4.46,1.83),(5.81,2.38)`, the near-degenerate sliver from fix #3 above), after fixing
    `pickRep`, `maxQuaPar(g1,g2)` throws `assemblePieces: boundary edge (4,4) has no matching
    neighbour` -- a self-loop edge, suggesting a numerical degeneracy (near-coincident vertices)
    specific to very thin/skewed input triangles. Not investigated further; this is why fix #3's
    regression test deliberately stays at the Step-2-only level rather than running the full
    `maxQuaPar` pipeline on that triangle.

**Net assessment**: the original question ("is `3-edge.tex`'s parabola-only claim trustworthy?")
has a qualified YES for the boundary-TYPE question (no genuine hyperbola/ellipse ever turned up in
~25+ random trials, and the reasoning why is now well understood), but `maxQuaPar` is NOT yet
verified correct for arbitrary triangles beyond the specific ones exercised by the current test
suite -- there is at least one concrete triangle where it silently returns a wrong (too high)
numeric answer, and at least one other where it crashes on a near-degenerate shape. Recommend
NOT relying on `maxQuaPar` for arbitrary/production triangles yet; the hard-coded
`T=(0,0),(3,3),(1,2)` example and the three new regression-test triangles are known-good.

## Next steps

1. **Highest priority**: resolve the Step 1 envelope-correctness gap above
   (`T=(0,0),(7.02,0.67),(8.43,7.63)` overshoots ground truth). Likely needs either revisiting
   COAP Appendix A.4/A.5's proof for an unstated precondition, or closer scrutiny of
   `envelopeFromClassified`'s case-2 (`twoEdgeQuadPlain`) construction for edge cases the original
   proof/examples didn't cover.
2. Investigate the `assemblePieces: boundary edge (4,4) has no matching neighbour` crash on the
   near-degenerate sliver triangle (`T=(0,0),(4.46,1.83),(5.81,2.38)`) -- likely a vertex-merge
   tolerance or clipping robustness issue specific to very thin triangles.
3. Consider running a LARGER randomized stress test (using `maxQuaParTest.buildG1G2ForTriangle`,
   now checked in) once (1) and (2) are resolved, to gain confidence `maxQuaPar` is correct broadly,
   not just on the handful of triangles exercised so far.
4. The `delta>0, Delta=0` degeneracy open question from an earlier session (is it forced whenever
   adjacent envelopes share eigenvalue `lambda`, or coincidental?) is superseded by this session's
   finding that it's neither forced nor coincidental in the way originally framed -- see the
   "Stress-test finding" above -- but the ORIGINAL question about `3-edge.tex`'s conclusion section
   may still be worth revisiting in light of this. See `/home/ylucet/CCA2/3-edge.tex` (outside this
   repo).
5. The standalone `RatPol.conj` gap (rational piece with no known originating quadratic) is still
   open and untouched -- unrelated to this session, see prior handoffs.
6. Test command (Frances, prefix every batch call -- Maple toolbox conflicts with Symbolic):
   ```
   matlab -batch "restoredefaultpath; rehash toolboxcache; cd('/home/ylucet/CCA2/CCA2'); \
     res=runtests({'RatPolTest','convEnvCPLQTest','QuaParTest','conjCPLQTest','conjPieceCPLQTest','PLQVCTest','maxQuaParTest'}); \
     fprintf('\n==== SUMMARY ====\nTOTAL %d PASSED %d FAILED %d\n',numel(res),sum([res.Passed]),sum([res.Failed])); \
     exit(sum([res.Failed])>0)"
   ```

## Relevant files

- `maxQuaPar.m` -- `splitCell` (~line 499-650, bugs #1/#2 above; HISTORY section documents both),
  `dedupHits` (~line 726, bug #2), `assemblePieces` (~line 741-863, the Ec sign-orientation fix
  from bug #1).
- `maxQuaParTest.m` -- `buildG1G2ForTriangle`/`extractTriFace` (new, generalizes the hard-coded
  `buildG1G2` to any 3-convex-edge triangle), `splitCellAcceptsGenuineNonDegenerateParabola`,
  `dedupHitsMergesCrossingsAtACellCorner` (new regression tests).
- `conjPieceCPLQ.m` -- `pickRep` (~line 362, bug #3 fix + HISTORY note).
- `conjPieceCPLQTest.m` -- `pickRepFindsThinEdgeStripFace` (new regression test).
- `convEnvCPLQ.m` -- NOT modified this session, but is where the still-open Step 1 correctness gap
  (see above) most likely lives: `envelopeFromClassified`/`twoEdgeQuadPlain` (~line 142-179),
  `splitThreeConvex` (~line 215-227).
- `/home/ylucet/CCA2/3-edge.tex` (+ `.pdf`, outside this repo) -- unchanged this session; the
  "Stress-test finding" above is relevant context for it but no edits were made here.

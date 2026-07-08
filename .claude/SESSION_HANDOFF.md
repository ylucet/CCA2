# Session Handoff

_Last updated: 2026-07-07 (later session)_

## What happened this session
Picked up from the prior handoff. First correction: the prior handoff said the
`conjPSDRank1QuadTriangle` work was left uncommitted -- it was in fact committed and pushed since
(commit `918b1f1`, author-approved, tip of `cplq-engine` at session start). Verified 101/101 PASS
on Frances before starting new work.

**This session's work: handled the "known remaining gap" from the prior handoff** -- the
degenerate sub-case of `conjPSDRank1QuadTriangle` where two triangle vertices tie in the rotated
t-coordinate (the edge between them is exactly perpendicular to `A`'s nonzero eigenvector `u`).
Previously this raised `notImplemented` and was worked around in tests by picking a triangle that
avoids it (e.g. `(1,1),(4,3),(3,5)` instead of the natural `(0,0),(2,1),(1,2)`).

**Key fact used:** since `A = lam*u*u'` (rank 1), `A*v = lam*(u'v)*u` depends only on `t=u'v`, so
the two tied vertices `va,vb` (same `t0`) have the SAME dual point `s_ab = A*va+b = A*vb+b`. The
usual three collinear dual "vertices" `s1,s2,s3` collapse to two (`s_ab`, `s_c=A*vc+b`), and the
"kink" that normally splits one boundary chain into two segments disappears (there's no third
vertex at an intermediate t to create it) -- both boundary chains `va-vc` and `vb-vc` become single
segments spanning the full `[t0,tc]` range. Result: **5 faces instead of 6** -- two parabolic edge
strips `Ea` (through `va`), `Eb` (through `vb`), and three linear vertex cones `La,Lb,Lc`, meeting
at `s_ab` (apex of `Ea,La,Lb,Eb`) and `s_c` (apex of `Ea,Lc,Eb`). Implemented as
`conjPSDRank1QuadTriangleTie` (new function in `conjPieceCPLQ.m`), dispatched from
`conjPSDRank1QuadTriangle` when `t2-t1` or `t3-t2` is below tolerance (replacing the old
`notImplemented` error).

**Two sign bugs found and fixed via numeric ring-scans (same "search, don't derive by hand"
lesson as the prior session), not caught by the first test case alone:**
1. The ray direction `dray(m)` used for the `Ea`/`La` boundary always has gamma-rate `+1`
   (`uperp'*dray(m) == 1` identically, same fact the non-degenerate code already documents via
   `sgn13/sgn12`) -- but `Ea` (built from `ra`, the smaller of the two tied r-values) is *always*
   the gamma<0 chain (proof: maximizing `gamma*r` prefers the smaller `r` iff `gamma<0`), so its
   ray needs `-dray(ma)`, not `+dray(ma)`. `Eb`'s ray (`+dray(mb)`) was already correct since it's
   the gamma>0 chain by construction.
2. The `La`/`Lb` boundary lies on the line `perp(va-vb)` through `s_ab` -- but that whole line is
   an exact algebraic tie (`La(s)==Lb(s))` for its ENTIRE length, not just at one point, so no
   bisector-of-two-rays direction reliably lands inside one wedge, and there is no clean, general
   sign rule the way there was for `dA`/`dB` (tested two hand-derived candidates that worked on one
   example each but not both). Fixed with a new `pickRepSweep` helper: a full angular sweep (every
   10 degrees, a few magnitudes) around the wedge apex instead of a single hand-derived bisector
   direction -- more robust, in the same spirit as the existing `pickRep`'s magnitude search.

**Validated:** the two known reference triangles (the `xy`-envelope's `(0,0),(2,1),(1,2)` case, and
a hand-built generic rank-1 quadratic with a manufactured tie) both match numeric sup-sampling to
~1e-5/1e-6 (grid-noise level), PLUS a 40-trial random stress test (random eigenvector angle, random
tied r-values, random third vertex, random dual query points) -- 40/40 passed, worst error 1.4e-4
(consistent with the coarse grid, not a real error). New tests: `psdRank1QuadraticTieConjugate`
(hand-built, direct) and `psdRank1QuadraticTieEndToEnd` (through the real `convEnvCPLQ` ->
`conjPieceCPLQ` pipeline, using the exact `(0,0),(2,1),(1,2)` triangle that used to be the
documented degenerate case). Updated the header docstring and the stale comment in
`psdRank1QuadraticEndToEnd` (which used to say the tie case was "not yet handled").

Full suite **103/103 PASS** on Frances (101 previous + 2 new).

## Where things stand
- Branch: `cplq-engine` @ `918b1f1` (this session's starting point), working tree has this
  session's changes UNCOMMITTED (author's standing instruction: ask before every commit/push, no
  exceptions -- always ask, never assume a prior "yes" carries forward).
- Files changed this session: `conjPieceCPLQ.m` (new `conjPSDRank1QuadTriangleTie` +
  `rank1TieWinner` + `pickRepSweep` helpers, dispatch branch inside `conjPSDRank1QuadTriangle`,
  updated header docstring), `conjPieceCPLQTest.m` (new `psdRank1QuadraticTieConjugate`,
  `psdRank1QuadraticTieEndToEnd`, updated stale comment in `psdRank1QuadraticEndToEnd`), this
  handoff file.
- Not pushed.

## Next steps
- Ask the author whether to commit/push this session's work (standing rule: always ask).
- Conjugate of RATIONAL pieces (the 1-convex-edge `RatPol` output of Step 1) is still the
  remaining BIG gap before the full `conjCPLQ` pipeline (Step1 `convEnvCPLQ` -> Step2
  `conjPieceCPLQ` per piece -> Step3 max) is wired end to end for nonconvex inputs. This needs
  genuinely curved (parabolic) edges in the OUTPUT too, not just the input -- `QuaPar`'s
  curved-edge `orderEdges` support (commit `ccbaba6`) is already there for this. Consider reading
  `biconjugate2/plq_1piece.m`'s type-1 branch (rational envelope conjugate, general polygon) in
  the reference repo (`github.com/tanmaya11/convex`, cloned to scratchpad in a prior session, not
  part of this repo) before re-deriving from scratch.
- Then Step 3 (max of per-piece conjugates) to complete `conj(QuaPoly)` for nonconvex inputs.
- Test command (run directly on Frances -- MATLAB now runs there; prefix every batch call with
  `restoredefaultpath; rehash toolboxcache;` because the Maple toolbox conflicts with Symbolic):
  ```
  matlab -batch "restoredefaultpath; rehash toolboxcache; cd('/home/ylucet/CCA2/CCA2'); \
    res=runtests({'RatPolTest','convEnvCPLQTest','QuaParTest','conjCPLQTest','conjPieceCPLQTest','PLQVCTest'}); \
    fprintf('\n==== SUMMARY ====\nTOTAL %d PASSED %d FAILED %d\n',numel(res),sum([res.Passed]),sum([res.Failed])); \
    exit(sum([res.Failed])>0)"
  ```

## Relevant files
- `conjPieceCPLQ.m` -- this session's main file; new `conjPSDRank1QuadTriangleTie` +
  `rank1TieWinner` + `pickRepSweep` near the end (right after `rank1Winner`), dispatch branch
  inside `conjPSDRank1QuadTriangle` (checks `t2-t1 < tol || t3-t2 < tol` before the old collinear-
  dual-vertex construction).
- `conjPieceCPLQTest.m` -- `psdRank1QuadraticTieConjugate`, `psdRank1QuadraticTieEndToEnd`.
- `/home/ylucet/CCA2/DESIGN.md` -- design proposal (outside this repo).
- `/home/ylucet/CCA2/s10589-026-00781-5.pdf` (COAP) -- Appendix A (convex envelope) and Appendix B
  (conjugate); the tied-vertex case isn't explicitly in the paper (it's a measure-zero geometric
  coincidence of the *triangle*, not a distinct case in the theory) -- this session's derivation
  was from first principles, following the same rotate-to-(t,r) approach as the non-degenerate
  `conjPSDRank1QuadTriangle`.
- Reference repo (cloned to scratchpad in a prior session, not part of this repo):
  `github.com/tanmaya11/convex`, branches `b2`/`b3`/`biconjugate`/`biconjugate2` -- worth reading
  before tackling the rational-piece conjugate next.

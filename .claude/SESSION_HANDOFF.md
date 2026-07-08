# Session Handoff

_Last updated: 2026-07-08_

## What happened this session
Picked up from the prior handoff. First correction: the prior handoff said the tied-vertex
`conjPSDRank1QuadTriangleTie` work was left uncommitted -- it was in fact committed and pushed
since (commit `8084a41`, author-approved, tip of `cplq-engine` at session start). Verified
103/103 PASS on Frances before starting new work.

**This session's work: generalized `conjPieceCPLQ`'s indefinite-quadratic branch** from the exact
pure-`x*y` bilinear special case to ANY genuinely indefinite quadratic `q(x)=1/2 x'Ax+b'x+c`
(`beta*x*y+linear`, arbitrary rotation and shift) -- one of the two TODOs the file header had
flagged (`beta*x*y+linear (scale/shift)`).

**Key realization that scoped the work:** the prior handoff framed the remaining big gap as
"conjugate of RATIONAL pieces (the 1-convex-edge `RatPol` output of Step 1)". Investigation (incl.
re-reading the reference repo's `plq_1piece.m` type-1 branch, and the COAP paper's actual Appendix
B, which turns out to cover ONLY the conjugate of a bilinear function's own envelope, not a general
rational quad/linear conjugate) showed that a literal "conjugate an arbitrary `RatPol` with nonzero
denominator" is a genuinely harder, different problem: there is no known reduction to an
already-solved case when the caller has no memory of an original pre-envelope piece (e.g. a direct
`RatPol.conj()` call). That remains UNSOLVED and is explicitly still rejected
(`conjPieceCPLQTest/rationalRejected`, unchanged).

However: since `f* = (conv f)*`, Step 2 does NOT need to consume Step 1's rational envelope output
at all when the ORIGINAL per-triangle piece is still available to the caller -- it can conjugate
that original piece directly. For the one-convex-edge bilinear case, the original piece is exactly
`beta*x*y+linear` in some rotated/shifted frame -- so generalizing the indefinite-quadratic branch
is what actually lets a future `conjCPLQ` orchestrator wire Step 1 -> Step 2 -> Step 3 for
nonconvex inputs, without needing the harder standalone-`RatPol.conj` capability.

**Implementation** (`conjIndefiniteQuadTriangle` + `pushforwardQuaParDual`, new functions in
`conjPieceCPLQ.m`, mechanical linear algebra, no new curved-region derivation):
1. Rotate to the bilinear frame `u=Mx` via `bilinearFrame(A)` (same construction as
   `convEnvCPLQ`'s, duplicated locally): `1/2 x'Ax = u1*u2` exactly, so `beta=1` always with this
   specific `M` -- no separate "scale" case needed.
2. Complete the square on the residual linear part: `u1*u2+d*u1+e*u2+f0 = (u1+e)*(u2+d)+f0-ed`, a
   PURE TRANSLATION `u'=u+[e,d]` (edge convexity depends only on slope, so classification is
   unaffected by the shift -- classify once on the shifted triangle).
3. Reuse `conjBilinearXYzeroCE`/`conjBilinearXYoneCE` UNMODIFIED on the shifted triangle for the
   pure-bilinear conjugate in `s_u`-space.
4. Undo the shift: subtract a LINEAR-IN-`s_u` correction from every face's VALUE only (the region
   partition -- `Vd`,`E`,`Ec`,`F` -- is unaffected, since the correction is constant w.r.t. the
   argmax variable).
5. Undo the rotation: push the whole `QuaPar` through the linear dual-variable change
   `s_u = M^-T s` (`pushforwardQuaParDual`) -- transforms dual vertices (`Vd*M`), each face's
   quadratic (Hessian/linear part via `Minv*H*Minv'`, `Minv*L`), each edge's conic the same way
   (using its own PLAIN, non-doubled coefficient convention), and swaps left/right (`F` columns,
   `Ec` sign) if `det(M)<0` (a reflection).

**One real bug found and fixed via the stress test** (not a math error, a numerical-conditioning
one): `pushforwardQuaParDual`'s conic transform can amplify coefficient magnitudes, and
`QuaPar.assertParabolicEdges`'s `b^2-4ac<=sqrt(eps)` check uses a fixed ABSOLUTE tolerance, so a
large-magnitude (but perfectly valid) conic row could spuriously fail it. Fixed by normalizing
each transformed conic row by its own max-abs entry before constructing (a conic equation is
scale-invariant, so this is free and doesn't touch the shared `QuaPar` validation code).

**Validated:** 200-trial random stress test (random rotation angle, random eigenvalues of both
signs, random shift, random triangle, random dual query points) against numeric sup-sampling over
a fine triangle grid, run ad hoc on Frances (not committed) -- 0 mismatches across 498 checked
dual points (100/200 trials skipped for landing in the documented two-convex-edge dead end), worst
error 2.4e-5 (grid-noise level). Two new PERMANENT tests committed:
`indefiniteQuadraticZeroConvexEdgeConjugate` (`A=[0 3;3 0]`, i.e. `beta=3` not 1, plus a linear
shift, over the unit triangle -- exercises BOTH the scale-into-eigenvalues path and the shift
correction) and `indefiniteQuadraticOneConvexEdgeConjugate` (`A=[2 1;1 -2]`, genuinely rotated, not
axis- or diagonal-aligned, plus a shift -- exercises a nontrivial `det(M)` and the full
pushforward). Both cross-checked against numeric sup over a fine grid, matching the existing style
(e.g. `psdRank1QuadraticConjugate`).

Full suite **105/105 PASS** on Frances (103 previous + 2 new).

## Where things stand
- Branch: `cplq-engine` @ `8084a41` (this session's starting point). This session's changes are
  committed as `<FILL IN AFTER COMMIT>` (author-approved) and pushed.
- Files changed this session: `conjPieceCPLQ.m` (new `conjIndefiniteQuadTriangle`,
  `bilinearFrame`, `pushforwardQuaParDual`; replaced the old narrow pure-`x*y` dispatch branch and
  the old catch-all error with the general indefinite dispatch + an explicit concave-quadratic
  error; updated header docstring), `conjPieceCPLQTest.m` (new
  `indefiniteQuadraticZeroConvexEdgeConjugate`, `indefiniteQuadraticOneConvexEdgeConjugate`), this
  handoff file.

## Next steps
- Standalone `RatPol.conj` (conjugate of a genuine rational quad/linear piece with NO known
  originating quadratic) is still the remaining BIG unsolved gap -- needs a genuinely new 2D
  critical-point / curved-region derivation (unlike the indefinite-quadratic case, there is no
  reduction to an already-solved case). Reading `biconjugate2/plq_1piece.m`'s type-1 branch
  (`github.com/tanmaya11/convex`, NOT cloned into this repo -- re-clone to scratchpad if needed,
  see this session's notes) is still the best starting point; expect it to need genuinely curved
  edges in the OUTPUT domain (`QuaPar`'s curved-edge `orderEdges` support, commit `ccbaba6`, is
  already there for this) and a nontrivial amount of new geometric case analysis (interior
  critical point in a 2D rotated frame that does NOT decouple the way the rank-1-quadratic case
  did, since the rational numerator's dependence on both rotated coordinates does not vanish).
- Separately/in parallel: wire an actual `conjCPLQ` orchestrator (Step1 `convEnvCPLQ` -> Step2
  `conjPieceCPLQ` per ORIGINAL triangle piece (not Step 1's RatPol output, using this session's
  generalized indefinite-quadratic support) -> Step3 max of conjugates) for nonconvex inputs; this
  is now realistic without waiting on the `RatPol.conj` gap above, since `conjPieceCPLQ` can
  already conjugate every triangle-piece shape convEnvCPLQ's classification produces (affine, PD
  quadratic, rank-1 PSD quadratic incl. the tie sub-case, and now general indefinite with 0 or 1
  convex edge -- only the 2-convex-edge indefinite case is a genuine dead end, and Step 1 already
  proves that never arises after its own rank-1-PSD convexification).
- Test command (run directly on Frances -- MATLAB now runs there; prefix every batch call with
  `restoredefaultpath; rehash toolboxcache;` because the Maple toolbox conflicts with Symbolic):
  ```
  matlab -batch "restoredefaultpath; rehash toolboxcache; cd('/home/ylucet/CCA2/CCA2'); \
    res=runtests({'RatPolTest','convEnvCPLQTest','QuaParTest','conjCPLQTest','conjPieceCPLQTest','PLQVCTest'}); \
    fprintf('\n==== SUMMARY ====\nTOTAL %d PASSED %d FAILED %d\n',numel(res),sum([res.Passed]),sum([res.Failed])); \
    exit(sum([res.Failed])>0)"
  ```

## Relevant files
- `conjPieceCPLQ.m` -- this session's main file; new `conjIndefiniteQuadTriangle`,
  `bilinearFrame`, `pushforwardQuaParDual` right after the bilinear-conjugate helpers; dispatch
  changed inside the main `conjPieceCPLQ` function (the `elseif` chain after the rank-1 PSD case).
- `conjPieceCPLQTest.m` -- `indefiniteQuadraticZeroConvexEdgeConjugate`,
  `indefiniteQuadraticOneConvexEdgeConjugate`.
- `/home/ylucet/CCA2/DESIGN.md` -- design proposal (outside this repo); II.5.1 describes the
  3-step `cplq` engine.
- `/home/ylucet/CCA2/s10589-026-00781-5.pdf` (COAP) -- Appendix B is SPECIFICALLY "Conjugate of
  convex envelope of a bilinear function over a triangle" (B.1/B.2/B.3, matching
  `conjBilinearXYzeroCE`/`OneCE`/the documented 2-CE dead end), NOT a general rational-piece
  conjugate formula -- that gap is genuinely open, not just unread.
- Reference repo `github.com/tanmaya11/convex` (branches `b2`/`b3`/`biconjugate`/`biconjugate2`,
  NOT part of this repo, NOT currently cloned anywhere persistent -- re-clone to scratchpad):
  `plq_1piece.m`'s `conjugateFunction`/type-1 branch (~line 2387) is the general symbolic
  (Maple/MATLAB symbolic-toolbox) rational-conjugate machinery for arbitrary polygons; worth
  reading again before attempting the standalone `RatPol.conj` case, though it computes normal
  cones/subdifferentials symbolically for a much more general polygon setting than this project's
  triangle-only, numeric approach, so expect to adapt rather than port directly.

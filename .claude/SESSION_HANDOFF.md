# Session Handoff

_Last updated: 2026-07-17T23:40:00Z (superseded below)_

## What happened this session (2026-07-17, third session)

Continued the Part 2 investigation (the "2-convex-edge tightness bug, `q1` not tight throughout
its own sub-triangle" gap left open by the prior session). Found and fixed a genuine, previously
undiagnosed root cause (Part 2a), verified it against ground truth and the full test suite, then
went looking for why it doesn't close the gap to 0% and found a second, deeper, still-open issue
(Part 2b). Full derivation trail is in `DESIGN.md`'s `convEnvCPLQ.m` entry -- search "Part 2a,
FIXED this session" and "Part 2b, STILL OPEN".

**Part 2a (fixed, verified, kept)**: the cevian LOCATION used by `splitTwoConvexEdges` was itself
wrong -- a separate bug from Part 1 (which only fixed the FORMULA for the "other" sub-triangle, not
where the cevian sits). Root cause: the old `seamPoint`/`buildEdgeAffinePiece` pair located the seam
by factoring `q1` against a placeholder that merely matches the weak edge's AFFINE CHORD -- a valid
boundary condition on the triangle's boundary, but not the actual condition that bounds `q1`'s own
region of tightness. Re-derived the correct criterion directly from [COAP] Section A.1's own dual
formalism (the constructive envelope theorem: the implied dual point at `(x,y)` is
`(a,b) = grad q1(x,y)`; `q1` is valid exactly while each classified edge's own tangency point
`s_j(a,b)` stays inside that edge's segment; since `q1` is quadratic, `s_j(x,y)` is affine, so its
boundary is a genuine straight line, matching Locatelli). Confirmed the old seam sat at almost
exactly HALF the correct distance from the shared vertex `P` on the paper's own worked example
(`V=(2,1),(0,0),(1,0)`: old seam at x=0.293 along the target edge, correct at x=2-sqrt(2)=0.586),
which also means **the paper's own Appendix A.4 claim ("the domain is the entire triangle") is false
even on the paper's own worked example** -- confirmed concretely: `q1(0.474343,0)=-0.042780` there
instead of the true value 0 (since `f=xy` is identically 0 along that whole weak edge, hence `>=0`
is a valid global lower bound there via the trivial constant-0 minorant). New `edgeClipCevian`
(replacing `seamPoint`+`buildEdgeAffinePiece` in `convEnvCPLQ.m`) builds the correct line directly
from `q1`'s own plain coefficients and intersects it with the other convex edge -- same final
"intersect with `y=m*x+q`" step as before, just fed the right line. Verified against ground truth
(numerically-maximized biconjugate) on the paper's own example including at the specific wedge point
the old seam misclassified; full suite 79/79 (two pre-existing tests had hardcoded expectations that
were themselves artifacts of the two bugs above -- both updated with explanatory comments, not
silently changed: `convEnvCPLQTest.bilinearTwoConvexEdgesQuadratic`,
`conjPieceCPLQTest.psdRank1QuadraticEndToEnd`). **Quantified impact**: a 1500-triangle `rng(2026)`
stress test (regenerated fresh this session, not the prior session's uncommitted script) shows the
aggregate gap rate drop from ~45% to **23.1%** (108/468 genuine 2-convex-edge splits still show a
detectable gap) -- a real, substantial, independently-verified improvement, but not a full fix.

**Part 2b (found this session, NOT fixed, deeper than 2a, but substantially narrowed)**:
investigating why Part 2a doesn't close the gap to 0% found a second, distinct phenomenon. Concrete
repro: `V=(-5.2645,3.4904),(3.1062,0.5450),(5.0430,7.7766)`,
`f6=[0.3963 0.8289 0.0284 -2.4834 2.1149 -5.5999]` (QuaPoly's stored `[H11 H12 H22 L1 L2 kappa]`
convention, i.e. `q=0.5*x'Hx+L'x+kappa`). With Part 2a's corrected cevian, `q1`'s own sub-triangle
`T1={P,A,Ra}` (bilinear frame) is confirmed exact along ALL THREE of `T1`'s own edges in their
entirety, yet undershoots truth by up to ~0.38 in `T1`'s interior.
**Key finding**: `max(q1, R_Aw)` (`R_Aw` = `T2`'s OWN Appendix-A.3 formula, extended past its home
triangle `T2` into `T1`) matches ground truth EXACTLY everywhere checked in `T1` (~70 random
points) -- so the missing piece is not a new, undiscovered formula, it is `R_Aw` itself, valid over
a LARGER region than the current single cevian gives it. However `R_Aw` is unconditionally `>= q1`
throughout `T1` (including near `P`, where it wrongly overshoots by as much as ~2.8), so this is
NOT simply "take the pointwise max" -- the correct split is a genuine geometric boundary, not a
value comparison.
**Precise shape, mapped via a full barycentric grid scan of `T1`** (66 points): the gap is exactly
zero for more than ~30% of the way from `A-Ra` back toward `P`, and zero along the whole `P-A` and
`P-Ra` edges -- i.e. NOT a broad interior "bubble" as first thought from coarser sampling, but a
SMALL, THIN notch hugging the MIDDLE of the `A-Ra` edge, vanishing at both the `A` and `Ra` ends and
fading out within ~20-30% of the depth back to `P`. Each row of the scan (fixed distance from `A-Ra`)
has TWO transitions, not one, meaning the correction is not a single half-plane cut -- more like a
small triangular notch with its own apex, closer to `A-Ra` than to `P`. Bisected several transition
points precisely; checked whether the "lower" transition set lies on the extended `P-Ra` line
itself -- close but with a real, growing residual (not solver noise), so it is a nearby but
DIFFERENT line, not yet identified in closed form.
Ruled out as explanations (see DESIGN.md for full detail): the `twoEdgeQuadPlain` +/- branch choice
(the other branch is a strictly looser minorant everywhere, not relevant); a naive recursive
re-application of `edgeClipCevian` to `T1` itself (one candidate is degenerate, reproducing the
same `Ra`; the other lands just outside its target segment, `t~1.10`); the two
vertex-anchored-at-`P` rational candidates `R_Ph`/`R_Pw` (wildly invalid in the notch, ~-16 vs a
true value in the single digits); comparing `q1` directly against `R_Aw` as an equation (`q1=R_Aw`)
-- a genuine conic (symbolically confirmed), not the right boundary condition, consistent with
`R_Aw` never actually crossing back below `q1` inside `T1`.
**Leading hypothesis, not yet derived**: `R_Aw` has its own Section A.1 dual-tangency validity
limit, exactly analogous to `q1`'s (`edgeClipCevian`'s own derivation) -- most likely tied to edge
`h` (`mh,qh`) becoming "visible" again once extended past `A` into `T1`, since `R_Aw`'s own
derivation (Appendix A.3's "V-E" case, `V=A`, edge `w`) implicitly assumed edge `P-A` is NON-convex,
which is false once inside `T1`. The concrete next step is to redo `edgeClipCevian`'s exact
derivation but for `R_Aw`'s own gradient (a rational function -- needs the quotient rule, more work
than `q1`'s polynomial gradient) instead of numerically chasing more transition points. This gap
was NOT visible on the paper's own (near-symmetric, `mh=1,mw=0.5`) example -- confirmed via a fresh
check there (no notch found) -- so it should degenerate to "outside the triangle" as `mh/mw -> 1`
and grow as the ratio moves away from 1, a useful property to check any proposed closed form
against.

## Where things stand

- Branch: `main`. Contains this session's `convEnvCPLQ.m` fix (Part 2a: `edgeClipCevian` replaces
  `seamPoint`+`buildEdgeAffinePiece`), two test-expectation updates (both with explanatory
  comments), and `DESIGN.md`'s full diagnosis for both Part 2a (fixed) and Part 2b (still open).
- Not yet committed/pushed this session -- see below.
- Full MATLAB suite: 79/79 passing (`convEnvCPLQTest`, `conjPieceCPLQTest`, `conjCPLQTest`,
  `maxQuaParTest`, `RatPolTest`, `QuaParTest`).
- Untracked, not part of this repo: `cPLQ/` (clone of the original reference implementation); not
  touched this session, but its `quadQuad.m`/`quadLinear.m` scratch derivation scripts (using dual
  variables `a,b` and explicit `dl,du` edge-clip bounds) were a useful hint pointing toward the
  Section A.1 dual formalism this session's Part 2a fix is built on -- worth re-reading early next
  session if picking Part 2b back up, they may already contain the missing piece in an
  unfinished/uncommented form.
- Scratch diagnostic scripts (MATLAB + Python) live under the OS-temp scratchpad (session-specific,
  not persisted) -- regenerable from DESIGN.md's description if needed again. Note: use
  `C:\Users\ylucet\AppData\Local\miniforge3\python.exe` for `numpy`/`scipy`/`sympy` (has `pypdf`
  too, useful for reading the reference PDFs in `reference/` when grep-ing DESIGN.md isn't enough)
  -- the bare `python3` on PATH is a different interpreter without these.
- MATLAB CLI note: `matlab -batch "run('script.m')"` intermittently fails with "Incorrect number or
  types of inputs or outputs for function QuaPoly" on scripts that call the 4-arg `QuaPoly`
  constructor, even though the same code works fine when passed inline via `-batch "<code>"`
  directly (not through `run`) or when called as a function file (`function foo(); ...; end`) via
  `addpath(...); foo`. Cause not fully diagnosed (looks like a `run`-specific scoping quirk, not a
  real constructor bug) -- prefer inlining short scripts directly into `-batch` or wrapping longer
  ones in a named function file when scripting MATLAB from the shell.

## Next steps

- **Top priority**: find Part 2b's missing split for `T1`'s own interior bubble. Per Locatelli one
  is guaranteed to exist and must be polyhedral. Leading hypothesis: an interior-anchored fan split
  (3+ pieces meeting at a new interior point), not another boundary cevian -- untested. Re-reading
  `cPLQ/quadQuad.m`'s `dl,du` dual-bound derivation closely (see above) may shortcut re-deriving
  this from scratch.
- Decide whether/how to unblock `conjPieceCPLQ` for a genuinely rational (quadratic/linear) piece --
  still open, unchanged from before, not attempted this session either.
- Fill in the exact `[LOCATELLI]` citation in `DESIGN.md`'s reference list -- still not done.
- Still open, unchanged from before: the 2/741 (0.3%) residual `maxQuaPar:internal` crashes (see
  earlier handoff history / DESIGN.md for repro triangles) -- not investigated this session.
- Lower priority / untouched: `QuaPar.orderEdges`/`createP`'s "Face k has ... expected 2 but got
  N" error on near-degenerate/thin triangles; `partialConj` for `'cplq'`/`'pqp'`; `add` for
  `RatPol`/`RatPar`; conjugate engines `'pqp'`/`'graph'`; the standalone `RatPol.conj` gap.
- Decide whether to commit/push this session's changes (pending as of this writing).

## Relevant files

- `convEnvCPLQ.m` -- `splitTwoConvexEdges` (Part 2a fix: `edgeClipCevian`, new, replaces
  `seamPoint`+`buildEdgeAffinePiece`, both deleted). Locates the cevian from `q1`'s own gradient via
  the Section A.1 dual-tangency criterion instead of the old, wrong line-factoring heuristic. Part
  2b (the interior "bubble" in `T1`) is NOT fixed here -- `T1` still uses `q1` unchanged, over its
  now-correctly-sized-but-still-not-fully-tight domain.
- `convEnvCPLQTest.m` -- `bilinearTwoConvexEdgesQuadratic` updated: the paper's own single-quadratic
  formula is only checked at points confirmed to be on the `q1`-valid (`P`) side of the corrected
  boundary; the point that used to be (wrongly) asserted equal to `q1` there is removed, with a
  comment pointing to `bilinearTwoConvexEdgesSplitIsTight` for the comprehensive check.
- `conjPieceCPLQTest.m` -- `psdRank1QuadraticEndToEnd` updated: `T1`'s shape changed (now correctly
  smaller) under Part 2a's fix, which happens to move this specific example into the tied-vertex
  dispatch (`nf=5`, not the generic `nf=6`) -- updated with a comment explaining why, values
  reverified against ground truth.
- `DESIGN.md` -- full derivation trail for Part 2a (fixed) and Part 2b (open); search "Part 2a,
  FIXED this session" and "Part 2b, STILL OPEN" under `convEnvCPLQ.m`'s entry.
- `doc/bug.svg`, `doc/bug_pipeline.svg` -- prior sessions' figures for the ORIGINAL counterexample
  triangle (`T=(-9.95506,...)`); still accurate background for Part 1/the original diagnosis, not
  updated this session (this session's Part 2a/2b repros use different, smaller/cleaner examples --
  the paper's own `V=(2,1),(0,0),(1,0)` for 2a, and a fresh `rng(2026)`-stress-test triangle for 2b).
- `cPLQ/` -- untracked clone of the original reference implementation; not modified, but see the
  "Where things stand" note above about its `quadQuad.m`/`quadLinear.m` being a possible shortcut.

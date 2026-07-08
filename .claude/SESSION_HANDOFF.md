# Session Handoff

_Last updated: 2026-07-07 (local session end)_

## What happened this session
Started from the prior handoff's plan to extend the raw-bilinear `conjPieceCPLQ` construction to
nCE==2 (two convex edges). That plan turned out to be based on a flawed premise -- see "Key
finding" below -- so the session pivoted to implementing the case that's actually needed.

**Key finding (dead end, documented, not implemented):** the raw indefinite-bilinear "conjugate
of `x*y` over a triangle with TWO convex edges" needs a boundary between the two edges' own COAP
B.2 quadratic formulas, and that boundary is a genuine HYPERBOLA whenever the two slopes differ
(discriminant = `(m1-m2)^2/(4*m1*m2) > 0`) -- verified both numerically (2M+ random dual points,
zero mismatches vs. true 2D sup) and symbolically (sympy `factor`). `QuaPar` only supports
parabolic/linear edges, so this cannot be represented. This case also **does not arise from the
wired pipeline**: Step 1 (`convEnvCPLQ`) always convexifies a 2-convex-edge piece into a rank-1
PSD quadratic (its harmonic-mean formula, `twoEdgeQuadPlain`, has discriminant identically zero
for *any* slopes -- proved via sympy, `s^2-1=0` since `s=+-1`) before Step 2 ever sees it, so
`conjPieceCPLQ` is never actually asked to conjugate the raw indefinite piece in that case. The
old pinned test `bilinearTwoConvexEdgesNotImplemented` was renamed
`bilinearTwoConvexEdgesDocumentedLimitation` with this explanation; it still pins the same
`verifyError`/`notImplemented` behavior, just now as an intentional limitation, not a TODO.

**What was actually needed and got implemented:** `conjPSDRank1QuadTriangle` (in
`conjPieceCPLQ.m`) -- the conjugate of a rank-1 PSD quadratic `q(x)=1/2 x'Ax+b'x+c` (one zero
eigenvalue) over a triangle, which *is* the real Step-2 case for a 2-convex-edge piece. Six-face
QuaPar: three PARABOLIC edge strips (one per triangle edge -- every edge, not just "convex-for-xy"
ones, since `q` is jointly convex) + three LINEAR vertex cones, meeting at three COLLINEAR dual
points `s_i = A*v_i+b` (collinear because rank-1 `A*v` only varies along its eigenvector `u`).
Derived by rotating to `t=u'x` (curved direction) / `r=uperp'x` (flat direction, where `q` is
convex-but-affine), reducing to a 1D "conjugate of a concave... er, convex quadratic over a
possibly-kinked piecewise-linear chain" problem per half-plane (sign of `gamma=s'*uperp-w`).
Wired into `conjPieceCPLQ`'s dispatch as a new `elseif` between the strict-PD branch
(`conjConvexQuadTriangle`) and the pure-bilinear-frame branch (order matters: PD check must run
first). Representative-point selection for `QuaPar`'s `orientFaces` turned out to be fragile by
hand (two bugs: wrong bisector sign for a vertex wedge, and a strip-interior offset large enough
to overshoot past the strip's own valid t*-range) -- fixed with an oracle-based `pickRep` helper
that evaluates the 6 gated candidate formulas directly and searches over sign x magnitude rather
than trusting a hand-derived geometric argument. New tests: `psdRank1QuadraticConjugate` (hand-
built rank-1 quadratic, numeric-sup-validated) and `psdRank1QuadraticEndToEnd` (through the real
`convEnvCPLQ` -> `conjPieceCPLQ` pipeline). **Known remaining gap**: `conjPSDRank1QuadTriangle`
errors with `notImplemented` on a degenerate input where two triangle vertices tie in the rotated
t-coordinate (an edge exactly perpendicular to `A`'s eigenvector) -- hit by the very triangle used
elsewhere in the test file, `(0,0),(2,1),(1,2)`, for that reason `psdRank1QuadraticEndToEnd` uses
`(1,1),(4,3),(3,5)` instead. Not yet handled; would need a tie-breaking rule (e.g. tiny
perturbation, or a direct 2-vertex/1-vertex special case) if it turns out to matter in practice.

Full suite **101/101 PASS** on Frances (this session ran directly on Frances, not round-tripped).
Author was asked how to proceed after the hyperbola finding; answered with two new resources not
previously known to me: the convex-envelope paper is also expected at `./doc/...` (not found in
this repo or the `github.com/tanmaya11/convex` clone -- may be elsewhere or not yet added) and the
original research code is at `github.com/tanmaya11/convex` (cloned to scratchpad, branches `b2`,
`b3`, `biconjugate`, `biconjugate2` are relevant -- `biconjugate2/plq_1piece.m`'s general
`conjugateFunction` treats domains as general polygons via a `psi0/psi1/zeta` parametrization
(`convexExpr.m`, types 1=rational and 3=affine-or-bilinear-product; types 2/4 are marked
"disable this" by the original author) rather than case-splitting by convex-edge count -- worth
reading fully if/when tackling the rational-piece conjugate below, but was not needed for this
session's rank-1-quadratic derivation, which was done from first principles instead.

## Where things stand
- Branch: `cplq-engine` @ `04988db` locally, working tree has the changes described above
  UNCOMMITTED (author's standing instruction: ask before every commit/push, no exceptions).
- Files changed this session: `conjPieceCPLQ.m` (new dispatch branch + `conjPSDRank1QuadTriangle`
  + 3 helpers: `rank1EdgeQuad`, `pickRep`, `rank1Winner`; updated header docstring),
  `conjPieceCPLQTest.m` (renamed + rewrote the old nCE==2 test's docstring, added
  `psdRank1QuadraticConjugate` and `psdRank1QuadraticEndToEnd`), this handoff file.
- Not pushed (working tree has real uncommitted changes now, unlike the previous handoff's
  "no code changes" state).

## Next steps
- Ask the author whether to commit/push this session's work (standing rule: always ask).
- The degenerate rank-1 case (two vertices tying in rotated-t) -- decide if it's worth handling or
  documenting as a permanent limitation like the hyperbola case.
- Conjugate of RATIONAL pieces (the 1-convex-edge `RatPol` output of Step 1) is still the
  remaining big gap before the full `conjCPLQ` pipeline (Step1 `convEnvCPLQ` -> Step2
  `conjPieceCPLQ` per piece -> Step3 max) is wired end to end for nonconvex inputs. This needs
  genuinely curved (parabolic) edges in the OUTPUT too, not just the input -- `QuaPar`'s
  curved-edge `orderEdges` support (commit `ccbaba6`) is already there for this. Consider reading
  `biconjugate2/plq_1piece.m`'s type-1 branch (rational envelope conjugate, general polygon) in
  the cloned reference repo before re-deriving from scratch.
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
- `conjPieceCPLQ.m` -- this session's main file; new `conjPSDRank1QuadTriangle` + helpers near the
  end, dispatch branch near the top of `conjPieceCPLQ` itself.
- `conjPieceCPLQTest.m` -- `bilinearTwoConvexEdgesDocumentedLimitation`, `psdRank1QuadraticConjugate`,
  `psdRank1QuadraticEndToEnd`.
- `convEnvCPLQ.m` (`twoEdgeQuadPlain`, `envelopeFromClassified` case 2) -- Step 1's source of the
  rank-1 PSD quadratic that motivated this session's real work.
- `/home/ylucet/CCA2/DESIGN.md` -- design proposal (outside this repo).
- `/home/ylucet/CCA2/s10589-026-00781-5.pdf` (COAP) -- Appendix A (convex envelope, incl. A.4
  two-convex-edge harmonic mean) and Appendix B (conjugate, B.1/B.2 -- B.3's "two convex edges"
  section describes the case this session's `conjPSDRank1QuadTriangle` actually implements, though
  the derivation here was done independently rather than by transcribing B.3's prose).
- Reference repo (cloned to scratchpad this session, not part of this repo):
  `github.com/tanmaya11/convex`, branches `b2`/`b3`/`biconjugate`/`biconjugate2` -- the original
  (more general, polygon-based) research code; worth reading before tackling the rational-piece
  conjugate next.

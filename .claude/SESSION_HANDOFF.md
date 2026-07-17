# Session Handoff

_Last updated: 2026-07-17T19:48:35Z_

## What happened this session

Diagnosed the top-priority item from last session: the 15/741 (2.0%) silent wrong-answer cases
found in the prior session's `maxQuaPar` stress test. Reproduced the exact `rng(2026)`
1500-triangle dataset (724 ok / 2 crash / 15 wrong -- byte-for-byte matching last session's counts,
including the same crash triangles). Traced the root cause rigorously, via three independent
methods: (1) cross-checked `maxQuaPar`'s output against `conj(r)` computed by brute-force grid
sampling of `r`'s own raw quadratic pieces (bypassing `conjPieceCPLQ`/`maxQuaPar` entirely) --
matched to grid precision, proving Steps 2/3 are computing the conjugate of `r` correctly; (2)
verified `r` is a genuine globally convex minorant of `f=xy` (correct kink direction at the internal
seam, checked via gradient-jump sign); (3) verified tightness directly via the biconjugate identity
`env(f)(x0) = sup_s[s.x0 - f*(s)]`, using only the independently-correct closed-form
`supBilinearOverPoly` -- found a genuine, large gap (~2.14, not a numerical artifact) between `r`
and the true envelope at a concrete point. Conclusion: **all 15 wrong-answer cases are a Step 1
(`convEnvCPLQ`) envelope-tightness bug, not a `maxQuaPar`/`conjPieceCPLQ` bug** -- `maxQuaPar`
faithfully computes `conj(r)`, but `r` itself (`convEnvCPLQ`'s 2-convex-edge `splitTwoConvexEdges`,
12/15 cases, and the 3-convex-edge `splitThreeConvex` path, 3/15 cases -- though per the user,
the 3-convex-edge case is not actually distinct, since triangulating it reduces to the 2-convex-edge
case) is not always the tightest envelope, despite `DESIGN.md`'s own record of a prior session
believing this was already fixed (verified there via a 60-triangle stress test that happened not to
catch this).

Also verified: the assembled dual arrangement (`g1`, `g2`, and `h=maxQuaPar(g1,g2)`) for the
counterexample triangle **is** a genuine polyhedral complex (every pairwise face intersection is
empty, a shared vertex, or a full shared edge -- checked computationally via `shapely`, using the
actual sign-aware `polyConstraints`/`facePoly` logic from `maxQuaPar.m`, box-clipped for unbounded
faces) -- so the wrong-answer bug is cleanly separated from the historical "fan"/assembly-crash bug
class; assembly is working correctly here, Step 1's input to it is just wrong.

Cross-checked against the actual original reference implementation: the user pointed to
`github.com/tanmaya11/cPLQ`, cloned locally into `cPLQ/` (untracked, has its own `.git` -- not
committed into this repo). Confirmed its 2-convex-edge envelope formula
(`plq_1p.convexEnvelope1`, case `nCE==2`) is algebraically identical to this codebase's
`buildTwoEdge`/`twoEdgeQuadPlain` (no split logic exists there at all -- it's literally the paper's
Appendix A.4 single-quadratic formula, applied over the whole triangle). Ran its full pipeline
(`.triangulate` / `.maximum`) on one of the 15 wrong triangles and got a THIRD distinct value,
different from both truth and this codebase's result, with a large (~2.03), clean (non-degenerate)
discrepancy at a test dual point -- confirming the tightness gap is a genuine, pre-existing bug in
the published algorithm/original reference code itself, not something this codebase's later
`splitTwoConvexEdges` "fix" introduced or regressed.

Built two documentation figures (both new, committed) and one interactive artifact (published,
not part of the repo) illustrating the concrete counterexample
`T=(-9.95506,3.70366),(-9.345,-5.34552),(1.29049,5.31738)`:
- `doc/bug.svg` -- 2-panel: primal split diagram + a 1D cross-section showing `r` dipping 2.14
  below the true envelope at the gap point `x0=(-7.2947,-0.3917)`.
- `doc/bug_pipeline.svg` -- full 7-row pipeline (matching `doc/full_pipeline.svg`'s style/structure):
  primal domain -> split into `T1`/`T2` -> convex envelope of each (`q1`/`q_other`, with formulas)
  -> conjugates `g1`/`g2` (their own dual-space face subdivisions) -> `h=max(g1,g2)` (assembled
  14-cell arrangement) -> biconjugate of each piece (slice, confirms Step 2 correctness: `biconj(g1)
  = q1` on `T1`, `biconj(g2) = q_other` on `T2`, machine precision) -> overall biconjugate (slice,
  since no closed form is known -- this is the panel that actually shows the bug, comparing
  `biconj(f)` [the true envelope, independent of `r`] against `r`, NOT `conj(h)` [which would just
  trivially recover `r` itself via Fenchel-Moreau, since `r` is already convex -- would show nothing]).
- Interactive 3D artifact (published at
  `https://claude.ai/code/artifact/960a5c56-958d-411c-ba17-5ebd6d933ed4`, not in this repo): a
  custom canvas-based rotatable 3D plot (drag=orbit, scroll=zoom, matching MATLAB's `rotate3d`; no
  external libraries, per Artifact CSP) of `f=xy`, `r`, and `env(f)` over the triangle, camera
  pivoted on the gap point (the 2.14 gap is small next to the surface's ~87-unit height range, so a
  whole-triangle default view would hide it). Source file, generator script, and screenshots are in
  the OS-temp scratchpad (see below), not committed.

Finally, per the user's request, verified the gap via a SECOND, completely independent, more
elementary method (no conjugate/duality machinery at all): sampled ~14000 points of the graph
`(x,y,xy)` over `T` (dense barycentric grid + local refinement near `x0`), computed the 3D convex
hull (`scipy.spatial.ConvexHull` -- the lower-hull facets are exactly the graph of the convex
envelope, a standard fact), and evaluated the lower-hull facet at `x0`: **-1.811980**, matching the
biconjugate-based value (-1.811990) to 5 decimals, both confirming the same 2.143 gap against
`r=-3.955080`. The hull facet at `x0` has vertices at `B` and (essentially) `R` -- a concrete hint
that the true envelope's supporting plane there runs along the seam `R-B` extended, not along
`q_other`'s actual formula.

## Where things stand

- Branch: `main` @ `1695311` -- "Add figures documenting the convEnvCPLQ 2-convex-edge tightness
  bug". Previous commit `2cdc876` (last session's handoff) unchanged otherwise.
- Pushed: yes -- `2cdc876..eb18e94 main -> main`.
- Files changed/added this session: `doc/bug.svg`, `doc/bug_pipeline.svg` (both committed). No
  source `.m` files were modified this session -- this was pure diagnosis, no fix attempted (the
  user explicitly steered toward "diagnose and visualize," not "fix," given `DESIGN.md`'s own
  assessment that a real fix needs a genuine mathematical extension of Appendix A.4, not a quick
  patch).
- Untracked, not part of this repo: `cPLQ/` (clone of `github.com/tanmaya11/cPLQ`, has its own
  `.git`, left untracked deliberately -- see "What happened" above). Also not committed: the
  published interactive 3D artifact's HTML source and all of this session's scratch MATLAB/Python
  diagnostic scripts (`convEnvCPLQDebug.m`, `exportpipeline2.m`, `gen_3d_data.py`,
  `make_bug_pipeline.py`, `verify_hull.py`, etc.) -- all under the OS-temp scratchpad dir
  (session-specific, not persisted); easily regenerable from this file's description if needed
  again, but not guaranteed to survive to a future session.
- Full test suite: not re-run this session (no source changed).

## Next steps

- **Top priority, unchanged in kind from last session but now precisely scoped**: derive a
  genuinely tight convex-envelope construction for the 2-convex-edge case in `convEnvCPLQ.m`
  (`splitTwoConvexEdges`/`buildEdgeAffinePiece`/`seamPoint`). This is confirmed to be a real,
  open mathematical gap -- not a transcription bug relative to the published paper (the original
  reference code has the identical issue) and not something the existing split "fix" resolves in
  general (it was only verified against a 60-triangle sample that missed this counterexample). The
  hull-facet hint above (true envelope's support at `x0` passes through `B` and `R`, not through
  `q_other`) is a concrete lead for the next attempt, but no fix was attempted this session.
- Regenerate reproducible artifacts as needed: the counterexample triangle and its dual-space
  objects (`T1`, `T2`, `q1`, `q_other`, `g1`, `g2`, `h`) are fully described in this file and in the
  git history's own commit message / `bug_pipeline.svg`'s captions; the exact `rng(2026)` stress
  test (1500 triangles, same seed) reproduces the same 15 wrong-answer / 2 crash triangles if a
  broader sweep is wanted again.
- Still open, unchanged from before: the 2/741 (0.3%) residual `maxQuaPar:internal` crashes
  (`T=[-7.71708 -0.856785;-8.22504 0.349933;-7.27519 1.40742]` and `T=[3.9615 8.03275;2.94803
  -9.30541;-7.05705 6.43129]`) -- lower priority than the tightness bug above, not investigated
  this session.
- Update `DESIGN.md`'s `convEnvCPLQ.m`/`maxQuaPar.m` Implementation-status write-ups with this
  session's diagnosis (the "Important: this is NOT a correctness gap in the underlying math"
  passage there is now known to be wrong in general -- it IS a correctness gap, just isolated to
  Step 1, not Step 2/3 as that passage's surrounding context might suggest to a future reader) --
  not done this session.
- Lower priority / untouched, unchanged from before: `QuaPar.orderEdges`/`createP`'s separate
  "Face k has ... expected 2 but got N" error on near-degenerate/thin triangles; `partialConj` for
  `'cplq'`/`'pqp'`; `add` for `RatPol`/`RatPar`; conjugate engines `'pqp'`/`'graph'`; the standalone
  `RatPol.conj` gap.

## Relevant files

- `convEnvCPLQ.m` -- `splitTwoConvexEdges`/`buildEdgeAffinePiece`/`seamPoint`: confirmed this
  session to still admit non-tight envelopes on some 2-convex-edge triangles (counterexample
  T=(-9.95506,3.70366),(-9.345,-5.34552),(1.29049,5.31738), gap ~2.14 at x0=(-7.2947,-0.3917));
  not modified this session.
- `maxQuaPar.m` -- confirmed (not modified) to correctly compute `conj(r)` and to assemble a
  genuine polyhedral complex for this counterexample; the wrong answer is entirely attributable to
  its input `r`, not to this file's own logic.
- `doc/bug.svg`, `doc/bug_pipeline.svg` -- new this session, both illustrate the counterexample
  above; see their own captions for the full pipeline breakdown.
- `doc/full_pipeline.svg`, `doc/cplq_split.svg`, `doc/cplq_unsplit.svg` -- pre-existing sibling
  figures `bug_pipeline.svg` was styled to match.
- `cPLQ/` -- untracked clone of the original `github.com/tanmaya11/cPLQ` reference implementation;
  useful for cross-checking whether a future fix attempt's target behavior matches (or improves on)
  the original algorithm, not just this codebase's own history.
- `DESIGN.md` -- has the fullest prior history of the 2-convex-edge tightness saga (search for
  "2-convex-edge" and "Open research question"); NOT updated this session with this session's
  findings (see "Next steps").

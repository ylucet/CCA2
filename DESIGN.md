# CCA2 — Code Design Proposal for a 2D PLQ Computational Convex Analysis Toolbox

**Date:** 2026-06-26  ·  **Scope:** pure MATLAB (toolboxes allowed: Symbolic Math,
Optimization). 2D primary, 1D as a thin degenerate case. General PLQ (occasionally PLC),
**no grid/LLT algorithms**.

**Authoritative references** (PDFs in `reference/`):
- **[COAP]** Karmarkar & Lucet, *Computing the convex envelope of bivariate PLQ functions
  in linear time*, Comput. Optim. Appl. 94 (2026) 747–780. (`KARMARKAR-26-convex-envelope.pdf`)
- **[JOGO]** Karmarkar & Lucet, *A linear-time algorithm to compute the conjugate of
  nonconvex bivariate PLQ functions*, J. Glob. Optim. 94 (2026) 3–34. (`KARMARKAR-26-conjugate.pdf`)

These two papers define the data types and the algorithms this design implements. **Both handle
nonconvex PLQ** (the `'cplq'` engine, II.5.1); they build on a lineage of predecessor work that
handles the **convex** case only:
- **[KUMAR-20]** Kumar & Lucet, *Towards the Biconjugate of Bivariate Piecewise Quadratic
  Functions*, in *Optimization of Complex Systems*, Adv. Intell. Syst. Comput., Springer (2020)
  257–266. **Deepak Kumar**'s numeric convex-envelope-of-a-quadratic-over-a-polytope method
  (`codeOld/deepak`, I.4) is the direct predecessor that [COAP]/[JOGO] extend — it is **not** a
  parametric-QP method (a previous version of this document conflated the two; see I.4/I.5 for
  the correct attribution below).
- **[JAKEE-13]** Khan Md. Kamall Jakee, *Computational Convex Analysis using Parametric
  Quadratic Programming*, M.Sc. thesis, University of British Columbia (Okanagan), 2013.
  (`JAKEE-09-pQP-MSc.pdf`) — **Jakee Khan**'s work; the actual source of the **parametric-QP**
  (`'pqp'`) conjugate engine (II.5.2). Convex PLQ only.
- **[GARDINER-09]** Gardiner & Lucet, *Numerical Computation of Fitzpatrick Functions*,
  J. Convex Anal. 16 (2009) 779–790. (`GARDINER-09-Fitzpatrick.pdf`) — source of the
  Fitzpatrick-function routine kept as an example generator (I.6).
- **[GARDINER-11]** Gardiner & Lucet, *Graph-Matrix Calculus for Computational Convex
  Analysis*, in *Fixed-Point Algorithms for Inverse Problems in Science and Engineering*,
  Springer Optim. Appl. (2011) 243–259. (`GARDINER-11-graph-matrix.pdf`)
- **[GARDINER-13]** Gardiner & Lucet, *Computing the Conjugate of Convex Piecewise
  Linear-Quadratic Bivariate Functions*, Math. Program. 139 (2013) 161–184.
  (`GARDINER-13-conjugate-PLQ.pdf`) — convex-only bivariate conjugate (`lft2.sci`, I.5).
- **[GARDINER-14]** Gardiner, Khan & Lucet, *Computing the Partial Conjugate of Convex
  Piecewise Linear-Quadratic Bivariate Functions*, Comput. Optim. Appl. 58 (2014) 249–272.
  (`GARDINER-14-partial-conjugate.pdf`) — source of `partialConj` for `'pqp'` (II.4/II.5.2).
- **[HAQUE-17]** Haque, Tasnuva, *Computation of Convex Conjugates in Linear Time using
  Graph-Matrix Calculus*, M.Sc. thesis, University of British Columbia (Okanagan), 2017.
  (`HAQUE-17-conjugate-convex-PLQ-MSc.pdf`)
- **[HAQUE-18]** Haque & Lucet, *A Linear-time Algorithm to Compute the Conjugate of Convex
  Piecewise Linear-Quadratic Bivariate Functions*, Comput. Optim. Appl. 70 (2018) 593–613.
  (`HAQUE-18-conjugate-convex-PLQ.pdf`) — **Tasnuva Haque**'s work; the algorithm the
  point-cloud + neighbour-graph (`'graph'`) engine implements (II.5.3). Convex PLQ only.
- **[HIRIART-URRUTY-07]** Hiriart-Urruty & Lucet, *Parametric Computation of the
  Legendre–Fenchel Conjugate with Application to the Computation of the Moreau Envelope*,
  J. Convex Anal. 14 (2007) 657–666. Source of the "expand-the-square" identity that computes
  the Moreau envelope as a **single conjugate**, not an inf-convolution (II.6) — valid for
  **any** function, convex or not.
- **[LOCATELLI]** Locatelli, M. — proved that the convex envelope of a quadratic function over a
  polyhedron always admits a POLYHEDRAL subdivision (straight-edge domain partition). Full
  citation not yet filled in here (recorded 2026-07-17 per the user, not independently looked up
  this session) — fill in before citing in the paper. Directly relevant to `convEnvCPLQ.m`'s
  2-convex-edge case (see its Implementation-status entry below): a session's numerical
  exploration of the `q1`/Appendix-A.3 pairing suggested a curved transition boundary, which per
  this theorem must be an artifact of an incomplete search (the wrong split/pairing), not a
  genuine counterexample — the correct fix is a DIFFERENT polyhedral split, not a curved-domain
  extension of `RatPol`.

---

## Implementation status (as of 2026-07-13)

This document is a **design proposal**; large parts of it describe intended, not yet built,
code. Cross-check the actual repo file list before assuming an operator or engine exists.

**Implemented** (MATLAB, in the repo root):
- Classes `QuaPoly` (leaf/input type, formerly `PLQVC`), `QuaPar` (conjugate `f*`), `RatPol`
  (convex envelope `f**`) — each its own `classdef`, **not yet organized under a common `RatPar`
  parent** (II.3's class hierarchy is proposed, not built).
- Conjugate engine **`'cplq'` only** (`conjCPLQ.m`, `conjPieceCPLQ.m`, `convEnvCPLQ.m`) —
  incremental; each file's own header STATUS block lists exactly which piece-classification
  cases are covered so far.
- `maxQuaPar.m` — pointwise max of two full-domain `QuaPar` objects (needed when Step 1 splits a
  nonconvex piece into more than one sub-piece). Now wired into `conjCPLQ`'s own Step 3
  (`conjMaxOfSubTriangles`, in `conjCPLQ.m`): a single bounded triangle whose indefinite quadratic
  has 3 convex edges — the one case Step 1 (`convEnvCPLQ`) splits into more than one sub-piece —
  conjugates each of the two sub-triangles (Step 2, always exactly 2-convex-edge/rank-1-PSD by
  COAP Appendix A.5, so always polyhedral-domain-safe for `maxQuaPar`) and combines them via
  `maxQuaPar`. Exercised end to end (not just the pre-existing manual `conjPieceCPLQ`+`maxQuaPar`
  wiring in `maxQuaParTest.m`) by `conjCPLQTest.m`'s
  `indefiniteTriangleThreeConvexEdgesUsesStep3`.
  **Correctness fix (2026-07-14)**: `convEnvCPLQ.m`'s `splitThreeConvex` used to split the
  3-convex-edge triangle by a HORIZONTAL line through the middle vertex — wrong in general (only
  correct when the two sub-envelopes happen to be mirror-symmetric). The two sub-envelopes q1,q2
  both touch `u1*u2` exactly along the entire third (shared) edge, forcing `q1-q2` to factor as
  (that edge's line) × (a second line through the middle vertex); THAT second line, not a
  horizontal one, is the actual smooth-fit split direction (closed form in `splitThreeConvex`'s
  header). The old horizontal split left `q1`,`q2` agreeing only at the two seam endpoints, not
  along the interior of the seam — a small, silent discrepancy for the one hand-worked example
  used everywhere in tests, but large enough for other triangles to make the glued 2-piece
  "envelope" genuinely non-convex (a real jump discontinuity across the interior seam) and its
  `maxQuaPar` conjugate numerically too high. Fixed; re-verified against ground truth to machine
  precision on the previously-wrong triangle and on a 600-triangle randomized stress test (0 wrong
  answers, only pre-existing/separate `maxQuaPar:internal` assembly-topology crashes on ~4% of
  valid samples — a different, still-open gap, not a correctness issue). See `splitThreeConvex`'s
  HISTORY comment and the session handoff for the derivation and stress-test details.
  **Reliability fix (2026-07-14, later session)**: the `maxQuaPar:internal` assembly-topology
  crash noted above (`assemblePieces: boundary edge (r,c) has no matching neighbour`) was rooted in
  `assemblePieces` merging all pieces' vertices via one coordinate-distance tolerance and then
  matching half-edges by resulting vertex-index equality — unsound for near-degenerate triangles,
  since a genuine cross-arithmetic-noise gap and a genuine distinct-nearby-vertex gap can be the
  same order of magnitude (~1e-5), so no single tolerance separates them. Fixed by matching
  half-edges directly by geometry instead (comparing only DIFFERENT pieces' edges, never a piece's
  own), then deriving global vertex identity afterwards via union-find restricted to exactly those
  confirmed matches — see `maxQuaPar.m`'s header HISTORY and `assemblePieces`' own HISTORY comment
  for the full derivation. Crash rate on a ~5000-triangle randomized stress test dropped from
  ~4/800 valid samples to ~1/800, with zero new wrong-answer regressions (verified: every
  wrong-answer case reproduces identically on the unmodified pre-fix code). The residual ~1/800
  case is a genuinely ambiguous 3-way vertex cluster that needs vertex PROVENANCE (tracking which
  original `g1`/`g2` face boundary a vertex came from) to resolve, not just tighter geometry —
  still open, see the session handoff.
  **Reliability fix (2026-07-14, later session)**: diagnosed the residual ~1/800 case above and
  found full vertex provenance turned out to be unnecessary. Root cause: a 3-way (or more) vertex
  cluster where several pieces meeting near one point each compute a slightly different position
  (order ~1e-4 -- too coarse for `matchHalfEdges`' cross-arithmetic-noise tolerance, too fine to be
  a distinct feature) for what is mathematically ONE vertex; the best-first greedy matcher can
  pair up only 2 of the 3-plus mutually-close half-edges, always orphaning one. Diagnosis showed
  the orphaned edge's own two endpoints ALWAYS resolve to the very same global vertex via the
  OTHER, already-confirmed matches on its own piece's boundary -- i.e. it is provably zero-length
  once the rest of the topology is resolved, so it is safe to just drop (emit nothing for it)
  rather than error. Fixed via the new `checkOrphanHalfEdges` (called after global vertex identity
  is built): drops an orphaned edge only when its two endpoints already coincide globally,
  otherwise still raises the original error (a genuinely unresolved gap, or an orphaned ray -- no
  evidence yet that rays can be legitimately degenerate this way). Verified on 5 of 6
  randomly-found repro triangles end-to-end against ground truth (machine precision); crash rate
  on the "no matching neighbour" error specifically dropped to 0/455 valid samples on a 3000-
  triangle randomized stress test (down from ~4-6/800 before this session). See
  `checkOrphanHalfEdges`'s own header in `maxQuaPar.m` and
  `maxQuaParTest.checkOrphanHalfEdgesDropsProvablyDegenerateOrphanEdges`.
  **New issue found while verifying the fix above (still open, higher priority than it looks)**:
  one of the 6 repro triangles (a 5-piece cluster, more complex than the simple 3-way case) no
  longer crashes but instead SILENTLY returns the wrong answer (`Inf`, i.e. a domain-coverage gap)
  over a substantial region -- confirmed via `polyConstraints` that the original, pre-assembly
  piece geometrically DOES cover that region, so the bug is in final assembly/`QuaPar.eval`, not
  in piece generation. This is NOT caused by this session's fix (the crash previously masked it
  for this exact triangle), and is NOT the same bug as the vertex-cluster crash above -- likely a
  sign/orientation issue in `QuaPar.eval`'s per-edge conic membership test (`edgeConics`/
  `evalConic`) specific to an unbounded face whose two rays point in nearly the SAME direction (a
  degenerate "parallel-strip" shape), not yet isolated further. A broader stress test (2000
  triangles) found this silent-wrong-answer pattern in ~2/300 valid samples -- roughly the same
  order of magnitude as the crash rate this fix addressed, and worse in kind (silent, not loud).
  **This should be the top priority for the next session** -- see the session handoff for the
  repro triangle and diagnostic trail.
  **Correctness fix (2026-07-14, later session)**: diagnosed and fixed the silent-wrong-answer
  issue above. Root cause: `insertPassthroughVertices` (which re-inserts an original `g1`/`g2` face
  vertex into a clipped cell when it lies in the open interior of one of the cell's own edges, so
  the true neighbouring face never silently changes mid-edge -- see `clipByFace`'s HISTORY) checked
  "is this point ALREADY a vertex of the cell" using the SAME tight tolerance (1e-7) as
  `onOpenSegment`/`onOpenRay`'s own edge-matching test. A genuine cross-arithmetic-noise gap between
  the original face vertex and the already-present cell vertex it geometrically coincides with can
  be ~3e-5 (two different formulas computing "the same" point, the same noise floor documented
  elsewhere in this file) -- too coarse for that shared 1e-7 -- so the point narrowly failed the
  "already present" check and was instead inserted as a brand-new, near-duplicate vertex, creating a
  near-zero-length sliver edge whose line equation was dominated by floating-point noise in the tiny
  direction vector. `QuaPar.eval`'s exact (no-tolerance) membership test then wrongly excluded a
  real region of the plane behind that noise-direction line. Fixed by giving the "already a vertex"
  pre-check its OWN, wider tolerance (`tolSnap=1e-4`), decoupled from `onOpenSegment`/`onOpenRay`'s
  matching tolerance (left at 1e-7): widening the shared tolerance instead (tried first) broke two
  OTHER regression tests, each at a different noise scale (a genuinely distinct ~7.9e-4-away vertex
  wrongly absorbed into the wrong edge in one; a genuine ~1e-4 feature of a near-degenerate triangle
  wrongly merged in the other) -- no single shared value separated all three cases, but `tolSnap`
  only ever recognizes p as coincident with a vertex the cell construction already produced, never
  changing whether a genuinely new point gets inserted, so it can't manufacture that kind of wrong
  topology. See `insertPassthroughVertices`'s own header in `maxQuaPar.m` and
  `maxQuaParTest.insertPassthroughVerticesDropsNearDuplicateCrossingPoint`. A 3000-triangle
  randomized correctness stress test found 0/491 silent-wrong-answer (`Inf`) cases after this fix
  (down from ~2/300 before it); a handful of small-magnitude (not `Inf`) wrong-answer cases also
  turned up in that same stress test and were confirmed to reproduce identically on the unmodified
  pre-fix code -- a separate, still-open issue, diagnosed (not fixed) in the next paragraph.
  **Open research question found while diagnosing the small-magnitude wrong-answer cases
  (2026-07-14, later session; higher priority than it looks -- affects the published paper's own
  worked example)**: traced ALL of the small-magnitude wrong-answer cases to `convEnvCPLQ.m`'s
  2-convex-edge quadratic (`envelopeFromClassified` case 2 / `twoEdgeQuadPlain` / `buildTwoEdge`,
  implementing [COAP] Appendix A.4) NOT always being the true tightest convex envelope, even though
  it is a valid minorant that correctly touches `f=u1u2` along both classified convex edges.
  Reproduced and confirmed via three independent methods on **the paper's own Appendix A.4.3
  example**, `V=(2,1),(0,0),(1,0)`  (the same triangle used verbatim in
  `convEnvCPLQTest.bilinearTwoConvexEdgesQuadratic`): (1) hand-derived exact closed-form boundary
  analysis gives `sup_{x in V}[s.x-f(x)] = 0` at `s=(-0.008727,-0.999962)`, achieved at vertex
  `(0,0)`; (2) the REAL pipeline (`convEnvCPLQ` then `conjPieceCPLQ`, i.e. `gTest.eval(s)`) gives
  `0.03864091` instead; (3) an independent biconjugate reconstruction (
  `conv(f)(x0) = sup_s[s.x0 - supBilinearOverPoly(s,V)]`, solved numerically) confirms the TRUE
  envelope value at the discrepancy's own argmax point `x0=(0.474343,0)` is `~0` (matching `f`
  exactly there), while the formula's `q1(x0) = -0.042780` -- a genuine, reproducible ~0.043 gap,
  not a numerical artifact (also reproduced on the session's original stress-test triangle
  `T=(3.8398,5.0413),(8.8152,7.2338),(8.7969,5.6447)`, gap ~0.11 at its own analogous point).
  Root cause: `x0` lies on the triangle's THIRD, non-convex ("weak") edge, where `twoEdgeQuadPlain`
  is only constrained to match `f` at the 2 shared endpoints (not the whole edge) -- being a single
  convex quadratic forced through those 2 endpoint values, it is mathematically forced to dip
  strictly below `f` in the open interior of that edge, more than the TRUE 2D envelope needs to.
  Established (via the same biconjugate check, at several points) that the TRUE envelope along the
  weak edge equals the AFFINE CHORD between its 2 endpoint values (matches to 6 decimals at every
  point checked) -- but naively using `max(q1, affine-interpolation-of-the-3-vertex-values)` as a
  replacement is NOT globally valid: the affine piece exceeds `f` (and `q1`) in a region nearer the
  common vertex (by up to 1.46 in the stress-test triangle), and a systematic re-check over many
  random dual points found `s` values where that combination's conjugate is wrong in the OPPOSITE
  direction too. A correct general fix likely needs a genuine sub-partition of the 2-convex-edge
  triangle (an actual second piece with its own derivation, analogous to how the 3-convex-edge case
  above splits into two sub-triangles) rather than a simple pointwise max of two candidates -- this
  is a real extension of the published Appendix A.4 derivation, not a simple code fix. Consistent
  with the paper's own text for its worked example -- *"by checking the bounds we get the domain as
  the entire triangle"* -- which reads as a check that can fail for other triangles, not a
  universal guarantee; the code never implements that check, so it silently assumes the single
  quadratic is always valid.
  **Fixed (2026-07-15, later session)**: derived and implemented the genuine sub-partition
  anticipated above (`convEnvCPLQ.m`'s `splitTwoConvexEdges`/`buildEdgeAffinePiece`/`seamPoint`).
  Let `P` be the two convex edges' shared vertex and `A,B` the weak edge's endpoints (edges `P-A`,
  `P-B` convex). The triangle is split by a cevian from ONE of `A,B` into the OPPOSITE convex edge,
  giving two sub-triangles: the one containing `P` keeps `q1=twoEdgeQuadPlain(...)` UNCHANGED (it
  is still exactly tight there, touching `f` along both of ITS two convex-edge sub-segments); the
  other uses a NEW quadratic (`buildEdgeAffinePiece`) touching `f` along its remaining convex-edge
  segment and the affine chord along the (now fully-contained) weak edge. Derivation: matching a
  quadratic to `f` along one full line and to a DIFFERENT affine target along another full line is
  a linear system in the quadratic's 6 coefficients that is rank-deficient by exactly 1 (a pencil,
  like `twoEdgeQuadPlain`'s own derivation) -- but unlike `twoEdgeQuadPlain`'s two distinct rank-1
  (Hessian-singular) solutions needing a +/- branch choice, this pencil's rank-1 condition is a
  DOUBLE root (tangent), always picking a UNIQUE convex quadratic (proved symbolically with
  `sympy` and confirmed numerically). The cevian's direction (which of `A,B` to cevian from, and
  where it lands on the opposite convex edge) is forced, not free -- exactly analogous to
  `splitThreeConvex`: `q1` and the new piece both touch `f` along the SAME original convex edge
  being split, so their difference vanishes identically along that edge's line, hence factors as
  (that line) * (the seam's line); intersecting the seam's line with the target edge gives the
  split point. Of the two candidate cevians (from `A` into `P-B`, or from `B` into `P-A`), exactly
  one lands strictly inside its target edge -- verified never both, never neither, across ~6500
  random 2-convex-edge triangles this session, and used as the runtime selection rule (mirroring
  how `twoEdgeQuadPlain` already picks its own branch by validity). A special (measure-zero) case
  needs no split at all: if `q1`'s curvature along the weak edge is already zero (e.g. the mirror-
  symmetric triangle `(0,0),(2,1),(1,2)`, where `mh*mw=1`), `q1` already equals the chord along the
  WHOLE weak edge, so it is already the tight envelope everywhere (both candidate cevians degenerate
  to parallel/at-infinity in exactly this case, giving a direct, cheap detection). Verified via: (1)
  the paper's own Appendix A.4.3 example (`V=(2,1),(0,0),(1,0)`) -- the weak-edge dip point
  `x0=(0.474343,0)` now gives `~0` (was `-0.042780`), and the flagged bad dual point
  `s=(-0.008727,-0.999962)` now gives the correct `sup=0` via each sub-triangle's own conjugate
  (`conjPieceCPLQ`), matching `maxQuaParTest`-style ground truth (`supBilinearOverPoly`) to machine
  precision; (2) a ~60-triangle MATLAB-side randomized stress test comparing `max` of the two
  sub-triangles' own `conjPieceCPLQ` conjugates against the same exact ground truth, at 5 random
  dual points each: max error `5e-14` (machine precision, not the `~1e-6`-`1e-8` optimizer noise of
  this session's earlier from Python-side numerical verification) across every case that produced a
  genuine split. Regression test: `convEnvCPLQTest.bilinearTwoConvexEdgesSplitIsTight`. Full suite:
  146/146 pass (145 prior + 1 new; two prior tests -- `conjCPLQTest.
  indefiniteTriangleTwoConvexEdgesSidestepsToEnvelope` (uses the no-split-needed symmetric
  triangle, so unaffected) and `conjPieceCPLQTest.psdRank1QuadraticEndToEnd` (updated to extract
  one sub-triangle face from Step 1's now-2-face output, since it specifically tests Step 2 on a
  single rank-1 PSD triangle) -- needed no or minor updates, both documented at the call site).
  **New downstream gap found as a direct consequence (partially diagnosed and partially fixed a
  later same session, NOT fully resolved -- see below)**: `conjCPLQ`'s existing Step-3 fallback
  (`conjMaxOfSubTriangles`, which already combines a multi-face Step-1 envelope via per-face
  `conjPieceCPLQ` + `maxQuaPar`, exactly as needed here) now gets exercised for the 2-convex-edge
  case too, and `maxQuaPar` frequently (about half of random split cases, ~52% in a 138-triangle
  stress test) errors combining the two sub-triangles' conjugates, with `maxQuaPar:internal` (an
  unmatched boundary half-edge in `assemblePieces`) or occasionally `:notDegenerate`/
  `:notImplemented`. Likely root cause: `splitThreeConvex`'s sub-triangles are explicitly
  constructed so their conjugates paste together SMOOTHLY (C1) along the shared seam (see that
  function's own HISTORY) -- an assumption `maxQuaPar`'s `splitCell`/`assemblePieces` machinery
  leans on (per its own comments, e.g. "this should never happen for two conjugates of adjacent
  sub-pieces of the same originally-nonconvex domain"). `splitTwoConvexEdges`'s cevian, by
  contrast, only guarantees the two PRIMAL pieces agree in VALUE along the seam (proved tight
  against exact ground truth above), with no claim of matching gradients -- confirmed directly (by
  evaluating both pieces' gradients at several points along the seam on the paper's own example):
  gradients match at one seam endpoint (where the two original convex edges happen to already be
  tangent) but diverge increasingly toward the other, a genuine, continuously-varying kink, unlike
  the always-smooth 3-edge case.
  **Important: this is NOT a correctness gap in the underlying math.** `h(s) = max(g1(s),g2(s))` is
  *always* exactly `f*(s)` for `f` split into two pieces on a domain partition, by the elementary
  identity `sup_{x in A union B} phi(x) = max(sup_{x in A} phi(x), sup_{x in B} phi(x))` -- true for
  ANY partition, kinked or smooth. Verified directly: extracting each sub-triangle's own
  `conjPieceCPLQ` conjugate and taking a plain `max` of the two (bypassing `maxQuaPar`'s assembly
  entirely) matches `supBilinearOverPoly` ground truth to `5e-14` (machine precision) across every
  split case in a 60-triangle stress test, with zero exceptions. The gap is *entirely* in
  `maxQuaPar`'s ASSEMBLY: building one well-formed `QuaPar` (a clean partition of the dual plane,
  every boundary edge with exactly one neighbour) when combining `g1`,`g2` whose seam is kinked, not
  in any value `maxQuaPar` DOES manage to produce -- confirmed by re-running the same stress test
  end-to-end through `conjCPLQ`: 0 mismatches among the ~48% of cases that succeeded.
  Diagnosed root cause (this session): a kink along a shared primal seam gives `g1` and `g2` a
  genuine **positive-area TIE region** in dual space (`g1(s)=g2(s)` over a whole 2D set, not just a
  curve) -- expected whenever a seam's two adjacent pieces share a VERTEX (e.g. the seam's own
  endpoints), since both `g1`'s and `g2`'s OWN face structure independently produce a valid
  vertex-cone face there, and (unlike the smooth 3-edge case, where `maxQuaPar` had only ever been
  tested against) these two vertex-cones can legitimately overlap or sit adjacent to several of the
  OTHER side's faces at once (a "fan" of faces meeting at the dual image of the shared vertex).
  `maxQuaPar`'s per-`(k,l)`-pair double loop clips each `g1` face against each `g2` face
  independently and expects the results to tile the plane cleanly; when two DIFFERENT `(k,l)` pairs
  produce the exact same geometric cell (one instance of the tie phenomenon above), `assemblePieces`
  sees two competing claims on the same boundary, orphaning one copy's edges.
  **Partial fix implemented**: `dedupPieces` (new, in `maxQuaPar.m`, called right before
  `assemblePieces`) detects groups of pieces with IDENTICAL geometry (same vertex count, same real
  vertices as a set, same ray directions, no curved edge, tolerance `1e-6`) produced by different
  `(k,l)` pairs, and collapses each group to one: if every member agrees on which row (`f1`/`f2`)
  wins, keeps one copy; if they disagree (each drawn from a DIFFERENT, only partially-applicable
  candidate face), reconciles by evaluating every candidate row at a point verifiably interior to
  the shared cell and keeping the largest -- sound because both `g1`'s and `g2`'s own rows are
  independently exact everywhere (per the machine-precision check above), so whichever is larger at
  an interior point is provably correct for the whole cell (the same premise `decideWinner` already
  relies on). Verified this correctly detects and reconciles the EXACT-duplicate case (e.g. the
  paper's own example, `V=(2,1),(0,0),(1,0)`, where two `(k,l)` pairs produce an identical cell with
  DIFFERING winners, one right and one wrong by construction) and causes zero regressions (146/146
  suite still passes). **However, empirically (138-triangle stress test, with and without
  `dedupPieces`) it never independently flips a single failing triangle to success** -- every
  triangle exhibiting the exact-duplicate pattern ALSO exhibits a second, more complex pattern:
  several small single- or two-vertex cone/strip faces near the dual image of the kink's vertex
  whose RAYS fail to pair with any counterpart (`matchHalfEdges` finds no candidate within
  tolerance), rather than pairing with an EXACT duplicate. This looks like the same underlying
  "positive-area tie / fan of faces" phenomenon but manifesting as adjacent-but-not-identical
  cells (not caught by `dedupPieces`'s exact-match criterion) rather than as literal duplicates --
  not yet understood well enough to fix. Concrete repros for continued work (both from the
  `rng(12345)` stress test in this session's scratch `stress2edge.m`, NOT committed -- regenerate
  by sampling random triangles and keeping ones where `classifyConvexEdges` gives exactly 2 rows):
  `V=[2 1;0 0;1 0]` (the paper's own example) and `V=[1.508518 2.818371; 2.687354 4.499057;
  -1.870671 3.524095]`; both fail with `maxQuaPar:internal`, "piece 1" always the specific piece
  with an unmatched edge in these traces. Recommended direction for a full fix: give the shared
  seam-endpoint vertex's dual-side "fan" first-class handling BEFORE the general `(k,l)` double
  loop (e.g. explicitly detect, for each vertex `g1`/`g2` inherit from a shared primal point,
  whether their own vertex-cone/edge-strip faces there overlap or complement, and merge/resolve
  before the generic clip loop runs) -- likely needs the vertex-PROVENANCE tracking this file's own
  HISTORY comments already anticipate as a tool of last resort for a different (now-resolved)
  ambiguity, tying each dual face back to the specific primal vertex/edge that produced it. This is
  a real, fairly common gap for `conjCPLQ`/`maxQuaPar` (~52% of random 2-convex-edge splits), but it
  fails LOUDLY and cleanly (never a silent wrong answer, confirmed above) -- a substantial follow-up
  in its own right, not a blocker for Step 1's correctness fix, which stands independently verified.
  **Follow-up session (2026-07-16): two more concrete bugs found and fixed within the "fan"
  phenomenon, both in `maxQuaPar.m`, neither sufficient on its own.** Diagnosed by reproducing the
  paper's own `V=[2 1;0 0;1 0]` example end to end and directly inspecting the (normally internal)
  `pieces` list `assemblePieces` builds just before failing, then testing pairs of pieces for
  literal point-set overlap (not just relying on the crash message) via a coverage/sampling check:
  1. **`matchHalfEdges`' ray matching accepted a pairing on (apex, direction) equality alone, with
     no check that the two candidate pieces are on OPPOSITE sides of that ray.** When a g1 (or g2)
     face is cut into several `(k,l)` sub-pieces by different opposing faces, EVERY sub-piece
     independently inherits the parent face's own boundary rays as part of its own boundary --
     not because each is adjacent to some other sub-piece there, but simply because they all
     descend from the same ancestor. The old code paired up the first two such inheritors it found
     as if they were mutual neighbours; confirmed via direct point-sampling that two of this
     triangle's pieces (built from the SAME `g1` face, different `g2` faces) genuinely overlapped
     over a positive-area region, not merely touched a shared boundary -- exactly the "fan" pattern
     described above, but now pinned to a specific, fixable root cause rather than only its symptom.
     **Fixed** by `oppositeSides` (new, called from `matchHalfEdges`): each ray candidate's own
     adjacent geometry (the next real vertex, or the OTHER ray's direction for a pure 2-ray cone) is
     tested via the 2D cross product against the shared ray direction; two candidates are accepted
     as genuine twins only when their representative points fall on opposite sides. Verified
     directly: the previously-overlapping pair is no longer paired, and the piece that used to
     wrongly claim that neighbour's territory now correctly finds ITS true opposite neighbour
     instead (confirmed by re-inspecting the post-fix `pieces`/half-edge list by hand).
  2. **A related but distinct redundancy**: two different `(k,l)` pairs sharing the identical
     winning row `f` can produce pieces that are NOT exact geometric duplicates (`dedupPieces` only
     collapses those) yet still overlap, one wholly containing the other -- pure redundant coverage
     of territory the larger piece already legitimately claims. **Fixed** by `dropSubsumedPieces`
     (new, called from `maxQuaPar` right after `dedupPieces`): for every pair of pieces agreeing on
     `f`, drops whichever's every real vertex (and, if unbounded, both recession directions) lies
     inside the other's own half-plane/ray constraints (reusing `polyConstraints`). Confirmed on the
     same example: piece count drops from 12 to 10 with zero change in any resolved value.
  Both fixes are individually correct (each addresses a genuine, hand-confirmed geometric defect,
  not a heuristic tolerance tweak) and cause zero regressions: full suite still 147/147 (146 prior +
  1 new, `maxQuaParTest.matchHalfEdgesRejectsSameSideRayPairingAndDropsSubsumedPieces`, which pins
  the CURRENT documented failure mode below since neither fix alone resolves it). **Important,
  initially-surprising finding, established via a properly-controlled experiment**: neither fix
  moves the observed failure rate on GENERIC random 2-convex-edge triangles at all. A same-triangle-
  set A/B comparison (pre-generating every triangle candidate AND every dual test point up front, so
  the two runs -- with vs. without both fixes -- are byte-for-byte reproducible regardless of which
  triangles happen to succeed or fail; an earlier, naive attempt at this comparison that drew random
  numbers lazily inside the success/failure branches silently desynchronised the two runs' RNG
  streams after the first behavioural difference, producing a spurious-looking improvement that
  evaporated once corrected) across 1494 randomly generated 2-convex-edge triangles found **zero**
  triangles whose pass/fail outcome changed in either direction (931/1494, 62.3%, succeed with or
  without the fixes). So: both fixes are real, confirmed, worth keeping (the overlap bug in
  particular is a genuine correctness hazard, not just a crash-avoidance nicety, even though it
  happened not to flip any outcome in this sample), but the ~38% aggregate crash rate on generic
  triangles is driven by a THIRD, still-undiagnosed pattern within the same "fan" phenomenon, not by
  either of these two. Confirmed the paper's own example (and, unchecked but likely, the session's
  other repro `V=[1.508518 2.818371; 2.687354 4.499057; -1.870671 3.524095]`) still throws
  `maxQuaPar:internal` after both fixes: with `dropSubsumedPieces` active it fails one piece later in
  the loop (piece 1's edge instead of piece 3's) but the underlying obstruction is unchanged -- a
  small cluster of pieces near the dual image of the shared primal seam-vertex (this example: primal
  vertex `A=(0,0)`) still fails to close up combinatorially. Recommended direction unchanged from the
  diagnosis above (vertex-provenance-aware fan resolution); one attempted-and-reverted approach worth
  recording so it isn't retried: broadening `insertPassthroughVertices`' candidate points from just
  the current cell's own two parent faces (`polyK`/`polyL`) to EVERY vertex of both full inputs
  `g1`/`g2` (reasoning: a missing "T-junction" split could in principle come from a third, unrelated
  face) neither fixed the repro triangle nor helped the aggregate rate, and caused 7 new regressions
  in the existing suite -- reverted; the missing structure is not a simple T-junction of this kind.
  The fully general case -- a multi-face original
  domain (`nf>1`), or a single non-triangular face — remains open: `convEnvCPLQ`'s own
  multi-face triangulation can produce a triangle piece with exactly ONE convex edge (a genuinely
  rational envelope), which `conjPieceCPLQ` cannot conjugate yet (its own header TODO). **Revised
  (Part 2c, below): the 3-convex-edge single-triangle case above is NO LONGER exempt from this
  gap either** — it was solved on its own only as long as it never produced a rational piece; Part
  2c's fix means it now correctly does (see Part 2c's own entry, and `conjCPLQ`'s "NOT implemented"
  bullet, for the reasoning and current status).
  **Diagnosis (2026-07-17, later session): the 2026-07-15 "fix" above is itself NOT tight in
  general -- confirmed via ground truth independent of the whole conjugate/duality pipeline.**
  Found via a concrete counterexample, `T=(-9.95506,3.70366),(-9.345,-5.34552),(1.29049,5.31738)`
  (from a stress test's silent-wrong-answer case; the SAME triangle as `doc/bug.svg`/
  `doc/bug_pipeline.svg`, which document the ORIGINAL diagnosis this continues), cross-checked
  against TWO independent ground-truth methods: (1) `sup_s[s.x0-supBilinearOverPoly(s,V)]`
  (numerically maximized, matching the session handoff's earlier machine-precision value at its
  own reference point to 5+ decimals); (2) the 3D convex hull of a dense (300-refinement
  barycentric grid, ~45k points) sample of the graph `(x,y,xy)` over `T` (the lower-hull facets
  are exactly the graph of `conv(f)`, a standard fact -- same method as the interactive artifact
  from the prior session, but now driving a systematic sweep instead of one point).
  - **Part 1, FIXED**: the "other" sub-triangle's own formula (`buildEdgeAffinePiece`, touching
    `u1*u2` along its remaining convex-edge segment and the AFFINE CHORD along the weak edge) is
    a valid boundary condition on that sub-triangle's OWN BOUNDARY but not sufficient to pin down
    the tight envelope in its INTERIOR: fitting a general quadratic to ground truth over that
    sub-triangle left a real (~0.36 max, ~0.11 rms) residual, proving the true envelope there is
    not even a single quadratic, let alone `buildEdgeAffinePiece`'s specific one. Diagnosed the
    correct replacement by RECLASSIFYING the sub-triangle's own 3 edges from scratch: the
    "remaining convex-edge segment" is still convex (same line, hence same slope, as the
    original triangle's edge it's part of) and the weak edge is unchanged, but the NEW internal
    seam (the cevian itself) is ALSO consistently classified non-convex across ~17 random
    2-convex-edge triangles checked -- i.e. the sub-triangle is a genuine, ordinary
    ONE-convex-edge triangle, not a special case needing a bespoke construction at all. Replaced
    `buildEdgeAffinePiece`'s result with a direct call to the EXISTING, already-validated Appendix
    A.3 rational envelope (`envelopeFromClassified` case 1) on that sub-triangle's own
    reclassified edges -- exactly the same "reclassify each sub-triangle and call
    `envelopeFromClassified` again" pattern the 3-convex-edge case (`splitThreeConvex`) already
    uses, now unified across both split cases (see `convEnvCPLQ.m`'s `nCE==2` branch and
    `splitTwoConvexEdges`'s updated header). Verified against ground truth to grid resolution
    (~1e-4, i.e. the sampling/hull discretization noise floor, not a real residual) on the
    counterexample AND on 16/17 fresh random 2-convex-edge triangles that needed a genuine split
    (the 17th was within 0.15% relative error, plausibly just coarser-grid noise for that
    particular triangle's scale -- not independently re-verified with a finer grid this session).
    `buildEdgeAffinePiece` itself is KEPT (renamed in purpose, not deleted) as a placeholder used
    only to locate the cevian via `seamPoint`'s existing line-factoring argument, which still
    applies unchanged since the placeholder and the real replacement formula both touch `u1*u2`
    identically along the same shared full line and therefore agree at the shared opposite
    vertex, the only two facts `seamPoint`'s derivation actually needs. **Consequence for Step 2**:
    this makes the "other" sub-triangle's own piece genuinely rational, something
    `conjPieceCPLQ` cannot conjugate yet (the SAME pre-existing gap already flagged above for the
    multi-face/non-triangular case) -- `maxQuaParTest.buildG1G2ForTriangle`/`extractTriFace` now
    detect this and error clearly (`extractTriFace:rationalFaceNotSupported`) instead of silently
    dropping the denominator and building the wrong `QuaPoly`, which is what `
    matchHalfEdgesRejectsSameSideRayPairingAndDropsSubsumedPieces` (the one test that exercises
    this exact triangle) was unwittingly doing before; that test now pins the new loud failure
    instead of the old silent one (its own three `maxQuaPar` assembly-bug fixes remain covered by
    this file's OTHER repro triangles, all still passing). Full suite: 79/79 (unchanged count;
    one test's assertions replaced, none newly broken).
  - **Part 2, STILL OPEN, more fundamental than Part 1**: fixing the "other" sub-triangle alone
    does NOT make the SAME-triangle counterexample's overall wrong-answer disappear, because
    `q1` -- the UNCHANGED formula kept for the sub-triangle containing `P` -- turns out to be
    **not tight throughout that sub-triangle either**, a gap that predates this session's fix and
    was previously masked by the (much larger) Part 1 gap dominating whichever point a spot-check
    happened to land on. Confirmed directly: sampling points along the shared internal seam
    between the two sub-triangles shows `q1` matching ground truth only AT the seam's two
    endpoints (as expected -- both are original triangle vertices, or lie on an original convex
    edge) but undershooting truth by up to ~2.55 at interior points of the seam, and by up to
    ~10.2 at other interior points of its own sub-triangle -- i.e. `q1`'s zone of exact tightness
    is a genuine, smaller-than-the-sub-triangle 2D neighbourhood of `P`, not the whole sub-triangle
    the current cevian split hands it.
    **Established fact (per Locatelli): the convex envelope of a quadratic over a polyhedron
    always admits a polyhedral subdivision.** This session's numerical exploration of where `q1`
    stops matching ground truth (rays from `P`, bisection; a general conic fit to the empirically-
    traced transition points) fit a CURVE far better than any single straight line -- but per
    Locatelli's theorem, a polyhedral (straight-edge) subdivision is GUARANTEED to exist; that
    numerical finding therefore reflects an INCOMPLETE SEARCH (the specific combinatorial splits
    and formula pairings tried this session -- `q1` vs. a single Appendix-A.3 rational piece,
    variously placed -- were simply the wrong pairing/split, not evidence of genuine
    non-polyhedrality), not a counterexample to the theorem. Recorded here so the next session
    does not re-derive this and does not pursue a curved-`RatPol`-domain extension, which is now
    known to be unnecessary. Ruled out this session (all checked computationally against ground
    truth, all fail -- kept as negative results to avoid retrying): (a) `max(q1, caseA1)` or
    `max(q1,caseA1,caseB1)` over the WHOLE original triangle (the rational candidates are only
    valid minorants -- i.e. `<=` truth -- within a bounded sub-region; extended past that they
    exceed truth, so a naive pointwise max wrongly picks the now-invalid candidate); (b) splitting
    via a cevian from `P` into the weak edge at a THIRD point `S` and using two case-1 formulas
    with opposite vertex `S` (the two formulas agree exactly along the WHOLE cevian `P-S` for
    ANY `S`, an exact identity, but neither is tight throughout its own region for essentially any
    choice of `S` tried); (c) fitting an unconstrained general rational family (matching ground
    truth exactly at 5 points while still touching `u1*u2` along both original convex edges) --
    solvable but the fit does not generalize away from the 5 points (blows up near an
    uncontrolled pole). The correct polyhedral split (a DIFFERENT combinatorial structure than any
    tried this session -- e.g. more than 2 pieces, a different choice of which vertex/edge anchors
    each piece, or a genuinely different split direction) has NOT yet been found.
    **Practical impact, quantified**: a 1500-triangle `rng(2026)` resample (matching the earlier
    session's seed) found 508 genuine-2-convex-edge-split triangles, 277 with no detectable gap
    at 5 random interior points each (`AbsTol` 1e-3 against a numerically-maximized biconjugate)
    and 231 (~45%) with a real, positive gap (`r` too loose, `truth-r` up to ~4.2) at at least one
    of those points -- both BEFORE and AFTER this session's Part 1 fix, identically (same 231
    triangles), confirming Part 1 is a genuine, isolated improvement to the "other" sub-triangle
    specifically, not a fix to the aggregate wrong-answer rate, which is still dominated by Part
    2. Separately confirmed `r` remains a VALID minorant throughout (no case of `r` exceeding
    truth survived a more carefully multi-started numerical check -- an earlier quick pass using a
    single coarse-grid-seeded `fminsearch` per point produced ~90 apparent "invalid" points that
    were a ground-truth-oracle artifact, not a real defect in `r`, confirmed by re-solving those
    same points with a much finer multi-start search). **This should be the top priority for the
    next session**: Part 2 is the actual remaining source of wrong answers; Part 1's fix, while
    correct and worth keeping, does not move that needle. Given Locatelli's theorem, the next
    session's search should stay within straight-edge (polyhedral) splits -- likely more than 2
    pieces, or a different combinatorial pairing than `q1` + a single Appendix A.3 piece -- not a
    curved-boundary extension. See the session handoff for the concrete counterexample, the
    ground-truth methodology (reusable), and the ruled-out approaches above.
  - **Part 2a, FIXED this session (2026-07-17, later session): the cevian LOCATION itself was
    wrong, independent of Part 1.** Root cause found by deriving, from [COAP] Section A.1's own
    dual formalism (eq. 11-14) applied directly to `q1` (rather than assumed away, as the paper's
    own Appendix A.4 does when it asserts "the domain is the entire triangle" -- confirmed FALSE in
    general, even on the paper's own worked example, `V=(2,1),(0,0),(1,0)`: q1(0.474343,0) =
    -0.042780 there instead of the true value 0, since f=xy is identically 0 -- hence `>=` any
    valid minorant, in particular `>=0` -- along the WHOLE weak edge there): the constructive
    envelope theorem says the implied dual point at `(x,y)` is `(a,b) = grad q1(x,y)`, and `q1` is
    valid exactly while each classified edge's own tangency point `s_j(a,b) = (a+b*m_j-q_j)/(2 m_j)`
    stays inside that edge's OWN segment. Since `q1` is quadratic, its gradient is affine, so
    `s_j(x,y)` is affine in `(x,y)` -- meaning "`s_j(x,y)` = far endpoint's x-coordinate" is a
    genuine STRAIGHT LINE (consistent with Locatelli), always passing exactly through that far
    endpoint by construction. The OLD `seamPoint`/`buildEdgeAffinePiece` pair (a "factor `q1` minus
    a placeholder that merely matches the weak edge's AFFINE CHORD" argument) computed a
    DIFFERENT, wrong line: on the paper's own example it placed the cevian at x=0.2929 along the
    target edge instead of the correct x=2-sqrt(2)=0.5858 (confirmed via a closed form,
    `xR = xP + sqrt(mh/mw)*(xA-xP)` for the beta=1,lin=0 case, verified symbolically) -- almost
    exactly HALF the correct distance from `P`. Concretely, the OLD seam made `T1` (the `q1` piece)
    too BIG, wrongly annexing a wedge near the weak edge where `q1` undershoots (confirmed: e.g. at
    the wedge centroid (0.626,0.146) in the paper's own example, `q1`=0.0616 vs true 0.0640).
    New `edgeClipCevian` (replacing `seamPoint`+`buildEdgeAffinePiece` in `convEnvCPLQ.m`) builds
    the line directly from `q1`'s own plain coefficients (`p=2A+m*B, q=B+2m*C,
    r=D+m*E-q_edge-2m*xFar`) and intersects it with the OTHER convex edge -- same final
    "intersect a line with `y=m*x+q`" step the old code used, just fed the correct line, so the
    surrounding candidate-selection structure (try cevian-from-A first, fall back to cevian-from-B
    if it does not land inside its target edge) is unchanged. Verified: matches ground truth
    exactly (to solver precision) on the paper's own example, including at the wedge point the old
    seam misclassified; full suite 79/79 (two pre-existing tests' hardcoded expectations were
    themselves artifacts of the two bugs above -- `convEnvCPLQTest.bilinearTwoConvexEdgesQuadratic`
    asserted the paper's own (disproven) "whole triangle" claim at a point now correctly assigned to
    the other piece, and `conjPieceCPLQTest.psdRank1QuadraticEndToEnd` hardcoded a `nf` count for
    `T1`'s shape that changed now that `T1` is correctly smaller; both updated with comments
    explaining why, not just silently changed).
    **Practical impact, quantified**: the SAME 1500-triangle `rng(2026)`-style resample (468
    genuine 2-convex-edge splits found this run -- close to, not identical to, the prior session's
    508, since `rng` state is consumed differently by this session's own regenerated stress script,
    not the original uncommitted one) drops from ~45% of triangles showing a detectable gap to
    **23.1%** (108/468), and the single worst gap seen dropped somewhat (2.07 vs ~4.2) but did NOT
    disappear -- so this is a real, substantial, independently-verified improvement, not a full fix.
  - **Part 2b, STILL OPEN: `T1` itself needs a further split for sufficiently asymmetric edge
    slopes, found while investigating why Part 2a's fix does not close the gap to 0%; narrowed
    considerably this session but not yet resolved.** Concrete repro (from a fresh `rng(2026)`
    stress script, its first failing triangle): `V=(-5.2645,3.4904),(3.1062,0.5450),(5.0430,7.7766)`,
    `f6=[0.3963 0.8289 0.0284 -2.4834 2.1149 -5.5999]` (i.e. `q=0.5*x'Hx+L'x+kappa` with
    `H=[[0.3963,0.8289],[0.8289,0.0284]]`, `L=[-2.4834,2.1149]`, `kappa=-5.5999` -- QuaPoly's stored
    convention, NOT plain `[A B C D E F]`; `mh=0.275376,qh=-2.551975,mw=1.578980,qw=7.945196` in the
    bilinear frame, a much bigger `mh`/`mw` ratio than the paper's own `1`/`0.5` example, where no
    such further split is needed at all -- see below). With Part 2a's corrected cevian, `T1` (the
    `q1` sub-triangle, vertices `P,A,Ra` in the bilinear frame) is confirmed exact along ALL THREE of
    its own edges in their ENTIRETY (the two original full rays from `P`, AND the new `A-Ra` seam,
    checked pointwise along each, not just at endpoints -- the seam agreement is additionally a
    proven algebraic IDENTITY, `q1 == ` the corresponding Appendix-A.3 piece anchored at `A` (called
    `R_Aw` below -- exactly `T2`'s own formula from Part 1), along the whole line `s_h(x,y)=x_A`,
    verified symbolically for general `mh,qh,mw,qw`). Yet `q1` still UNDERSHOOTS truth by up to
    ~0.38 at points in `T1`'s STRICT INTERIOR (tight on the entire boundary of `T1`, not tight
    inside it).
    **Key finding this session (CORRECTED -- an earlier draft of this note overclaimed "matches
    everywhere"; see the actual data below)**: at the SPECIFIC points identified by the `q1`-gap
    scan (i.e. points where `q1` itself undershoots truth), `R_Aw` matches ground truth EXACTLY (6/6
    checked, to solver precision) -- i.e. **`T2`'s own formula, `R_Aw`, extended past its home
    triangle `T2` into `T1`, is the right formula for those specific failing points**, not a new,
    undiscovered third formula. HOWEVER, a separate, broader check (~60 UNIFORM random points across
    ALL of `T1`, not just the known-failing ones) found `R_Aw` OVERSHOOTS truth in the clear MAJORITY
    of them (worst case ~2.8, many in the 0.1-1.8 range), essentially everywhere except a small
    region -- so `R_Aw` is NOT simply valid over some large sub-region of `T1` glued to `q1` along
    one boundary; it is only correct in the SAME small notch the barycentric grid scan (below)
    independently identifies, and invalid (an overshooting, unusable candidate) essentially
    everywhere else in `T1`. `R_Aw` is `>= q1` EVERYWHERE checked in `T1` (including deep in `q1`'s
    correct territory near `P`), so the correct split is a genuine geometric boundary, not a value
    comparison. This resolves the earlier confusion about what the "bubble" even
    is: it isn't that `q1` is wrong in some mysterious interior region for no visible reason -- it's
    that `T2`'s own rightful territory (where `R_Aw` is the tight envelope) extends further into
    what Part 2a's single cevian currently calls `T1` than that cevian accounts for. The real
    question still open is where the true `R_Aw`-vs-`q1` boundary is and why it isn't simply the
    `A-Ra` cevian already found.
    Ruled out this session: (a) `twoEdgeQuadPlain`'s +/- branch choice -- both branches are
    legitimate rank-1-PSD, edge-touching solutions, but the other branch (`s=-1`) is a much LOOSER
    minorant everywhere checked, not relevant here; (b) a naive RECURSIVE re-application of
    `edgeClipCevian` to `T1` itself (treating it as its own fresh 2-convex-edge triangle with weak
    edge `A-Ra`): the "cevian from `A`" candidate is degenerate (reproduces the SAME `Ra`, since
    edge `h`'s far endpoint is still `A`, unchanged -- a general fact, not specific to this repro),
    and the "cevian from `Ra`" candidate (using edge `w`'s far endpoint now `Ra` instead of `B`)
    lands just OUTSIDE its target segment (`t~1.10`, not `(0,1)`) on this repro -- confirms plain
    recursion of the SAME two-candidate logic is not the mechanism; (c) the two
    vertex-anchored-at-`P` rational candidates (`R_Ph`, `R_Pw`, i.e. Appendix A.3 with anchor `P`
    instead of the far vertex) are wildly invalid in the bubble (~-16 vs a true value in the single
    digits), not relevant; (d) comparing `q1` directly against `R_Aw` as an equation (`q1=R_Aw`) to
    find a crossing curve was tried in the ORIGINAL (non-recursive) frame earlier this session and
    is a genuine CONIC, not a line (symbolically confirmed) -- consistent with `R_Aw` never actually
    crossing back below `q1` inside `T1` (it stays `>= q1` throughout), so that equation is not the
    right boundary condition to chase either.
    **Leading hypothesis for next session, not yet derived**: `R_Aw` (built for Appendix A.3's
    "V-E" case with `V=A`, edge `w`) has its OWN Section A.1 dual-tangency validity limit, exactly
    analogous to `q1`'s (`edgeClipCevian`'s derivation) -- i.e. `R_Aw` is valid only while its own
    implied dual point `(a,b) = grad R_Aw(x,y)` keeps edge `w`'s (or possibly edge `h`'s, since `A`
    is adjacent to edge `h` too) own tangency point inside a segment appropriate to `R_Aw`'s actual
    home triangle `T2 = {A,Ra,B}`. That home triangle's OWN third vertex is `B` (not `Ra` or `P`),
    so the relevant boundary is likely tied to `B` specifically, or to edge `h` becoming "visible"
    again once we cross into `T1` (where, unlike in `T2`, edge `P-A` is genuinely convex, a fact
    `R_Aw`'s derivation implicitly assumed false). This has NOT been derived or checked this session
    -- it is a concrete, well-scoped next step: redo `edgeClipCevian`'s exact derivation but for the
    Appendix A.3 rational formula's own gradient instead of `q1`'s, most likely against edge `h`
    (using `mh,qh`) since that is the edge whose "convexity" `R_Aw` implicitly ignores once extended
    past `A`. This was NOT visible in the paper's own (near-symmetric, `mh=1,mw=0.5`) worked example
    -- a fresh check there (25-40 random interior samples of the corrected `T1`, all matching `q1`
    exactly, no bubble) -- so whatever this boundary is, it must degenerate to "outside the
    triangle" (irrelevant) for `mh`/`mw` close to `1`, and move inside `T1` as the ratio grows -- a
    useful property to check any proposed closed form against.
    **Precise shape of the gap region, mapped this session via a full barycentric grid scan of
    `T1`** (66 points, `n=10` subdivisions, `(i,j,k)` = barycentric weights toward `(P,A,Ra)`): gap
    is EXACTLY zero for `i>=3` (i.e. everywhere more than ~30% of the way from the `A-Ra` edge
    toward `P`) and also exactly zero at `j=0` (the whole `P-Ra` edge) and `k=0` (the whole `P-A`
    edge), confirming those are genuinely fully tight, as found earlier. The gap is strictly
    positive only for `i=1` (`j=1..6` out of `k=9-j`, peaking at `j=3`, gap~0.11) and `i=2`
    (`j=1,2` only, gap~0.005-0.01) -- i.e. the failing region is a SMALL, THIN notch hugging the
    MIDDLE of the `A-Ra` edge, vanishing at both the `A` and `Ra` ends (consistent with the proven
    identity there) and fading out within about 20-30% of the distance back to `P`, NOT a broad
    bubble filling much of `T1`'s interior as first thought from the coarser earlier sampling. Each
    row (fixed `i`) has TWO transitions (gap turns on somewhere after `j=0`, then off again before
    reaching `k=0`), not one -- meaning whatever the true correction is, it is not a single
    half-plane cutting `T1` from one side (that would give only one transition per row); it is
    something closer to a small triangular "notch" carved out of `T1` near the middle of `A-Ra`,
    which by Locatelli must itself have straight edges, likely meeting at (or near) an interior
    apex point closer to `A-Ra` than to `P`. This is a materially sharper target for the next
    session than "an interior bubble" -- concretely, look for a THIRD formula (probably `R_Aw`
    again, given it already provably matches truth throughout the region checked) valid only in
    that small notch, bounded by two new short edges rather than by extending `T2`'s existing
    single cevian.
    **Precise transition points located this session (bisection, `i`=`P`-weight fixed, `j`=`A`-weight
    varied)**: at `i=0.10,0.15,0.20` the LOWER transition (near the `P-Ra` edge side, `j` near 0)
    sits at `j=0.0011,0.0021,0.0054` respectively (essentially hugging the `P-Ra` edge, but a real,
    growing offset from it, not exactly on it), while the UPPER transition (`j` further from 0)
    sits at `j=0.634,0.451,0.264` -- shrinking fast, extrapolating to 0 (i.e. the notch closing up)
    around `i~0.3`. In `(x,y)`: lower transition points `(-1.1275,5.2034), (-0.7795,5.3418),
    (-0.4201,5.4699)`; upper transition points `(1.9876,2.4345), (1.4325,3.3757), (0.8532,4.3382)`.
    Checked whether the lower set simply lies on the extended `P-Ra` LINE itself: close but not
    exact, with a growing residual (`~0.008,0.014,0.035`) inconsistent with solver noise alone --
    so the lower boundary is a DIFFERENT, nearby line (or a genuine curve, in which case the
    combinatorial split found so far is still not the Locatelli-guaranteed correct one). Neither
    the lower nor upper transition set was confirmed exactly colinear with the numerical precision
    available this session -- worth re-deriving symbolically (analogous to `q1`'s own
    `edgeClipCevian` derivation, i.e. from a genuine closed-form gradient/dual-tangency condition)
    rather than continuing to fit noisy numerical transition points, which is the likely next step.
  - **Part 2c, found later the SAME session -- Part 2b's notch closed, via a cleaner criterion,
    which also turned out to apply, forced and unconditionally, to the pre-existing `nCE==3` code
    path -- see "Status" below.** Re-derived the `T1`/`T2` boundary from scratch as a genuine
    C1 tangency condition rather than the earlier value-only match: writing `q1 - R_Aw` as one
    rational expression and clearing denominators, the numerator factors (using that the anchor
    vertex lies on its own edge's line) as `(that edge's line) * (a PERFECT SQUARE)`, i.e. a double
    root -- `q1` and `R_Aw` agree in both VALUE and GRADIENT along the entire line through the
    anchor vertex with slope `-sqrt(mh*mw)` (verified to ~7 significant figures numerically, and
    exactly via the symbolic factorization). This slope depends ONLY on `mh,mw` -- `beta` and the
    affine shift `lin` cancel out of the derivation completely (confirmed symbolically: the
    `linDen` contribution to the numerator is exactly `(lin1 x+lin2 y+lin3)*D`, canceling the same
    term from `qfull*D`). Replaced `edgeClipCevian` with `tangentCevian` (line through the anchor
    vertex with this slope, intersected with the other convex edge) in `convEnvCPLQ.m`. Verified
    exactly (ground truth, both resulting sub-regions) on the Part 2b repro triangle, and the
    aggregate stress-test gap rate (same rng(2026) methodology) dropped further, from 23.1% to
    9.3% (108/468 -> 68/731; the genuine-2CE-split COUNT itself changed between runs, 468 vs 731 --
    not yet reconciled, possibly a difference in how many triangles error out vs split cleanly
    between the two code versions, not independently confirmed as harmless).
    **Related discovery, NOT part of the original Part 2 scope**: the residual ~9% turned out to
    be entirely (at least in the one repro checked) attributable to a DIFFERENT, pre-existing code
    path -- `nCE==3` (three convex edges, `splitThreeConvex`) -- which calls
    `envelopeFromClassified`'s plain-quadratic case DIRECTLY on each of its 2 sub-triangles,
    without ever checking whether THAT sub-triangle (itself exactly 2-convex-edge) needs the same
    kind of further split `splitTwoConvexEdges` now checks for. Generalized `convEnvCPLQ.m` with a
    new recursive `solveTriangleBF`, used uniformly by both the `nCE==2` and `nCE==3` branches (via
    a new `assemblePiecesBF`, built on the already-existing, arbitrary-piece-count
    `assembleTriangles`), so either branch's sub-triangles can themselves recurse into further
    pieces. Checked against 9 sample "3 convex edge" triangles (the 8 used by currently-failing
    tests, described below, plus the paper's own famous `T=(0,0),(3,3),(1,2))`: **all 9 of 9** now
    produce 4 pieces (2 plain quadratic + 2 genuinely rational), not the previously-assumed-always
    2. Spot-checked `T=(0,0),(3,3),(1,2)` specifically against the PRE-session code at 6 interior
    points: 5 of 6 values are numerically IDENTICAL to machine precision, one differs by ~0.0018 --
    i.e. the correction is real (confirmed by the perfect-square/tangency derivation, not a
    hand-wave) but can be extremely small for a near-symmetric triangle, and apparently large
    enough in `curv`'s existing tolerance check (`1e-9*scale`) to trigger a further split anyway.
    **The "how often does this trigger" QUESTION is presented as a strong pattern from a SMALL,
    non-random sample (9 triangles, all originally hand-picked by past sessions for unrelated
    reasons -- assembly-bug reproduction, not randomly sampled) -- treat "essentially all
    3-convex-edge triangles need 4 pieces" as a plausible working hypothesis, NOT an established
    fact about FREQUENCY.** It has not been checked against a proper random sample the way the
    `nCE==2` stress test was. This is a genuinely open, but low-stakes, item for a future session
    (a random-sample stress test, analogous to the `nCE==2` one).
    **Whether to APPLY the fix at all is, by contrast, not open.** The tangency criterion itself
    (`tangentCevian`'s slope `-sqrt(mh*mw)`) is proven to depend only on a 2-convex-edge triangle's
    own `mh,mw` -- nothing about `beta`, the affine shift, or how that triangle arose (standalone
    `nCE==2`, or a `splitThreeConvex` sub-triangle) enters the derivation. A 3CE sub-triangle is an
    ordinary 2CE triangle; there is no basis for exempting it from a criterion already proven to
    hold for every 2CE triangle. So `solveTriangleBF`'s recursion is a forced consequence of Part
    2c, not a scope expansion up for debate -- reverting it would mean knowingly reinstating a
    formula already proven untight in some `nCE==3` cases. **Status: KEPT, decision resolved.**
    This still breaks 6 existing tests (`conjCPLQTest.indefiniteTriangleThreeConvexEdgesUsesStep3`
    plus 5 `maxQuaParTest` cases, all going through the shared `maxQuaParTest.buildG1G2ForTriangle`
    helper) because `conjPieceCPLQ` cannot yet conjugate a rational piece (the same pre-existing,
    separately-tracked gap Part 1 already hit -- see `conjPieceCPLQ.m`'s own entry below: the
    closed-form recipe for this already exists in the reference `cPLQ/plq_1piece.m`, it just hasn't
    been derived into CCA2's closed-form-coefficient style yet), and EVERY one of the 5 distinct
    hand-picked "3 convex edge" triangles those tests use now produces 2 rational sub-pieces --
    meaning this is not 6 independent assertion tweaks but the entire test-fixture STRATEGY those
    tests share (pick a 3-convex-edge triangle, derive `g1,g2` via the live Step1+Step2 pipeline)
    becoming structurally incompatible with Step 1's own (now more) correct behavior. Those 5 tests
    are regression tests for `maxQuaPar`'s OWN assembly bugs (dedup, orphan half-edges, near-
    duplicate vertices, etc.), not for Step 1's tightness, so the fix is to freeze each test's
    existing `g1,g2` (computed once from the pre-session code) as hardcoded `QuaPar` fixtures,
    decoupling them from the live pipeline. **DONE** (same session): all 6 fixed this way; full
    suite verified back to passing.
    **Two FURTHER tests broke for the identical reason** (not part of the original "6", found only
    once the full suite was rerun after keeping the fix): `conjPieceCPLQTest.
    pickRepFindsThinEdgeStripFace` and `convEnvCPLQTest.threeConvexEdgesSplit`, both hardcoding
    `nf==2` for a 3-convex-edge triangle. Fixed the same way as `conjCPLQTest.
    indefiniteTriangleThreeConvexEdgesUsesStep3` above (update the split-count assertion to 4;
    `pickRepFindsThinEdgeStripFace`'s extracted face-1 sub-triangle is now smaller — a further-split
    piece of the original first sub-triangle — landing in `conjPieceCPLQ`'s "tied vertex" 5-face
    case (`conjPSDRank1QuadTriangleTie`) instead of the generic 6-face one; both reverified against
    ground truth after updating).
    **Also found and fixed, while diagnosing `threeConvexEdgesSplit`'s remaining failure after the
    `nf` fix**: `RatPol.eval` (the main polytope-membership loop, not the edge-chain branch) blindly
    overwrote `fVal` with whichever matching face was processed last, with no check at all —
    unlike the documented policy ("`fval(i)=nan` if PC is discontinuous at `x(i)`"), which the
    edge-chain branch actually implements but this one never did. This was silently harmless as
    long as every piece was a plain quadratic (no matching face could ever evaluate to NaN on its
    own), but a genuinely rational (Appendix A.3) piece's denominator is *exactly* zero at its own
    anchor vertex by construction — a REMOVABLE singularity (the true limit exists and is finite),
    not a real discontinuity — so it can report NaN at a vertex where OTHER matching faces
    correctly report the true finite value. Confirmed exactly this on `threeConvexEdgesSplit`'s own
    apex vertex `(1,1)`: 2 of the 4 matching faces (the plain-quadratic ones) correctly evaluate to
    `1`, the other 2 (the rational ones, each anchored at that very vertex) evaluate to `NaN`, and
    since faces are visited in order, the LAST one processed (a NaN) clobbered the correct answer.
    Fixed by tracking, per point, whether a finite value has already been found: a later face's NaN
    no longer overwrites an existing finite value, while a later face's genuinely DIFFERENT finite
    value still correctly triggers NaN (preserving the documented real-discontinuity behavior).
    Verified: `(1,1)` now evaluates to exactly `1`; full suite re-run clean (147/147) after this fix
    plus the two `nf` fixes above.
- `scalarMul`/`negate` — an instance method on each of `QuaPoly`/`QuaPar`/`RatPol` (trivial:
  scales `f`, the numerator for `RatPol`; domain/mesh untouched). No `RatPar`, so no single
  shared implementation; each class has its own copy.
- `add` — `QuaPoly` and `QuaPar` (`QuaPoly.add`/`QuaPar.add`, in `addQuaPoly.m`/`addQuaPar.m`):
  overlays the two domain subdivisions by pairwise convex-polygon clipping and sums the two
  quadratics on each overlap cell (always exactly one quadratic per cell — no case-split needed,
  unlike `maxQuaPar`'s pointwise max). `QuaPoly.add` (straight edges only) is the base case,
  adapted from `maxQuaPar.m`'s facePoly/clipByFace machinery, generalized to allow an unmatched
  edge to become a genuine domain-boundary edge rather than an error, since unlike `maxQuaPar`'s
  inputs, `add`'s inputs need not be full-domain. `QuaPar.add` generalizes this further to allow
  one curved (parabolic-arc) edge per face: every boundary edge and every clipping constraint —
  straight or curved — is represented uniformly as a 6-coefficient conic row (`a=b=c=0` for the
  linear case); a straight edge is parametrized affinely in a scalar `t`, the one curved edge via
  an axis-rotation of its parabola (`x(u),y(u)`, each quadratic in `u`); either substituted into
  any constraint gives a plain polynomial (degree ≤2 for a straight edge, degree ≤4 for the curved
  edge) whose real roots are exactly the crossing points — found directly via full root-solving,
  not an endpoint-sign-disagreement shortcut, since a conic constraint (unlike a linear one) can
  cross even a straight edge twice while agreeing in sign at both endpoints. Handles bounded and
  unbounded domains (including replacing an entirely-discarded ray with a new one along a linear
  cutting constraint, mirroring `clipPolyHalfPlane`); a result needing two curved edges, a curved
  ray, a non-single-branch (double-line-degenerate) parabola, a new *unbounded curved* ray, or
  more than one disjoint surviving piece (other than the two-rays-survive-with-a-middle-gap
  shape, also not implemented) all error clearly rather than being silently mishandled — see
  `addQuaPar.m`'s header STATUS block. See `addQuaPolyTest.m`/`addQuaParTest.m`. **Not yet
  extended to `RatPol`** (would need a common-denominator sum).
  - Found and fixed while building this: `polyConstraints` (shared by `addQuaPoly.m` and
    `maxQuaPar.m`) looped `1:nv-1` for a BOUNDED poly, silently dropping its closing edge
    `(nv,1)` — under-constraining `clipByFace` by exactly one of `polyL`'s half-planes. Fixed in
    both files (loop now wraps for bounded polys, unaffected for unbounded ones); regression test
    added to `addQuaPolyTest.m`.
- `infConv.m` — `conj(add(conj(f,engine),conj(g,engine)),engine)` (II.6), valid for `f,g` both
  convex. Thin composition, no new geometry: since `conj(f,engine)` may come back as either
  `QuaPoly` (conjCPLQ's full-domain-quadratic shortcut) or `QuaPar` (its general single-piece
  case), and `QuaPoly.add`/`QuaPar.add` only accept same-class operands, both conjugates are
  first promoted to `QuaPar` (a lossless relabeling — `QuaPoly` and `QuaPar` share the same
  `V/E/f/F` layout, `QuaPar` just adds an all-zero `Ec`) before calling `add`. Exercised
  end-to-end (both `conj` calls and the final `conj` all round-trip) only on full-domain
  quadratics so far — a bounded-triangle `f,g` pair currently can't complete the final `conj`
  step, since that would land the summed `QuaPar` in `conjCPLQ`'s still-unimplemented Step 3
  (max of conjugates over a multi-face domain). See `infConvTest.m`.
- `addQuadratic`/`addScaledEnergy` — instance methods on both `QuaPoly` and `QuaPar`. Simpler
  than `add`: the added term is a FULL-DOMAIN quadratic, so it never restricts the domain and
  needs no overlay/clipping — just bump columns 5:10 (`[x^2 xy y^2 x y const]`) of `f` by
  `[A(1,1) A(1,2) A(2,2) b(1) b(2) c]` identically on every face/edge row. `addScaledEnergy(f,
  alpha)` adds `alpha*(x1^2+x2^2)` (note: no leading 1/2 — it's `addQuadratic(f, 2*alpha*eye(2),
  0, 0)`), matching exactly how `moreau`'s own formula uses it (see below).
- `moreau.m` — `e_mu*f` via the "expand the square" single-conjugate identity
  [HIRIART-URRUTY-07] (II.6): `g = addQuadratic(scalarMul(f,mu), eye(2), 0, 0)` (= `mu*f +
  1/2||.||^2`), `h = addScaledEnergy(scalarMul(conj(g,engine), -1/mu), 1/(2*mu))`. Pure algebraic
  regrouping, no biconjugation, so — unlike `infConv` — **no convexity assumption on `f` at
  all**; tested on a genuinely nonconvex full-domain quadratic (indefinite Hessian) as well as
  convex ones, cross-checked against an independently-solved stationarity condition for a
  shifted (linear+constant terms) case. See `moreauTest.m`.
- `lasryLions.m` — `h = negate(moreau(negate(moreau(f,lambda,engine)),mu,engine))`, a pure
  composition of `moreau`/`negate` (II.6): no new geometry, and — inherited from `moreau` — no
  convexity assumption on `f` (tested on a mixed-convexity diagonal quadratic, one convex
  component and one concave). See `lasryLionsTest.m`.
- `proxAverage.m` — the "sandwich of two conjugations around a weighted `add`" derivation (II.6):
  `gf=addQuadratic(scalarMul(f,mu),I,0,0)`, `gg=` likewise for `g`; `s = lambda*conj(gf,engine) +
  (1-lambda)*conj(gg,engine)` (both conjugates promoted to `QuaPar` via `toQuaPar.m` before the
  weighted `add`, same plumbing `infConv` needs and for the same reason); `h =
  addScaledEnergy(scalarMul(conj(s,engine),1/mu), -1/(2*mu))`. Valid (returns the actual
  proximal average, not just its convex envelope) only for `f,g` both convex, like `infConv` —
  the final step is a genuine biconjugation, unlike `moreau`'s pure algebra. Tested against `P`'s
  own defining characterization `e_mu*P = lambda*e_mu*f+(1-lambda)*e_mu*g`, cross-checked via the
  already-tested `moreau.m` rather than a hand-derived closed form for `P` itself. See
  `proxAverageTest.m`.
- `toQuaPar.m` — extracted shared helper (originally local to `infConv.m`, now also used by
  `proxAverage.m`): promotes a `QuaPoly` conjugate to the equivalent `QuaPar` (lossless
  relabeling, all-zero `Ec`) so two conjugates of possibly-different returned type can be
  `add`-ed together.

This completes the **nonconvex-PLQ operator pipeline** that has been this project's stated focus
(`conj`→`infConv`/`moreau`→`lasryLions`/`proxAverage`) — every operator in the II.4 method table
now has *some* working implementation, modulo the scope limits noted throughout (in particular
`conjCPLQ`'s own still-open Step 3, which is what currently confines `infConv`/`proxAverage` to
full-domain-quadratic `f,g` end to end — see their own bullets above).

**NOT implemented** (II.5.2/II.5.3 and parts of II.4/II.6 describe the intended design; no code
exists yet for):
- Conjugate engine **`'pqp'`** (parametric-QP, Jakee Khan [JAKEE-13]) — `conj(f,'pqp')` already
  errors explicitly in code (`QuaPoly.conj`, `PLQ:conj:engine`, "not implemented yet; use
  'cplq'"), it does not silently fall back or give a wrong answer.
- Conjugate engine **`'graph'`** (point-cloud + neighbour-graph, Tasnuva Haque
  [HAQUE-17]/[HAQUE-18]) — same: `conj(f,'graph')` errors explicitly, not ported.
- `RatPar` (the abstract storage-umbrella parent class, II.3).
- `add` for `RatPol` (`QuaPoly`/`QuaPar` are both done — see above).
- `partialConj` — not implemented for any engine.
- `convEnvDirect` (the `'direct'` envelope method built on Kumar's/Karmarkar's per-piece
  method, II.6) — not ported; `convEnv(f,'direct')` does not exist.
- `conjCPLQ`'s own Step 3 for a genuinely multi-face original domain (`nf>1`) or a single
  non-triangular face — needs `conjPieceCPLQ` to also handle a 1-convex-edge (genuinely rational)
  envelope piece, which `convEnvCPLQ`'s multi-face triangulation can produce and which
  `conjPieceCPLQ` cannot conjugate yet — this is what still confines `infConv`/`moreau`/
  `lasryLions`/`proxAverage` to full-domain-quadratic `f,g` for an *exact* end-to-end round trip on
  a genuinely multi-face `f,g`. **Revised (Part 2c, this session): the single-bounded-triangle/
  3-convex-edge case is NO LONGER exempt from this gap.** It was previously believed to work end to
  end because it "never needs a rational piece" (see `maxQuaPar.m`'s bullet above) — Part 2c's fix
  (recursively applying the proven `nCE==2` tightness split inside `splitThreeConvex`'s own
  sub-triangles, since that criterion depends only on a triangle's own `mh,mw`, not on how it
  arose) means a 3-convex-edge triangle now correctly produces 2 rational sub-pieces (alongside 2
  plain-quadratic ones), so it hits the SAME `conjPieceCPLQ` rational-piece gap as the general
  multi-face case. A bounded-triangle pair with 3 convex edges no longer works end to end for
  `conj`'s Step 3 until that gap is closed via the NUMERIC path; **it now DOES work via the Phase 1
  `cPLQ` integration below** (`quaPolyToPlq`/`evalFunctionNDomain`), same as the general multi-face
  case.
  **RESOLVED for the general multi-face case, via Phase 1's `cPLQ` integration (this session)**:
  `quaPolyToPlq.m` (CCA2 `QuaPoly` → cPLQ `plq`) + `evalFunctionNDomain.m` (numeric eval of a
  `functionNDomain` result) close this gap end to end for the symbolic engine — see
  `cplqAdapterTest.m`: a genuinely multi-triangle nonconvex PLQ (`f=xy` over a square split into 2
  independent triangles — exactly the case `maxQuaPar`'s own curved-edge restriction blocks
  numerically) now conjugates correctly through `quaPolyToPlq` → `.triangulate` → `.maximum` →
  `evalFunctionNDomain`, checked against numeric sup-sampling ground truth at 9 points. One
  narrow, documented exception: an exact symmetric TIE point (`s=(0.5,0.5)`, where both original
  triangles' own vertex cones meet) is not covered by the assembled region partition — a limitation
  of the vendored `functionNDomain.mergeL`/`region.removeTangent`, not of the adapter — flagged,
  not fixed, per `.claude/SESSION_HANDOFF.md`. This does not yet wire into `conjCPLQ.m`'s own
  `conj(f,'cplq')` call path (the adapter is currently invoked directly, not from `conjCPLQ`) —
  see "Next planned" below.

**Next planned — REVISED (this session), superseding the "closed-form derivation" plan below**:
prior sessions tried to close `conjCPLQ`'s Step 3 gap by deriving a closed-form numeric conjugate
for a single genuinely-rational `RatPol` piece (item 1, struck through below). Two findings changed
the plan:
- That single-triangle rational-piece conjugate turns out to be **dead code from the wired
  pipeline's point of view**: `conjSingleTriangle` always succeeds by conjugating the *original*
  quadratic piece directly, or an envelope that is affine/rank-1-PSD-quadratic — never a genuinely
  rational one. The *actual* remaining blocker for Step 3 on a genuinely multi-triangle input is
  `maxQuaPar`'s own documented restriction to purely polyhedral (non-curved) domains — confirmed by
  a direct test (`f=xy` on a 2-triangle square: each triangle's conjugate is a correct 6-face curved
  `QuaPar`, but `maxQuaPar(g1,g2)` errors immediately on the curved edges).
- Per **user decision (2026-07-18)**: rather than deriving new closed-form numeric machinery for
  this (curved-edge `maxQuaPar`, or the rational-piece conjugate) from scratch, **integrate the
  reference `cPLQ` package (almost as-is, symbolic) first** — this is actually what engine `'cplq'`
  was always specified to be (see §0 point 3: *"`'cplq'` — symbolic per-piece conjugate (`cPLQ`
  Lagrange-multiplier method)"*) — the numeric closed-form code built so far
  (`conjCPLQ.m`/`conjPieceCPLQ.m`/`convEnvCPLQ.m`/`maxQuaPar.m`) is a partial, faster reimplementation
  of a SUBSET of what `cPLQ` already does symbolically, not a replacement for it. **Two-phase plan:**

1. **Phase 1 — integrate `cPLQ` (almost as-is) and get it fully tested. IN PROGRESS, core adapter
   DONE (this session).** Copied the 9-file runtime dependency closure (`plq.m`, `plq_1p.m`,
   `plq_1piece.m`, `functionNDomain.m`, `region.m`, `domain.m`, `symbolicFunction.m`,
   `conjugateExpr.m`, `yIntercept.m` — NOT `quadQuad`/`qq_conj`, which turned out to be unused
   offline derivation scripts, correcting this section's own earlier claim) plus their 5 test
   suites into the CCA2 repo root (from the untracked `cPLQ/` reference clone). Fixed 2 toolbox-
   version-compat bugs found while establishing a passing baseline (`isequal(sym,double)` no
   longer works in the current Symbolic Math Toolbox; `plq.m`'s `pieces` property was mistyped,
   silently breaking all of `testMaxMultiRegion.m`'s 24 tests) — see `.claude/SESSION_HANDOFF.md`
   for full diagnosis. Built the thin adapter: `quaPolyToPlq.m` (CCA2 `QuaPoly` → cPLQ `plq`, per-
   face via `QuaPoly.matrixForm` → a `sym` expression → `domain`+`symbolicFunction`+`plq_1p`) and
   `evalFunctionNDomain.m` (numeric eval of a `functionNDomain` result at a dual point, for this
   codebase's standard numeric-sup-sampling test convention — no conversion back to `QuaPar` yet,
   deliberately, since that would reintroduce the exact curved-edge assembly problem `cPLQ` already
   solves symbolically). Validated in `cplqAdapterTest.m`: a single one-convex-edge triangle
   through `.conjugate` matches CCA2's own existing numeric `conjPieceCPLQ` exactly; a genuinely
   multi-triangle nonconvex PLQ (`f=xy` over a square split into 2 independent triangles — the case
   `maxQuaPar` itself cannot do numerically) through `.triangulate`→`.maximum` matches the true
   numeric sup at 9 test points (one exact symmetric tie point excluded and documented — a
   narrow `mergeL`/`removeTangent` limitation, not an adapter bug).
   **Confirmed timing/quality**: envelope+conjugate+max are correct and complete for every case
   tried; `cPLQ`'s own code needed no math fixes for these steps (only the 2 toolbox-compat bugs
   above). **Biconjugate (`plq.biconjugateF`) is NOT yet reliable — 2 distinct bugs found so far**:
   (1) `region.poly2orderUnbounded` had a reproducible unhandled-loop-fallthrough bug (index
   overflow) — **FIXED** (new session, 2026-07-18): the caller's `nv~=size(ineqs,2)` dispatch
   signal for "call `poly2orderUnbounded`" is not exclusive to genuinely unbounded regions (a
   BOUNDED region with redundant/tangent inequality rows triggers the same mismatch, same root
   cause as the `testRegion/testCreation` flakiness above); fixed by searching only up to `obj.nv`,
   returning early when no ray-adjacent vertex is found, and wrapping the "found at the last
   angular position" case via `mod` (same wraparound-fix pattern already used once in this
   codebase's own `maxQuaPar.m`). Verified: no regressions (`testPCE0`/`testPCE3`/`testPCE1`,
   `testcPLQ`'s 6 non-biconjugate tests, `cplqAdapterTest`'s 2 tests all still pass), and the fix
   demonstrably lets `testMaxMultiRegion/testMax` proceed further than before (confirmed by hitting
   a LATER, different error next). (2) `region.getNormalConeVertexQ` throws
   `MATLAB:catenate:dimensionMismatch` (horzcat of inconsistently-sized candidate-point arrays) —
   **NOT fixed yet**, found immediately after fixing (1), not deeply diagnosed. See
   `.claude/SESSION_HANDOFF.md` for the full trace and the open prioritization decision (keep
   fixing biconjugate bugs one at a time vs. pivot to the `mergeL` tie-point issue vs. wire
   Steps 1-3 into `conjCPLQ.m` now and defer both).
   **Wired into `conjCPLQ.m` itself (this session, user's choice of 3 options)**: `conjCPLQ.m` now
   has a Case C (general bounded domain, `nf>1` and/or non-triangular) that calls `quaPolyToPlq`
   -> `.triangulate` -> `.maximum` instead of erroring -- generalized `quaPolyToPlq.m` to accept
   non-triangular faces too (relies on `plq.triangulate`'s own fan-splitting). `g` for Case C is a
   `functionNDomain` array (NOT `QuaPoly`/`QuaPar`) -- documented in `conjCPLQ.m`'s own header;
   composition (`biconj`/`infConv`/`moreau`/...) is not supported for Case C yet for that reason.
   Verified via `conjCPLQTest.m`'s new `multiFaceBoundedDomainViaCPLQIntegration` (through the
   actual public `conj('cplq')` entry point) and a full-suite regression check (18/18 pass).
   **Still open for Phase 1**: fix `getNormalConeVertexQ` (and whatever else surfaces after it) so
   biconjugate is reliable, and the `mergeL`/`removeTangent` exact-tie-point gap -- both deferred
   per the user's choice this session, unaffected by the Case C wiring (Case C never calls
   biconjugate); convert `evalFunctionNDomain`'s output back into a `QuaPar` if/when a caller needs
   the structured (not just numerically-evaluable) result, so Case C can support composition too;
   run the remaining ~14 untested `testMaxMultiRegion` cases. `cPLQ`'s own code being slow and
   noisy (`maximumP`/`mergeL`/`removeTangent`'s repeated symbolic `isAlways` "truth unknown"
   warnings/retries) remains
   expected and fine for Phase 1 — it is exactly Phase 2's target, not a Phase 1 blocker.
2. **Phase 2 (later) — improve performance**: once Phase 1 is integrated and fully tested, replace
   the symbolic computation with direct closed-form numeric formulas incrementally, one case/step
   at a time, validating each replacement against the Phase-1 symbolic result on the same inputs
   before moving to the next (the same "derive once, implement in plain arithmetic" pattern already
   used for e.g. `conjBilinearXYoneCE`/`conjPSDRank1QuadTriangle`) — go carefully and slowly,
   comparing against the working (Phase 1) code at each step, rather than re-deriving the whole
   pipeline's math from scratch in one jump (the failure mode of the prior two sessions' abandoned
   `conjPieceCPLQ` rational-piece attempt).
3. ~~**`conjPieceCPLQ`'s own rational-piece TODO**~~ — superseded by the plan above; do not resume
   the single-RatPol-piece closed-form derivation as a standalone task (see `conjPieceCPLQ.m`'s
   header and `.claude/SESSION_HANDOFF.md` for the full diagnosis of why it was mis-scoped).
4. **`partialConj`** for the `'cplq'`/`'pqp'` engines (II.4) — not started; `cPLQ`'s own
   `conjugateExpr`/Lagrange-multiplier machinery already supports a single-variable restriction, so
   this may fall out of the Phase 1 integration for the `'cplq'` engine specifically.
5. **`add` for `RatPol`** (common-denominator sum) and the **`RatPar`** parent class (II.3) —
   deprioritized, nothing in the operator pipeline calls for either.

---

## 0. TL;DR — what to build

1. **Anchor data structure:** keep the released `PLQVC` mesh (`V/E/f/F` + adjacency `P`
   + `dom`). It already represents quadratic-on-polyhedral functions and all degenerate
   1D cases, and ships the linear-time convexity test.
2. **Three function families must be storable** (per [COAP]/[JOGO]):
   - **`QuaPoly`** — quadratic on a *polyhedral* subdivision (the input; the released type).
   - **`RatPol`** — quadratic÷linear on a *polyhedral* subdivision (the convex envelope `f**`).
   - **`QuaPar`** — quadratic on a *parabolic* subdivision (the conjugate `f*`).
   Organize these as a **class hierarchy** whose most general parent **`RatPar`** is *rational
   (cubic÷linear) on parabolic* — a pure storage umbrella the user may populate with **any** such
   function; the three *operated* families (**`RatPol`**, **`QuaPar`**, **`QuaPoly`**) are its
   specializations, ending at the released quadratic-on-polyhedral class.
3. **Three conjugate engines** behind one `conj(f, engine)`:
   - **`'pqp'`** — exact parametric-QP conjugate (Scilab `pQP` port; Jakee Khan's M.Sc. thesis
     [JAKEE-13]). **Convex `f` only.**
   - **`'graph'`** — point-cloud + neighbour-graph conjugate (`PLQVG` primal↔dual swap;
     Tasnuva Haque [HAQUE-17]/[HAQUE-18]). **Convex `f` only.**
   - **`'cplq'`** — symbolic per-piece conjugate (`cPLQ` Lagrange-multiplier method) — this is
     the engine that follows the [JOGO] 3-step algorithm exactly, and the only one of the three
     that handles **nonconvex** `f`.
4. **Convex envelope = biconjugate**, built *on top of the `'cplq'` conjugate engine*
   (`convEnv = biconj = conj∘conj`, `engine='cplq'`). This is how `cPLQ` does it natively; it does
   **not** work over `'pqp'`/`'graph'` for a genuinely nonconvex `f` (they only accept convex
   input, so the first conjugation is undefined — II.5/II.6). A **direct** envelope path (Deepak
   Kumar's per-piece method [KUMAR-20], `codeOld/deepak`, extended by Karmarkar's [COAP] Step 4
   assembly — **not** parametric-QP) is kept as a faster alternative and an independent
   cross-check.
5. **Everything else** (`add`, `lasryLions`, `partialConj`) is a short composition over `conj` +
   `add` + scalar ops. `infConv(f,g)=conj(conj(f)+conj(g))` is this same kind of composition but
   is only valid (returns the true inf-convolution) for `f,g` convex (II.6). `moreau` looks
   similar but is **not** built on `infConv`: it is one direct conjugate call via the "expand the
   square" identity of Hiriart-Urruty & Lucet [HIRIART-URRUTY-07], which needs no convexity
   assumption at all. `proxAverage` is **not** a composition of `moreau` calls either: it reduces
   directly to `conj` — a sandwich of two conjugations around a weighted `add` — and, like
   `infConv`, is only valid for `f,g` convex (II.6).

> **The pillars:** conjugate = {`pQP` (Jakee Khan, convex), `PLQVG` (Tasnuva Haque, convex),
> `cPLQ` (Karmarkar, convex + nonconvex)} · convex envelope = `biconj` with `engine='cplq'`
> (+ Kumar/Karmarkar `direct`). All other operators compose on these — except `infConv` and
> `proxAverage` (both convex-only, both reduce to `conj`+`add`) and `moreau` (a single direct
> conjugate, not `infConv`; convexity-free) — see II.6.

---

## PART I — Inventory of the existing code

### I.1 Released code (CCA2, the last release)

| File | Type | Description |
|------|------|-------------|
| `CCA2/PLQVC.m` | `classdef` (value) | **The anchor.** 2D piecewise linear-quadratic/cubic function. Domain = vertex/edge/face mesh (`V`,`E`,`F`) + face adjacency list `P` + `dom` struct. Function = `nf×10` cubic coefficient matrix `f` (quadratic uses 6). Methods: constructor, `eval`, `createP`/`createDom`/`orderEdges` (topology), `isConvex`/`isFaceConvex`/`isEdgeConvex`/`isEdgeContinuous`/`isQuadratic`/`isDomBounded` (tests), `disp`/`plot`/`plotDomain`/`plotDomainGraph`/`plotVertexes`. Static helpers: `evalHessian`, `evalPoly`, `matrixForm`, `lineEquation`, `isPositiveSemidefinite`, `isCubicConvexOnEdge`, `belongToEdge`, `isCollinear`. Static factories: `examples`, `examples2`, `examplesNonconvex`, `examplesDiscontinuous`, plus `oneNorm`, `oneNormConjugate`, `energy`, `cubic1`. |
| `CCA2/PLQVCTest.m` | test class | Unit tests (`matlab.unittest`). |

**Naming convention** (class header): **P**iecewise **L**inear-**C**ubic, **V** = vertex
representation of the polyhedral pieces, **C** = coefficient representation. Variants:
**H** = half-space, **P** = pointwise, **G** = graph (pointwise) function rep, **Q** =
quadratic vs cubic. The release ships the data structure + convexity tests; **no operators**.

### I.2 `codeOld/plq2-dataStructure` — competing MATLAB data structures

| Class / file | Domain rep | Function rep | Key methods / operations |
|--------------|-----------|--------------|--------------------------|
| **`PLCVC`** | vertex (V/E/F/P) | coefficient (cubic) | `eval`, `isConvex`, `isFaceConvex`, `isEdgeConvex`, `isDD`, `toPLCHC`, `toPLQVG`. **Direct ancestor of `PLQVC`.** |
| **`PLCHC`** | half-space per piece | coefficient | `lineEquation`, `toPLCVC`, **`add`** (polyshape intersection + coeff sum — working pointwise addition), `isFaceBounded`. |
| **`PLQVG`** | vertex via `Entity` | **graph/pointwise** (values `y` + subgrads `s`) | **`conjugate`** (primal↔dual entity swap: Face↔Vertex, Ray↔Ray, Segment↔Segment), `toPLCVC`, `toPLCVCDual`. **= the point-cloud + neighbour-graph conjugate engine**, implementing **Tasnuva Haque**'s entity-graph algorithm [HAQUE-17]/[HAQUE-18] (this specific code is an unpublished refinement by Hatton/Kagdiwala/Patodia, not Haque's own implementation). **Convex PLQ only** (per [HAQUE-17]/[HAQUE-18]'s own scope). |
| **`PLQList`** | inequality list | quadratic coeffs | **`computeDual`** (symbolic conjugate), **`computeMoreauEnvelope`**, `scalarMultiply`, `additionOfPLQAndQuadraticFunction`. |
| `Entity`/`Entitytype` | — | — | Primal/dual graph cell (type ∈ {Vertex,Face,Ray,Segment,Line}); enables the primal-dual conjugate. |
| `lcon2vert`,`vert2lcon`,`qlcon2vert` | — | — | **Reusable:** half-space ↔ vertex conversion. |
| `unionHull`,`intersectionHull`,`addBounds`,`separateBounds` | — | — | Convex hull of union / intersection of polyhedra (used by `add`). |
| `Hausdorff` | — | — | Polygon Hausdorff distance (validation metric). |

### I.3 `codeOld/cPLQ` — symbolic per-piece conjugate & convex envelope (→ Engine `'cplq'`)

| File | Description |
|------|-------------|
| `plq.m` | Container of `plq_1p` pieces; orchestrates `convexEnvelope`, `maximumConjugate`, `biconjugateF`. |
| `plq_1p.m` | One PLQ piece. `convexEnvelope` (closed form ≤2 convex edges incl. harmonic-mean curvature `a=(mₕmᵥ)/(mₕ+mᵥ+2√(mₕmᵥ))`; iterative merge otherwise), `conjugateFunction`, `maximumConjugate`. |
| `conjugateExpr.m` | **Conjugate of a piece on an edge** via Lagrange multipliers `sup{ s·x − f(x) : g(x)=0 }`. |
| `functionNDomain.m` | (function, region) pair; `maximumP` (pointwise max via region split), subdifferential at vertices/edges, domain merge. |
| `region.m` | Region as `≤0` inequalities (linear **or quadratic/parabolic**); vertex enumeration, normal cones, `splitmax3`, tangent removal. |
| `quadQuad.m`,`qq_conj.m`,`quadLinear.m`,`conjQuad.m` | Closed-form conjugate building blocks for the **bilinear `x*y` / quad-on-quad** case. |
| `domain.m`,`symbolicFunction.m`,`checkConvex.m` | Triangular domain metadata; symbolic wrapper; convexity check. |

Implements **conjugate, convex envelope, biconjugate, subdifferential, max** symbolically.
This is the engine that mirrors [JOGO]'s 3-step pipeline and handles parabolic regions and
rational pieces directly.

### I.4 `codeOld/deepak` — Deepak Kumar's convex-envelope front-end ([KUMAR-20]; predecessor to [COAP]/[JOGO])

**Not a parametric-QP method** — a previous version of this document mislabeled it as one,
confusing it with the unrelated `'pqp'` engine (I.5/[JAKEE-13]). Deepak Kumar's code
([KUMAR-20]) is a **numeric, eigenvalue/classification-based** method for the convex envelope of
*one* quadratic over *one* polytope; **Tanmaya Karmarkar**'s [COAP]/[JOGO] algorithms extend it
into the full symbolic, multi-piece pipeline (`codeOld/cPLQ`, I.3) that this design is built on.

| Entry point | Description |
|-------------|-------------|
| `cvxEnv2d` | **Convex envelope of a quadratic over a polytope** ([COAP] Step 1+4): indefinite→bilinear `xy`, classify convex edges vs saddle vertices, build `η(a,b)` supporting families, subdivide dual space, closed-form `η(a,b)`-family solve per subregion, map back. Output = rational(quad÷linear)/polyhedral. |
| `cvxHull2d_solver`, `OP_Solver`, `compute_OP*`, `solve_OP_*` | Parametric 1D optimization solver core (closed-form QQ/QL/LL cases) — a per-edge 1D parametrization, unrelated to the multivariate KKT `'pqp'` engine of I.5. |
| `Conjugate/compute_conjugate*` | Conjugate of quadratic/rational pieces; domain division into parabolic subdivision. |
| `transform_quadxp_to_xy` | Bilinear handling via Hessian eigendecomposition. |

Used as the **direct convex-envelope engine** and a numeric cross-check.

### I.5 `codeOld/MP_plq2-lft` + `COAP_plq2-lftpartial` — Scilab operator algebra (the blueprint)

The Scilab toolbox where **one conjugate generates the whole operator algebra**:
`plq2_lft` (full conjugate), `plq2_lft_partial` (partial), `plq2_me`/`plq2_me2`/`plq2_me_partial`
(Moreau, two formulas), `plq2_pa`/`plq2_pa_partial` (proximal average), `plq2_add`/`plq2_subtract`/
`plq2_scalar`, `plq2_s` (self-smoothing). PLQ2 = arrangement (DCEL) + `6×n` coeff matrix.
`macros/lft2.sci` (the main conjugate) and the partial-conjugate algorithm are **Bryan
Gardiner**'s work ([GARDINER-13], [GARDINER-14] — the latter co-authored with **Jakee Khan**).
Nested inside this same tree, `plq2_lft_pQP/parametric_ME.sci` + `…/partial_conjugate.sci` are a
**different, separate** method — **Jakee Khan**'s **parametric-QP** engine (KKT solve, symbolic
in the dual parameter, [JAKEE-13]) → Engine `'pqp'`. Both are **convex-only** ([GARDINER-13]/
[GARDINER-14]/[JAKEE-13] all state "convex" in their titles/scope). The example scripts
(`l1-norm-lft+me`, `l1-linf-norm-pa`, `house-lft`, `diamond-lft`, `four-faces-lft+me`) define
the intended public API.

### I.6 `codeOld/JCA_fitz`, `JCA_plt-scripts`, `NA_FMEI-scipts` — Scilab

`JCA_fitz/gph.sci`: **Fitzpatrick function** `F_{m,A}` as a PLQ, per **Bryan Gardiner**'s
[GARDINER-09]; `plq_pds` (primal-dual symmetric
antiderivative) via `plq_lft`+`plq_pa`; `subdiff_plq`, `diff_plq`. `JCA_plt-scripts` & `NA_FMEI`:
1D grid/LLT Moreau-envelope demos/timing (`me_llt/pe/nep/plt/direct`) on `abs`,`sqr`,indicators,
nonconvex `|‖x‖−1|`,`(x²−1)²`. **The 1D grid scripts are out of scope** (no grid algorithms);
keep only the **Fitzpatrick** routine and the **nonconvex test functions** as example generators.

### I.7 `codeOld/polyPak` — Fortran-90 SOS / polynomial-positivity prover — **OUT OF SCOPE.**

---

## PART II — Proposed code design

### II.1 Design principles

1. **One data-structure family, organized as a hierarchy.** The released `PLQVC` mesh is the
   substrate. Generalize it along two independent axes — *boundary type* (linear ⊂ parabolic)
   and *function type* (polynomial ⊂ rational) — and expose the combinations the math needs as
   a class hierarchy (II.3). The most general parent stores everything (incl. cubics, usable by
   `isConvex`); the leaf is today's quadratic-on-polyhedral class, unchanged.
2. **Operators compose.** Implement two primitives well — **`conj`** and **`add`** (+ trivial
   `scalarMul`/`negate`/`addQuadratic`) — then *derive* `infConv`, `moreau`, `proxAverage`,
   `lasryLions`, and `convEnv`(=`biconj`). Mirrors the Scilab `MP_plq2-lft` design.
3. **Three conjugate engines, one interface.** `conj(f,'pqp'|'graph'|'cplq')`, same in/out type.
   Every derived operator inherits all three engines for free.
4. **Storage is general; operators are restricted.** The user may store **any `RatPar`**
   function (rational cubic÷linear on parabolic, including cubic numerators). **Operators apply
   only to the specific subclasses** `QuaPoly`/`QuaPar`/`RatPol`: the three operated families
   round-trip under conjugation (II.2), but a bare `RatPar` (rational-on-parabolic, or any cubic
   numerator) is **rejected** by `conj`/`add`/`moreau`/… The **one exception is `isConvex`** (and
   the convexity sub-tests), which **is allowed on cubic polynomials** — the released
   `isCubicConvexOnEdge`/`evalHessian` machinery already supports it.
5. **1D = degenerate 2D.** Needle/segment/ray/chain already first-class in `dom`; one code path.
6. **Match CCA2 conventions:** value `classdef`; `% objective/[input]/[output]` doc blocks;
   `methods(Static)` helpers + `examples*` factories; camelCase. **Operators use exact
   symbolic + rational arithmetic by default** (per [COAP]/[JOGO], to avoid degenerate
   floating-point subdivisions); `sqrt(eps)` tolerances only in geometric predicates/tests.

### II.2 Function families and the type cycle

Per [COAP] §2–3 and [JOGO] §2–3, a piece is defined by a **function type** and a **domain
boundary type**:

| Class | Function on a piece | Domain boundary | Role | Operator input? |
|-------|---------------------|-----------------|------|-----------------|
| **`QuaPoly`** | quadratic `½xᵀQx+qᵀx+κ` | polyhedral (lines) | input `f` (released `PLQVC`) | ✅ yes (`isConvex` also accepts cubic) |
| **`RatPol`** | `(quad)÷(linear)` `r=(ax²+bxy+cy²+dx+ey+f)/(gx+hy+k)` | polyhedral | convex envelope `f**` = `conv f` | ✅ yes |
| **`QuaPar`** | quadratic, often parabolic `qₚ` (`a>0,c>0,b²−4ac=0`) | parabolic (`b²−4ac=0` conics) | conjugate `f*` | ✅ yes |
| **`RatPar`** (parent) | `(cubic)÷(linear)` | parabolic | storage umbrella for **any** of the above | ❌ storable, **not** an operator input |

**Parabola** ([COAP] Def.): `{ax²+bxy+cy²+dx+ey+f=0}`, coeffs not all zero, **b²−4ac=0**
(lines are the case `a=b=c=0`). **Parabolic region**: intersection of such inequalities `≤0`.
**Parabolic subdivision**: union of parabolic regions whose pairwise intersections lie in a
parabola (polyhedral subdivision is the special case of all-linear).

**Conjugate of a PLQ has exactly four piece kinds** ([JOGO]/[COAP] §3): (1) linear/polyhedral,
(2) linear on the convex side of a parabola, (3) parabolic function on the nonconvex side of a
parabola, (4) parabolic function on an unbounded 3-edge polyhedral region. All are covered by
the *parabolic* family.

**The type cycle** (closed under conjugation — this is why three families suffice):

```
              conj                          conj
  QuaPoly  ───────────▶  QuaPar (f*)  ───────────▶  RatPol (f**)
 (quad/polyhedral)            │                          │
        ▲                     │ conj                     │ conj
        │                     ▼                          ▼
        └──────────────  RatPol  ◀───────────────  QuaPar
            biconj = convEnv     (f** = conv f, quad÷linear / polyhedral)
```

- `conj(QuaPoly)` → `QuaPar`  ([JOGO]: `f*` is quadratic on parabolic).
- `conj(QuaPar)` = `f**` → `RatPol`  ([COAP] Step 4: `conv f` is quad÷linear on polyhedral).
- `conj(RatPol)` = `conj(conv f)` = `f*` → `QuaPar` (conjugate of the envelope = conjugate).
- `convEnv(f)` = `biconj(f)` = `conj(conj(f))` → `RatPol`.

### II.3 Class hierarchy

```
              ┌──────────────────────────────────────────────────┐
              │  RatPar   (abstract storage umbrella)              │
              │  rational (cubic÷linear) on parabolic              │
              │  data: V, E, Ec(conic), num(:,10), den(:,3), F,P,dom│
              │  capabilities: eval, plot, convert, isConvex; NO ops│
              └───────────────┬───────────────────┬────────────────┘
        specialize: den=1     │                   │  specialize: linear
        (polynomial)          ▼                   ▼  edges (polyhedral)
        ┌──────────────────────────┐   ┌──────────────────────────┐
        │ QuaPar                    │   │ RatPol                    │
        │ quadratic on parabolic    │   │ quad÷linear on polyhedral │
        │ = conjugate f*            │   │ = convex envelope f**     │
        └────────────┬─────────────┘   └──────────────────────────┘
                     │ specialize: linear edges (polyhedral)
                     ▼
        ┌────────────────────────────────────┐
        │ QuaPoly   (= released PLQVC)         │
        │ quadratic on polyhedral             │
        │ input type; linear-time isConvex    │
        └────────────────────────────────────┘
```

- **Inheritance = mathematical specialization (is-a).** A quadratic-on-polyhedral function *is a*
  quadratic-on-parabolic function (a line is a degenerate parabola); a quadratic *is a* rational
  with `den=1`. Each subclass tightens an invariant and may add specialized methods.
- **`RatPar`** (parent): holds the union of all fields — a `num(:,10)` cubic-numerator matrix,
  a `den(:,3)` linear-denominator matrix (`den=[0 0 1]` ⇒ polynomial), and edges that may be
  **conic arcs**. Provides `eval`/`plot`/converters **and `isConvex`** (which accepts cubics);
  **no convex-analysis operators**. The user may instantiate and store any `RatPar`.
- **`RatPol`** (`f**`), **`QuaPar`** (`f*`), **`QuaPoly`** (input): the three *operated*
  families. `QuaPoly` is the released `PLQVC` essentially verbatim (the leaf class).
- **Reject rule:** every operator (`conj`,`add`,`moreau`,…) calls `obj.assertOperable`, which
  errors unless `obj` is one of `QuaPoly`/`QuaPar`/`RatPol` with a **quadratic** numerator;
  a bare `RatPar` or any cubic numerator ⇒ `error('PLQ:op:unsupportedType', …)`. **`isConvex`
  is exempt** — it accepts cubic polynomials (released convexity machinery).

#### Generalized properties (added to the released `PLQVC` fields)

| Property | Meaning |
|----------|---------|
| `V (:,2)` | vertices (unchanged) |
| `E (:,3)` | edges `[v1 v2 isSegment]` (unchanged for **linear** edges) |
| `Ec (:,6)` *(new)* | per-edge conic coefficients `[a b c d e f]` of the boundary `ax²+bxy+cy²+dx+ey+f=0` with `b²−4ac=0`; **all-zero row ⇒ the edge is the straight line through its endpoints** (polyhedral case). Curved edges are arcs of this parabola between `V(v1)` and `V(v2)`. |
| `num (:,10)` | per-face numerator coefficients in the released cubic basis (quadratic uses last 6; cubic allowed for `isConvex` only) |
| `den (:,3)` *(new)* | per-face denominator `[g h k]` for `g x + h y + k`; `[0 0 1]` ⇒ polynomial (PLQ/parabolic), nonzero `g,h` ⇒ rational |
| `F (:,2)`, `P {cell}`, `dom` | left/right face per edge, ordered adjacency, domain descriptor (unchanged) |
| `degree` *(cached)* | 1/2/3 numerator degree; lets operators reject cubics |

This is a **strict superset** of `PLQVC`: with `Ec` all-zero, `den=[0 0 1]`, and `num` quadratic,
the object *is* a released `PLQVC` and every existing method works unchanged.

### II.4 Public method catalogue

| Method | Signature | Operation | Built on |
|--------|-----------|-----------|----------|
| `conj` | `g = conj(f, engine)` | Fenchel conjugate `f*` | engine ∈ {`'pqp'`,`'graph'`,`'cplq'`} (default `'cplq'`); **`'pqp'`/`'graph'` require `f` convex, `'cplq'` handles both** (II.5) |
| `partialConj` | `g = partialConj(f, idx, engine)` | conjugate w.r.t. variable `idx` | engine ∈ {`'pqp'`,`'cplq'`} only |
| `biconj` | `g = biconj(f, engine)` | `f** = conj(conj(f))` | `conj∘conj` — nonconvex `f` needs `engine='cplq'` for the first `conj` (II.5/II.6) |
| `convEnv` | `h = convEnv(f, method)` | convex envelope; `method`=`'biconj'`(default `'cplq'`) \| `'direct'`(Kumar+Karmarkar, [KUMAR-20]/[COAP]) | `biconj` / `convEnvDirect` — see II.6 note on engine choice |
| `add` | `h = add(f, g)` | pointwise `f+g` | domain overlay + coeff/rational add — **implemented for `QuaPoly` and `QuaPar`** (`addQuaPoly.m`/`addQuaPar.m`); `RatPol` still open, see Implementation status above |
| `sub` | `h = sub(f, g)` | `f − g` | `add(f, negate(g))` |
| `scalarMul`,`negate` | `c·f`, `−f` | coeff scaling | — implemented on all three classes |
| `addQuadratic` | `addQuadratic(f, A,b,c)` | `f + (½xᵀAx+bᵀx+c)` | per-face coeff update |
| `infConv` | `infConv(f,g,engine)` | `(f□g)(x)=inf_z f(z)+g(x−z)` | `conj((conj f)+(conj g))` — **valid (returns the true inf-conv) only for `f,g` convex**; see II.6 |
| `moreau` | `moreau(f,mu,engine)` | Moreau envelope `e_μ f` | **single conjugate** (`conj`, once) via the expand-the-square identity [HIRIART-URRUTY-07] — **NOT** `infConv(f,½μ‖·‖²)`; see II.6 |
| `lasryLions` | `lasryLions(f,lambda,mu,engine)` | double envelope | `−moreau(−moreau(f,λ),μ)` |
| `proxAverage` | `proxAverage(f,g,lambda,mu,engine)` | proximal average | **reduces directly to `conj`** (two conjugations sandwiching a weighted `add`, II.6) — **not** a composition of `moreau` calls; convex `f,g` only, like `infConv` |
| `eval`,`isConvex`,`plot`,… | (unchanged) | evaluation / tests / display | existing `PLQVC` code |

Static factories: keep `QuaPoly.examples*`,`oneNorm`,`energy`; add `QuaPoly.l1Norm`,`linfNorm`,
`indicator`, `nonconvexW`, and the [JOGO]/[COAP] worked examples (Fig. 2 functions) as a
golden test set.

### II.5 The three conjugate engines

```matlab
function g = conj(obj, engine)
% objective: Fenchel conjugate of an allowed 2D PLQ-family function
% [input]  obj    : QuaPoly | RatPol | QuaPar  (bare RatPar & cubic rejected)
%          engine : 'cplq' (default, symbolic, exact) | 'pqp' | 'graph'
% [output] g      : the conjugate (type per the II.2 cycle)
    obj.assertOperable();                 % reject cubic / rational-on-parabolic
    if nargin<2, engine='cplq'; end
    switch lower(engine)
        case 'cplq',  g = conjCPLQ(obj);   % II.5.1 — symbolic per-piece ([JOGO] 3 steps)
        case 'pqp',   g = conjPQP(obj);    % II.5.2 — parametric QP / KKT
        case 'graph', g = conjGraph(obj);  % II.5.3 — point cloud + neighbour graph
        otherwise, error('PLQ:conj:engine','Unknown engine "%s"',engine);
    end
end
```

#### II.5.1 Engine `'cplq'` — symbolic per-piece conjugate (default; the [JOGO] algorithm)

*Source:* integrated (2026-07-18, Phase 1 — see "Next planned" above) from the reference `cPLQ`
package, almost as-is: `plq.m`, `plq_1p.m`, `functionNDomain.m`, `region.m`, `domain.m`,
`symbolicFunction.m`, `conjugateExpr.m`, `yIntercept.m` now live at the CCA2 repo root (copied
from the untracked `cPLQ/` reference clone, which stays untracked/unused going forward). These 8
files are `plq.m`'s own runtime dependency closure, confirmed by grep. `plq_1piece.m` (2603 lines)
was copied in too, separately: it is an OLDER, largely-parallel per-piece implementation with the
same method names as `plq_1p.m` (`triangulate`/`convexEnvelope`/`conjugateFunction`/
`maximumConjugate`/`biconjugateP`) but is never called by `plq.m` itself (which only ever
constructs `plq_1p` objects) — it is however the class `testMaxMultiRegion.m`'s 24 tests
(including 8 `testBiconjugate*` cases) exercise directly, so it earns its place as "relevant"
despite being dead code from `plq.m`'s own point of view; treat it as a second, independently-
tested implementation of the same per-piece algorithm, not a stale leftover to delete.
`quadQuad.m`/`qq_conj.m`/`quadLinear.m`/`conjQuad.m`/`conjR.m`/`checkConvex.m`/
`subdiffParabolicEdge.m`/`vertexNan.m`/`quadquad1.m`/`plots.m`/`plot2.m` remain standalone offline
derivation *scripts* with no live call sites from any of the 9 files above or their test suites —
confirmed by grep across all of them, not just the production files — and were NOT copied in
(correcting this section's earlier "quadQuad/qq_conj" source claim below, which was aspirational/
inaccurate). The existing test suites (`testcPLQ.m`, `testfunctionNDomain.m`, `testRegion.m`,
`testSymbolicFunction.m`, `testMaxMultiRegion.m`) were copied alongside the code and mostly pass
as-is (see `.claude/SESSION_HANDOFF.md` for the one known exception and timings).

**Method** — exactly the three steps of [JOGO]/[COAP]:
1. **Convex envelope of each quadratic piece** `conv(qᵢ+I_{Pᵢ})` → rational(quad÷linear) on a
   polyhedral subdivision (`plq_1p.convexEnvelope`).
2. **Conjugate of each rational piece** via Lagrange multipliers `sup{s·x − r(x) : g(x)=0}`
   (`conjugateExpr`, called from `plq_1p.conjugateFunction`) → quadratic on a parabolic
   subdivision.
3. **Maximum of the conjugates** (`functionNDomain.maximumP`, region splitting) → `f*` as a
   quadratic-on-parabolic function. (Computing the conjugate of *that* again gives `f**`.)

**Why default:** it is exact (symbolic + rational arithmetic, avoiding the degenerate-subdivision
floating-point issues [COAP]/[JOGO] explicitly warn about), and it naturally produces and
consumes all three families — so `biconj`, `moreau`, etc. compose through it without restriction.
`partialConj` is supported here (conjugate restricted to one variable, symbolic Lagrange solve).

#### II.5.2 Engine `'pqp'` — parametric quadratic programming (exact)

*Source:* **Jakee Khan**'s M.Sc. thesis [JAKEE-13] (Scilab `parametric_ME.sci`); partial variant
co-authored with Bryan Gardiner [GARDINER-14] (`partial_conjugate.sci`); closed forms adapted
from Gardiner's `lft2.sci` [GARDINER-13].

**Convex PLQ only** — [JAKEE-13] computes the KKT stationarity conditions of
`v(s)=sup_x sᵀx − r(x) s.t. Ax≤β` per primal entity; this characterizes the true supremum only
when the maximization is concave, i.e. `r` convex, so `conj(f,'pqp')` requires `f` convex
(`assertOperable` should reject/flag a nonconvex `f` here rather than silently returning a
KKT critical point that isn't the actual conjugate).

For each primal entity collect `(Q,q,κ)` and active constraints `Ax≤β`, then solve, symbolically
in the dual parameter `s`, via KKT
`[Q Aᵀ; A 0][x;λ]=[s−q; β]`. Optimizer `x*(s)` is affine in `s`; back-substitution gives the
piece value; the active-set critical region is the dual cell. Assemble cells into the dual mesh.
`partialConj` = the same solve with `s` on selected coordinates only
(`_plq2_quad_face_partconj`/`_plq2_lin_face_partconj`). Reuse `PLQVC.matrixForm`/`lineEquation`/
`orderEdges` for `(Q,q)` and adjacency; build the result through the constructor.

#### II.5.3 Engine `'graph'` — point cloud + neighbour graph

*Source:* `PLQVG` + `Entity`/`Entitytype`, implementing **Tasnuva Haque**'s entity-graph
algorithm [HAQUE-17]/[HAQUE-18] (this codebase's specific `PLQVG` implementation is an
unpublished refinement by other students, not Haque's own code — see I.2).

**Convex PLQ only** — [HAQUE-17]/[HAQUE-18] are both titled "...of **convex** piecewise
linear-quadratic bivariate functions": the primal↔dual entity swap and the well-poised
subgradient sampling below both presume a genuine (single-valued-per-point, convex) subgradient
at each entity, which is not well defined for a nonconvex `f`. So `conj(f,'graph')` requires `f`
convex, same restriction as `'pqp'`.

Represent the graph of `f` as a **point cloud with a neighbour graph**: each entity carries a
point `x`, value `y=f(x)`, and subgradient(s) `s∈∂f(x)`; adjacency = the mesh neighbour graph.
The Legendre transform is a **primal↔dual swap**: dual point = subgradient `s`, dual value =
`sᵀx−f(x)`, entity types swap (Face↔Vertex, Ray↔Ray, Segment↔Segment); reconnect via the same
adjacency. Convert `QuaPoly→PLQVG` (gradients from `evalHessian`; vertex normal-cone fans from
`region`/`getSubdiffVertex*` give subgradients at kinks), `PLQVG.conjugate`, then back to the
class. Fast, geometric, no symbolic solve; great for many pieces, visualization, and as an
independent cross-check of the exact engines. **`partialConj` is not implemented for `'graph'`**
(use `'pqp'` or `'cplq'`).

> **Cross-validation oracle:** the three engines must agree (to `sqrt(eps)` on `eval`) on the
> convex test suite, and against the analytic pairs in II.8 — an internal correctness check
> analogous to FLTW's exact-backend validation.

### II.6 Convex envelope, and derived operators (thin compositions)

```matlab
function h = convEnv(f, method)             % convex envelope = closed convex hull = f**
    if nargin<2, method='biconj'; end
    switch lower(method)
        case 'biconj', h = biconj(f);           % conj∘conj — see NOTE below on engine choice
        case 'direct', h = convEnvDirect(f);    % Kumar's per-piece envelope [KUMAR-20], assembled
                                                 % across pieces per [COAP] Step 4
    end
end
% Output is RatPol (quad÷linear on polyhedral), per [COAP] Step 4 / [JOGO] Prop 1.
```

**Which engine actually computes the envelope of a nonconvex `f`?** Only `'cplq'` (II.5.1) — the
`'pqp'` and `'graph'` engines (II.5.2/II.5.3) are convex-only by construction of the algorithms
they implement ([JAKEE-13], [HAQUE-17]/[HAQUE-18]), so `conj(f,engine)` for a nonconvex `f` is
not a meaningful operation for those two. Consequently:
- `convEnv(f,'biconj')` (default `engine='cplq'`) is the general nonconvex→envelope path: `f*`
  (via `conjCPLQ`'s own internal per-piece envelope + conjugate + max, II.5.1 steps 1–3) is
  already convex, so the second conjugation is an ordinary (convex-input) conjugate in any engine.
- `convEnv(f,'biconj', 'pqp')` / `('graph')` only make sense when `f` is **already convex**
  (where `biconj(f)=f` is a round-trip identity check, not an actual envelope computation) —
  using them on a genuinely nonconvex `f` is invalid, since the *first* conjugation is undefined.
- `convEnv(f,'direct')` never conjugates at all: it is Deepak Kumar's numeric per-piece method
  ([KUMAR-20], `codeOld/deepak`, I.4), with cross-piece assembly per [COAP] Step 4 — a genuinely
  different, faster route, valid for nonconvex `f` by construction (that is the whole point of
  [COAP]). Kept as an independent cross-check against `'biconj'`+`'cplq'`.

```matlab
function h = infConv(f,g,engine)            % (f □ g)
    h = conj( add( conj(f,engine), conj(g,engine) ), engine );
end
```

**Only valid (returns the true `f□g`) when `f` and `g` are both convex.** The identity
`(f□g)* = f*+g*` holds for *any* `f,g` (a pure sup-interchange, no convexity needed), so this
code actually computes `conj(f*+g*) = (f□g)**`, the **biconjugate** of the inf-convolution. By
Fenchel–Moreau, `h**=h` iff `h` is convex, proper and lsc — guaranteed when `f,g` are both convex
(inf-convolution of two convex functions is convex). For nonconvex `f` or `g`, `infConv` as coded
instead returns `conv(f□g)`, the convex envelope of the true inf-convolution, **not** `f□g` itself.

```matlab
function h = moreau(f,mu,engine)            % Moreau envelope e_μ f   (cf. plq2_me)
    %   e_μ f = (1/(2μ))‖·‖² − (1/μ)·conj( μ·f + ½‖·‖² )       [HIRIART-URRUTY-07]
    g = addQuadratic( scalarMul(f,mu), eye(2), [0;0], 0 );
    h = addScaledEnergy( scalarMul( conj(g,engine), -1/mu ), 1/(2*mu) );
end
```

**Deliberately NOT `infConv(f, ½μ‖·‖²)`**, even though `e_μf = f □ (1/(2μ))‖·‖²` by definition —
composing through `infConv` would (per the note above) silently require `f` convex. Instead this
expands the square in the Moreau-envelope definition and regroups to make a **single** conjugate
appear: `e_μf(x) = inf_z[f(z)+(1/(2μ))‖x−z‖²] = (1/(2μ))‖x‖² − sup_z[x·z−(μf(z)+½‖z‖²)]/μ =
(1/(2μ))‖x‖² − (1/μ)·(μf+½‖·‖²)*(x)`, exactly the formula coded above. This is pure algebraic
regrouping of the definition of conjugate — no biconjugation, no Fenchel–Moreau argument, hence
**no convexity requirement on `f`** ([HIRIART-URRUTY-07] gives this identity in full). This is
also why `moreau` only ever calls `conj` **once**, unlike `infConv`'s two calls.

```matlab
function h = lasryLions(f,lambda,mu,engine)      % −e_μ( −e_λ f )
    h = negate( moreau( negate( moreau(f,lambda,engine) ), mu, engine ) );
end

function h = proxAverage(f,g,lambda,mu,engine)   % P_μ(f,g;λ,1−λ), λ∈[0,1]
    gf = addQuadratic( scalarMul(f,mu), eye(2), [0;0], 0 );  % T_μf := μf + ½‖·‖²
    gg = addQuadratic( scalarMul(g,mu), eye(2), [0;0], 0 );  % T_μg := μg + ½‖·‖²
    s  = add( scalarMul(conj(gf,engine),lambda), scalarMul(conj(gg,engine),1-lambda) );
    h  = addScaledEnergy( scalarMul(conj(s,engine),1/mu), -1/(2*mu) );
end
```

**Derivation of the `proxAverage` formula (reduces to `conj`, not to `moreau`).** Write
`T_μf := μf+½‖·‖²`, so `moreau`'s inner term is `T_μf` and `e_μf = (1/(2μ))‖·‖² − (1/μ)(T_μf)^*`
(II.6 above). The proximal average `P=P_μ(f,g;λ)` is *characterized* by its Moreau envelope being
the convex combination of the two envelopes: `e_μP = λe_μf+(1-λ)e_μg`. Substituting the identity
above into both sides and cancelling the `(1/(2μ))‖·‖²` term (since `λ+(1-λ)=1`) gives
`(T_μP)^* = λ(T_μf)^*+(1-λ)(T_μg)^*`. The right-hand side is **always** convex, proper and lsc
(a nonnegative combination of conjugates, and every conjugate is convex/lsc by construction) —
so, **provided `T_μP=μP+½‖·‖²` is itself convex** (guaranteed when `f,g` are convex and `μ>0`,
the paper's standing assumption), Fenchel–Moreau lets us take the conjugate of both sides and
recover `T_μP` exactly: `T_μP = conj(λ(T_μf)^*+(1-λ)(T_μg)^*)`. Solving for `P` gives the code
above: two conjugations (`conj(gf)`,`conj(gg)`) feed a weighted `add`, then a third conjugation
recovers `T_μP`, and `P` falls out after subtracting back `½‖·‖²` and rescaling by `1/μ`. Unlike
`moreau` (one conjugate call, no convexity needed — the identity is pure algebra with no
biconjugation step), this derivation's last step **is** a biconjugation, so — like `infConv` —
**`proxAverage` is only valid (returns the actual proximal average) for `f,g` both convex**; for
nonconvex input the formula instead computes the biconjugate of the intended object, same
caveat as `infConv` (II.6 above). This also means `proxAverage`'s only two new-code
dependencies are `add` (on `QuaPar`, the type `conj(gf)`/`conj(gg)` land in — now implemented,
see Implementation status) and `addQuadratic`/`addScaledEnergy`; it needs no call to the
`moreau` function at all, despite computing something moreau-envelope-flavored.

`add` is the other primitive needing real geometry: **overlay the two domain subdivisions and
add functions per overlapping piece** (sum quadratics, or add rationals over a common
denominator) — implemented for `QuaPoly`/`QuaPar` (`addQuaPoly.m`/`addQuaPar.m`, curved
boundaries handled by intersecting parabolic arcs via an axis-rotation parametrization, see
Implementation status); `RatPol` (common-denominator sum) remains open.

### II.7 File / module layout (pure MATLAB)

```
CCA2/
  RatPar.m            % abstract parent: rational(cubic÷linear) on parabolic; storage + eval/plot/convert/isConvex
  RatPol.m            % quad÷linear on polyhedral    (convex envelope f**)
  QuaPar.m            % quadratic on parabolic        (conjugate f*)
  QuaPoly.m           % quadratic on polyhedral       (= released PLQVC; input type)
  PLQVC.m             % thin alias of QuaPoly (backward compatibility)
  conjCPLQ.m          % Engine 'cplq': symbolic per-piece conjugate (default; [JOGO] 3 steps)
  conjPQP.m           % Engine 'pqp' : parametric-QP / KKT conjugate (+ partial variant); convex only [JAKEE-13]
  conjGraph.m         % Engine 'graph': point-cloud + neighbour-graph conjugate (PLQVG); convex only [HAQUE-17]/[HAQUE-18]
  convEnvDirect.m     % direct convex-envelope path (Kumar cvxEnv2d [KUMAR-20] + Karmarkar assembly [COAP])
  addQuaPoly.m        % IMPLEMENTED: QuaPoly.add's domain-overlay sum (polygon clipping adapted
                      % from maxQuaPar.m); RatPol add not yet extended -- see II.4
  addQuaPar.m         % IMPLEMENTED: add extended to QuaPar (curved/Ec-edge clipping via an
                      % axis-rotation parabola parametrization) -- prerequisite for infConv and
                      % proxAverage; see Implementation status for its STATUS/scope notes
  % addQuadratic/addScaledEnergy: IMPLEMENTED as instance methods directly on QuaPoly/QuaPar
  % (no separate file -- trivial per-face coefficient bump, no domain overlay/clipping needed,
  % unlike add's addQuaPoly.m/addQuaPar.m); prerequisite for moreau and proxAverage
  toQuaPar.m          % IMPLEMENTED: promote a QuaPoly conjugate to QuaPar (lossless relabeling)
                      % so infConv/proxAverage can add two conjugates of possibly-different type
  infConv.m           % IMPLEMENTED: conj(add(conj f,conj g)); convex f,g only (II.6); see
                      % Implementation status for scope (full-domain quadratics round-trip end to
                      % end; a bounded-triangle pair hits conjCPLQ's Step 3 TODO on the final conj)
  moreau.m            % IMPLEMENTED: single conjugate via expand-the-square [HIRIART-URRUTY-07];
                      % no convexity assumption on f -- see moreauTest.m's nonconvex case
  lasryLions.m        % IMPLEMENTED: composition of moreau (II.6); no convexity assumption on f
  proxAverage.m       % IMPLEMENTED: reduces to conj directly, NOT to moreau (II.6 derivation);
                      % convex f,g only, same Step-3 scope caveat as infConv (see above)
  +internal/
    Entity.m  Entitytype.m            % graph engine
    kktConjFace.m                     % QQ/QL/LL closed-form face conjugates (lft2.sci)
    quadQuadEnv.m                     % cPLQ symbolic bilinear x*y / quad-quad envelope
    overlayAdd.m                      % domain-overlay addition (PLCHC.add + unionHull) -- for
                                      % QuaPoly this role is filled by the real addQuaPoly.m above
    parabolicGeom.m                   % conic-arc edge intersection / parabolic-region ops
    lcon2vert.m vert2lcon.m unionHull.m intersectionHull.m
  PLQTest.m  PLQConjugateTest.m  PLQOperatorTest.m  PLQTypeCycleTest.m
```

### II.8 Validation strategy

- **Engine agreement:** `conjCPLQ` vs `conjPQP` vs `conjGraph` match on every convex example
  (`conjPQP`/`conjGraph` have no nonconvex examples to agree on — they reject those inputs, II.5).
- **Type cycle:** `conj`/`convEnv` outputs land in the correct family (`QuaPar` /
  `RatPol`); `biconj(convex f)=f` (any engine); `biconj` of nonconvex `f` `=convEnv(f)`, computed
  with `engine='cplq'` (II.6).
- **Analytic pairs:** `energy*=energy`; `oneNorm*=ι_{‖·‖∞≤1}`; `l1Norm*=ι_{‖·‖∞≤1}`,
  `linfNorm*=ι_{‖·‖1≤1}`; `moreau(energy,μ)` closed form (nonconvex test functions welcome here
  too, since `moreau` needs no convexity, II.6); `proxAverage(l1,linf;λ)` from `l1-linf-norm-pa`.
- **Paper replication:** reproduce [JOGO] Fig. 2 / Table 1 and [COAP] Fig. 2 / Table 1 exactly
  (rational arithmetic) — the strongest correctness test, since those are the published outputs.
- **Direct vs biconj envelope:** `convEnv(f,'direct')` (Kumar/Karmarkar) ≡ `convEnv(f,'biconj')`
  with `engine='cplq'` (the only engine valid for nonconvex `f` on the `'biconj'` side, II.6).

---

## PART III — Mapping summary (what to port from where)

| Target capability | Primary source(s) | Notes |
|-------------------|-------------------|-------|
| Data structure `V/E/f/F/P/dom` | released `PLQVC` | becomes leaf class `QuaPoly`; superset adds `Ec`,`den` |
| Class hierarchy (rational/parabolic types) | [COAP]/[JOGO] defs; `cPLQ.region` (parabolic ineqs) | `RatPar`→`RatPol`/`QuaPar`→`QuaPoly` |
| Convexity tests | released `PLQVC` | unchanged |
| **conjugate Engine `'cplq'`** (default) | **`cPLQ`** (`conjugateExpr`,`convexEnvelope`,`maximumP`), Tanmaya Karmarkar | the symbolic [JOGO] 3-step pipeline; bilinear `x*y` via `quadQuad`; convex **and** nonconvex |
| **conjugate Engine `'pqp'`** | Scilab `parametric_ME.sci` (Jakee Khan, [JAKEE-13]) + `lft2.sci` closed forms (Bryan Gardiner, [GARDINER-13]) | KKT parametric-QP + QQ/QL/LL closed forms; **convex only** |
| `partialConj` | Scilab `plq2_lft_partial`, COAP `partial_conjugate.sci` (Gardiner, Khan & Lucet, [GARDINER-14]) | partial branch of `'pqp'` |
| **conjugate Engine `'graph'`** | `PLQVG` + `Entity`/`Entitytype`, implementing Tasnuva Haque's algorithm ([HAQUE-17]/[HAQUE-18]) | primal-dual swap (point cloud + neighbour graph); **convex only** |
| `convEnv` `'direct'` | `deepak/cvxEnv2d` (Deepak Kumar, [KUMAR-20]) + [COAP] Step 4 assembly (Karmarkar) | faster envelope + numeric cross-check; nonconvex-capable |
| `convEnv` `'biconj'` | `conj∘conj`, `engine='cplq'` for nonconvex `f` | default; how `cPLQ` does it natively; `'pqp'`/`'graph'` only valid when `f` already convex |
| `add`/`sub`/`scalarMul` | `PLCHC.add`, Scilab `plq2_add`/`plq2_scalar` | overlay + coeff/rational arithmetic |
| `infConv` | Scilab `plq2_pa`-adjacent; `PLQList` | `conj(add(conj f,conj g))`; valid only for convex `f,g` (II.6); `add` on `QuaPar` (its prerequisite) is now done |
| `moreau` | Scilab `plq2_me` per the Hiriart-Urruty–Lucet identity [HIRIART-URRUTY-07] | **single** conjugate (expand-the-square), not `infConv`; valid for nonconvex `f` too (II.6) |
| `lasryLions` | Scilab `plq2_pa`-adjacent | `−moreau(−moreau(f,λ),μ)`; pure composition, no new geometry (II.6) |
| `proxAverage` | Scilab `plq2_pa`; `PLQList` (concept only — this design's formula is a fresh derivation, II.6) | **reduces directly to `conj`**: two conjugations sandwiching a weighted `add`, **not** three `moreau` calls; valid only for convex `f,g` (II.6) |
| polytope/parabolic utilities | `lcon2vert`/`vert2lcon`/`unionHull`/`intersectionHull`; `cPLQ.region` | reuse; extend intersection to conic arcs |
| Fitzpatrick (Gardiner, [GARDINER-09]); 1D grid-ME; SOS prover | `JCA_fitz`; `*_FMEI`; `polyPak` | Fitzpatrick + nonconvex fns kept as examples; rest **out of scope** |

---

## PART IV — Resolved decisions & remaining notes

**Resolved:**
1. **Class names & hierarchy** → `RatPar` (rational cubic÷linear on parabolic, storage umbrella)
   ⊃ `RatPol` (quad÷linear / polyhedral, `f**`), `QuaPar` (quadratic / parabolic, `f*`) ⊃
   `QuaPoly` (quadratic / polyhedral, released `PLQVC` leaf, input).
2. **Curved boundaries** → quadratic level-set arcs with **b²−4ac=0** (parabolas; lines
   degenerate). Ellipse/hyperbola arcs do **not** occur ([COAP]/[JOGO]).
3. **Storage vs operators** → the user may store **any `RatPar`**; operators apply only to
   `QuaPoly`/`QuaPar`/`RatPol` with quadratic numerator (validated by `assertOperable`). A bare
   `RatPar` / rational-on-parabolic is rejected as an operator input.
4. **Cubic** → quadratic numerator for all operators; **cubic polynomials are allowed for
   `isConvex` only** (released convexity machinery). Cubic otherwise storable, not operated on.
5. **Default arithmetic** → **symbolic + rational** everywhere in the operators (exact; avoids
   the degenerate-subdivision floating-point issues [COAP]/[JOGO] warn about).
6. **Conjugate engines** → three: `'cplq'` (default, symbolic, handles convex **and** nonconvex
   `f`, [JOGO]/[COAP]/`cPLQ`, Tanmaya Karmarkar); `'pqp'` (parametric-QP/KKT, **convex `f`
   only**, Jakee Khan's M.Sc. thesis [JAKEE-13], with closed forms reused from Bryan Gardiner's
   `lft2.sci` [GARDINER-13]); `'graph'` (point-cloud + neighbour-graph, **convex `f` only**,
   Tasnuva Haque [HAQUE-17]/[HAQUE-18]).
7. **Convex envelope** → `biconj` (default `engine='cplq'`, the only engine valid for a
   genuinely nonconvex `f` — `'pqp'`/`'graph'` are convex-only, II.5/II.6), or Deepak Kumar's
   per-piece method [KUMAR-20] (`codeOld/deepak`) as extended by Karmarkar's [COAP] Step 4
   assembly, `'direct'`, as a faster alternative/cross-check. **Not** parametric-QP (a previous
   version of this document mislabeled Kumar's method as such; see I.4).
8. **Backward compatibility** → keep **`PLQVC` as a thin alias of `QuaPoly`** (so existing code
   and tests calling `PLQVC(...)` keep working).
9. **`partialConj`** → exposed for the **`'pqp'` and `'cplq'`** engines only (the `'graph'`
   partial conjugate is left unimplemented for now).
10. **`infConv` vs `moreau`** → `infConv(f,g)=conj(conj(f)+conj(g))` computes `(f□g)**`
    (Fenchel–Moreau biconjugate), which equals the true inf-convolution `f□g` only when `f,g`
    are both convex; it is **not** used to implement `moreau`. `moreau` is instead computed as a
    **single** conjugate via the "expand the square" identity of Hiriart-Urruty & Lucet
    [HIRIART-URRUTY-07], which needs no convexity assumption on `f` at all (II.6).
11. **`proxAverage` (revised 2026-07-13)** → **not** implemented as `−moreau(−λ₁M_μf−λ₂M_μg)`
    (a previous draft of this document). Instead it reduces directly to `conj`, mirroring
    `moreau`'s "expand the square" style but with a biconjugation: `T_μP = conj(λ(T_μf)^*+
    (1-λ)(T_μg)^*)` where `T_μh:=μh+½‖·‖²`, so `P` needs three `conj` calls (two in, one out)
    around a weighted `add`, and **no call to `moreau` itself** (II.6 derivation above). Like
    `infConv`, valid only for `f,g` both convex. Build-order consequence: `proxAverage` and
    `infConv` both require `add` to work on `QuaPar` (not just `QuaPoly`) — that extension is now
    implemented (`addQuaPar.m`, see Implementation status) — **not** the `RatPol` extension of
    `add`, which nothing in the `conj`/`infConv`/`moreau`/`lasryLions`/`proxAverage` pipeline
    currently calls for (see Implementation status's "Next planned").

*(No open questions remain — the design is fully specified.)*

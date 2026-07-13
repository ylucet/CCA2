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
  `indefiniteTriangleThreeConvexEdgesUsesStep3`. The fully general case — a multi-face original
  domain (`nf>1`), or a single non-triangular face — remains open: `convEnvCPLQ`'s own
  multi-face triangulation can produce a triangle piece with exactly ONE convex edge (a genuinely
  rational envelope), which `conjPieceCPLQ` cannot conjugate yet (its own header TODO); the
  3-convex-edge single-triangle case above never hits that gap, which is why it alone is solved.
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
  non-triangular face — the single-bounded-triangle/3-convex-edge case IS now implemented (see
  `maxQuaPar.m`'s bullet above), but the fully general case needs `conjPieceCPLQ` to also handle a
  1-convex-edge (genuinely rational) envelope piece, which `convEnvCPLQ`'s multi-face
  triangulation can produce and which `conjPieceCPLQ` cannot conjugate yet — this is what still
  confines `infConv`/`moreau`/`lasryLions`/`proxAverage` to full-domain-quadratic `f,g` for an
  *exact* end-to-end round trip on a genuinely multi-face `f,g` (a bounded-triangle pair with 3
  convex edges now works end to end; a multi-face pair still errors clearly at the final `conj`
  call rather than silently giving a wrong answer).

**Next planned**: the **nonconvex-PLQ operator pipeline** (`conj`→`infConv`/`moreau`→
`lasryLions`/`proxAverage`) that has been this project's stated focus is now code-complete (see
above) — remaining work is either widening its *scope* or growing the toolbox in directions the
pipeline itself doesn't need:
1. **`conjPieceCPLQ`'s own rational-piece TODO** (conjugate of a 1-convex-edge, genuinely
   rational envelope piece) is the highest-value next step: it is the single remaining gap in
   `conjCPLQ`'s Step 3, needed for `infConv`/`moreau`/`lasryLions`/`proxAverage` to work end to end
   on a genuinely multi-face (not just single-triangle) `f,g`. The single-triangle/3-convex-edge
   case is already done (`maxQuaPar.m`'s bullet above), since it never needs a rational piece.
2. **`partialConj`** for the `'cplq'`/`'pqp'` engines (II.4) — not started.
3. **`add` for `RatPol`** (common-denominator sum) and the **`RatPar`** parent class (II.3) —
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

*Source:* `cPLQ` (`conjugateExpr`, `plq_1p`, `functionNDomain`, `region`, `quadQuad`/`qq_conj`).

**Method** — exactly the three steps of [JOGO]/[COAP]:
1. **Convex envelope of each quadratic piece** `conv(qᵢ+I_{Pᵢ})` → rational(quad÷linear) on a
   polyhedral subdivision (`plq_1p.convexEnvelope`; bilinear `x*y` / quad-on-quad handled by the
   closed-form `quadQuad`/`qq_conj`, the symbolic complement to deepak's numeric subdivision).
2. **Conjugate of each rational piece** via Lagrange multipliers `sup{s·x − r(x) : g(x)=0}`
   (`conjugateExpr`) → quadratic on a parabolic subdivision.
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

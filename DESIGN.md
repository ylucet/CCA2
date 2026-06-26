# CCA2 — Code Design Proposal for a 2D PLQ Computational Convex Analysis Toolbox

**Date:** 2026-06-26  ·  **Scope:** pure MATLAB (toolboxes allowed: Symbolic Math,
Optimization). 2D primary, 1D as a thin degenerate case. General PLQ (occasionally PLC),
**no grid/LLT algorithms**.

**Authoritative references** (in `/home/ylucet/CCA2`):
- **[COAP]** Karmarkar & Lucet, *Computing the convex envelope of bivariate PLQ functions
  in linear time*, Comput. Optim. Appl. 94 (2026) 747–780. (`s10589-026-00781-5.pdf`)
- **[JOGO]** Karmarkar & Lucet, *A linear-time algorithm to compute the conjugate of
  nonconvex bivariate PLQ functions*, J. Glob. Optim. 94 (2026) 3–34. (`s10898-025-01503-7.pdf`)

These two papers define the data types and the algorithms this design implements.

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
   - **`'pqp'`** — exact parametric-QP conjugate (Scilab `pQP` port).
   - **`'graph'`** — point-cloud + neighbour-graph conjugate (`PLQVG` primal↔dual swap).
   - **`'cplq'`** — symbolic per-piece conjugate (`cPLQ` Lagrange-multiplier method) — this is
     the engine that follows the [JOGO] 3-step algorithm exactly.
4. **Convex envelope = biconjugate**, built *on top of any conjugate engine*
   (`convEnv = biconj = conj∘conj`). This is how cPLQ does it natively, and it works equally
   over the `'pqp'` and `'graph'` engines. A **direct** envelope path (`deepak`, [COAP]) is
   kept as a faster alternative and an independent cross-check.
5. **Everything else** (`add`, `infConv`, `moreau`, `proxAverage`, `lasryLions`, `partialConj`)
   is a short composition over `conj` + `add` + scalar ops.

> **The pillars:** conjugate = {`pQP`, `PLQVG`, `cPLQ`} · convex envelope = `biconj` over any
> engine (+ `deepak` direct). All other operators compose on these.

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
| **`PLQVG`** | vertex via `Entity` | **graph/pointwise** (values `y` + subgrads `s`) | **`conjugate`** (primal↔dual entity swap: Face↔Vertex, Ray↔Ray, Segment↔Segment), `toPLCVC`, `toPLCVCDual`. **= the point-cloud + neighbour-graph conjugate engine.** |
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

### I.4 `codeOld/deepak` — parametric-QP convex envelope / conjugate ([COAP]/[JOGO] reference impl)

| Entry point | Description |
|-------------|-------------|
| `cvxEnv2d` | **Convex envelope of a quadratic over a polytope** ([COAP] Step 1+4): indefinite→bilinear `xy`, classify convex edges vs saddle vertices, build `η(a,b)` supporting families, subdivide dual space, parametric-QP per subregion, map back. Output = rational(quad÷linear)/polyhedral. |
| `cvxHull2d_solver`, `OP_Solver`, `compute_OP*`, `solve_OP_*` | Parametric 1D/QP solver core (QQ/QL/LL). |
| `Conjugate/compute_conjugate*` | Conjugate of quadratic/rational pieces; domain division into parabolic subdivision. |
| `transform_quadxp_to_xy` | Bilinear handling via Hessian eigendecomposition. |

Used as the **direct convex-envelope engine** and a numeric cross-check.

### I.5 `codeOld/MP_plq2-lft` + `COAP_plq2-lftpartial` — Scilab operator algebra (the blueprint)

The Scilab toolbox where **one conjugate generates the whole operator algebra**:
`plq2_lft` (full conjugate), `plq2_lft_partial` (partial), `plq2_me`/`plq2_me2`/`plq2_me_partial`
(Moreau, two formulas), `plq2_pa`/`plq2_pa_partial` (proximal average), `plq2_add`/`plq2_subtract`/
`plq2_scalar`, `plq2_s` (self-smoothing). PLQ2 = arrangement (DCEL) + `6×n` coeff matrix.
The `plq2_lft_pQP/parametric_ME.sci` + `…/partial_conjugate.sci` are the **parametric-QP**
engine (KKT solve, symbolic in the dual parameter) → Engine `'pqp'`. The example scripts
(`l1-norm-lft+me`, `l1-linf-norm-pa`, `house-lft`, `diamond-lft`, `four-faces-lft+me`) define
the intended public API.

### I.6 `codeOld/JCA_fitz`, `JCA_plt-scripts`, `NA_FMEI-scipts` — Scilab

`JCA_fitz/gph.sci`: **Fitzpatrick function** `F_{m,A}` as a PLQ; `plq_pds` (primal-dual symmetric
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
| `conj` | `g = conj(f, engine)` | Fenchel conjugate `f*` | engine ∈ {`'pqp'`,`'graph'`,`'cplq'`} (default `'cplq'`) |
| `partialConj` | `g = partialConj(f, idx, engine)` | conjugate w.r.t. variable `idx` | engine ∈ {`'pqp'`,`'cplq'`} only |
| `biconj` | `g = biconj(f, engine)` | `f** = conj(conj(f))` | `conj∘conj` |
| `convEnv` | `h = convEnv(f, method)` | convex envelope; `method`=`'biconj'`(any engine) \| `'direct'`(deepak) | `biconj` / `convEnvDirect` |
| `add` | `h = add(f, g)` | pointwise `f+g` | domain overlay + coeff/rational add |
| `sub` | `h = sub(f, g)` | `f − g` | `add(f, negate(g))` |
| `scalarMul`,`negate` | `c·f`, `−f` | coeff scaling | — |
| `addQuadratic` | `addQuadratic(f, A,b,c)` | `f + (½xᵀAx+bᵀx+c)` | per-face coeff update |
| `infConv` | `infConv(f,g,engine)` | `(f□g)(x)=inf_z f(z)+g(x−z)` | `conj((conj f)+(conj g))` |
| `moreau` | `moreau(f,mu,engine)` | Moreau envelope `e_μ f` | `infConv(f, ½μ‖·‖²)` (II.6) |
| `proxAverage` | `proxAverage(f,g,lambda,mu,engine)` | proximal average | compose `moreau`,`add`,`negate` |
| `lasryLions` | `lasryLions(f,lambda,mu,engine)` | double envelope | `−moreau(−moreau(f,λ),μ)` |
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

*Source:* Scilab `parametric_ME.sci`, COAP `partial_conjugate.sci`; closed forms from `lft2.sci`.

For each primal entity collect `(Q,q,κ)` and active constraints `Ax≤β`, then solve, symbolically
in the dual parameter `s`, `v(s)=sup_x sᵀx − r(x) s.t. Ax≤β` via KKT
`[Q Aᵀ; A 0][x;λ]=[s−q; β]`. Optimizer `x*(s)` is affine in `s`; back-substitution gives the
piece value; the active-set critical region is the dual cell. Assemble cells into the dual mesh.
`partialConj` = the same solve with `s` on selected coordinates only
(`_plq2_quad_face_partconj`/`_plq2_lin_face_partconj`). Reuse `PLQVC.matrixForm`/`lineEquation`/
`orderEdges` for `(Q,q)` and adjacency; build the result through the constructor.

#### II.5.3 Engine `'graph'` — point cloud + neighbour graph

*Source:* `PLQVG` + `Entity`/`Entitytype`.

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
        case 'biconj', h = biconj(f);           % conj∘conj over ANY engine (default)
        case 'direct', h = convEnvDirect(f);    % deepak cvxEnv2d pipeline ([COAP], faster)
    end
end
% Output is RatPol (quad÷linear on polyhedral), per [COAP] Step 4 / [JOGO] Prop 1.

function h = infConv(f,g,engine)            % (f □ g)
    h = conj( add( conj(f,engine), conj(g,engine) ), engine );
end

function h = moreau(f,mu,engine)            % Moreau envelope e_μ f   (cf. plq2_me)
    %   e_μ f = (1/(2μ))‖·‖² − (1/μ)·conj( μ·f + ½‖·‖² )
    g = addQuadratic( scalarMul(f,mu), eye(2), [0;0], 0 );
    h = addScaledEnergy( scalarMul( conj(g,engine), -1/mu ), 1/(2*mu) );
end

function h = proxAverage(f,g,lambda,mu,engine)   % −M_μ(−λ₁M_μf − λ₂M_μg), λ₁+λ₂=1
    Mf=moreau(f,mu,engine); Mg=moreau(g,mu,engine);
    h = negate( moreau( add(scalarMul(Mf,-lambda),scalarMul(Mg,-(1-lambda))), mu, engine ) );
end

function h = lasryLions(f,lambda,mu,engine)      % −e_μ( −e_λ f )
    h = negate( moreau( negate( moreau(f,lambda,engine) ), mu, engine ) );
end
```

`add` is the only other primitive needing real geometry: **overlay the two domain subdivisions
and add functions per overlapping piece** (sum quadratics, or add rationals over a common
denominator). Reuse `PLCHC.add` (polyshape intersection + coeff add) and the
`unionHull`/`intersectionHull` utilities; for curved boundaries, intersect parabolic regions.

### II.7 File / module layout (pure MATLAB)

```
CCA2/
  RatPar.m            % abstract parent: rational(cubic÷linear) on parabolic; storage + eval/plot/convert/isConvex
  RatPol.m            % quad÷linear on polyhedral    (convex envelope f**)
  QuaPar.m            % quadratic on parabolic        (conjugate f*)
  QuaPoly.m           % quadratic on polyhedral       (= released PLQVC; input type)
  PLQVC.m             % thin alias of QuaPoly (backward compatibility)
  conjCPLQ.m          % Engine 'cplq': symbolic per-piece conjugate (default; [JOGO] 3 steps)
  conjPQP.m           % Engine 'pqp' : parametric-QP / KKT conjugate (+ partial variant)
  conjGraph.m         % Engine 'graph': point-cloud + neighbour-graph conjugate (PLQVG)
  convEnvDirect.m     % direct convex-envelope path (deepak cvxEnv2d)
  +internal/
    Entity.m  Entitytype.m            % graph engine
    kktConjFace.m                     % QQ/QL/LL closed-form face conjugates (lft2.sci)
    quadQuadEnv.m                     % cPLQ symbolic bilinear x*y / quad-quad envelope
    overlayAdd.m                      % domain-overlay addition (PLCHC.add + unionHull)
    parabolicGeom.m                   % conic-arc edge intersection / parabolic-region ops
    lcon2vert.m vert2lcon.m unionHull.m intersectionHull.m
  PLQTest.m  PLQConjugateTest.m  PLQOperatorTest.m  PLQTypeCycleTest.m
```

### II.8 Validation strategy

- **Engine agreement:** `conjCPLQ` vs `conjPQP` vs `conjGraph` match on every convex example.
- **Type cycle:** `conj`/`convEnv` outputs land in the correct family (`QuaPar` /
  `RatPol`); `biconj(convex f)=f`; `biconj` of nonconvex `=convEnv`.
- **Analytic pairs:** `energy*=energy`; `oneNorm*=ι_{‖·‖∞≤1}`; `l1Norm*=ι_{‖·‖∞≤1}`,
  `linfNorm*=ι_{‖·‖1≤1}`; `moreau(energy,μ)` closed form; `proxAverage(l1,linf;λ)` from
  `l1-linf-norm-pa`.
- **Paper replication:** reproduce [JOGO] Fig. 2 / Table 1 and [COAP] Fig. 2 / Table 1 exactly
  (rational arithmetic) — the strongest correctness test, since those are the published outputs.
- **Direct vs biconj envelope:** `convEnv(f,'direct')` (deepak) ≡ `convEnv(f,'biconj')`.

---

## PART III — Mapping summary (what to port from where)

| Target capability | Primary source(s) | Notes |
|-------------------|-------------------|-------|
| Data structure `V/E/f/F/P/dom` | released `PLQVC` | becomes leaf class `QuaPoly`; superset adds `Ec`,`den` |
| Class hierarchy (rational/parabolic types) | [COAP]/[JOGO] defs; `cPLQ.region` (parabolic ineqs) | `RatPar`→`RatPol`/`QuaPar`→`QuaPoly` |
| Convexity tests | released `PLQVC` | unchanged |
| **conjugate Engine `'cplq'`** (default) | **`cPLQ`** (`conjugateExpr`,`convexEnvelope`,`maximumP`) | the symbolic [JOGO] 3-step pipeline; bilinear `x*y` via `quadQuad` |
| **conjugate Engine `'pqp'`** | Scilab `parametric_ME.sci` + `lft2.sci` | KKT parametric-QP + QQ/QL/LL closed forms |
| `partialConj` | Scilab `plq2_lft_partial`, COAP `partial_conjugate.sci` | partial branch of `'pqp'` |
| **conjugate Engine `'graph'`** | `PLQVG` + `Entity`/`Entitytype` | primal-dual swap (point cloud + neighbour graph) |
| `convEnv` `'direct'` | `deepak/cvxEnv2d` | faster envelope + numeric cross-check |
| `convEnv` `'biconj'` | `conj∘conj` (any engine) | default; how `cPLQ` does it natively |
| `add`/`sub`/`scalarMul` | `PLCHC.add`, Scilab `plq2_add`/`plq2_scalar` | overlay + coeff/rational arithmetic |
| `moreau`/`infConv`/`proxAverage`/`lasryLions` | Scilab `plq2_me`/`plq2_pa`; `PLQList` | compose on `conj`+`add` |
| polytope/parabolic utilities | `lcon2vert`/`vert2lcon`/`unionHull`/`intersectionHull`; `cPLQ.region` | reuse; extend intersection to conic arcs |
| Fitzpatrick; 1D grid-ME; SOS prover | `JCA_fitz`; `*_FMEI`; `polyPak` | Fitzpatrick + nonconvex fns kept as examples; rest **out of scope** |

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
6. **Conjugate engines** → three: `'cplq'` (default, symbolic), `'pqp'`, `'graph'`.
7. **Convex envelope** → `biconj` over any engine (default), `deepak` `'direct'` as alternative.
8. **Backward compatibility** → keep **`PLQVC` as a thin alias of `QuaPoly`** (so existing code
   and tests calling `PLQVC(...)` keep working).
9. **`partialConj`** → exposed for the **`'pqp'` and `'cplq'`** engines only (the `'graph'`
   partial conjugate is left unimplemented for now).

*(No open questions remain — the design is fully specified.)*

# CCA2

**Computational Convex Analysis for bivariate piecewise linear-quadratic (PLQ) functions.**

CCA2 computes Fenchel conjugates, convex envelopes and the operators derived from them for
bivariate PLQ functions — **including nonconvex ones** — in exact/closed form rather than on a
grid. No LLT or other grid-based algorithms are used.

Copyright (C) 2008–2026, Yves Lucet. Licensed under **GPL-3** (see `LICENSE`); contributors listed
in `copyright.txt`.

---

## Scope

**CCA2's goal is the conjugate and biconjugate of a `QuaPol`.** That is the design target, and it
is what determines which cases are implemented.

Everything downstream of a `QuaPol` is a *special case*, and the special-ness is load-bearing:
the convex envelope of a `QuaPol` triangle is a very special `RatPol`, and the conjugate of those
triangles is a very special `QuaPar` (for instance, a parabolic edge only ever occurs surrounded by
two parallel rays). **Hyperbolic edges therefore never arise** and are not representable — not as a
limitation, but because no `QuaPol` conjugate/biconjugate or intermediate computation can produce
one.

This matters when reading the source: several `notImplemented` guards protect *unreachable*
branches. They are assertions, not missing features. `SUPPORT_MATRIX.md` classifies every guard as
**OK** / **GAP** / **not reachable** / **internal invariant**, generated from the code itself.

---

## Status

**Alpha — not yet released.** The API is still changing; there is no tagged version.

Test suite: **238 passed / 1 failed** across 20 suites (the single failure, `testRegion/testCreation`,
is a longstanding toolbox-compatibility issue unrelated to the conjugate pipeline).

Before depending on CCA2, read `SUPPORT_MATRIX.md` §8 — the summary of what actually blocks a
general release. The largest gaps today: `partialConj` is unimplemented for every engine; the
`'pqp'` and `'graph'` conjugate engines do not exist; unbounded multi-face domains error.

---

## Types

A piecewise function is built from **two independent choices**: the *function* type and the
*subdivision* type. Each is an axis with a general case and one specialization:

| Axis | General | Specialization | Pinned by the specialization |
|---|---|---|---|
| Function | `Rat` — rational | `Qua` — polynomial | denominator ≡ 1 |
| Subdivision | `Par` — parabolic | `Pol` — polyhedral | edge conics ≡ 0 |

The four named types are the four cells of that grid, and inheritance is the grid's own partial
order — a genuine diamond, because that is what a product lattice is:

```
             RatPar        rational on parabolic     (abstract; nothing produces one yet)
             /    \
        RatPol    QuaPar   rational on polyhedral / quadratic on parabolic
             \    /
             QuaPol       quadratic on polyhedral  — the input type

        QuaParCPLQ         same maths as QuaPar, still in symbolic form (no mesh)
```

```matlab
g = f.conj();          % always a RatPar
g.kind()               % 'QuaPol' | 'RatPol' | 'QuaPar' | 'QuaParCPLQ'
g.isMeshed()           % false only for the symbolic QuaParCPLQ form

isa(g, 'Qua')          % ask about ONE axis: is the function a polynomial?
isa(g, 'Pol')          % ...is the subdivision polyhedral?
```

The axis markers are what let you query one axis without enumerating combinations, and they keep
working when a type is added. All data lives on `RatPar` alone — MATLAB makes a property defined in
two superclasses fatal *and* unresolvable, so `QuaPol < RatPol & QuaPar` requires it. Each type's
pinned values are enforced instead by validators that read the object's own traits, so no type can
be made to lie about itself. `kind()` is derived from the real class, never stored, so it cannot
drift. See `RETURN_TYPE.md` and `RatPar.m`.

> **Renamed 2026-07-27:** `QuaPoly` → `QuaPol` (and with it `addQuaPoly` → `addQuaPol`,
> `quaPolyToPlq` → `quaPolToPlq`), so all four type names are the uniform 3+3 combinations of the
> two axes. No shim was left for the old name — CCA2 has no tagged release, so nothing external can
> depend on it. `PLQVC`, which *was* released, keeps its alias.

**Not modelled as a type: the numerator's degree.** `f` is stored in the 10-wide cubic basis, but
nothing *dispatches* on degree — every consumer either rejects a cubic outright or ignores the
distinction. A type should track what changes the algorithm, not what changes the data's shape, so
degree stays a runtime check rather than becoming a `Cub` class (which would double the lattice for
no dispatch benefit).

---

## Quick start

```matlab
addpath('/path/to/CCA2');

% f(x,y) = x*y over the triangle (0,0),(1,0),(0,1) — nonconvex on that domain
V = [0 0; 1 0; 0 1];
E = [1 2 1; 2 3 1; 3 1 1];      % [from to isSegment]
F = [1 0; 1 0; 1 0];            % [leftFace rightFace], 0 = +infinity
f = QuaPol(V, E, [0 1 0 0 0 0], F);   % coeffs [x^2 xy y^2 x y 1]

g = f.conj();                   % Fenchel conjugate  -> RatPar (kind 'QuaPar')
g.eval([0.3 0.7])               % ans = 0.7

% full-domain strictly convex quadratic
q = QuaPol([1 0 1 0 0 0]);     % (x^2 + y^2)/2
q.conj().kind()                 % 'QuaPol'
q.biconj().kind()               % 'QuaPol' -- back to itself
```

**`biconj` does not yet work for every input.** It is `conj∘conj`, so it needs the *conjugate* to
be conjugable too, and the conjugate of a bounded-domain function is finite everywhere — an
unbounded multi-face domain, which `conjCPLQ` does not handle. Current coverage:

| Input | `conj` | `biconj` |
|---|---|---|
| Full-domain strictly convex quadratic | ✅ `QuaPol` | ✅ |
| **Single bounded triangle** | ✅ `QuaPar` | ❌ `PLQ:conjCPLQ:notImplemented` |
| General bounded multi-face domain | ✅ `QuaParCPLQ` | ✅ (symbolic) |

Note the shape of this: the **symbolic** path supports the biconjugate, while the faster **numeric**
single-triangle path does not. Closing it needs conjugation of an unbounded multi-face `QuaPar`.

---

## Operators

| Operator | Meaning | Notes |
|---|---|---|
| `conj(f, engine)` | Fenchel conjugate `f*` | `engine='cplq'` (default) is the **only** one implemented |
| `biconj(f)` | `f** = conj(conj(f))` | closed convex envelope |
| `convEnv(f)` | convex envelope | via `biconj`, or `convEnvCPLQ` directly |
| `add`, `sub` | pointwise `f±g` | `QuaPol` and `QuaPar`; **not** `RatPol` |
| `scalarMul`, `negate` | `c·f`, `−f` | all types |
| `addQuadratic(f,A,b,c)` | `f + (½xᵀAx+bᵀx+c)` | all types |
| `infConv(f,g)` | inf-convolution | **convex `f,g` only** |
| `moreau(f,mu)` | Moreau envelope | a **single** conjugate (expand-the-square); valid for nonconvex `f` |
| `lasryLions(f,lam,mu)` | double envelope | as `moreau` |
| `proxAverage(f,g,lam,mu)` | proximal average | **convex `f,g` only** |
| `eval`, `isConvex`, `plot` | evaluation / tests / display | `isConvex` is untested |
| `partialConj` | conjugate in one variable | **not implemented** (any engine, any type) |

---

## Algorithms and references

The `'cplq'` engine implements the linear-time algorithms of:

- **[COAP]** Karmarkar & Lucet, *Computing the convex envelope of bivariate PLQ functions in
  linear time*, Comput. Optim. Appl. **94** (2026) 747–780.
- **[JOGO]** Karmarkar & Lucet, *A linear-time algorithm to compute the conjugate of nonconvex
  bivariate PLQ functions*, J. Glob. Optim. **94** (2026) 3–34.

Both handle **nonconvex** PLQ. `DESIGN.md` lists the full lineage of predecessor work (convex-only
methods: parametric-QP, graph-matrix, and Kumar's convex-envelope front-end) and maps each to the
engine it would supply.

**On the original `cPLQ` package.** CCA2 contains a substantially rewritten version of the symbolic
per-piece machinery (`plq.m`, `region.m`, `functionNDomain.m`, `symbolicFunction.m`, and
companions). The original `cPLQ` package is **deliberately not bundled**: it is GPL-3 too, so there
is no licensing obstacle — the reason is that CCA2's rewrite supersedes it, and shipping both would
leave users wondering which to use.

---

## Repository map

| File | Purpose |
|---|---|
| `Rat.m`, `Qua.m` | function-type axis markers (rational ⊃ polynomial) |
| `Par.m`, `Pol.m` | subdivision-type axis markers (parabolic ⊃ polyhedral) |
| `RatPar.m` | abstract parent of all function types; holds all data; `kind()`, `isMeshed()` |
| `QuaPol.m` | quadratic on polyhedral — the input type |
| `QuaPar.m` | quadratic on parabolic — the conjugate's type |
| `RatPol.m` | rational on polyhedral — the convex envelope's type |
| `QuaParCPLQ.m` | conjugate still in symbolic form |
| `conjCPLQ.m` | conjugate entry point (Cases A/B/C) |
| `convEnvCPLQ.m` | Step 1 — convex envelope per piece |
| `conjPieceCPLQ.m` | Step 2 — conjugate of one envelope piece |
| `maxQuaPar.m` | Step 3 — pointwise max of two full-domain conjugates |
| `clipArcByHalfPlane.m`, `parabolaArcFrame.m` | parabola-arc geometry primitives |
| `addQuaPol.m`, `addQuaPar.m` | pointwise addition |
| `infConv.m`, `moreau.m`, `lasryLions.m`, `proxAverage.m` | derived operators |
| `SUPPORT_MATRIX.md` | what works, what does not, generated from the error guards |
| `RETURN_TYPE.md` | why `RatPar` exists |
| `DESIGN.md` | full design proposal and implementation history |

---

## Tests

```matlab
runtests('RatParTest')        % the type hierarchy contract
runtests('conjCPLQTest')      % the conjugate pipeline end to end
runtests('maxQuaParTest')     % Step 3, including curved edges
```

Run everything with `runtests` over the `*Test.m` / `test*.m` files in the repository root.

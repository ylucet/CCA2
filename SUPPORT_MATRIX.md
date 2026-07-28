# CCA2 — Supported / unsupported matrix

**Generated from the error guards in the source**, not from `DESIGN.md` (which is a design
*proposal* and describes intended as well as built code). Every "not supported" row below
corresponds to a specific `error(...)` call, cited by file and identifier, so this file can be
re-derived and re-checked mechanically.

Last regenerated: 2026-07-27, against a full test census of **238 passed / 1 failed** across 20
suites (the single failure is `testRegion/testCreation`, a longstanding toolbox-compatibility
issue unrelated to the conjugate pipeline).

---

## 0. Scope — read this first

**CCA2's goal is `QuaPoly` conjugate and biconjugate.** It is *not* to cover every possible
`RatPol` or `QuaPar`.

Everything downstream of a `QuaPoly` is a **special case**, and that special-ness is load-bearing:

- the convex envelope of any triangle coming from a `QuaPoly` is a very special case of `RatPol`;
- the conjugate of those triangles is a very special case of `QuaPar` — e.g. a parabolic edge only
  ever occurs **surrounded by two parallel rays**.

Consequently **hyperbolic edges never need to be stored**: they never arise from a `QuaPoly`
conjugate/biconjugate, nor from any intermediate computation ([COAP]/[JOGO], Karmarkar & Lucet
2026 — see `DESIGN.md`'s reference list).

This is why the table below separates **"unreachable by design"** from **"gap"**. A guard that
fires only on an input the wired pipeline cannot produce is an *assertion protecting an invariant*,
not a missing feature — there is nothing there to implement. Counting such guards as gaps
overstates what is missing, and has caused repeated mis-reporting of this project's status.

Legend:

| Mark | Meaning |
|---|---|
| **OK** | Implemented and covered by tests |
| **GAP** | Genuinely not implemented; a caller can reach it with a legitimate `QuaPoly` |
| **N/R** | Not reachable from a `QuaPoly` conjugate/biconjugate — out of scope by design, not a gap |
| **INV** | Internal invariant; firing indicates a bug, not an unsupported input |

---

## 1. Conjugate — `conj(f, engine)`

### 1.1 By engine

| Engine | Status | Guard |
|---|---|---|
| `'cplq'` (default) | **OK** — the only implemented engine | — |
| `'pqp'` | **GAP** | `PLQ:conj:engine` — `QuaPoly.m:610`, `QuaPar.m:661` |
| `'graph'` | **GAP** | same guard |
| unknown name | **INV** (input validation) | `PLQ:conj:engine` — `QuaPoly.m:613`, `QuaPar.m:664` |

Two of the three engines named in the design are **not implemented at all**. This is the single
largest functional gap.

### 1.2 `'cplq'` engine, by input shape (`conjCPLQ.m`)

| Input | Status | Returns | Guard |
|---|---|---|---|
| Full-domain **strictly convex** quadratic (`nv==0, nf==1`) | **OK** | `QuaPoly` | — |
| Full-domain **non-strictly-convex** quadratic | **GAP** | — | `PLQ:conjCPLQ:notImplemented` — `conjCPLQ.m:78` |
| Single **bounded triangle** (`nf==1, nv==3, ne==3`) | **OK** | `QuaPar` | — |
| General **bounded** domain (multi-face and/or non-triangular) | **OK** (symbolic; slow) | `QuaParCPLQ` | — |
| **Unbounded** multi-face domain | **GAP** | — | `PLQ:conjCPLQ:notImplemented` — `conjCPLQ.m:103` |
| Step 3 with a **non-triangular** envelope piece | **GAP** | — | `PLQ:conjCPLQ:notImplemented` — `conjCPLQ.m:161` |
| Cubic (`PLC`) input | **N/R** — cubic is for `isConvex` only | — | `assertOperable`; `quaPolyToPlq:cubic` |

### 1.3 Per-piece conjugate (`conjPieceCPLQ.m`, Step 2)

| Piece | Status | Guard |
|---|---|---|
| Affine | **OK** | — |
| Strictly convex (PD) quadratic | **OK** | — |
| Convex rank-1 PSD quadratic (COAP A.4 envelope output) | **OK** | — |
| Indefinite quadratic, **0 or 1** convex edge (any frame/shift) | **OK** | — |
| Indefinite quadratic, **2** convex edges | **N/R** | `conjPieceCPLQ:notImplemented` — `conjPieceCPLQ.m:586` |
| Concave (negative semidefinite) | **N/R** — call `convEnvCPLQ` first | `conjPieceCPLQ.m:128` |
| Genuinely **rational** piece (nonzero `RatPol` denominator) | **N/R** | `conjPieceCPLQ.m:107` |
| Non-triangle / unbounded piece | **N/R** — Step 1 only emits bounded triangles | `conjPieceCPLQ.m:103` |
| Cubic numerator | **N/R** | `conjPieceCPLQ.m:112` |

Notes on the **N/R** rows — each is protected by an upstream invariant, not by luck:

- **2 convex edges**: Step 1 (`convEnvCPLQ`) always convexifies such a piece into a rank-1 PSD
  quadratic before Step 2 sees it. Reaching this branch requires hand-feeding a raw indefinite
  piece, bypassing Step 1. See `conjPieceCPLQTest.bilinearTwoConvexEdgesIsUnreachableFromTheWiredPipeline`.
- **Rational piece**: `conjCPLQ`'s `conjSingleTriangle` conjugates the *original* quadratic piece
  directly, never a materialized rational envelope — so a genuinely rational piece is never handed
  to Step 2 by the wired pipeline.

---

## 2. Partial conjugate — `partialConj(f, idx, engine)`

| Engine | Status | Guard |
|---|---|---|
| `'cplq'` | **GAP** | `PLQ:conjCPLQ:partialNotImplemented` — `conjCPLQ.m:53` |
| `'pqp'` | **GAP** | `PLQ:partialConj:engine` — `QuaPoly.m:628`, `QuaPar.m:679` |
| `RatPol` input | **GAP** | `RatPol:partialConj:notImplemented` — `RatPol.m:659` |

**`partialConj` is not implemented for any engine or any type.** The dispatch exists; every path
errors.

---

## 3. Biconjugate and convex envelope

`biconj` is `conj∘conj`, so it needs the **conjugate** to be conjugable too. The conjugate of a
bounded-domain function is finite everywhere — an unbounded multi-face domain — which §1.2 shows
`conjCPLQ` rejects. Measured coverage, by the input's `conjCPLQ` case:

| Input | `conj` | `biconj` | Why |
|---|---|---|---|
| **A** full-domain strictly convex quadratic | **OK** → `QuaPoly` | **OK** | conjugate is again a full-domain quadratic (Case A twice) |
| **B** single bounded triangle | **OK** → `QuaPar` | **GAP** | conjugate is a mesh `QuaPar` on an unbounded multi-face domain → `PLQ:conjCPLQ:notImplemented` (`conjCPLQ.m:103`) |
| **C** general bounded multi-face domain | **OK** → `QuaParCPLQ` | **OK** | `QuaParCPLQ.conj` routes through cPLQ's own symbolic machinery, which handles it |

Note the shape of this gap: the **symbolic** path supports the biconjugate while the faster
**numeric** single-triangle path does not. Since conjugate *and biconjugate* of a `QuaPoly` is the
project's stated goal (§0), this is a first-class gap, not a corner case. Closing it means
conjugating an unbounded multi-face `QuaPar` — either by extending `conjCPLQ`, or by routing Case B
through the cPLQ path when a second conjugation is wanted.
Regression tests: `conjCPLQTest.biconjCoverageByInputCase`.

| Operation | Status | Notes |
|---|---|---|
| `QuaPoly.biconj` | **partial** | see the table above |
| `QuaPar.biconj` | **partial** | same root cause — a full-domain `QuaPar` is not conjugable |
| `RatPol.biconj` | **GAP** | `RatPol:biconj:notImplemented` — `RatPol.m:664` |
| `RatPol.conj` | **GAP** | `RatPol:conj:notImplemented` — `RatPol.m:653` |
| `convEnvCPLQ`, nonconvex quadratic over a **bounded triangle** | **OK** | 0/1/2/3-convex-edge splits all implemented |
| `convEnvCPLQ`, nonconvex quadratic over a non-triangle/unbounded domain | **GAP** | `convEnvCPLQ:notImplemented` — `convEnvCPLQ.m:65` |
| `convEnvCPLQ`, multi-piece with an **unbounded** face | **GAP** | `convEnvCPLQ:notImplemented` — `convEnvCPLQ.m:457` |
| `convEnvCPLQ` of a nonconvex **full-domain** quadratic | **N/R** — envelope is not finite | `convEnvCPLQ.m:60` |

---

## 4. Pointwise max of conjugates — `maxQuaPar` (Step 3)

| Case | Status | Guard |
|---|---|---|
| Both inputs purely polyhedral | **OK** | — |
| Exactly **one** input with parabolic edges | **OK** (added 2026-07-27) | — |
| **Both** inputs curved | **GAP** — needs conic–conic intersection | `maxQuaPar:notImplemented` — `maxQuaPar.m:192` |
| Splitting a cell that **already carries an arc** | **GAP** — ~26% of sampled valid splits | `maxQuaPar.m:1264` |
| Clip line cutting one arc **twice** (arc bulging across) | **GAP** (defensive; 0 occurrences observed) | `maxQuaPar.m:802` |
| Face vertex inside the **open interior of an arc** (arc must be split) | **GAP** (defensive; 0 occurrences observed) | `maxQuaPar.m:617` |
| Curved **ray** edge | **N/R** — `QuaPar` has no unbounded curved edge | `maxQuaPar.m:400` |
| Input not finite everywhere | **N/R** — Step 3 combines full-domain conjugates | `maxQuaPar:notFullDomain` — `maxQuaPar.m:409` |
| Difference of the two candidate quadratics is an irreducible ellipse/hyperbola | **N/R** | `maxQuaPar:notDegenerate` — `maxQuaPar.m:1306` |
| `maxQuaPar:internal` (8 sites) | **INV** | assembly-topology invariants |

Accuracy of the curved path, measured on a randomized sweep of convex quadrilaterals split by a
diagonal (115 splits, 85 assembled): exact to ~1e-14 at **all 340** straight-edge midpoints and
**all 5100** interior sample points, with **zero** `maxQuaPar:internal` crashes and **zero**
arrangement-validity violations. See `maxQuaPar.m`'s header VALIDATION block.

---

## 5. Arithmetic and derived operators

| Operation | `QuaPoly` | `QuaPar` | `RatPol` |
|---|---|---|---|
| `add` / `sub` | **OK** | **OK** (partial, see below) | **GAP** — not implemented |
| `scalarMul`, `negate` | **OK** | **OK** | **OK** |
| `addQuadratic` | **OK** | **OK** | **OK** |
| `eval` | **OK** | **OK** (see §7) | **OK** |
| `isConvex` | untested (`%TO BE TESTED` in source) | untested | untested |
| `plotDomain` | **OK** | linear edges only — `QuaPar:plotDomain:curved` (`QuaPar.m:570`) | **OK** |

`addQuaPar` carries **11 distinct `addQuaPar:notImplemented` guards** (curved rays, >1 curved edge
per face, degenerate conics, disconnected clip results, new unbounded curved rays, …). These are
**GAP**s for arbitrary `QuaPar` input, but most are **N/R** for `QuaPar`s that arise from a
`QuaPoly` conjugate. They have not been individually classified — doing so requires the same
"can a QuaPoly conjugate reach this?" analysis applied per guard, and is listed as follow-up work.

Derived operators — all thin compositions, all **OK** within the limits of the `conj` they call:

| Operator | Built on | Validity precondition |
|---|---|---|
| `infConv(f,g)` | `conj((conj f)+(conj g))` | `f,g` **convex** |
| `moreau(f,mu)` | a **single** `conj` (expand-the-square, [HIRIART-URRUTY-07]) | any `f`, convex or not |
| `lasryLions(f,lam,mu)` | `−moreau(−moreau(f,λ),μ)` | as `moreau` |
| `proxAverage(f,g,lam,mu)` | two `conj` around a weighted `add` | `f,g` **convex** |

---

## 6. Class hierarchy

| Item | Status |
|---|---|
| `QuaPoly`, `QuaPar`, `RatPol`, `QuaParCPLQ` | **OK** — each a standalone `classdef` |
| Common `RatPar` parent (`DESIGN.md` II.3) | **GAP** — proposed, not built |

The missing parent is what forces `conj`'s return type to vary by input shape (`QuaPoly` /
`QuaPar` / `QuaParCPLQ`). See `RETURN_TYPE.md`.

---

## 7. Known defects (not "unimplemented" — these are wrong answers)

| Defect | Impact | Where |
|---|---|---|
| `mergeL` / `removeTangent` exact-tie-point gap | Wrong result at exact symmetric tie points | vendored cPLQ `functionNDomain` / `plq.biconjugateF`; inherited by `QuaParCPLQ.conj` |
| `QuaPar.eval` exactly **at a vertex** | ~1.4% of result vertices return `Inf` or a wrong value; correct to ~1e-15 at radius 1e-8 | `QuaPar.eval`'s exact, no-tolerance point location — affects the polyhedral path too, not just curved |
| `testRegion/testCreation` | 1 failing test | toolbox-compatibility, unrelated to the conjugate pipeline |

---

## 8. Summary — what actually blocks a general release

Ordered by how likely a downstream caller is to hit it:

1. **`biconj` fails for a single bounded triangle** (§3) — half the project's stated goal, broken
   on the numeric path. Root cause is the same as item 3 below.
2. **`partialConj` is entirely unimplemented** (§2).
3. **Unbounded multi-face domains error** (§1.2) — not an exotic input, and the cause of item 1.
4. **`'pqp'` and `'graph'` engines missing** (§1.1).
5. **`RatPol.conj`/`biconj`/`add` missing** (§3, §5).
6. **Two known wrong-answer defects** (§7).
7. Performance: general bounded domains route through the symbolic pipeline (Phase 2).

**Resolved 2026-07-27:** `conj`'s return type no longer varies by input shape — every type now
inherits from `RatPar` and `kind()` reports the concrete one. See `RETURN_TYPE.md` and `RatPar.m`.

Items **not** on this list, deliberately: every **N/R** row above. They are unreachable from a
`QuaPoly` conjugate/biconjugate and require no work.

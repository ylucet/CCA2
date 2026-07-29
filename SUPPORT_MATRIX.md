# CCA2 — Supported / unsupported matrix

**Generated from the error guards in the source**, not from `DESIGN.md` (which is a design
*proposal* and describes intended as well as built code). Every "not supported" row below
corresponds to a specific `error(...)` call, cited by file and identifier, so this file can be
re-derived and re-checked mechanically.

Last regenerated: 2026-07-28, against a full census of **262 passed / 1 failed** across 22 suites.
Broken out: the 17 CCA2 suites (`*Test.m`) are **187 / 0**, up from 176 — the 11 new tests are
`biconjCPLQTest` (§3) and `conjCPLQTest.indefiniteTriangleTwoConvexEdgesSplitViaCPLQStep2` (§1.2).
The 5 vendored cPLQ suites (`test*.m`) are 75 / 1; the single failure, `testRegion/testCreation`,
is a longstanding toolbox-compatibility issue unrelated to the conjugate pipeline.

---

## 0. Scope — read this first

**CCA2's goal is `QuaPol` conjugate and biconjugate.** It is *not* to cover every possible
`RatPol` or `QuaPar`.

Everything downstream of a `QuaPol` is a **special case**, and that special-ness is load-bearing:

- the convex envelope of any triangle coming from a `QuaPol` is a very special case of `RatPol`;
- the conjugate of those triangles is a very special case of `QuaPar` — e.g. a parabolic edge only
  ever occurs **surrounded by two parallel rays**.

Consequently **hyperbolic edges never need to be stored**: they never arise from a `QuaPol`
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
| **GAP** | Genuinely not implemented; a caller can reach it with a legitimate `QuaPol` |
| **N/R** | Not reachable from a `QuaPol` conjugate/biconjugate — out of scope by design, not a gap |
| **INV** | Internal invariant; firing indicates a bug, not an unsupported input |

---

## 1. Conjugate — `conj(f, engine)`

### 1.1 By engine

| Engine | Status | Guard |
|---|---|---|
| `'cplq'` (default) | **OK** — the only implemented engine | — |
| `'pqp'` | **GAP** | `PLQ:conj:engine` — `QuaPol.m:589`, `QuaPar.m:633` |
| `'graph'` | **GAP** | same guard |
| unknown name | **INV** (input validation) | `PLQ:conj:engine` — `QuaPol.m:592`, `QuaPar.m:636` |

Two of the three engines named in the design are **not implemented at all**. This is the single
largest functional gap.

### 1.2 `'cplq'` engine, by input shape (`conjCPLQ.m`)

| Input | Status | Returns | Guard |
|---|---|---|---|
| Full-domain **strictly convex** quadratic (`nv==0, nf==1`) | **OK** | `QuaPol` | — |
| Full-domain **non-strictly-convex** quadratic | **GAP** | — | `PLQ:conjCPLQ:notImplemented` — `conjCPLQ.m:78` |
| Single **bounded triangle**, Step 1 envelope has **one** face | **OK** (numeric; fast) | `QuaPar` | — |
| …envelope split with a **rational** face, 2-face split | **OK** (falls back to cPLQ's Step 2/3; slow, ~100 s) | `QuaParCPLQ` | — |
| …envelope split with a **rational** face, 4-face split | **GAP** — cPLQ's **Step 3** (cross-piece max), see below | — | `PLQ:conjCPLQ:cplqFailed` — `conjCPLQ.m`'s `conjEnvelopeViaCPLQ` |
| General **bounded** domain (multi-face and/or non-triangular) | **OK** (symbolic; slow) | `QuaParCPLQ` | — |
| **Unbounded** multi-face domain | **GAP** | — | `PLQ:conjCPLQ:notImplemented` — `conjCPLQ.m:103` |
| Step 3 with a **non-triangular** envelope piece | **GAP** | — | `PLQ:conjCPLQ:notImplemented` — `conjCPLQ.m:161` |
| Cubic (`PLC`) input | **N/R** — cubic is for `isConvex` only | — | `assertOperable`; `quaPolToPlq:cubic` |

Those three rows replace a single previously-claimed **OK**. Since Step 1 became recursive
(`convEnvCPLQ`'s `solveTriangleBF`/`splitTwoConvexEdges`, 2026-07-17/18), the envelope of an
indefinite triangle can contain genuinely rational faces, and `conjMaxOfSubTriangles` hands one to
Step 2, which has no rational branch (§1.3).

**Step 2 now falls back to cPLQ's own symbolic Step 2/3 for exactly those envelopes**
(`conjCPLQ.m`'s `conjEnvelopeViaCPLQ`, via the new `ratPolToPlq.m` adapter), which closes the
2-convex-edge case: on `f=xy` over `(2,1),(0,0),(1,0)` the result matches `sup_x ⟨s,x⟩ − f(x)` to
≤ 8.9e-16 at all 10 sampled dual points. Pinned by
`conjCPLQTest.indefiniteTriangleTwoConvexEdgesSplitViaCPLQStep2`.

Note **which half** of cPLQ is reused — its Step 2/3, on *CCA2's* Step 1 output. Running cPLQ end
to end instead does not work, and it is worth recording why, because "cPLQ already does this" is
the natural assumption:

- for a **2-convex-edge** triangle cPLQ's own envelope is the single Appendix A.4 formula, which
  CCA2 established is not always tight (`convEnvCPLQTest.bilinearTwoConvexEdgesSplitIsTight`). Its
  conjugate inherits that: on the triangle above it leaves the paper's own flagged dual point
  `s=(-0.008727,-0.999962)` covered by **no region at all** (evaluates `NaN`).
- for a **3-convex-edge** triangle cPLQ has **no branch whatsoever** — neither `convexEnvelope1`
  nor `conjugateFunction` in `plq_1p.m` has an `nCE==3` case — so it silently returns an empty
  envelope and an empty conjugate.

CCA2's Step 1 is ahead of cPLQ's on both counts (the tightness split, and [COAP] Appendix A.5's
3-convex-edge split, are CCA2's own), so the working combination is CCA2 Step 1 + cPLQ Step 2/3.

There is **no 3-convex-edge case to implement**, in CCA2 or in cPLQ. [COAP] Appendix A.5's split
reduces such a triangle to 2-convex-edge sub-triangles and Step 1 already applies it, so every
face reaching Step 2 has come through it. The edge count describes the *input*; it is not what
fails.

What still fails is the **4-face** split, and the failure is now in **Step 3** — cPLQ's
cross-piece maximum (`plq.maximumConjugate` → `functionNDomain.maximumP` → `region.maximum`).
Step 2 completes on all four faces. Getting there took resolving a chain of removable `0/0`s (see
the note below) plus two robustness guards in `region.m`.

**Update (2026-07-29) — Step 3 now runs; only its ASSEMBLY is wrong.** Four further blockers were
fixed and the 4-face pipeline completes end to end (~15 min). Steps 1 and 2 are **exact** on the
reference case `f = xy` over `conv{(0,0),(3,3),(1,2)}`: the pointwise max of the four per-piece
conjugates matches `sup_{x∈T} ⟨s,x⟩ − xy` at **all 289 points** of a 17×17 dual grid, and Step 1's
envelope reproduces `f*` to grid resolution. The four fixes were

| Site | Defect |
|---|---|
| `region.maxArray` (`region.probeAlong`/`probePerp`) | tie-break probe directions came from constraint *slopes*; a symbolically-zero-but-not-syntactically-zero slope made `d == 0` false, so `-1/d` survived as an unevaluated expression and detonated later in `subsF`'s `simplifyFraction` as `symbolic:kernel:DivisionByZero` |
| `functionNDomain.mergeL` | second accumulation loop passed the **first** region's vertices to `removeTangent`, and none at all (`Unrecognized function or variable 'nP'`) when that first region was empty |
| `region.linear3pt` | translated an edge-local vertex index into a region index **in place**, so the next iteration indexed a 3-element array out of range |
| `plq_1p.conjugateFunction`, `nCE==2` | the `grad` half-plane per edge was the hard-coded pattern `+,−,+` — a statement about one edge *ordering*. For any other ordering two edges got each other's half-plane, so `f*` came out as the **min** of the two edge candidates. Now derived from geometry: the rank-1 Hessian's kernel `k = [b,−2a]` gives `grad = ⟨s−L,k⟩`, so edge *j* owns `{grad ≥ 0}` exactly when its outward normal has `⟨k,n_j⟩ > 0`; an edge parallel to `k` is a level edge and contributes no region |

The remaining error is **entirely in the assembly**, and has two named causes, both of which let a
region claim territory that was never its own — carrying the wrong value there (~12% of that grid):

1. **`region.merge` is unsound.** It unions `A ∩ {g≤0}` and `B ∩ {−g≤0}` by deleting the shared
   facet and intersecting the rest, giving `A ∩ B`. That equals the union only when `A` and `B`
   are the *same* constraint set. Measured: three same-valued Step 3 regions merge into one
   covering `s=(1,1)`, where none of them does.
2. **`region.simplifyUnboundedRegion` drops non-redundant constraints.** Its rule is "delete any
   constraint not passing through a finite vertex", which is not a redundancy test for unbounded
   or curved regions. Measured: a region keeping `(s₁+s₂)²/4` loses `−s₁−s₂ ≤ 0` and reports
   `f*(−3,−3) = 9` where the true value is `0`.

Guarding (1) alone with a constraint-set-equality test is provably sound but makes things **worse**
(36 → 125 wrong of 289), because refusing merges leaves more regions for (2) to damage. **Fix them
together.** Until then `conjCPLQ`'s `assertStep3MatchesPieces` cross-checks the assembled maximum
against the per-piece max — the same `f*` by a different route — over an 11×11 dual grid and raises
`PLQ:conjCPLQ:cplqFailed` on disagreement, so this remains a loud failure rather than a silently
wrong answer.

**Correction (2026-07-28).** An earlier edition of this section blamed the numeric/symbolic
boundary and called for exact arithmetic in Step 1. That was **wrong**, and the measurement that
refutes it is worth keeping: at the vertex where a face's denominator vanishes, the numerator
vanishes too, with ∇N parallel to ∇den — residual 0 to 1.3e-16 — so the `0/0` is cleanly
**removable** and the envelope coefficients are accurate to ~1e-16. Nothing needed to be made
exact. What was missing was *resolving* the `0/0` by a limit instead of substituting into it.
`symbolicFunction.subsF` already flagged it (returning `NaN`) and `region.funcVertices` already
used a `NaN`-then-limit idiom, but `symbolicFunction.limit` takes **iterated univariate** limits,
which return `NaN` again for a bivariate `0/0`. The new `symbolicFunction.limitDirectional`
restricts to a line through the point — one univariate limit, which a rational function always
answers — and requires several directions to agree, which is itself the check that the limit
exists. Pinned by `conjCPLQTest.indefiniteTriangleThreeConvexEdgesUsesStep3`. `biconj` is
unaffected throughout — see §3.

### 1.3 Per-piece conjugate (`conjPieceCPLQ.m`, Step 2)

| Piece | Status | Guard |
|---|---|---|
| Affine | **OK** | — |
| Strictly convex (PD) quadratic | **OK** | — |
| Convex rank-1 PSD quadratic (COAP A.4 envelope output) | **OK** | — |
| Indefinite quadratic, **0 or 1** convex edge (any frame/shift) | **OK** | — |
| Indefinite quadratic, **2** convex edges | **N/R** | `conjPieceCPLQ:notImplemented` — `conjPieceCPLQ.m:586` |
| Concave (negative semidefinite) | **N/R** — call `convEnvCPLQ` first | `conjPieceCPLQ.m:128` |
| Genuinely **rational** piece (nonzero `RatPol` denominator) | **GAP** (was listed N/R — see below) | `conjPieceCPLQ.m:107` |
| Non-triangle / unbounded piece | **N/R** — Step 1 only emits bounded triangles | `conjPieceCPLQ.m:103` |
| Cubic numerator | **N/R** | `conjPieceCPLQ.m:112` |

Notes on these rows:

- **2 convex edges** (**N/R**): Step 1 (`convEnvCPLQ`) always convexifies such a piece into a
  rank-1 PSD quadratic before Step 2 sees it. Reaching this branch requires hand-feeding a raw
  indefinite piece, bypassing Step 1. See
  `conjPieceCPLQTest.bilinearTwoConvexEdgesIsUnreachableFromTheWiredPipeline`.
- **Rational piece** (**GAP**, reclassified 2026-07-28): the previous edition claimed
  `conjSingleTriangle` conjugates the *original* quadratic piece directly and so never
  materializes a rational envelope. That reasoning covers only the branch where Step 2 accepts
  the raw piece. When it does not (concave, or ≥2 convex edges), `conjSingleTriangle` falls back
  to Step 1 — and since Step 1 became recursive, that envelope can carry rational faces, which
  `conjMaxOfSubTriangles` then feeds straight into Step 2. A legitimate single-triangle `QuaPol`
  reaches this guard, so it is a gap in `conjPieceCPLQ`, not an unreachable assertion.
  **It is no longer a gap in the pipeline**: `conjCPLQ` catches it and routes the whole envelope
  through cPLQ's Step 2/3 instead (§1.2). Implementing a native numeric rational branch here would
  buy speed — ~100 s drops to milliseconds — not new coverage.

---

## 2. Partial conjugate — `partialConj(f, idx, engine)`

| Engine | Status | Guard |
|---|---|---|
| `'cplq'` | **GAP** | `PLQ:conjCPLQ:partialNotImplemented` — `conjCPLQ.m:53` |
| `'pqp'` | **GAP** | `PLQ:partialConj:engine` — `QuaPol.m:607`, `QuaPar.m:651` |
| `RatPol` input | **GAP** | `RatPol:partialConj:notImplemented` — `RatPol.m:631` |

**`partialConj` is not implemented for any engine or any type.** The dispatch exists; every path
errors.

---

## 3. Biconjugate and convex envelope

`biconj` dispatches by input shape in `biconjCPLQ.m`, exactly as `conjCPLQ` does. Cases A and C
are the literal `conj∘conj`; **Case B is not**, and that is what closed this section's former
blocker. For a compact domain `T` and `q` continuous on it, `f = q + I_T` gives

    f** = cl conv f = conv f = conv(q + I_T)

which is precisely what **Step 1** (`convEnvCPLQ`) computes. So the biconjugate of a bounded
triangle is available in closed form with no second conjugation, and therefore without needing the
unbounded multi-face conjugate of §1.2 — which stays open, and is still the right long-term fix
for `conj` itself.

| Input | `conj` | `biconj` | Route |
|---|---|---|---|
| **A** full-domain strictly convex quadratic | **OK** → `QuaPol` | **OK** → `QuaPol` | `conj∘conj` (conjugate is again a full-domain quadratic, so Case A twice) |
| **B** single bounded triangle | **OK** → `QuaPar` (except §1.2's rational-sub-piece row) | **OK** → `RatPol` | Step 1's convex envelope, `convEnvCPLQ` |
| **C** general bounded multi-face domain | **OK** → `QuaParCPLQ` | **OK** → `QuaParCPLQ` | `conj∘conj`; `QuaParCPLQ.conj` routes through cPLQ's own symbolic machinery |
| **unbounded** domain | **GAP** | **GAP** | falls through to `conj∘conj`, which errors at `conjCPLQ.m:103` |

Case B's route is **strictly stronger** than `conj∘conj` would have been, not merely equal: it
succeeds on the very triangles whose *conjugate* Step 2 cannot compute (§1.2's rational-sub-piece
row) — Step 1 produces the envelope regardless, so `f**` is available even where `f*` is not.

Correctness is checked against a ground truth that uses **none** of the conjugate pipeline:
`f*(s)` by exact maximization of `⟨s,x⟩ − q(x)` over the triangle, then `f**(x) = sup_s ⟨s,x⟩ −
f*(s)` numerically. Agreement ≤ 8.9e-16 at interior points of 7 triangles covering every Step 1
branch (affine, convex, concave, indefinite with 0/1, 2 and 3 convex edges — i.e. 1-, 2- and
4-face envelopes). Regression tests: `biconjCPLQTest` (10 tests) and
`conjCPLQTest.biconjCoverageByInputCase`.

| Operation | Status | Notes |
|---|---|---|
| `QuaPol.biconj` | **OK** for every bounded domain | see the table above; unbounded domains still error |
| `QuaPar.biconj` | **OK** for a polyhedral bounded triangle | a triangle with a genuinely **parabolic** side is not a `convEnvCPLQ` input, so it falls through to `conj∘conj` and still errors |
| `RatPol.biconj` | **GAP** | `RatPol:biconj:notImplemented` — `RatPol.m:636` |
| `RatPol.conj` | **GAP** | `RatPol:conj:notImplemented` — `RatPol.m:625` |
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

| Operation | `QuaPol` | `QuaPar` | `RatPol` |
|---|---|---|---|
| `add` / `sub` | **OK** | **OK** (partial, see below) | **GAP** — not implemented |
| `scalarMul`, `negate` | **OK** | **OK** | **OK** |
| `addQuadratic` | **OK** | **OK** | **OK** |
| `eval` | **OK** | **OK** (see §7) | **OK** |
| `isConvex` | untested (`%TO BE TESTED` in source) | untested | untested |
| `plotDomain` | **OK** | linear edges only — `QuaPar:plotDomain:curved` (`QuaPar.m:542`) | **OK** |

`addQuaPar` carries **11 distinct `addQuaPar:notImplemented` guards** (curved rays, >1 curved edge
per face, degenerate conics, disconnected clip results, new unbounded curved rays, …). These are
**GAP**s for arbitrary `QuaPar` input, but most are **N/R** for `QuaPar`s that arise from a
`QuaPol` conjugate. They have not been individually classified — doing so requires the same
"can a QuaPol conjugate reach this?" analysis applied per guard, and is listed as follow-up work.

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
| `QuaPol`, `QuaPar`, `RatPol`, `QuaParCPLQ` | **OK** |
| Common `RatPar` parent (`DESIGN.md` II.3) | **OK** — built 2026-07-27 |
| Axis markers `Rat`/`Qua` × `Par`/`Pol` | **OK** — property-less abstract traits, so either axis can be queried alone |
| `RatPar` itself as a *result* type (rational-on-parabolic) | **N/R** — abstract; no operator produces one |

Every type now inherits from `RatPar`, so `conj`/`biconj`/`convEnv` are statically typed: they
return a `RatPar`, and `kind()` reports the concrete one (`'QuaPol'` / `'RatPol'` / `'QuaPar'` /
`'QuaParCPLQ'`). Use `isMeshed()` to tell a `V/E/Ec/F/P` mesh from the still-symbolic
`QuaParCPLQ`. See `RETURN_TYPE.md` and `RatPar.m`; the contract is pinned entry by entry in
`RatParTest`.

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

1. **`partialConj` is entirely unimplemented** (§2).
2. **`conj` of a triangle whose Step 1 envelope splits into 4 faces** (§1.2) — Steps 1 and 2 are
   now exact on the reference case; cPLQ's Step 3 *assembly* is not. Two named causes,
   `region.merge`'s unsound facet-deletion union and `region.simplifyUnboundedRegion`'s
   not-a-redundancy-test constraint dropping, must be fixed together. Guarded by
   `assertStep3MatchesPieces`, so it fails loudly rather than returning wrong numbers. The 2-face
   split works end to end; `biconj` works for both.
3. **Unbounded multi-face domains error** (§1.2) — the remaining reason `conj` is not closed under
   itself.
4. **`'pqp'` and `'graph'` engines missing** (§1.1).
5. **`RatPol.conj`/`biconj`/`add` missing** (§3, §5).
6. **Two known wrong-answer defects** (§7).
7. **`maxQuaPar` cannot split a cell that already carries an arc** (§4) — ~26% of sampled splits.
8. Performance: general bounded domains route through the symbolic pipeline (Phase 2).

**Resolved 2026-07-28:** `biconj` works for every bounded domain, including the single bounded
triangle that used to be blocker 1 — via Step 1's convex envelope rather than a second
conjugation. See §3 and `biconjCPLQ.m`. Also `conj` of a 2-convex-edge split triangle, by falling
back to cPLQ's Step 2/3 (§1.2, `ratPolToPlq.m`). That work also fixed four latent bugs in the
vendored symbolic layer, all previously unreachable because the `0/0` errored first:
`plq_1p.m`'s `nCE==2` branch left the coefficient `d` unassigned when the envelope had no linear
term; `region.vertexOfEdge` and `region.simplifyUnboundedRegion` substituted into a removable
`0/0` instead of taking its limit; and `region.plus` indexed into an empty region array.

**Resolved 2026-07-27:** `conj`'s return type no longer varies by input shape — every type now
inherits from `RatPar` and `kind()` reports the concrete one. See `RETURN_TYPE.md` and `RatPar.m`.

Items **not** on this list, deliberately: every **N/R** row above. They are unreachable from a
`QuaPol` conjugate/biconjugate and require no work.

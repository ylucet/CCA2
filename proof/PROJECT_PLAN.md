# Project Plan

Phase 0 is the specification and does not change without a `DECISIONS.md` entry.
Phases 1-8 are the roadmap. Do not skip ahead.

> **Status, 2026-08-22: Phases 1-5 are complete.** `conj_isQuaCon` is proved,
> `sorry`-free, with `#print axioms` clean, for the Stage 1 input class. Where the
> Lean development deviates from Phase 0 below -- the ambient type, the encoding of
> a piece's quadratic, the shape of the conic classification, the emptiness of the
> `top` cell -- each deviation has a `DECISIONS.md` entry. Phase 0 is left as
> written so the deviations stay visible.

---

# Phase 0 — The target statement

## 0.1 Setting

`E := EuclideanSpace ℝ (Fin 2)`. Chosen over `ℝ × ℝ` because mathlib puts an
`InnerProductSpace` instance on it directly. Lemmas that do not mention conics
are stated over a general finite-dimensional real inner product space where that
is free, and specialized to `E` only in the conic part.

Coefficients are **real**. Rationality is a separate, already-answered question
(`../CONJ_FIELD_PROOF.md`) and is orthogonal to shape.

## 0.2 Input — `QuaPol`

```lean
structure QuaPiece (E) where
  verts : Finset E
  verts_nonempty : verts.Nonempty
  Q : E →L[ℝ] E
  Q_symm : IsSelfAdjoint Q
  beta : E
  gamma : ℝ

def QuaPiece.T (p) : Set E := convexHull ℝ ↑p.verts
def QuaPiece.q (p) (x) : ℝ := (1/2) * ⟪x, p.Q x⟫ + ⟪p.beta, x⟫ + p.gamma

structure QuaPol (E) where
  pieces : Finset (QuaPiece E)
  pieces_nonempty : pieces.Nonempty

def QuaPol.eval (f) (x) : EReal :=
  ⨅ p ∈ f.pieces, (if x ∈ p.T then (p.q x : EReal) else ⊤)
```

Four deliberate choices, each with a reason:

* **V-representation.** Pieces are convex hulls of finite vertex sets, not
  intersections of half-planes. Recovering vertices from half-planes needs
  Minkowski–Weyl, which mathlib does not have.
* **No consistency hypothesis.** The pieces are not required to have disjoint
  interiors, nor is `f` required to be continuous or single-valued piece to
  piece. `f` is defined as the pointwise `min` over pieces, and the conjugate
  depends only on that `min`, so no such hypothesis is needed. This is a genuine
  simplification over the informal proof, which discusses continuity at length
  for the *field* question only.
* **Stage 1 = bounded pieces.** `convexHull` of a `Finset` is compact for free,
  which gives attainment of the per-piece supremum from
  `IsCompact.exists_isMaxOn` and avoids Frank–Wolfe entirely. Unbounded pieces
  are Phase 7.
* **`Q` arbitrary.** Indefinite and singular Hessians are allowed. This matters:
  the singular case is the one hard branch of the selection lemma (0.6, step S8).

## 0.3 Output — the conjugate

```lean
def conj (f : QuaPol E) (s : E) : EReal := ⨆ x : E, ((⟪s, x⟫ : ℝ) : EReal) - f.eval x
```

`EReal`-valued, per the design decision to carry the `+∞` region rather than
restrict to `dom f*`. With `pieces_nonempty` and `verts_nonempty`,
`conj f s > ⊥` always; at Stage 1 it is finite everywhere, so the `⊤` region is
empty until Phase 7 makes it real.

## 0.4 The candidate family

A `Finset` of quadratics on `E` (deduplicated — see 0.7). For each piece `p` with
vertex set `V`:

| branch | included when | value |
|---|---|---|
| `g_v`, `v ∈ V` | always | `⟪s, v⟫ - q v` — affine |
| `g_vw`, `{v,w} ⊆ V`, `d = w - v` | `α := ⟪d, Q d⟫ > 0` | `g_v s + L(s)² / (2α)` with `L(s) = ⟪s, d⟫ - ⟪∇q(v), d⟫` — rank-one quadratic |
| `g_int` | `Q` invertible | `(1/2)⟪s - β, Q⁻¹ (s - β)⟫ - γ` — full-rank quadratic |

Note the pairs range over **all** pairs of vertices, not over edges of the
polygon. That is what the Carathéodory route buys: no polytope face theory is
needed, and the argument generalizes to `ℝⁿ` unchanged.

## 0.5 Cells

```lean
def active (f) (s) : Finset (Quad E) :=
  (cand f).filter (fun q => (q.eval s : EReal) = conj f s)

def cell (f) (S : Finset (Quad E)) : Set E := {s | active f s = S}
```

`cell ∅` is exactly `{s | conj f s = ⊤}`, so the distinguished `+∞` cell needs no
special construction — it falls out of the same definition. Cells are pairwise
disjoint and cover `E` **by construction**; the content of the theorem is that
`cell ∅` is the complement of `dom (conj f)` (that is the selection lemma) and
that multi-active cells sit inside conics.

## 0.6 The key lemma

> **Selection.** For every `s` with `conj f s ≠ ⊤`, there is `q ∈ cand f` with
> `(q.eval s : EReal) = conj f s`.

This is Lemma 1 + Theorem 1 of `../CONJ_FIELD_PROOF.md`, restated. **It is not
the statement `conj f = max over cand f`, which is false** — an edge branch
`g_vw` overshoots the true supremum when its stationary point falls outside the
segment. Selection is what the cell construction actually needs, and it is what
the informal Theorem 1 means by "a max selects; it never manufactures".

Proof route. The whole risk of the project lives here.

* **S1** `conj f s = max over pieces p of (⨆ x ∈ p.T, ⟪s,x⟫ - p.q x)`. A supremum
  over a finite union splits as a finite max of suprema.
* **S2** Each per-piece sup is **attained**: `p.T` is compact (convex hull of a
  finite set) and nonempty, and `ψ := fun x => ⟪s,x⟫ - p.q x` is continuous.
* **S3** Carathéodory: `x* ∈ convexHull ↑V` gives an affinely independent
  `W ⊆ V` with `x* ∈ convexHull ↑W`, so `|W| ≤ 3`. Take `W` of **minimal**
  cardinality; then all barycentric coordinates of `x*` in `W` are strictly
  positive, i.e. `x*` lies in the relative interior of the simplex `conv W`.
* **S4** `x*` maximizes `ψ` over `conv W` as well, since it maximizes over the
  larger set `p.T ⊇ conv W`.
* **S5** First-order condition: for every `d` in the direction space of
  `affineSpan ℝ W`, `⟪s - Q x* - β, d⟫ = 0`. (Relative-interior maximum.)
* **S6** `|W| = 1`: `x*` is a vertex, value `= g_v s`.
* **S7** `|W| = 2`, `d = w - v`, `α = ⟪d, Q d⟫`. The second-order condition gives
  `α ≥ 0`. If `α > 0`, S5 forces the stationary point and the value is `g_vw s`.
  If `α = 0`, `ψ` is affine along `d` with zero slope, so the value is `g_v s`.
* **S8** `|W| = 3`, so `affineSpan ℝ W = ⊤` in dimension 2, and S5 gives
  `Q x* = s - β`. If `Q` is invertible the value is `g_int s`. If `Q` is
  singular, pick `0 ≠ d ∈ ker Q`; then `ψ (x* + t • d) = ψ x*` for every `t`,
  because `∇q(x*) = Q x* + β = s` and `Q d = 0`. Barycentric coordinates are
  affine in `t`, so at the largest admissible `t` one coordinate vanishes: the
  maximum is also attained on a proper face, contradicting the minimality of `W`
  and descending to `|W| ≤ 2`.
* **S9** Assemble by strong induction on `|W|`.

S8 is the riskiest single step. Write it first, as a standalone file with its own
sanity `example`.

## 0.7 Conics

```lean
def IsConic (C : Set E) : Prop := ∃ a b c d e g : ℝ,
  ¬(a = 0 ∧ b = 0 ∧ c = 0 ∧ d = 0 ∧ e = 0 ∧ g = 0) ∧
  C = {s | a*(s 0)^2 + b*(s 0)*(s 1) + c*(s 1)^2 + d*(s 0) + e*(s 1) + g = 0}
```

Candidates are stored as a `Finset` of quadratic **data**, so two distinct
elements are distinct quadratics and their difference is a nonzero quadratic
polynomial — which is what makes `{q_i = q_j}` a conic rather than the whole
plane. Deduplication is load-bearing, not cosmetic.

Classification to prove, by `Δ := b² - 4ac` together with the 3x3 determinant:

| `Δ` | nondegenerate | degenerate |
|---|---|---|
| `< 0` | ellipse (incl. circle) | a single point, or empty |
| `= 0` | parabola | one line, two parallel lines, or empty |
| `> 0` | hyperbola | two crossing lines |

plus the case `a = b = c = 0`: a line, or empty. The whole plane is excluded by
the nonzero-coefficient condition together with dedup.

Degenerate conics are **legal edges**: Theorem 3 of the informal proof shows that
adjacent pieces produce exactly line pairs, and those are real edges of `f*`.

## 0.8 The theorem

```lean
theorem conj_isQuaCon (f : QuaPol E) :
    -- cells are pairwise disjoint and cover the plane
    (∀ S T, S ≠ T → Disjoint (cell f S) (cell f T))
  ∧ (⋃ S, cell f S) = Set.univ
    -- on each cell, f* is one fixed quadratic
  ∧ (∀ S q, q ∈ S → ∀ s ∈ cell f S, conj f s = q.eval s)
    -- only finitely many cells are nonempty
  ∧ {S | (cell f S).Nonempty}.Finite
    -- the +infinity region is exactly the empty-activity cell
  ∧ cell f ∅ = {s | conj f s = ⊤}
    -- every multi-active cell sits inside a conic
  ∧ (∀ S, 2 ≤ S.card → ∀ q₁ ∈ S, ∀ q₂ ∈ S, q₁ ≠ q₂ →
        cell f S ⊆ {s | q₁.eval s = q₂.eval s} ∧ IsConic {s | q₁.eval s = q₂.eval s})
```

The first three conjuncts are near-free once the cells are defined by activity;
the fourth is a cardinality count; the fifth **is** the selection lemma; the
sixth is 0.7.

What the theorem deliberately does **not** claim, per the agreed regularity
level: nothing about dimension, connectedness, arcs, or a face-to-face CW
structure. `DECISIONS.md` records why, and what it would take to add.

---

# Phase 1 — Scaffold  ✅ DONE
- [x] `lakefile.toml` pinned to mathlib `v4.33.0`, `lean-toolchain` `v4.33.0`
- [x] `lake exe cache get` succeeds from `~/.cache/mathlib`; `lake build` green
- [x] `.gitignore` for `.lake/`; first commit

# Phase 2 — Quadratics and conics (no conjugate yet)  — core done, normal forms open
- [x] `Quad` structure, `eval`, `DecidableEq`, extensionality
- [x] `IsConic`; the zero set of any nonzero real quadratic polynomial is one
- [ ] the discriminant trichotomy of 0.7, degenerate branches included
- [x] sanity `example`s: the unit circle, `y = x²`, `xy = 1`, a line pair

# Phase 3 — QuaPol, conjugate, candidates, full statement  ✅ DONE
- [x] `QuaPiece`, `QuaPol`, `eval`, `conj`
- [x] sanity `example`s: conjugate of one piece with `Q = I` on a point; on a segment
- [x] `cand`, `active`, `cell`
- [x] `conj_isQuaCon` stated in full with `sorry`; `SORRY_LEDGER.md` seeded

# Phase 4 — Selection lemma (the work)
- [x] S8 `|W| = 3`, singular-`Q` descent — **highest risk, write first**
- [x] S1 supremum over a finite union splits
- [x] S2 attainment by compactness
- [x] S3 (done differently: induction on the face, no Carathéodory) gives positive barycentric coordinates
- [x] S5 first-order condition at a relative-interior maximum
- [x] S6 `|W| = 1`
- [ ] S7 `|W| = 2`, both `α` branches
- [x] S9 induction; selection lemma `sorry`-free

# Phase 5 — Assembly  ✅ DONE
- [x] disjoint, cover, quadratic-on-cell, finiteness
- [x] `cell ∅ = {s | conj f s = ⊤}`
- [x] conic containment; `conj_isQuaCon` `sorry`-free
- [x] `#print axioms conj_isQuaCon` clean

# Phase 6 — The witnesses that the theorem is not vacuous
- [ ] formalize the three-piece example of `../doc/QuaConExample.md` §2
- [ ] check in Lean that its `{g₁ = g₃}` really has `Δ < 0` (`norm_num`)
- [ ] this is what shows the theorem is not vacuous, and that `QuaPar` is too narrow
- [ ] add the two `det₃ ≠ 0` witnesses below: they are the reason `IsConic`'s
      classification must be decided algebraically and not by appearance

Concrete `norm_num` targets, all from the five-piece example's `f*` (faces `U₁..U₇`
of `../doc/QuaCon.svg` row 3, re-derived from the primal data — `DECISIONS.md`,
2026-08-22). Every one of these is a pure integer arithmetic check:

| pair | `Δ = B² - 4AC` | `det₃` | must classify as |
|---|---|---|---|
| `U₁\|U₆` | `+45054240` | `+260148962304` | hyperbola, **non**-degenerate |
| `U₁\|U₂` | `+49275072` | `+2474882726976` | hyperbola, **non**-degenerate |
| `U₃\|U₆` | `-13064` | `-8650208` | ellipse, non-degenerate |
| `U₂\|U₄` | `+14140` | `0` | **degenerate**: crossing line pair |
| `U₃\|U₆` (again) | `-13064` | `-8650208`, and `(a+c)·det₃ < 0` | a **real** ellipse, not an empty one |
| `U₁\|U₅` | `0` | `0` | **degenerate**: parallel/repeated lines |

`det₃` is cubic in the coefficients, so its sign depends on the overall sign
chosen for the equation; only its **vanishing** is an invariant, and that is what
the degeneracy test uses. The single criterion that reads a sign is real-versus-
empty ellipse, in the scale-invariant form `(a+c)·det₃ < 0`.

The first two are the trap: both render as straight segments in the figure, and
both are genuine hyperbolas. A classification predicate that got these wrong would
still typecheck, so they belong in the sanity `example`s that `CLAUDE.md` →
Verification point 3 requires. The full 21-pair census is 7 hyperbolas, 5 ellipses,
5 line pairs, 4 parallel/repeated lines, and **no parabola** — so a parabola
witness has to come from panel (3b), the indefinite piece `q = x₁x₂` on
`(0,0),(1,1),(2,0)`, whose conjugate has the parabolic edge `¼(s₁+s₂)² = 2s₁`.

# Phase 7 — Stage 2, unbounded pieces
- [ ] extend `QuaPiece` with a recession cone (vertices plus rays)
- [ ] Frank–Wolfe style attainment, or a direct argument for quadratics
- [ ] `dom (conj f)` proper; the `⊤` cell becomes nonempty
- [ ] re-prove selection in the extended setting

# Phase 8 — Write-up
- [ ] `PROOF_NOTES.md` mapping each Lean lemma to its informal counterpart
- [ ] decide with Yves whether this becomes a paper or a note inside CCA2

## Decision log

Decisions, dead ends and rejected approaches go in **`DECISIONS.md`**, not here.
This file is the roadmap; that one is the record of what was ruled out and why.

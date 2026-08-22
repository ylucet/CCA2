/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.QuaCon

/-!
# Which conics arise: the shape of a tie set, from the shapes of the branches

`conj_isQuaCon` says every multi-active cell lies inside a conic. This file says
**which** conic, as a function of what kind of branch the two active quadratics
are. The answer is sharp, and it explains the gap in [JOGO] Theorem 6.

## The three shapes of a branch

A candidate quadratic in `s` has a quadratic part of rank 0, 1 or 2:

| branch | quadratic part | rank |
|---|---|---|
| `vertexBranch q v` | zero — the branch is affine | 0 |
| `edgeBranch q v w` | `(1/2α)·δδᵀ`, positive semidefinite | 1 |
| `interiorBranch q` | `(1/2)·H⁻¹` | 2 |

Since the discriminant of a tie set depends only on the two quadratic *parts*,
the type of the conic is decided by those ranks — and, crucially, this holds
**across pieces**, not just within one.

## The type theorems

With `Δ = b² - 4ac`, so `Δ < 0` elliptic, `Δ = 0` parabolic, `Δ > 0` hyperbolic:

| pair | `Δ` of the difference | consequence |
|---|---|---|
| vertex, vertex | `0`, and the quadratic part vanishes | a line |
| vertex, edge | `0` | parabolic |
| edge, edge | `cross(δ₁,δ₂)² / (α₁α₂) ≥ 0` | **never elliptic** |
| interior, vertex | `-1 / hessDet` | elliptic iff the Hessian is definite |
| interior, edge / interior, interior | unconstrained | any type |

## The headline

`disc_nonneg_of_flat` and `not_flat_of_disc_neg`:

> **An elliptical tie conic requires an interior branch.**

That is exactly the gap in [JOGO] Theorem 6 (`../doc/QuaConExample.md` §1). Its
proof argues that "when we compare two functions we always get one of them as
linear", which forces `Δ ≥ 0` and hence a parabolic subdivision. An elliptical
edge is precisely a comparison in which *neither* function is flat — which
happens as soon as two interior branches meet, i.e. as soon as two non-adjacent
pieces are compared.

## The degeneracy theorems

Type is not the whole story: a conic with `Δ > 0` may be a hyperbola or a
crossing pair of lines, and `det₃ = 0` is what tells them apart. Degeneracy,
unlike type, is **not** determined by the quadratic parts, so these theorems are
sensitive to whether the branches come from the *same* piece:

* `det3_interiorBranch_sub_vertexBranch_self` — interior against a vertex branch
  of the **same** piece is always degenerate, and elliptic when the Hessian is
  definite, hence a **single point**. This is why
  `../doc/QuaConExample.md` §3.3 finds three pairs that "touch at a single point
  and are not edges". Across pieces this fails, and §3.3's `I4|V3` is a genuine
  ellipse — the two branches there come from pieces 4 and 1.
* `det3_vertexBranch_sub_edgeBranch_self` — vertex against edge of the same
  piece is a genuine parabola exactly when the vertex is off the edge's line.
-/

namespace QuaConProof

open scoped Classical

/-! ### Naming the seven cases -/

/-- The seven possible shapes of the zero set of a real quadratic in two
variables, as decided by `Quad.disc` and `Quad.det3`.

`linesParallelOrRepeated` also covers the case of a single line, which is what a
quadratic with vanishing quadratic part cuts out. -/
inductive ConicKind
  /-- `Δ < 0`, `det₃ ≠ 0` -/
  | ellipse
  /-- `Δ = 0`, `det₃ ≠ 0` -/
  | parabola
  /-- `Δ > 0`, `det₃ ≠ 0` -/
  | hyperbola
  /-- `Δ > 0`, `det₃ = 0`: two lines meeting at a point -/
  | crossingLines
  /-- `Δ = 0`, `det₃ = 0`: one line, two parallel lines, or empty -/
  | linesParallelOrRepeated
  /-- `Δ < 0`, `det₃ = 0`: a single point, or empty -/
  | pointOrEmpty
  /-- the zero quadratic, whose zero set is the whole plane -/
  | wholePlane
  deriving DecidableEq, Repr

/-- The shape of `{q = 0}`, read off from the two invariants. -/
noncomputable def Quad.kind (q : Quad) : ConicKind :=
  if q = 0 then .wholePlane
  else if q.det3 = 0 then
    if q.disc < 0 then .pointOrEmpty
    else if q.disc = 0 then .linesParallelOrRepeated
    else .crossingLines
  else
    if q.disc < 0 then .ellipse
    else if q.disc = 0 then .parabola
    else .hyperbola

lemma Quad.kind_eq_ellipse_iff (q : Quad) :
    q.kind = .ellipse ↔ q ≠ 0 ∧ q.det3 ≠ 0 ∧ q.disc < 0 := by
  unfold Quad.kind
  split_ifs <;> simp_all

/-! ### The quadratic parts of the three branches -/

@[simp] lemma edgeBranch_a (q : Quad) (v w : Plane) :
    (edgeBranch q v w).a = (w.1 - v.1) ^ 2 / (2 * edgeCurv q v w) := rfl
@[simp] lemma edgeBranch_b (q : Quad) (v w : Plane) :
    (edgeBranch q v w).b = (w.1 - v.1) * (w.2 - v.2) / edgeCurv q v w := rfl
@[simp] lemma edgeBranch_c (q : Quad) (v w : Plane) :
    (edgeBranch q v w).c = (w.2 - v.2) ^ 2 / (2 * edgeCurv q v w) := rfl

@[simp] lemma interiorBranch_a (q : Quad) : (interiorBranch q).a = q.c / q.hessDet := rfl
@[simp] lemma interiorBranch_b (q : Quad) : (interiorBranch q).b = -q.b / q.hessDet := rfl
@[simp] lemma interiorBranch_c (q : Quad) : (interiorBranch q).c = q.a / q.hessDet := rfl

/-- **The quadratic part of an edge branch is degenerate.** Its determinant
`ac - b²/4` vanishes: the form has rank one. -/
theorem edgeBranch_quadratic_part_singular (q : Quad) (v w : Plane)
    (hα : edgeCurv q v w ≠ 0) :
    (edgeBranch q v w).b ^ 2 - 4 * (edgeBranch q v w).a * (edgeBranch q v w).c = 0 := by
  simp only [edgeBranch_a, edgeBranch_b, edgeBranch_c]
  field_simp
  ring

/-! ### Type theorems: what the discriminant of a tie set must be

All of these hold **across pieces**: the discriminant of a difference depends
only on the two quadratic parts, and those are determined by which kind of
branch each side is. -/

/-- **Vertex against vertex: no quadratic part at all.** The tie set is a line,
or empty, or the whole plane — never a curve. -/
@[simp] theorem disc_vertexBranch_sub_vertexBranch (q₁ q₂ : Quad) (v w : Plane) :
    (vertexBranch q₁ v - vertexBranch q₂ w).disc = 0 := by
  simp only [Quad.disc, vertexBranch, Quad.sub_a, Quad.sub_b, Quad.sub_c]
  ring

@[simp] theorem det3_vertexBranch_sub_vertexBranch (q₁ q₂ : Quad) (v w : Plane) :
    (vertexBranch q₁ v - vertexBranch q₂ w).det3 = 0 := by
  simp only [Quad.det3, vertexBranch, Quad.sub_a, Quad.sub_b, Quad.sub_c, Quad.sub_d,
    Quad.sub_e, Quad.sub_f]
  ring

/-- **Vertex against edge is parabolic.** Subtracting an affine function from a
rank-one quadratic leaves a rank-one quadratic, so `Δ = 0`. -/
@[simp] theorem disc_vertexBranch_sub_edgeBranch (q₁ q₂ : Quad) (z v w : Plane)
    (hα : edgeCurv q₂ v w ≠ 0) :
    (vertexBranch q₁ z - edgeBranch q₂ v w).disc = 0 := by
  simp only [Quad.disc, vertexBranch, edgeBranch_a, edgeBranch_b, edgeBranch_c,
    Quad.sub_a, Quad.sub_b, Quad.sub_c]
  field_simp
  ring

/-- **Edge against edge, in closed form.** The discriminant is a square over the
product of the two curvatures. -/
theorem disc_edgeBranch_sub_edgeBranch (q₁ q₂ : Quad) (v w z u : Plane)
    (h₁ : edgeCurv q₁ v w ≠ 0) (h₂ : edgeCurv q₂ z u ≠ 0) :
    (edgeBranch q₁ v w - edgeBranch q₂ z u).disc
      = cross (w - v) (u - z) ^ 2 / (edgeCurv q₁ v w * edgeCurv q₂ z u) := by
  simp only [Quad.disc, edgeBranch_a, edgeBranch_b, edgeBranch_c, Quad.sub_a, Quad.sub_b,
    Quad.sub_c, cross, Prod.fst_sub, Prod.snd_sub]
  field_simp
  ring

/-- **Edge against edge is never elliptic.** -/
theorem disc_edgeBranch_sub_edgeBranch_nonneg (q₁ q₂ : Quad) (v w z u : Plane)
    (h₁ : 0 < edgeCurv q₁ v w) (h₂ : 0 < edgeCurv q₂ z u) :
    0 ≤ (edgeBranch q₁ v w - edgeBranch q₂ z u).disc := by
  rw [disc_edgeBranch_sub_edgeBranch q₁ q₂ v w z u (ne_of_gt h₁) (ne_of_gt h₂)]
  positivity

/-- **Interior against vertex: the discriminant is `-1/hessDet`.**

So the tie set is elliptic exactly when the Hessian is definite, and hyperbolic
exactly when it is indefinite. The vertex branch may come from any piece — it
contributes no quadratic part. -/
theorem disc_interiorBranch_sub_vertexBranch (q₁ q₂ : Quad) (v : Plane)
    (hD : q₁.hessDet ≠ 0) :
    (interiorBranch q₁ - vertexBranch q₂ v).disc = -(1 / q₁.hessDet) := by
  simp only [Quad.disc, interiorBranch_a, interiorBranch_b, interiorBranch_c, vertexBranch,
    Quad.sub_a, Quad.sub_b, Quad.sub_c, Quad.hessDet] at *
  field_simp
  ring

/-! ### The headline: an ellipse needs an interior branch -/

/-- A **flat** branch is one whose quadratic part is degenerate: a vertex branch
(rank zero) or an edge branch (rank one). Equivalently, one of the two branch
types that [JOGO] Theorem 6's proof assumes it is always comparing. -/
inductive IsFlatBranch : Quad → Prop
  | vertex (q : Quad) (v : Plane) : IsFlatBranch (vertexBranch q v)
  | edge (q : Quad) (v w : Plane) (h : 0 < edgeCurv q v w) : IsFlatBranch (edgeBranch q v w)

/-- **Two flat branches never tie along an ellipse.** -/
theorem disc_nonneg_of_flat {b₁ b₂ : Quad} (h₁ : IsFlatBranch b₁) (h₂ : IsFlatBranch b₂) :
    0 ≤ (b₁ - b₂).disc := by
  cases h₁ with
  | vertex q₁ v =>
      cases h₂ with
      | vertex q₂ w => rw [disc_vertexBranch_sub_vertexBranch]
      | edge q₂ z u h => rw [disc_vertexBranch_sub_edgeBranch _ _ _ _ _ (ne_of_gt h)]
  | edge q₁ v w h =>
      cases h₂ with
      | vertex q₂ z =>
          have hkey : (edgeBranch q₁ v w - vertexBranch q₂ z).disc = 0 := by
            simp only [Quad.disc, vertexBranch, edgeBranch_a, edgeBranch_b, edgeBranch_c,
              Quad.sub_a, Quad.sub_b, Quad.sub_c]
            have : edgeCurv q₁ v w ≠ 0 := ne_of_gt h
            field_simp
            ring
          rw [hkey]
      | edge q₂ z u h' => exact disc_edgeBranch_sub_edgeBranch_nonneg q₁ q₂ v w z u h h'

/-- **The [JOGO] gap, in one line.**

If the conic on which two branches agree is an **ellipse**, then at least one of
them is an interior branch — neither can be flat.

[JOGO] Theorem 6 argues that "when we compare two functions we always get one of
them as linear". Under that assumption `disc_nonneg_of_flat` applies and the
subdivision really is parabolic. The assumption fails as soon as Step 3b compares
two non-adjacent pieces, because then both compared functions can be interior
branches. See `../doc/QuaConExample.md` §1. -/
theorem not_flat_of_disc_neg {b₁ b₂ : Quad} (h : (b₁ - b₂).disc < 0) :
    ¬ IsFlatBranch b₁ ∨ ¬ IsFlatBranch b₂ := by
  by_contra hc
  push Not at hc
  exact absurd (disc_nonneg_of_flat hc.1 hc.2) (not_le.2 h)

/-- Stated the other way round: a positive definite piece contributes an interior
branch whose tie with **any** vertex branch is elliptic. So a single positive
definite piece already puts an elliptic conic in the candidate list. -/
theorem disc_neg_of_hessDet_pos {q₁ q₂ : Quad} (v : Plane) (hD : 0 < q₁.hessDet) :
    (interiorBranch q₁ - vertexBranch q₂ v).disc < 0 := by
  rw [disc_interiorBranch_sub_vertexBranch q₁ q₂ v (ne_of_gt hD)]
  simp only [neg_lt_zero]
  positivity

/-! ### Degeneracy theorems: same piece versus across pieces

Type is decided by the quadratic parts alone; **degeneracy is not**. These two
are sensitive to whether the branches come from the same piece. -/

/-- **Interior against a vertex branch of the SAME piece is degenerate.**

Combined with `disc_interiorBranch_sub_vertexBranch`, for a definite Hessian this
makes the tie set a **single point** — the vertex branch's plane is tangent to
the interior branch's paraboloid. Those pairs are therefore not edges, which is
what `../doc/QuaConExample.md` §3.3 observes for `I1|V3`, `I2|V3` and `I5|V7`.

The hypothesis that both branches come from the same `q` is essential: §3.3's
`I4|V3` compares pieces 4 and 1 and is a genuine, non-degenerate ellipse. -/
theorem det3_interiorBranch_sub_vertexBranch_self (q : Quad) (v : Plane)
    (hD : q.hessDet ≠ 0) :
    (interiorBranch q - vertexBranch q v).det3 = 0 := by
  have h2 : q.hessDet ^ 2 ≠ 0 := pow_ne_zero 2 hD
  have h3 : q.hessDet ^ 3 ≠ 0 := pow_ne_zero 3 hD
  simp only [Quad.det3, interiorBranch, vertexBranch, Quad.eval,
    Quad.sub_a, Quad.sub_b, Quad.sub_c, Quad.sub_d, Quad.sub_e, Quad.sub_f]
  field_simp
  ring_nf
  rw [Quad.hessDet] at *
  ring

/-- **Vertex against edge of the same piece: a genuine parabola exactly when the
vertex is off the edge's line.**

`det₃` is `cross(w-v, z-v)² / (2α)`, so it vanishes precisely when `z` lies on
the line through `v` and `w`. This is the "edge versus *far* vertex branch"
labelled in `../doc/QuaCon.svg`: far means off the line. -/
theorem det3_vertexBranch_sub_edgeBranch_self (q : Quad) (z v w : Plane)
    (_hα : edgeCurv q v w ≠ 0) :
    (vertexBranch q z - edgeBranch q v w).det3
      = cross (w - v) (z - v) ^ 2 / (2 * edgeCurv q v w) := by
  simp only [Quad.det3, vertexBranch, edgeBranch, Quad.eval, edgeCurv, Quad.alongCurv,
    Quad.gradAt, dot, cross, Quad.sub_a, Quad.sub_b, Quad.sub_c, Quad.sub_d, Quad.sub_e,
    Quad.sub_f, Prod.fst_sub, Prod.snd_sub] at *
  field_simp
  ring

/-! ### Theorem 3: adjacency forces degeneracy

`../CONJ_FIELD_PROOF.md` Theorem 3. Two pieces sharing an edge, with `f`
continuous across it, have `q₂ - q₁ = l·m` with `l` the affine form vanishing on
the shared edge. The hypothesis below is exactly that factorisation, and it is
purely **algebraic** — continuity is what produces it, but the theorem does not
mention continuity, so the input class needs no extra hypothesis. -/

set_option maxHeartbeats 2000000 in
/-- **Theorem 3, in coefficients.** If `q₂ - q₁` is the product of two affine
forms, the two interior branches tie along a **degenerate** conic. -/
theorem det3_interiorBranch_sub_interiorBranch_of_coeffs {q₁ q₂ : Quad}
    (hD₁ : q₁.hessDet ≠ 0) (hD₂ : q₂.hessDet ≠ 0) {u₁ u₂ c₀ v₁ v₂ d₀ : ℝ}
    (ha : q₂.a = q₁.a + u₁ * v₁)
    (hb : q₂.b = q₁.b + (u₁ * v₂ + u₂ * v₁))
    (hc : q₂.c = q₁.c + u₂ * v₂)
    (hd : q₂.d = q₁.d + (u₁ * d₀ + c₀ * v₁))
    (he : q₂.e = q₁.e + (u₂ * d₀ + c₀ * v₂))
    (hf : q₂.f = q₁.f + c₀ * d₀) :
    (interiorBranch q₂ - interiorBranch q₁).det3 = 0 := by
  have h1 : q₁.hessDet ≠ 0 := hD₁
  have h2 : q₂.hessDet ≠ 0 := hD₂
  simp only [Quad.det3, interiorBranch, Quad.sub_a, Quad.sub_b, Quad.sub_c, Quad.sub_d,
    Quad.sub_e, Quad.sub_f]
  field_simp
  rw [Quad.hessDet, Quad.hessDet] at *
  rw [ha, hb, hc, hd, he, hf] at *
  ring

/-- **Theorem 3.** If the difference of two pieces' quadratics factors as a
product of affine forms -- which is what continuity across a shared edge forces --
then their interior branches tie along a degenerate conic: a pair of lines, never
an ellipse, parabola or hyperbola.

This is the case [JOGO] Theorem 6's proof does cover, and it is sound. What its
proof does not cover is the *non*-adjacent comparison, where no such
factorisation exists and `not_flat_of_disc_neg` allows an ellipse. -/
theorem det3_interiorBranch_sub_of_factorisation {q₁ q₂ : Quad}
    (hD₁ : q₁.hessDet ≠ 0) (hD₂ : q₂.hessDet ≠ 0) (u₁ u₂ c₀ v₁ v₂ d₀ : ℝ)
    (hfac : ∀ x : Plane, q₂.eval x - q₁.eval x
      = (u₁ * x.1 + u₂ * x.2 + c₀) * (v₁ * x.1 + v₂ * x.2 + d₀)) :
    (interiorBranch q₂ - interiorBranch q₁).det3 = 0 := by
  have hsub : q₂ - q₁ = (⟨u₁ * v₁, u₁ * v₂ + u₂ * v₁, u₂ * v₂,
      u₁ * d₀ + c₀ * v₁, u₂ * d₀ + c₀ * v₂, c₀ * d₀⟩ : Quad) := by
    refine Quad.eq_of_eval_eq fun x => ?_
    rw [Quad.eval_sub, hfac x]
    simp only [Quad.eval]
    ring
  have ha : q₂.a = q₁.a + u₁ * v₁ := by
    have := congrArg Quad.a hsub; simp only [Quad.sub_a] at this; linarith
  have hb : q₂.b = q₁.b + (u₁ * v₂ + u₂ * v₁) := by
    have := congrArg Quad.b hsub; simp only [Quad.sub_b] at this; linarith
  have hc : q₂.c = q₁.c + u₂ * v₂ := by
    have := congrArg Quad.c hsub; simp only [Quad.sub_c] at this; linarith
  have hd : q₂.d = q₁.d + (u₁ * d₀ + c₀ * v₁) := by
    have := congrArg Quad.d hsub; simp only [Quad.sub_d] at this; linarith
  have he : q₂.e = q₁.e + (u₂ * d₀ + c₀ * v₂) := by
    have := congrArg Quad.e hsub; simp only [Quad.sub_e] at this; linarith
  have hf : q₂.f = q₁.f + c₀ * d₀ := by
    have := congrArg Quad.f hsub; simp only [Quad.sub_f] at this; linarith
  exact det3_interiorBranch_sub_interiorBranch_of_coeffs hD₁ hD₂ ha hb hc hd he hf

end QuaConProof

/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Linarith

/-!
# Real quadratics in two variables

A `Quad` is the **coefficient vector** of a real quadratic polynomial in two
variables,

`a s₁² + b s₁s₂ + c s₂² + d s₁ + e s₂ + f`.

Storing coefficients rather than the function they define is deliberate and
load-bearing for the project. `PROJECT_PLAN.md` §0.7: the candidate family is a
`Finset Quad`, so two distinct members are distinct *coefficient vectors*, their
difference is a nonzero `Quad`, and the set where they agree is therefore the
zero set of a genuinely nonzero polynomial — a conic, not the whole plane.
`eq_zero_of_eval_eq_zero` below is what turns that data-level distinctness into
the function-level statement, and it is the reason the deduplication matters.

## Main results

* `Quad.eval` — the polynomial function a coefficient vector defines.
* `Quad.eq_zero_of_eval_eq_zero` — a quadratic vanishing at every point of the
  plane has all six coefficients zero. This is what stops a "conic" from being
  the whole plane.
* `Quad.translate` — precomposition with a translation, together with the fact
  that it does not disturb the quadratic part.

## Implementation notes

The plane is `ℝ × ℝ`, not `EuclideanSpace ℝ (Fin 2)`. Nothing in this file or in
`Conic.lean` needs an inner product, and `ℝ × ℝ` makes the coefficient algebra
directly accessible to `ring`. See `DECISIONS.md`, 2026-08-22.
-/

namespace QuaConProof

/-- The plane. Carries its usual topology and real vector space structure; no
inner product is needed for the conic theory. -/
abbrev Plane := ℝ × ℝ

/-- The coefficient vector of a real quadratic in two variables,
`a s₁² + b s₁s₂ + c s₂² + d s₁ + e s₂ + f`. -/
@[ext]
structure Quad where
  /-- coefficient of `s₁²` -/
  a : ℝ
  /-- coefficient of `s₁s₂` -/
  b : ℝ
  /-- coefficient of `s₂²` -/
  c : ℝ
  /-- coefficient of `s₁` -/
  d : ℝ
  /-- coefficient of `s₂` -/
  e : ℝ
  /-- constant coefficient -/
  f : ℝ

namespace Quad

/-- `Quad` has no decidable equality (its fields are real), so one is supplied
classically. It is needed only to form `Finset Quad`. -/
noncomputable instance : DecidableEq Quad := Classical.decEq _

instance : Zero Quad := ⟨⟨0, 0, 0, 0, 0, 0⟩⟩

@[simp] lemma zero_a : (0 : Quad).a = 0 := rfl
@[simp] lemma zero_b : (0 : Quad).b = 0 := rfl
@[simp] lemma zero_c : (0 : Quad).c = 0 := rfl
@[simp] lemma zero_d : (0 : Quad).d = 0 := rfl
@[simp] lemma zero_e : (0 : Quad).e = 0 := rfl
@[simp] lemma zero_f : (0 : Quad).f = 0 := rfl

instance : Sub Quad :=
  ⟨fun p q => ⟨p.a - q.a, p.b - q.b, p.c - q.c, p.d - q.d, p.e - q.e, p.f - q.f⟩⟩

@[simp] lemma sub_a (p q : Quad) : (p - q).a = p.a - q.a := rfl
@[simp] lemma sub_b (p q : Quad) : (p - q).b = p.b - q.b := rfl
@[simp] lemma sub_c (p q : Quad) : (p - q).c = p.c - q.c := rfl
@[simp] lemma sub_d (p q : Quad) : (p - q).d = p.d - q.d := rfl
@[simp] lemma sub_e (p q : Quad) : (p - q).e = p.e - q.e := rfl
@[simp] lemma sub_f (p q : Quad) : (p - q).f = p.f - q.f := rfl

/-- The polynomial function defined by a coefficient vector. -/
def eval (q : Quad) (s : Plane) : ℝ :=
  q.a * s.1 ^ 2 + q.b * s.1 * s.2 + q.c * s.2 ^ 2 + q.d * s.1 + q.e * s.2 + q.f

lemma eval_apply (q : Quad) (x y : ℝ) :
    q.eval (x, y) = q.a * x ^ 2 + q.b * x * y + q.c * y ^ 2 + q.d * x + q.e * y + q.f := rfl

@[simp] lemma eval_zero (s : Plane) : (0 : Quad).eval s = 0 := by
  simp [eval]

@[simp] lemma eval_sub (p q : Quad) (s : Plane) : (p - q).eval s = p.eval s - q.eval s := by
  simp only [eval, sub_a, sub_b, sub_c, sub_d, sub_e, sub_f]; ring

/-- A `Quad` is zero exactly when all six coefficients are. -/
lemma eq_zero_iff (q : Quad) :
    q = 0 ↔ q.a = 0 ∧ q.b = 0 ∧ q.c = 0 ∧ q.d = 0 ∧ q.e = 0 ∧ q.f = 0 := by
  constructor
  · rintro rfl; exact ⟨rfl, rfl, rfl, rfl, rfl, rfl⟩
  · rintro ⟨h1, h2, h3, h4, h5, h6⟩; ext <;> simp [h1, h2, h3, h4, h5, h6]

lemma sub_eq_zero_iff (p q : Quad) : p - q = 0 ↔ p = q := by
  rw [eq_zero_iff]
  constructor
  · rintro ⟨h1, h2, h3, h4, h5, h6⟩
    simp only [sub_a, sub_b, sub_c, sub_d, sub_e, sub_f, sub_eq_zero] at h1 h2 h3 h4 h5 h6
    ext <;> assumption
  · rintro rfl; simp

/-- **A quadratic that vanishes everywhere has all coefficients zero.**

Six evaluations suffice: the origin pins `f`, the four points `(±1, 0)` and
`(0, ±1)` pin `a, d` and `c, e` in pairs, and `(1,1)` then pins `b`.

This is the lemma that makes `Conic.IsConic` non-vacuous: a nonzero coefficient
vector cannot cut out the whole plane. -/
lemma eq_zero_of_eval_eq_zero {q : Quad} (h : ∀ s : Plane, q.eval s = 0) : q = 0 := by
  have h00 := h (0, 0)
  have h10 := h (1, 0)
  have hm10 := h (-1, 0)
  have h01 := h (0, 1)
  have h0m1 := h (0, -1)
  have h11 := h (1, 1)
  simp only [eval_apply] at h00 h10 hm10 h01 h0m1 h11
  norm_num at h00 h10 hm10 h01 h0m1 h11
  rw [eq_zero_iff]
  refine ⟨by linarith, by linarith, by linarith, by linarith, by linarith, by linarith⟩

/-- Contrapositive of `eq_zero_of_eval_eq_zero`: a nonzero quadratic is nonzero
somewhere. -/
lemma exists_eval_ne_zero {q : Quad} (h : q ≠ 0) : ∃ s : Plane, q.eval s ≠ 0 := by
  by_contra hc
  push Not at hc
  exact h (eq_zero_of_eval_eq_zero hc)

/-- Two coefficient vectors defining the same function are equal. -/
lemma eq_of_eval_eq {p q : Quad} (h : ∀ s : Plane, p.eval s = q.eval s) : p = q := by
  rw [← sub_eq_zero_iff]
  refine eq_zero_of_eval_eq_zero fun s => ?_
  rw [eval_sub, h s, sub_self]

/-! ### Translation

Precomposition with `s ↦ s + v`. The quadratic part is untouched, which is why
`Conic.disc` is a translation invariant; `Conic.det3` is one too, and that is
proved in `Conic.lean`. -/

/-- The quadratic `s ↦ q.eval (s + v)`. -/
def translate (q : Quad) (v : Plane) : Quad where
  a := q.a
  b := q.b
  c := q.c
  d := 2 * q.a * v.1 + q.b * v.2 + q.d
  e := q.b * v.1 + 2 * q.c * v.2 + q.e
  f := q.eval v

@[simp] lemma translate_a (q : Quad) (v : Plane) : (q.translate v).a = q.a := rfl
@[simp] lemma translate_b (q : Quad) (v : Plane) : (q.translate v).b = q.b := rfl
@[simp] lemma translate_c (q : Quad) (v : Plane) : (q.translate v).c = q.c := rfl

@[simp] lemma eval_translate (q : Quad) (v s : Plane) :
    (q.translate v).eval s = q.eval (s + v) := by
  simp only [eval, translate, Prod.fst_add, Prod.snd_add]; ring

/-! ### Sanity checks

`CLAUDE.md` → Verification, point 3: every definition meant to model something
is checked against a hand computation, so that a typo cannot make a downstream
theorem vacuously true. -/

section Sanity

/-- `s₁² + s₂² - 1` evaluated at `(3, 4)` is `24`. -/
example : (⟨1, 0, 1, 0, 0, -1⟩ : Quad).eval (3, 4) = 24 := by norm_num [eval_apply]

/-- The cross term is read correctly: `s₁s₂` at `(2, 5)` is `10`. -/
example : (⟨0, 1, 0, 0, 0, 0⟩ : Quad).eval (2, 5) = 10 := by norm_num [eval_apply]

/-- The linear terms are read correctly, and are not swapped:
`3s₁ + 7s₂` at `(1, 0)` is `3`, not `7`. -/
example : (⟨0, 0, 0, 3, 7, 0⟩ : Quad).eval (1, 0) = 3 := by norm_num [eval_apply]

/-- A hand-computed general value: `2s₁² - 3s₁s₂ + 4s₂² + 5s₁ - 6s₂ + 7` at
`(1, 2)` is `2 - 6 + 16 + 5 - 12 + 7 = 12`. -/
example : (⟨2, -3, 4, 5, -6, 7⟩ : Quad).eval (1, 2) = 12 := by norm_num [eval_apply]

/-- Translation really is precomposition: `s₁²` shifted by `v = (1, 0)` is
`(s₁+1)²`, which at `s = (2, 0)` is `9`. -/
example : ((⟨1, 0, 0, 0, 0, 0⟩ : Quad).translate (1, 0)).eval (2, 0) = 9 := by
  norm_num [eval_apply, translate, eval]

end Sanity

end Quad

end QuaConProof

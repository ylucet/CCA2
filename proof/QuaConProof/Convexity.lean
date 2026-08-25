/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.FrankWolfe
import Mathlib.Topology.Semicontinuity.Basic
import Mathlib.Topology.Instances.EReal.Lemmas

/-!
# The conjugate is convex, lower semicontinuous, and has a convex domain

`TODO.md` items A1, A2, A3. None of this is special to a `QuaPol` — it is the
standard fact that a supremum of affine functions is convex and lower
semicontinuous — but it has to exist in this development before anything can be
said about the *bi*conjugate, and stating it for an `EReal`-valued function
needs a little care.

## How convexity is stated, and why

The obvious statement,

    f*(a·s₁ + b·s₂) ≤ a * f*(s₁) + b * f*(s₂),

is an inequality in `EReal`, which drags in `EReal` multiplication and its
`0 * ⊤ = 0` convention. The epigraph form avoids all of it: the bound `t₁`, `t₂`
and the combination `a * t₁ + b * t₂` are **real numbers**, and the only `EReal`
reasoning left is that a point where `f = ⊤` contributes `⊥` to the supremum.

`conj_le_of_combo` is that statement, and `convex_epigraph_conj` packages it as
convexity of the epigraph. `conj_add_le` recovers the direct `EReal` inequality
for readers who want it.

## Main results

* `conj_le_of_combo` — the epigraph inequality, and the workhorse.
* `convex_epigraph_conj` — the epigraph of `f*` is convex.
* `convex_dom_conj` — `dom f*` is convex, so the `⊤` region is the complement of
  a convex set. That sharpens the fifth conjunct of `conj_isQuaCon`, which so far
  only knows the region can be inhabited.
* `lowerSemicontinuous_conj` — `f*` is lsc.
-/

namespace QuaConProof

namespace QuaPol

/-! ### Convexity -/

lemma dot_combo (a b : ℝ) (s₁ s₂ x : Plane) :
    dot (a • s₁ + b • s₂) x = a * dot s₁ x + b * dot s₂ x := by
  simp only [dot, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
  ring

/-- A single term of the supremum defining `f*`, bounded at one point. -/
lemma term_le_conj (f : QuaPol) (s x : Plane) :
    ((dot s x : ℝ) : EReal) - f.eval x ≤ f.conj s :=
  le_iSup (fun y => ((dot s y : ℝ) : EReal) - f.eval y) x

/-- **The epigraph inequality.** If `t₁` and `t₂` bound `f*` at `s₁` and `s₂`,
then `a * t₁ + b * t₂` bounds it at the convex combination.

This is convexity of `f*`, in the form that keeps every arithmetic step in `ℝ`. -/
theorem conj_le_of_combo (f : QuaPol) {s₁ s₂ : Plane} {a b t₁ t₂ : ℝ}
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : a + b = 1)
    (h₁ : f.conj s₁ ≤ (t₁ : EReal)) (h₂ : f.conj s₂ ≤ (t₂ : EReal)) :
    f.conj (a • s₁ + b • s₂) ≤ ((a * t₁ + b * t₂ : ℝ) : EReal) := by
  rw [conj_def]
  refine iSup_le fun x => ?_
  by_cases hx : f.eval x = ⊤
  · rw [hx, EReal.sub_top]
    exact bot_le
  · -- at a point where `f` is finite the whole inequality is about real numbers
    set E : ℝ := (f.eval x).toReal with hEdef
    have hE : f.eval x = (E : EReal) := (EReal.coe_toReal hx (eval_ne_bot f x)).symm
    have key : ∀ (s : Plane) (t : ℝ), f.conj s ≤ (t : EReal) → dot s x - E ≤ t := by
      intro s t hst
      have hle : ((dot s x : ℝ) : EReal) - ((E : ℝ) : EReal) ≤ (t : EReal) := by
        rw [← hE]
        exact le_trans (term_le_conj f s x) hst
      rw [← EReal.coe_sub, EReal.coe_le_coe_iff] at hle
      exact hle
    have k₁ := key s₁ t₁ h₁
    have k₂ := key s₂ t₂ h₂
    have e1 : a * (dot s₁ x - E) ≤ a * t₁ := mul_le_mul_of_nonneg_left k₁ ha
    have e2 : b * (dot s₂ x - E) ≤ b * t₂ := mul_le_mul_of_nonneg_left k₂ hb
    have expand : a * (dot s₁ x - E) + b * (dot s₂ x - E)
        = a * dot s₁ x + b * dot s₂ x - E := by
      linear_combination (-E) * hab
    rw [hE, ← EReal.coe_sub, EReal.coe_le_coe_iff, dot_combo]
    linarith

/-- **The epigraph of `f*` is convex.** -/
theorem convex_epigraph_conj (f : QuaPol) :
    Convex ℝ {p : Plane × ℝ | f.conj p.1 ≤ ((p.2 : ℝ) : EReal)} := by
  rintro ⟨s₁, t₁⟩ h₁ ⟨s₂, t₂⟩ h₂ a b ha hb hab
  exact conj_le_of_combo f ha hb hab h₁ h₂

/-- **The domain of `f*` is convex**, so the region where `f*` is `⊤` is the
complement of a convex set.

`conj_isQuaCon`'s fifth conjunct identifies that region with the empty-activity
cell, and `Sanity.cell_empty_rayPol_nonempty` shows it can be inhabited; this
says what shape it has. -/
theorem convex_dom_conj (f : QuaPol) : Convex ℝ {s : Plane | f.conj s ≠ ⊤} := by
  intro s₁ h₁ s₂ h₂ a b ha hb hab
  have hc : ∀ s : Plane, f.conj s ≠ ⊤ → f.conj s ≤ (((f.conj s).toReal : ℝ) : EReal) :=
    fun s h => le_of_eq (EReal.coe_toReal h (conj_ne_bot f s)).symm
  have hle := conj_le_of_combo f ha hb hab (hc s₁ h₁) (hc s₂ h₂)
  show f.conj (a • s₁ + b • s₂) ≠ ⊤
  intro htop
  rw [htop, top_le_iff] at hle
  exact EReal.coe_ne_top _ hle

/-- The direct `EReal` form of convexity, for the record. Both sides are `≠ ⊥`,
so the sum on the right is unambiguous. -/
theorem conj_combo_le (f : QuaPol) {s₁ s₂ : Plane} {a b : ℝ}
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : a + b = 1)
    (h₁ : f.conj s₁ ≠ ⊤) (h₂ : f.conj s₂ ≠ ⊤) :
    f.conj (a • s₁ + b • s₂)
      ≤ ((a * (f.conj s₁).toReal + b * (f.conj s₂).toReal : ℝ) : EReal) :=
  conj_le_of_combo f ha hb hab
    (le_of_eq (EReal.coe_toReal h₁ (conj_ne_bot f s₁)).symm)
    (le_of_eq (EReal.coe_toReal h₂ (conj_ne_bot f s₂)).symm)

/-! ### Lower semicontinuity -/

/-- Each term of the supremum is continuous in `s`: it is affine where `f` is
finite, and the constant `⊥` where `f` is `⊤`. -/
lemma continuous_conj_term (f : QuaPol) (x : Plane) :
    Continuous (fun s : Plane => ((dot s x : ℝ) : EReal) - f.eval x) := by
  by_cases hx : f.eval x = ⊤
  · simp only [hx, EReal.sub_top]
    exact continuous_const
  · have hE : f.eval x = (((f.eval x).toReal : ℝ) : EReal) :=
      (EReal.coe_toReal hx (eval_ne_bot f x)).symm
    rw [hE]
    have hcont : Continuous (fun s : Plane => (dot s x - (f.eval x).toReal : ℝ)) := by
      unfold dot
      fun_prop
    simpa only [Function.comp_def, EReal.coe_sub] using
      continuous_coe_real_ereal.comp hcont

/-- **`f*` is lower semicontinuous**, being a pointwise supremum of continuous
functions. -/
theorem lowerSemicontinuous_conj (f : QuaPol) : LowerSemicontinuous f.conj := by
  have : LowerSemicontinuous (fun s : Plane =>
      ⨆ x : Plane, ((dot s x : ℝ) : EReal) - f.eval x) :=
    lowerSemicontinuous_iSup fun x => (continuous_conj_term f x).lowerSemicontinuous
  exact this

end QuaPol

end QuaConProof

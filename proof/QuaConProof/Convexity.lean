/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.FrankWolfe
import Mathlib.Topology.Semicontinuity.Basic
import Mathlib.Topology.Instances.EReal.Lemmas

/-!
# Suprema of affine functions: convexity, lower semicontinuity, convex domains

`TODO.md` items A1, A2, A3, and the shared half of C1. None of this is special to
a `QuaPol` — it is the standard fact that a supremum of affine functions is
convex and lower semicontinuous — but it has to exist in this development before
anything can be said about the biconjugate, and stating it for an `EReal`-valued
function needs a little care.

## One definition for both conjugates

`f*(s) = sup_x ⟨s,x⟩ - f(x)` and `f**(x) = sup_s ⟨s,x⟩ - f*(s)` are the same
shape: a pointwise supremum over `y` of `z ↦ ⟨y,z⟩ - g y`. `supAffine` is that
shape, and everything below is proved once for it, then specialised twice. The
symmetry of `dot` is what lets the two orientations share a definition.

## How convexity is stated, and why

The obvious statement,

    F(a·z₁ + b·z₂) ≤ a * F(z₁) + b * F(z₂),

is an inequality in `EReal`, which drags in `EReal` multiplication and its
`0 * ⊤ = 0` convention. The epigraph form avoids all of it: the bounds `t₁`, `t₂`
and the combination `a * t₁ + b * t₂` are **real numbers**, and the only `EReal`
reasoning left is that a point where `g = ⊤` contributes `⊥` to the supremum.

## Main results

* `supAffine_le_of_combo` — the epigraph inequality, and the workhorse.
* `convex_epigraph_supAffine`, `convex_dom_supAffine`,
  `lowerSemicontinuous_supAffine`.
* `QuaPol.convex_epigraph_conj`, `QuaPol.convex_dom_conj`,
  `QuaPol.lowerSemicontinuous_conj` — the specialisations to `f*`. The domain
  statement sharpens the fifth conjunct of `conj_isQuaCon`, which so far only
  knows the `⊤` region can be inhabited.
-/

namespace QuaConProof

/-! ### The shape shared by both conjugates -/

/-- The pointwise supremum of the affine family `z ↦ ⟨y,z⟩ - g y`, indexed by
`y`. Both `conj` and `biconj` are of this shape. -/
noncomputable def supAffine (g : Plane → EReal) (z : Plane) : EReal :=
  ⨆ y : Plane, ((dot y z : ℝ) : EReal) - g y

lemma term_le_supAffine (g : Plane → EReal) (y z : Plane) :
    ((dot y z : ℝ) : EReal) - g y ≤ supAffine g z :=
  le_iSup (fun w => ((dot w z : ℝ) : EReal) - g w) y

lemma dot_combo (a b : ℝ) (s₁ s₂ x : Plane) :
    dot (a • s₁ + b • s₂) x = a * dot s₁ x + b * dot s₂ x := by
  simp only [dot, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
  ring

lemma dot_combo_right (a b : ℝ) (y z₁ z₂ : Plane) :
    dot y (a • z₁ + b • z₂) = a * dot y z₁ + b * dot y z₂ := by
  rw [dot_comm, dot_combo, dot_comm y z₁, dot_comm y z₂]

/-- Anything that is not `⊤` is bounded by a real number — with `⊥` bounded by
`0`, which is what keeps the domain statement free of a `≠ ⊥` hypothesis. -/
lemma exists_real_bound {v : EReal} (h : v ≠ ⊤) : ∃ t : ℝ, v ≤ (t : EReal) := by
  rcases eq_or_ne v ⊥ with rfl | hb
  · exact ⟨0, bot_le⟩
  · exact ⟨v.toReal, le_of_eq (EReal.coe_toReal h hb).symm⟩

/-- **The epigraph inequality.** If `t₁` and `t₂` bound the supremum at `z₁` and
`z₂`, then `a·t₁ + b·t₂` bounds it at the convex combination.

This is convexity, in the form that keeps every arithmetic step in `ℝ`. -/
theorem supAffine_le_of_combo {g : Plane → EReal} (hg : ∀ y, g y ≠ ⊥)
    {z₁ z₂ : Plane} {a b t₁ t₂ : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : a + b = 1)
    (h₁ : supAffine g z₁ ≤ (t₁ : EReal)) (h₂ : supAffine g z₂ ≤ (t₂ : EReal)) :
    supAffine g (a • z₁ + b • z₂) ≤ ((a * t₁ + b * t₂ : ℝ) : EReal) := by
  refine iSup_le fun y => ?_
  by_cases hy : g y = ⊤
  · rw [hy, EReal.sub_top]
    exact bot_le
  · -- where `g` is finite the whole inequality is about real numbers
    set G : ℝ := (g y).toReal with hGdef
    have hG : g y = (G : EReal) := (EReal.coe_toReal hy (hg y)).symm
    have key : ∀ (z : Plane) (t : ℝ), supAffine g z ≤ (t : EReal) → dot y z - G ≤ t := by
      intro z t hzt
      have hle : ((dot y z : ℝ) : EReal) - ((G : ℝ) : EReal) ≤ (t : EReal) := by
        rw [← hG]
        exact le_trans (term_le_supAffine g y z) hzt
      rw [← EReal.coe_sub, EReal.coe_le_coe_iff] at hle
      exact hle
    have k₁ := key z₁ t₁ h₁
    have k₂ := key z₂ t₂ h₂
    have e1 : a * (dot y z₁ - G) ≤ a * t₁ := mul_le_mul_of_nonneg_left k₁ ha
    have e2 : b * (dot y z₂ - G) ≤ b * t₂ := mul_le_mul_of_nonneg_left k₂ hb
    have expand : a * (dot y z₁ - G) + b * (dot y z₂ - G)
        = a * dot y z₁ + b * dot y z₂ - G := by
      linear_combination (-G) * hab
    rw [hG, ← EReal.coe_sub, EReal.coe_le_coe_iff, dot_combo_right]
    linarith

/-- **The epigraph is convex.** -/
theorem convex_epigraph_supAffine {g : Plane → EReal} (hg : ∀ y, g y ≠ ⊥) :
    Convex ℝ {p : Plane × ℝ | supAffine g p.1 ≤ ((p.2 : ℝ) : EReal)} := by
  rintro ⟨z₁, t₁⟩ h₁ ⟨z₂, t₂⟩ h₂ a b ha hb hab
  exact supAffine_le_of_combo hg ha hb hab h₁ h₂

/-- **The domain is convex.** No `≠ ⊥` hypothesis is needed: a value of `⊥` is
bounded by `0` like any other. -/
theorem convex_dom_supAffine {g : Plane → EReal} (hg : ∀ y, g y ≠ ⊥) :
    Convex ℝ {z : Plane | supAffine g z ≠ ⊤} := by
  intro z₁ h₁ z₂ h₂ a b ha hb hab
  obtain ⟨t₁, ht₁⟩ := exists_real_bound h₁
  obtain ⟨t₂, ht₂⟩ := exists_real_bound h₂
  have hle := supAffine_le_of_combo hg ha hb hab ht₁ ht₂
  show supAffine g (a • z₁ + b • z₂) ≠ ⊤
  intro htop
  rw [htop, top_le_iff] at hle
  exact EReal.coe_ne_top _ hle

/-- Each term of the supremum is continuous in `z`: affine where `g` is finite,
and the constant `⊥` where `g` is `⊤`. -/
lemma continuous_supAffine_term {g : Plane → EReal} (hg : ∀ y, g y ≠ ⊥) (y : Plane) :
    Continuous (fun z : Plane => ((dot y z : ℝ) : EReal) - g y) := by
  by_cases hy : g y = ⊤
  · simp only [hy, EReal.sub_top]
    exact continuous_const
  · have hG : g y = (((g y).toReal : ℝ) : EReal) := (EReal.coe_toReal hy (hg y)).symm
    rw [hG]
    have hcont : Continuous (fun z : Plane => (dot y z - (g y).toReal : ℝ)) := by
      unfold dot
      fun_prop
    simpa only [Function.comp_def, EReal.coe_sub] using
      continuous_coe_real_ereal.comp hcont

/-- **A supremum of affine functions is lower semicontinuous.** -/
theorem lowerSemicontinuous_supAffine {g : Plane → EReal} (hg : ∀ y, g y ≠ ⊥) :
    LowerSemicontinuous (supAffine g) :=
  lowerSemicontinuous_iSup fun y => (continuous_supAffine_term hg y).lowerSemicontinuous

namespace QuaPol

/-! ### The conjugate -/

/-- `f*` is a supremum of affine functions in the sense of `supAffine`. The two
orientations of `dot` agree because `dot` is symmetric. -/
lemma conj_eq_supAffine (f : QuaPol) : f.conj = supAffine f.eval := by
  funext s
  rw [conj_def, supAffine]
  exact iSup_congr fun x => by rw [dot_comm]

/-- A single term of the supremum defining `f*`, bounded at one point. -/
lemma term_le_conj (f : QuaPol) (s x : Plane) :
    ((dot s x : ℝ) : EReal) - f.eval x ≤ f.conj s :=
  le_iSup (fun y => ((dot s y : ℝ) : EReal) - f.eval y) x

/-- **The epigraph inequality for `f*`.** -/
theorem conj_le_of_combo (f : QuaPol) {s₁ s₂ : Plane} {a b t₁ t₂ : ℝ}
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : a + b = 1)
    (h₁ : f.conj s₁ ≤ (t₁ : EReal)) (h₂ : f.conj s₂ ≤ (t₂ : EReal)) :
    f.conj (a • s₁ + b • s₂) ≤ ((a * t₁ + b * t₂ : ℝ) : EReal) := by
  rw [conj_eq_supAffine] at h₁ h₂ ⊢
  exact supAffine_le_of_combo (eval_ne_bot f) ha hb hab h₁ h₂

/-- **The epigraph of `f*` is convex.** -/
theorem convex_epigraph_conj (f : QuaPol) :
    Convex ℝ {p : Plane × ℝ | f.conj p.1 ≤ ((p.2 : ℝ) : EReal)} := by
  simpa only [conj_eq_supAffine] using convex_epigraph_supAffine (eval_ne_bot f)

/-- **The domain of `f*` is convex**, so the region where `f*` is `⊤` is the
complement of a convex set.

`conj_isQuaCon`'s fifth conjunct identifies that region with the empty-activity
cell, and `Sanity.cell_empty_rayPol_nonempty` shows it can be inhabited; this
says what shape it has. -/
theorem convex_dom_conj (f : QuaPol) : Convex ℝ {s : Plane | f.conj s ≠ ⊤} := by
  simpa only [conj_eq_supAffine] using convex_dom_supAffine (eval_ne_bot f)

/-- **`f*` is lower semicontinuous.** -/
theorem lowerSemicontinuous_conj (f : QuaPol) : LowerSemicontinuous f.conj := by
  simpa only [conj_eq_supAffine] using lowerSemicontinuous_supAffine (eval_ne_bot f)

/-- The direct `EReal` form of convexity, for the record. Both values are `≠ ⊥`,
so the sum on the right is unambiguous. -/
theorem conj_combo_le (f : QuaPol) {s₁ s₂ : Plane} {a b : ℝ}
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : a + b = 1)
    (h₁ : f.conj s₁ ≠ ⊤) (h₂ : f.conj s₂ ≠ ⊤) :
    f.conj (a • s₁ + b • s₂)
      ≤ ((a * (f.conj s₁).toReal + b * (f.conj s₂).toReal : ℝ) : EReal) :=
  conj_le_of_combo f ha hb hab
    (le_of_eq (EReal.coe_toReal h₁ (conj_ne_bot f s₁)).symm)
    (le_of_eq (EReal.coe_toReal h₂ (conj_ne_bot f s₂)).symm)

end QuaPol

end QuaConProof

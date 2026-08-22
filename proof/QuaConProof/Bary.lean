/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.QuaPol
import Mathlib.Analysis.Convex.Combination

/-!
# Barycentric bookkeeping for the convex hull of a finite set

The three facts about `convexHull ℝ ↑W` that the selection lemma needs, and
nothing more. `PROJECT_PLAN.md` §0.6, steps S3, S5 and the descent in S8.

* `mem_convexHull_erase` — a point whose weight at `v₀` vanishes already lies in
  the hull of `W.erase v₀`. This is what makes the induction on `W.card` work,
  and it is why "all weights positive" may be assumed in the main case.
* `mem_convexHull_perturb` — with all weights positive one may step a little way
  along `w - v` in **either** direction and stay inside. This is what turns
  maximality into a two-sided inequality, and hence into the first-order
  condition.
* `exists_descent` — given a direction written with weights summing to zero, one
  may travel along it until some weight first vanishes, landing in the hull of a
  strictly smaller set. This is the bookkeeping half of S8.

Everything here is about convexity alone; no quadratic appears.
-/

namespace QuaConProof

open scoped Classical
open Finset

/-! ### Passing between a point of the hull and its weights -/

/-- Explicit convex weights for a point of the hull of a finite set. -/
lemma exists_weights {W : Finset Plane} {x : Plane}
    (hx : x ∈ convexHull ℝ (↑W : Set Plane)) :
    ∃ lam : Plane → ℝ, (∀ v ∈ W, 0 ≤ lam v) ∧ (∑ v ∈ W, lam v) = 1 ∧
      (∑ v ∈ W, lam v • v) = x := by
  rw [Finset.convexHull_eq] at hx
  obtain ⟨lam, h0, h1, h2⟩ := hx
  refine ⟨lam, h0, h1, ?_⟩
  rwa [Finset.centerMass_eq_of_sum_1 _ _ h1] at h2

/-- Conversely, any convex combination lies in the hull. -/
lemma mem_convexHull_of_weights {W : Finset Plane} {x : Plane} {lam : Plane → ℝ}
    (h0 : ∀ v ∈ W, 0 ≤ lam v) (h1 : (∑ v ∈ W, lam v) = 1)
    (h2 : (∑ v ∈ W, lam v • v) = x) :
    x ∈ convexHull ℝ (↑W : Set Plane) := by
  rw [Finset.convexHull_eq]
  refine ⟨lam, h0, h1, ?_⟩
  rw [Finset.centerMass_eq_of_sum_1 _ _ h1]
  exact h2

/-! ### A vanishing weight drops a vertex -/

/-- **If a weight vanishes, the point already lies in the hull of the smaller
set.** This is the step that lets the selection lemma induct on `W.card`. -/
lemma mem_convexHull_erase {W : Finset Plane} {x : Plane} {lam : Plane → ℝ}
    (h0 : ∀ v ∈ W, 0 ≤ lam v) (h1 : (∑ v ∈ W, lam v) = 1)
    (h2 : (∑ v ∈ W, lam v • v) = x) {v₀ : Plane} (hzero : lam v₀ = 0) :
    x ∈ convexHull ℝ (↑(W.erase v₀) : Set Plane) := by
  refine mem_convexHull_of_weights (fun v hv => h0 v (Finset.mem_of_mem_erase hv)) ?_ ?_
  · rw [Finset.sum_erase _ hzero]; exact h1
  · rw [Finset.sum_erase _ (by rw [hzero, zero_smul])]; exact h2

/-! ### Stepping along `w - v` -/

/-- **Two-sided room to move.** If every weight is positive, one may step from
`x` along `w - v` by any `t` with `-lam w ≤ t ≤ lam v` and stay in the hull.

Both signs are allowed, which is exactly what upgrades "`x` is a maximum" to a
*first-order condition* rather than a mere inequality. -/
lemma mem_convexHull_perturb {W : Finset Plane} {x : Plane} {lam : Plane → ℝ}
    (h0 : ∀ u ∈ W, 0 ≤ lam u) (h1 : (∑ u ∈ W, lam u) = 1)
    (h2 : (∑ u ∈ W, lam u • u) = x)
    {v w : Plane} (hv : v ∈ W) (hw : w ∈ W) (hvw : v ≠ w)
    {t : ℝ} (htv : t ≤ lam v) (htw : -t ≤ lam w) :
    x + t • (w - v) ∈ convexHull ℝ (↑W : Set Plane) := by
  classical
  have hsum_w : (∑ u ∈ W, (if u = w then (1 : ℝ) else 0)) = 1 := by
    simp [Finset.sum_ite_eq', hw]
  have hsum_v : (∑ u ∈ W, (if u = v then (1 : ℝ) else 0)) = 1 := by
    simp [Finset.sum_ite_eq', hv]
  have hsmul_w : (∑ u ∈ W, (if u = w then (1 : ℝ) else 0) • u) = w := by
    simp [ite_smul, Finset.sum_ite_eq', hw]
  have hsmul_v : (∑ u ∈ W, (if u = v then (1 : ℝ) else 0) • u) = v := by
    simp [ite_smul, Finset.sum_ite_eq', hv]
  refine mem_convexHull_of_weights
    (lam := fun u => lam u + t * ((if u = w then (1 : ℝ) else 0) - (if u = v then (1 : ℝ) else 0)))
    ?_ ?_ ?_
  · intro u hu
    by_cases huw : u = w
    · subst huw
      rw [if_pos rfl, if_neg (Ne.symm hvw)]
      linarith [htw]
    · by_cases huv : u = v
      · subst huv
        rw [if_neg huw, if_pos rfl]
        linarith [htv]
      · rw [if_neg huw, if_neg huv]
        simpa using h0 u hu
  · rw [Finset.sum_add_distrib, h1, ← Finset.mul_sum, Finset.sum_sub_distrib,
      hsum_w, hsum_v]
    ring
  · have hexp : ∀ u : Plane,
        (lam u + t * ((if u = w then (1 : ℝ) else 0) - (if u = v then (1 : ℝ) else 0))) • u
        = lam u • u + t • ((if u = w then (1 : ℝ) else 0) • u)
          - t • ((if u = v then (1 : ℝ) else 0) • u) := by
      intro u; module
    simp only [hexp]
    rw [Finset.sum_sub_distrib, Finset.sum_add_distrib, ← Finset.smul_sum, ← Finset.smul_sum,
      h2, hsmul_w, hsmul_v]
    module

/-! ### Travelling until a weight first vanishes -/

/-- **The descent.** Suppose every weight of `x` over `W` is positive, and `d` is
written as `∑ gam v • v` with the `gam` summing to zero and not all zero. Then
there is a step `t` and a vertex `v₀` with `x + t·d` in the hull of
`W.erase v₀`.

The weights are affine in `t`, so one simply travels until the first of them
reaches zero. This is the bookkeeping half of step S8: combined with
`psi_const_along_kernel`, which says travelling along `ker H` costs nothing, it
moves a maximiser onto a proper face. -/
lemma exists_descent {W : Finset Plane} {x d : Plane} {lam gam : Plane → ℝ}
    (hpos : ∀ v ∈ W, 0 < lam v) (h1 : (∑ v ∈ W, lam v) = 1)
    (h2 : (∑ v ∈ W, lam v • v) = x)
    (hg0 : (∑ v ∈ W, gam v) = 0) (hgd : (∑ v ∈ W, gam v • v) = d)
    (hne : ∃ v ∈ W, gam v ≠ 0) :
    ∃ (t : ℝ) (v₀ : Plane), v₀ ∈ W ∧
      x + t • d ∈ convexHull ℝ (↑(W.erase v₀) : Set Plane) := by
  classical
  -- some weight of the direction is strictly negative
  set N := W.filter (fun v => gam v < 0) with hN
  have hNne : N.Nonempty := by
    by_contra hc
    rw [Finset.not_nonempty_iff_eq_empty, hN, Finset.filter_eq_empty_iff] at hc
    have hall : ∀ v ∈ W, 0 ≤ gam v := fun v hv => not_lt.1 (hc hv)
    obtain ⟨v, hv, hvne⟩ := hne
    exact hvne ((Finset.sum_eq_zero_iff_of_nonneg hall).1 hg0 v hv)
  -- travel until the first weight hits zero
  obtain ⟨v₀, hv₀N, hv₀min⟩ :=
    Finset.exists_min_image N (fun v => lam v / (-gam v)) hNne
  have hv₀W : v₀ ∈ W := (Finset.mem_filter.1 hv₀N).1
  have hv₀neg : gam v₀ < 0 := (Finset.mem_filter.1 hv₀N).2
  set t : ℝ := lam v₀ / (-gam v₀) with ht
  set mu : Plane → ℝ := fun u => lam u + t * gam u with hmu
  have hgne : gam v₀ ≠ 0 := ne_of_lt hv₀neg
  have hmuv₀ : mu v₀ = 0 := by
    have hneg : (-gam v₀) ≠ 0 := neg_ne_zero.2 hgne
    have hrw : mu v₀ = lam v₀ + lam v₀ / (-gam v₀) * gam v₀ := rfl
    have hkey : lam v₀ / (-gam v₀) * gam v₀ = -lam v₀ := by
      rw [div_mul_eq_mul_div, div_eq_iff hneg]; ring
    rw [hrw, hkey]; ring
  have hmu0 : ∀ u ∈ W, 0 ≤ mu u := by
    intro u hu
    simp only [hmu]
    by_cases hgu : gam u < 0
    · have huN : u ∈ N := Finset.mem_filter.2 ⟨hu, hgu⟩
      have hle : t ≤ lam u / (-gam u) := hv₀min u huN
      have hneg : 0 < -gam u := by linarith
      rw [le_div_iff₀ hneg] at hle
      linarith
    · have ht0 : 0 ≤ t := by
        rw [ht]
        exact div_nonneg (hpos v₀ hv₀W).le (by linarith)
      have : 0 ≤ t * gam u := mul_nonneg ht0 (not_lt.1 hgu)
      linarith [hpos u hu]
  have hmusum : (∑ u ∈ W, mu u) = 1 := by
    simp only [hmu]
    rw [Finset.sum_add_distrib, h1, ← Finset.mul_sum, hg0, mul_zero, add_zero]
  have hmucombo : (∑ u ∈ W, mu u • u) = x + t • d := by
    simp only [hmu]
    have : ∀ u : Plane, (lam u + t * gam u) • u = lam u • u + t • (gam u • u) := by
      intro u; module
    simp only [this]
    rw [Finset.sum_add_distrib, ← Finset.smul_sum, h2, hgd]
  exact ⟨t, v₀, hv₀W, mem_convexHull_erase hmu0 hmusum hmucombo hmuv₀⟩


/-! ### The plane: parallelism and perpendicularity

Two-dimensional linear algebra, done by hand with the scalar cross product. This
is cheaper than invoking `AffineIndependent` and `Collinear`, and it is all the
selection lemma needs: whether a vertex set is contained in a line, and what
follows when it is not. -/

/-- The scalar cross product `a₁b₂ - a₂b₁`. It vanishes exactly on parallel
pairs, and it is the determinant deciding whether two directions span. -/
def cross (a b : Plane) : ℝ := a.1 * b.2 - a.2 * b.1

lemma cross_self (a : Plane) : cross a a = 0 := by simp only [cross]; ring

/-- **Parallel vectors.** If `cross d e = 0` and `d ≠ 0` then `e` is a multiple
of `d`. -/
lemma exists_smul_of_cross_eq_zero {d e : Plane} (hd : d ≠ 0) (h : cross d e = 0) :
    ∃ t : ℝ, e = t • d := by
  simp only [cross, sub_eq_zero] at h
  by_cases h1 : d.1 = 0
  · have h2 : d.2 ≠ 0 := by
      intro h2; exact hd (Prod.ext_iff.2 ⟨by simpa using h1, by simpa using h2⟩)
    refine ⟨e.2 / d.2, Prod.ext_iff.2 ⟨?_, ?_⟩⟩
    · simp only [Prod.smul_fst, smul_eq_mul, h1, mul_zero]
      rw [h1] at h
      have hz : d.2 * e.1 = 0 := by linarith
      rcases mul_eq_zero.1 hz with h' | h'
      · exact absurd h' h2
      · exact h'
    · simp only [Prod.smul_snd, smul_eq_mul]; field_simp
  · refine ⟨e.1 / d.1, Prod.ext_iff.2 ⟨?_, ?_⟩⟩
    · simp only [Prod.smul_fst, smul_eq_mul]; field_simp
    · simp only [Prod.smul_snd, smul_eq_mul]
      field_simp
      linarith [h]

/-- **Two independent constraints kill a vector.** If `r` is perpendicular to
both `d` and `e`, and `d`, `e` are not parallel, then `r = 0`. This is what turns
the first-order conditions at a two-dimensional face into `∇q(x) = s`. -/
lemma eq_zero_of_dot_eq_zero_of_cross_ne_zero {r d e : Plane}
    (hd : dot r d = 0) (he : dot r e = 0) (h : cross d e ≠ 0) : r = 0 := by
  simp only [dot] at hd he
  refine Prod.ext_iff.2 ⟨?_, ?_⟩
  · simp only [Prod.fst_zero]
    have hkey : r.1 * cross d e = 0 := by
      simp only [cross]; linear_combination e.2 * hd - d.2 * he
    exact (mul_eq_zero.1 hkey).resolve_right h
  · simp only [Prod.snd_zero]
    have hkey : r.2 * cross d e = 0 := by
      simp only [cross]; linear_combination -e.1 * hd + d.1 * he
    exact (mul_eq_zero.1 hkey).resolve_right h

/-- **Spanning.** If `d` and `e` are not parallel, every vector of the plane is a
combination of them, explicitly by Cramer's rule. -/
lemma exists_combo_of_cross_ne_zero {d e : Plane} (h : cross d e ≠ 0) (z : Plane) :
    ∃ b₁ b₂ : ℝ, z = b₁ • d + b₂ • e := by
  refine ⟨cross z e / cross d e, cross d z / cross d e, ?_⟩
  refine Prod.ext_iff.2 ⟨?_, ?_⟩ <;>
  · simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
    rw [div_mul_eq_mul_div, div_mul_eq_mul_div, ← add_div, eq_div_iff h]
    simp only [cross]; ring

/-! ### Collinearity of the weights -/

/-- Components of a weighted sum. -/
lemma fst_weighted_sum (W : Finset Plane) (lam : Plane → ℝ) :
    (∑ u ∈ W, lam u • u).1 = ∑ u ∈ W, lam u * u.1 := by
  classical
  induction W using Finset.induction_on with
  | empty => simp
  | insert a t ha ih => rw [Finset.sum_insert ha, Finset.sum_insert ha]; simp [ih]

lemma snd_weighted_sum (W : Finset Plane) (lam : Plane → ℝ) :
    (∑ u ∈ W, lam u • u).2 = ∑ u ∈ W, lam u * u.2 := by
  classical
  induction W using Finset.induction_on with
  | empty => simp
  | insert a t ha ih => rw [Finset.sum_insert ha, Finset.sum_insert ha]; simp [ih]

/-- **A convex combination of points parallel to `d` is itself parallel to `d`.**

If every vertex of `W` lies on the line through `v` with direction `D`, so does
every point of the hull — in particular `x`. -/
lemma cross_eq_zero_of_forall {W : Finset Plane} {x v D : Plane} {lam : Plane → ℝ}
    (h1 : (∑ u ∈ W, lam u) = 1) (h2 : (∑ u ∈ W, lam u • u) = x)
    (hall : ∀ u ∈ W, cross D (u - v) = 0) : cross D (x - v) = 0 := by
  have hx1 : x.1 = ∑ u ∈ W, lam u * u.1 := by rw [← h2, fst_weighted_sum]
  have hx2 : x.2 = ∑ u ∈ W, lam u * u.2 := by rw [← h2, snd_weighted_sum]
  have hv1 : v.1 = ∑ u ∈ W, lam u * v.1 := by rw [← Finset.sum_mul, h1, one_mul]
  have hv2 : v.2 = ∑ u ∈ W, lam u * v.2 := by rw [← Finset.sum_mul, h1, one_mul]
  simp only [cross, Prod.fst_sub, Prod.snd_sub]
  rw [hx1, hx2]
  nth_rewrite 1 [hv2]
  nth_rewrite 1 [hv1]
  rw [← Finset.sum_sub_distrib, ← Finset.sum_sub_distrib, Finset.mul_sum, Finset.mul_sum,
    ← Finset.sum_sub_distrib]
  refine Finset.sum_eq_zero fun u hu => ?_
  have hu0 := hall u hu
  simp only [cross, Prod.fst_sub, Prod.snd_sub] at hu0
  linear_combination lam u * hu0

end QuaConProof

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

/-! ### Directions inside a piece: weight representations

Phase 7. A point of a piece is `∑ λ_v v + ∑ μ_r r` over a vertex support `W` and
a ray support `S`. A *direction* along which the point may be moved is recorded
the same way, by how it changes those weights: the vertex weights must change by
a vector summing to zero, since they are normalised; the ray weights may change
freely, since they are not.

Recording directions this way is what lets one induction cover both kinds. A
difference of two vertices and a recession direction differ only in which of the
two weight families they touch, and every lemma below is stated for the record
rather than for either kind. -/

/-- `gam` and `nu` move the weights of a point of `conv W + cone S` along `d`. -/
def IsDirRep (W S : Finset Plane) (gam nu : Plane → ℝ) (d : Plane) : Prop :=
  (∑ v ∈ W, gam v) = 0 ∧ (∑ v ∈ W, gam v • v) + (∑ r ∈ S, nu r • r) = d

namespace IsDirRep

variable {W S : Finset Plane} {gam nu gam' nu' : Plane → ℝ} {d d' : Plane}

lemma add (h : IsDirRep W S gam nu d) (h' : IsDirRep W S gam' nu' d') :
    IsDirRep W S (fun z => gam z + gam' z) (fun z => nu z + nu' z) (d + d') := by
  refine ⟨?_, ?_⟩
  · rw [Finset.sum_add_distrib, h.1, h'.1, add_zero]
  · have hv : (∑ v ∈ W, (gam v + gam' v) • v)
        = (∑ v ∈ W, gam v • v) + (∑ v ∈ W, gam' v • v) := by
      rw [← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun v _ => by module
    have hr : (∑ r ∈ S, (nu r + nu' r) • r)
        = (∑ r ∈ S, nu r • r) + (∑ r ∈ S, nu' r • r) := by
      rw [← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun r _ => by module
    rw [hv, hr, ← h.2, ← h'.2]
    abel

lemma smul (h : IsDirRep W S gam nu d) (c : ℝ) :
    IsDirRep W S (fun z => c * gam z) (fun z => c * nu z) (c • d) := by
  refine ⟨?_, ?_⟩
  · rw [← Finset.mul_sum, h.1, mul_zero]
  · have hv : (∑ v ∈ W, (c * gam v) • v) = c • (∑ v ∈ W, gam v • v) := by
      rw [Finset.smul_sum]
      exact Finset.sum_congr rfl fun v _ => by rw [mul_smul]
    have hr : (∑ r ∈ S, (c * nu r) • r) = c • (∑ r ∈ S, nu r • r) := by
      rw [Finset.smul_sum]
      exact Finset.sum_congr rfl fun r _ => by rw [mul_smul]
    rw [hv, hr, ← h.2, smul_add]

end IsDirRep

/-- A difference of two vertices of the support is a direction. -/
lemma isDirRep_vert_sub {W S : Finset Plane} {v w : Plane} (hv : v ∈ W) (hw : w ∈ W) :
    IsDirRep W S (fun z => (if z = w then (1 : ℝ) else 0) - (if z = v then (1 : ℝ) else 0))
      (fun _ => 0) (w - v) := by
  classical
  constructor
  · rw [Finset.sum_sub_distrib]
    simp [Finset.sum_ite_eq', hv, hw]
  · have hexp : ∀ z : Plane,
        ((if z = w then (1 : ℝ) else 0) - (if z = v then (1 : ℝ) else 0)) • z
        = (if z = w then (1 : ℝ) else 0) • z - (if z = v then (1 : ℝ) else 0) • z := by
      intro z; module
    simp only [hexp]
    rw [Finset.sum_sub_distrib]
    simp [ite_smul, Finset.sum_ite_eq', hv, hw]

/-- A recession direction of the support is a direction. -/
lemma isDirRep_ray {W S : Finset Plane} {r : Plane} (hr : r ∈ S) :
    IsDirRep W S (fun _ => 0) (fun z => if z = r then (1 : ℝ) else 0) r := by
  classical
  refine ⟨by simp, ?_⟩
  simp [ite_smul, Finset.sum_ite_eq', hr]

/-! ### Membership of a piece, from weights -/

/-- A nonnegative normalised combination of a vertex subset, plus a nonnegative
combination of a ray subset, lies in the piece. -/
lemma mem_T_of_weights {p : QuaPiece} {W S : Finset Plane} {lam mu : Plane → ℝ} {x : Plane}
    (hWV : W ⊆ p.verts) (hSR : S ⊆ p.rays)
    (h0 : ∀ v ∈ W, 0 ≤ lam v) (h1 : (∑ v ∈ W, lam v) = 1) (hm0 : ∀ r ∈ S, 0 ≤ mu r)
    (hx : (∑ v ∈ W, lam v • v) + (∑ r ∈ S, mu r • r) = x) : x ∈ p.T := by
  classical
  refine ⟨∑ v ∈ W, lam v • v, ?_, ∑ r ∈ S, mu r • r, ?_, hx⟩
  · exact convexHull_mono (Finset.coe_subset.2 hWV) (mem_convexHull_of_weights h0 h1 rfl)
  · refine ⟨fun z => if z ∈ S then mu z else 0, ?_, ?_⟩
    · intro z _
      by_cases hz : z ∈ S
      · simpa [hz] using hm0 z hz
      · simp [hz]
    · rw [← Finset.sum_subset hSR (fun z _ hz => by simp [hz])]
      exact Finset.sum_congr rfl fun z hz => by simp [hz]

/-- **Stepping along a direction.** If the perturbed weights stay nonnegative,
the perturbed point stays in the piece. -/
lemma mem_T_of_perturb {p : QuaPiece} {W S : Finset Plane} {lam mu gam nu : Plane → ℝ}
    {x d : Plane} {t : ℝ}
    (hWV : W ⊆ p.verts) (hSR : S ⊆ p.rays)
    (h1 : (∑ v ∈ W, lam v) = 1)
    (hx : (∑ v ∈ W, lam v • v) + (∑ r ∈ S, mu r • r) = x)
    (hrep : IsDirRep W S gam nu d)
    (hl : ∀ v ∈ W, 0 ≤ lam v + t * gam v) (hm : ∀ r ∈ S, 0 ≤ mu r + t * nu r) :
    x + t • d ∈ p.T := by
  refine mem_T_of_weights hWV hSR hl ?_ hm ?_
  · rw [Finset.sum_add_distrib, h1, ← Finset.mul_sum, hrep.1, mul_zero, add_zero]
  · have hv : (∑ v ∈ W, (lam v + t * gam v) • v)
        = (∑ v ∈ W, lam v • v) + t • (∑ v ∈ W, gam v • v) := by
      rw [Finset.smul_sum, ← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun v _ => by rw [← mul_smul]; module
    have hr : (∑ r ∈ S, (mu r + t * nu r) • r)
        = (∑ r ∈ S, mu r • r) + t • (∑ r ∈ S, nu r • r) := by
      rw [Finset.smul_sum, ← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun r _ => by rw [← mul_smul]; module
    rw [hv, hr, ← hx, ← hrep.2, smul_add]
    abel

/-- **Two-sided room to move.** With every weight strictly positive there is a
neighbourhood of `0` in which the step may be taken in either direction. This is
what upgrades "`x` is a maximum" to a first-order condition. -/
lemma exists_perturb_radius {W S : Finset Plane} {lam mu : Plane → ℝ} (gam nu : Plane → ℝ)
    (hlam : ∀ v ∈ W, 0 < lam v) (hmu : ∀ r ∈ S, 0 < mu r) :
    ∃ ε > (0 : ℝ), ∀ t : ℝ, |t| < ε →
      (∀ v ∈ W, 0 ≤ lam v + t * gam v) ∧ (∀ r ∈ S, 0 ≤ mu r + t * nu r) := by
  classical
  set A : Finset ℝ :=
    insert (1 : ℝ) ((W.image fun v => lam v / (|gam v| + 1))
      ∪ (S.image fun r => mu r / (|nu r| + 1))) with hA
  have hAne : A.Nonempty := ⟨1, by simp [hA]⟩
  set ε : ℝ := A.min' hAne with hε
  have hmem : ε ∈ A := A.min'_mem hAne
  have hpos : 0 < ε := by
    rw [hA] at hmem
    simp only [Finset.mem_insert, Finset.mem_union, Finset.mem_image] at hmem
    rcases hmem with h | h | h
    · rw [h]; norm_num
    · obtain ⟨v, hv, hveq⟩ := h
      rw [← hveq]; exact div_pos (hlam v hv) (by positivity)
    · obtain ⟨r, hr, hreq⟩ := h
      rw [← hreq]; exact div_pos (hmu r hr) (by positivity)
  refine ⟨ε, hpos, fun t ht => ⟨?_, ?_⟩⟩
  · intro v hv
    have hin : lam v / (|gam v| + 1) ∈ A := by
      rw [hA]
      simp only [Finset.mem_insert, Finset.mem_union, Finset.mem_image]
      exact Or.inr (Or.inl ⟨v, hv, rfl⟩)
    have hle : ε ≤ lam v / (|gam v| + 1) := A.min'_le _ hin
    have hden : (0 : ℝ) < |gam v| + 1 := by positivity
    rw [le_div_iff₀ hden] at hle
    have h1 : |t * gam v| ≤ |t| * (|gam v| + 1) := by
      rw [abs_mul]
      exact mul_le_mul_of_nonneg_left (by linarith [le_abs_self (gam v)]) (abs_nonneg t)
    have h2 : |t| * (|gam v| + 1) < ε * (|gam v| + 1) := mul_lt_mul_of_pos_right ht hden
    have h3 : |t * gam v| < lam v := by linarith
    have h4 : -(lam v) < t * gam v := by
      have := neg_abs_le (t * gam v); linarith
    linarith
  · intro r hr
    have hin : mu r / (|nu r| + 1) ∈ A := by
      rw [hA]
      simp only [Finset.mem_insert, Finset.mem_union, Finset.mem_image]
      exact Or.inr (Or.inr ⟨r, hr, rfl⟩)
    have hle : ε ≤ mu r / (|nu r| + 1) := A.min'_le _ hin
    have hden : (0 : ℝ) < |nu r| + 1 := by positivity
    rw [le_div_iff₀ hden] at hle
    have h1 : |t * nu r| ≤ |t| * (|nu r| + 1) := by
      rw [abs_mul]
      exact mul_le_mul_of_nonneg_left (by linarith [le_abs_self (nu r)]) (abs_nonneg t)
    have h2 : |t| * (|nu r| + 1) < ε * (|nu r| + 1) := mul_lt_mul_of_pos_right ht hden
    have h3 : |t * nu r| < mu r := by linarith
    have h4 : -(mu r) < t * nu r := by
      have := neg_abs_le (t * nu r); linarith
    linarith

/-! ### Travelling until a weight first vanishes, with rays -/

/-- **Forward descent.** If some weight is strictly decreasing along the
direction, travel to the first zero. -/
lemma exists_descent_fwd {W S : Finset Plane} {lam mu gam nu : Plane → ℝ}
    (hlam : ∀ v ∈ W, 0 < lam v) (hmu : ∀ r ∈ S, 0 < mu r)
    (hneg : (∃ v ∈ W, gam v < 0) ∨ (∃ r ∈ S, nu r < 0)) :
    ∃ t : ℝ, 0 < t ∧ (∀ v ∈ W, 0 ≤ lam v + t * gam v) ∧ (∀ r ∈ S, 0 ≤ mu r + t * nu r) ∧
      ((∃ v₀ ∈ W, lam v₀ + t * gam v₀ = 0) ∨ (∃ r₀ ∈ S, mu r₀ + t * nu r₀ = 0)) := by
  classical
  set A : Finset ℝ :=
    ((W.filter fun v => gam v < 0).image fun v => lam v / (-gam v))
      ∪ ((S.filter fun r => nu r < 0).image fun r => mu r / (-nu r)) with hA
  have hAne : A.Nonempty := by
    rcases hneg with ⟨v, hv, hvneg⟩ | ⟨r, hr, hrneg⟩
    · refine ⟨lam v / (-gam v), ?_⟩
      rw [hA]
      refine Finset.mem_union_left _ (Finset.mem_image.2 ⟨v, ?_, rfl⟩)
      exact Finset.mem_filter.2 ⟨hv, hvneg⟩
    · refine ⟨mu r / (-nu r), ?_⟩
      rw [hA]
      refine Finset.mem_union_right _ (Finset.mem_image.2 ⟨r, ?_, rfl⟩)
      exact Finset.mem_filter.2 ⟨hr, hrneg⟩
  set t : ℝ := A.min' hAne with ht
  have hmem : t ∈ A := A.min'_mem hAne
  have htpos : 0 < t := by
    rw [hA] at hmem
    simp only [Finset.mem_union, Finset.mem_image, Finset.mem_filter] at hmem
    rcases hmem with ⟨v, ⟨hv, hvneg⟩, hveq⟩ | ⟨r, ⟨hr, hrneg⟩, hreq⟩
    · rw [← hveq]; exact div_pos (hlam v hv) (by linarith)
    · rw [← hreq]; exact div_pos (hmu r hr) (by linarith)
  refine ⟨t, htpos, ?_, ?_, ?_⟩
  · intro v hv
    by_cases hvneg : gam v < 0
    · have hin : lam v / (-gam v) ∈ A := by
        rw [hA]
        refine Finset.mem_union_left _ (Finset.mem_image.2 ⟨v, ?_, rfl⟩)
        exact Finset.mem_filter.2 ⟨hv, hvneg⟩
      have hle : t ≤ lam v / (-gam v) := A.min'_le _ hin
      rw [le_div_iff₀ (by linarith)] at hle
      linarith
    · have : 0 ≤ t * gam v := mul_nonneg htpos.le (not_lt.1 hvneg)
      linarith [hlam v hv]
  · intro r hr
    by_cases hrneg : nu r < 0
    · have hin : mu r / (-nu r) ∈ A := by
        rw [hA]
        refine Finset.mem_union_right _ (Finset.mem_image.2 ⟨r, ?_, rfl⟩)
        exact Finset.mem_filter.2 ⟨hr, hrneg⟩
      have hle : t ≤ mu r / (-nu r) := A.min'_le _ hin
      rw [le_div_iff₀ (by linarith)] at hle
      linarith
    · have : 0 ≤ t * nu r := mul_nonneg htpos.le (not_lt.1 hrneg)
      linarith [hmu r hr]
  · rw [hA] at hmem
    simp only [Finset.mem_union, Finset.mem_image, Finset.mem_filter] at hmem
    rcases hmem with ⟨v, ⟨hv, hvneg⟩, hveq⟩ | ⟨r, ⟨hr, hrneg⟩, hreq⟩
    · refine Or.inl ⟨v, hv, ?_⟩
      have hne : (-gam v) ≠ 0 := by linarith
      have hkey : lam v / (-gam v) * gam v = -lam v := by
        rw [div_mul_eq_mul_div, div_eq_iff hne]; ring
      rw [← hveq, hkey]; ring
    · refine Or.inr ⟨r, hr, ?_⟩
      have hne : (-nu r) ≠ 0 := by linarith
      have hkey : mu r / (-nu r) * nu r = -mu r := by
        rw [div_mul_eq_mul_div, div_eq_iff hne]; ring
      rw [← hreq, hkey]; ring

/-- **The descent.** A nonzero weight change moves the point onto a proper face,
in one direction or the other. Which direction is free, and that is exactly what
the singular case of the selection lemma needs: travelling along the kernel of
the Hessian costs nothing either way. -/
lemma exists_descent_gen {W S : Finset Plane} {lam mu gam nu : Plane → ℝ}
    (hlam : ∀ v ∈ W, 0 < lam v) (hmu : ∀ r ∈ S, 0 < mu r)
    (hne : (∃ v ∈ W, gam v ≠ 0) ∨ (∃ r ∈ S, nu r ≠ 0)) :
    ∃ t : ℝ, (∀ v ∈ W, 0 ≤ lam v + t * gam v) ∧ (∀ r ∈ S, 0 ≤ mu r + t * nu r) ∧
      ((∃ v₀ ∈ W, lam v₀ + t * gam v₀ = 0) ∨ (∃ r₀ ∈ S, mu r₀ + t * nu r₀ = 0)) := by
  by_cases hneg : (∃ v ∈ W, gam v < 0) ∨ (∃ r ∈ S, nu r < 0)
  · obtain ⟨t, _, h1, h2, h3⟩ := exists_descent_fwd hlam hmu hneg
    exact ⟨t, h1, h2, h3⟩
  · push Not at hneg
    have hneg' : (∃ v ∈ W, (-gam v) < 0) ∨ (∃ r ∈ S, (-nu r) < 0) := by
      rcases hne with ⟨v, hv, hvne⟩ | ⟨r, hr, hrne⟩
      · exact Or.inl ⟨v, hv, by
          have := hneg.1 v hv; rcases lt_trichotomy (gam v) 0 with h | h | h
          · exact absurd h (not_lt.2 this)
          · exact absurd h hvne
          · linarith⟩
      · exact Or.inr ⟨r, hr, by
          have := hneg.2 r hr; rcases lt_trichotomy (nu r) 0 with h | h | h
          · exact absurd h (not_lt.2 this)
          · exact absurd h hrne
          · linarith⟩
    obtain ⟨t, _, h1, h2, h3⟩ :=
      exists_descent_fwd (gam := fun z => -gam z) (nu := fun z => -nu z) hlam hmu hneg'
    refine ⟨-t, fun v hv => by have := h1 v hv; linarith [this], fun r hr => by
      have := h2 r hr; linarith [this], ?_⟩
    rcases h3 with ⟨v₀, hv₀, he⟩ | ⟨r₀, hr₀, he⟩
    · exact Or.inl ⟨v₀, hv₀, by linarith [he]⟩
    · exact Or.inr ⟨r₀, hr₀, by linarith [he]⟩

end QuaConProof

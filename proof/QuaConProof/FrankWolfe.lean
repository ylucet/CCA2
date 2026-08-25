/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.Selection

/-!
# Attainment on an unbounded piece

`PROJECT_PLAN.md` Phase 7. At Stage 1 every piece is a polytope, so `psi` attains
its supremum there by compactness. With recession directions the piece is
unbounded and compactness is gone, but the supremum is still attained **whenever
it is finite**. That is the Frank–Wolfe theorem for a quadratic on a polyhedron,
and this file proves it in the V-representation, by induction on the number of
rays.

## The argument

Write the piece as `T = conv V + cone R` and let `Δ = conv R`, a compact set with
`cone R = ℝ₊ · Δ ∪ {0}`. Everything turns on the curvature `α(δ) = δᵀHδ` of the
piece's quadratic along the directions `δ ∈ Δ`.

* **Some `δ ∈ Δ` has `α(δ) ≤ 0`.** If `α(δ) < 0` then `ψ` grows without bound
  along `x + ℝ₊δ`, contradicting finiteness. So `α(δ) = 0`, and then finiteness
  forces the slope `⟨s - ∇q(y), δ⟩` to be `≤ 0` at every `y`, so stepping along
  `δ` never increases `ψ`. Every point of `T` can be pushed *back* along `δ`
  until a ray weight vanishes, so the supremum is already the supremum over the
  finitely many pieces `conv V + cone (R \ {r})`, and the induction applies.

* **Every `δ ∈ Δ` has `α(δ) > 0`.** Then `α` has a positive minimum `m` on the
  compact `Δ`, and for `x = c + λδ` the value `ψ(x) ≤ P + λB - λ²m/2` runs off to
  `-∞` in `λ`, uniformly in `c ∈ conv V` and `δ ∈ Δ`. So beyond some `λ₀` no
  point can be optimal, and the supremum is attained on the **compact** truncation
  `conv V + [0, λ₀]·Δ`.

Note that the dichotomy is over the whole of `Δ`, not over the generators: three
rays can each have positive curvature while a combination of two of them has
negative curvature, so a test on `R` alone would be unsound.

## Main results

* `exists_isMaxOn_of_bddAbove` — Frank–Wolfe: bounded above implies attained.
* `attained_of_conj_ne_top` — the form the selection lemma consumes: wherever the
  conjugate is finite, every piece attains its supremum.
-/

namespace QuaConProof

open scoped Classical

/-! ### Generalised polyhedra -/

/-- A generalised polyhedron in V-representation: `conv V + cone R`. This is
`QuaPiece.T`, freed from the piece so that the induction can shrink `R`. -/
def poly (V R : Finset Plane) : Set Plane :=
  {x | ∃ c ∈ convexHull ℝ (↑V : Set Plane), ∃ k ∈ coneHull R, c + k = x}

lemma QuaPiece.T_eq_poly (p : QuaPiece) : p.T = poly p.verts p.rays := rfl

lemma mem_poly_of_mem_convexHull {V R : Finset Plane} {c : Plane}
    (hc : c ∈ convexHull ℝ (↑V : Set Plane)) : c ∈ poly V R :=
  ⟨c, hc, 0, zero_mem_coneHull _, by simp⟩

@[simp] lemma poly_empty (V : Finset Plane) : poly V ∅ = convexHull ℝ (↑V : Set Plane) := by
  ext x
  simp only [poly, coneHull_empty, Set.mem_ofPred_eq, Set.mem_singleton_iff]
  constructor
  · rintro ⟨c, hc, k, hk, rfl⟩
    rw [hk, add_zero]
    exact hc
  · intro hx
    exact ⟨x, hx, 0, rfl, by simp⟩

lemma poly_nonempty {V R : Finset Plane} (hV : V.Nonempty) : (poly V R).Nonempty := by
  obtain ⟨v, hv⟩ := hV
  exact ⟨v, mem_poly_of_mem_convexHull (subset_convexHull ℝ _ hv)⟩

/-- Dropping a ray shrinks the polyhedron. -/
lemma coneHull_mono {R R' : Finset Plane} (h : R' ⊆ R) : coneHull R' ⊆ coneHull R := by
  rintro k ⟨mu, hmu, rfl⟩
  refine ⟨fun z => if z ∈ R' then mu z else 0, fun z _ => by
    by_cases hz : z ∈ R'
    · simpa [hz] using hmu z hz
    · simp [hz], ?_⟩
  rw [← Finset.sum_subset h (fun z _ hz => by simp [hz])]
  exact Finset.sum_congr rfl fun z hz => by simp [hz]

lemma poly_mono_rays {V R R' : Finset Plane} (h : R' ⊆ R) : poly V R' ⊆ poly V R := by
  rintro x ⟨c, hc, k, hk, rfl⟩
  exact ⟨c, hc, k, coneHull_mono h hk, rfl⟩

/-- `Δ = conv R` sits inside the cone it generates. -/
lemma convexHull_subset_coneHull (R : Finset Plane) :
    convexHull ℝ (↑R : Set Plane) ⊆ coneHull R := by
  intro δ hδ
  obtain ⟨nu, h0, _, h2⟩ := exists_weights hδ
  exact ⟨nu, h0, h2⟩

lemma smul_mem_coneHull_of_mem_convexHull {R : Finset Plane} {δ : Plane}
    (hδ : δ ∈ convexHull ℝ (↑R : Set Plane)) {l : ℝ} (hl : 0 ≤ l) :
    l • δ ∈ coneHull R :=
  smul_mem_coneHull (convexHull_subset_coneHull R hδ) hl

lemma add_smul_mem_poly {V R : Finset Plane} {x δ : Plane} (hx : x ∈ poly V R)
    (hδ : δ ∈ convexHull ℝ (↑R : Set Plane)) {t : ℝ} (ht : 0 ≤ t) :
    x + t • δ ∈ poly V R := by
  obtain ⟨c, hc, k, hk, rfl⟩ := hx
  exact ⟨c, hc, k + t • δ, add_mem_coneHull hk (smul_mem_coneHull_of_mem_convexHull hδ ht),
    by abel⟩

/-- **A nonzero cone vector is a positive multiple of a point of `Δ`.** -/
lemma exists_scaling_of_mem_coneHull {R : Finset Plane} {k : Plane}
    (hk : k ∈ coneHull R) (hk0 : k ≠ 0) :
    ∃ l : ℝ, 0 < l ∧ ∃ δ ∈ convexHull ℝ (↑R : Set Plane), k = l • δ := by
  obtain ⟨mu, hmu, hsum⟩ := hk
  set l : ℝ := ∑ r ∈ R, mu r with hl
  have hlpos : 0 < l := by
    rcases lt_or_eq_of_le (Finset.sum_nonneg hmu) with h | h
    · exact h
    · exfalso
      refine hk0 ?_
      rw [← hsum]
      refine Finset.sum_eq_zero fun r hr => ?_
      have : mu r = 0 := by
        have hz := (Finset.sum_eq_zero_iff_of_nonneg hmu).1 h.symm
        exact hz r hr
      rw [this, zero_smul]
  refine ⟨l, hlpos, l⁻¹ • k, ?_, by rw [smul_smul, mul_inv_cancel₀ (ne_of_gt hlpos), one_smul]⟩
  refine mem_convexHull_of_weights (lam := fun z => l⁻¹ * mu z)
    (fun z hz => mul_nonneg (le_of_lt (inv_pos.2 hlpos)) (hmu z hz)) ?_ ?_
  · rw [← Finset.mul_sum, ← hl, inv_mul_cancel₀ (ne_of_gt hlpos)]
  · rw [← hsum, Finset.smul_sum]
    exact Finset.sum_congr rfl fun z _ => by rw [mul_smul]

/-! ### Pushing back along a direction until a ray weight vanishes -/

/-- **Ray descent.** Given a point of the cone and a convex combination `δ` of
the rays, one may subtract a nonnegative multiple of `δ` and land in the cone of
a strictly smaller ray set. -/
lemma exists_ray_descent {R : Finset Plane} {k δ : Plane} {mu nu : Plane → ℝ}
    (hmu : ∀ r ∈ R, 0 ≤ mu r) (hk : (∑ r ∈ R, mu r • r) = k)
    (hnu : ∀ r ∈ R, 0 ≤ nu r) (hnusum : (∑ r ∈ R, nu r) = 1)
    (hδ : (∑ r ∈ R, nu r • r) = δ) :
    ∃ (t : ℝ) (r₀ : Plane), 0 ≤ t ∧ r₀ ∈ R ∧ k - t • δ ∈ coneHull (R.erase r₀) := by
  classical
  set N : Finset Plane := R.filter (fun r => 0 < nu r) with hN
  have hNne : N.Nonempty := by
    by_contra hc
    rw [Finset.not_nonempty_iff_eq_empty, hN, Finset.filter_eq_empty_iff] at hc
    have hall : ∀ r ∈ R, nu r = 0 := fun r hr => le_antisymm (not_lt.1 (hc hr)) (hnu r hr)
    rw [Finset.sum_congr rfl (fun r hr => hall r hr), Finset.sum_const_zero] at hnusum
    exact zero_ne_one hnusum
  obtain ⟨r₀, hr₀N, hr₀min⟩ := Finset.exists_min_image N (fun r => mu r / nu r) hNne
  have hr₀R : r₀ ∈ R := (Finset.mem_filter.1 hr₀N).1
  have hr₀pos : 0 < nu r₀ := (Finset.mem_filter.1 hr₀N).2
  set t : ℝ := mu r₀ / nu r₀ with ht
  have htnn : 0 ≤ t := div_nonneg (hmu r₀ hr₀R) hr₀pos.le
  have hnn : ∀ r ∈ R, 0 ≤ mu r - t * nu r := by
    intro r hr
    by_cases hp : 0 < nu r
    · have hle : t ≤ mu r / nu r := hr₀min r (Finset.mem_filter.2 ⟨hr, hp⟩)
      rw [le_div_iff₀ hp] at hle
      linarith
    · have hz : nu r = 0 := le_antisymm (not_lt.1 hp) (hnu r hr)
      rw [hz, mul_zero, sub_zero]
      exact hmu r hr
  have hnu0 : nu r₀ ≠ 0 := ne_of_gt hr₀pos
  have hzero : mu r₀ - t * nu r₀ = 0 := by
    have hcan : t * nu r₀ = mu r₀ := by rw [ht, div_mul_cancel₀ _ hnu0]
    rw [hcan, sub_self]
  have hfz : (fun z : Plane => (mu z - t * nu z) • z) r₀ = 0 := by
    show (mu r₀ - t * nu r₀) • r₀ = 0
    rw [hzero, zero_smul]
  refine ⟨t, r₀, htnn, hr₀R, ?_⟩
  refine ⟨fun z => mu z - t * nu z, fun z hz => hnn z (Finset.mem_of_mem_erase hz), ?_⟩
  show (∑ z ∈ R.erase r₀, (mu z - t * nu z) • z) = k - t • δ
  rw [Finset.sum_erase (f := fun z : Plane => (mu z - t * nu z) • z) R hfz]
  have hexp : ∀ z : Plane, (mu z - t * nu z) • z = mu z • z - t • (nu z • z) := by
    intro z; module
  simp only [hexp]
  rw [Finset.sum_sub_distrib, ← Finset.smul_sum, hk, hδ]

/-! ### Continuity of the two functions the compactness arguments use -/

lemma continuous_alongCurv (q : Quad) : Continuous q.alongCurv := by
  unfold Quad.alongCurv
  fun_prop

lemma continuous_absSlope (q : Quad) (s : Plane) :
    Continuous (fun z : Plane × Plane => |dot s z.2 - dot (q.gradAt z.1) z.2|) := by
  refine Continuous.abs ?_
  unfold dot Quad.gradAt
  fun_prop

/-! ### Frank–Wolfe -/

private theorem attain_aux (q : Quad) (s : Plane) :
    ∀ n : ℕ, ∀ (V R : Finset Plane), R.card = n → V.Nonempty → ∀ M : ℝ,
      (∀ x ∈ poly V R, psi q s x ≤ M) →
      ∃ x ∈ poly V R, IsMaxOn (psi q s) (poly V R) x := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro V R hcard hVne M hM
    rcases Finset.eq_empty_or_nonempty R with hRe | hRne
    · -- **No rays**: the piece is compact.
      subst hRe
      rw [poly_empty]
      exact (V.finite_toSet.isCompact_convexHull ℝ).exists_isMaxOn
        (by obtain ⟨v, hv⟩ := hVne; exact ⟨v, subset_convexHull ℝ _ hv⟩)
        (continuous_psi q s).continuousOn
    · have hDne : (convexHull ℝ (↑R : Set Plane)).Nonempty := by
        obtain ⟨r, hr⟩ := hRne
        exact ⟨r, subset_convexHull ℝ _ hr⟩
      by_cases hA : ∃ δ ∈ convexHull ℝ (↑R : Set Plane), q.alongCurv δ ≤ 0
      · -- **A flat or negatively curved direction.**
        obtain ⟨δ, hδ, hδc⟩ := hA
        obtain ⟨x₀, hx₀⟩ := poly_nonempty (R := R) hVne
        -- negative curvature is impossible: `psi` would be unbounded along `δ`
        have hδ0 : q.alongCurv δ = 0 := by
          rcases lt_or_eq_of_le hδc with hlt | heq
          · exfalso
            set a : ℝ := -(q.alongCurv δ) / 2 with ha
            have hapos : 0 < a := by rw [ha]; linarith
            set L : ℝ := dot s δ - dot (q.gradAt x₀) δ with hL
            set E : ℝ := M - psi q s x₀ with hE
            have hEnn : 0 ≤ E := by rw [hE]; linarith [hM x₀ hx₀]
            set t : ℝ := max ((|L| + 1) / a) (E + 1) with ht
            have ht1 : (|L| + 1) / a ≤ t := le_max_left _ _
            have ht2 : E + 1 ≤ t := le_max_right _ _
            have htnn : 0 ≤ t := le_trans (by linarith [abs_nonneg L]) ht2
            have hat : |L| + 1 ≤ a * t := by
              rw [div_le_iff₀ hapos] at ht1
              linarith
            have hmem : x₀ + t • δ ∈ poly V R := add_smul_mem_poly hx₀ hδ htnn
            have hval := hM _ hmem
            rw [psi_along_dir] at hval
            have hLb : -|L| ≤ L := neg_abs_le L
            have hcurv : q.alongCurv δ = -2 * a := by rw [ha]; ring
            rw [hcurv] at hval
            nlinarith [abs_nonneg L]
          · exact heq
        -- the slope along `δ` is nonpositive everywhere on the piece
        have hslope : ∀ y ∈ poly V R, dot s δ - dot (q.gradAt y) δ ≤ 0 := by
          intro y hy
          by_contra hc
          push Not at hc
          set L : ℝ := dot s δ - dot (q.gradAt y) δ with hL
          set t : ℝ := (M - psi q s y + 1) / L with ht
          have hEnn : 0 ≤ M - psi q s y := by linarith [hM y hy]
          have htnn : 0 ≤ t := div_nonneg (by linarith) hc.le
          have hmem : y + t • δ ∈ poly V R := add_smul_mem_poly hy hδ htnn
          have hval := hM _ hmem
          rw [psi_along_dir, hδ0] at hval
          have htL : t * L = M - psi q s y + 1 := by
            rw [ht, div_mul_cancel₀ _ (ne_of_gt hc)]
          rw [htL] at hval
          linarith
        -- every point is dominated by a point of a smaller piece
        have hdrop : ∀ x ∈ poly V R, ∃ r₀ ∈ R, ∃ y ∈ poly V (R.erase r₀), psi q s x ≤ psi q s y := by
          intro x hx
          obtain ⟨c, hc, k, ⟨mu, hmu, hksum⟩, rfl⟩ := hx
          obtain ⟨nu, hnu0, hnu1, hnu2⟩ := exists_weights hδ
          obtain ⟨t, r₀, htnn, hr₀, hmem⟩ := exists_ray_descent hmu hksum hnu0 hnu1 hnu2
          refine ⟨r₀, hr₀, c + (k - t • δ), ⟨c, hc, k - t • δ, hmem, rfl⟩, ?_⟩
          have hy : c + (k - t • δ) ∈ poly V R :=
            poly_mono_rays (Finset.erase_subset _ _) ⟨c, hc, k - t • δ, hmem, rfl⟩
          have hback : c + k = (c + (k - t • δ)) + t • δ := by abel
          rw [hback, psi_along_dir, hδ0]
          have := hslope _ hy
          nlinarith [this, htnn]
        -- so the supremum is the best of finitely many smaller problems
        have hsmaller : ∀ r₀ ∈ R, ∃ y ∈ poly V (R.erase r₀),
            IsMaxOn (psi q s) (poly V (R.erase r₀)) y := by
          intro r₀ hr₀
          refine ih (R.erase r₀).card ?_ V (R.erase r₀) rfl hVne M ?_
          · rw [← hcard]
            exact Finset.card_erase_lt_of_mem hr₀
          · exact fun z hz => hM z (poly_mono_rays (Finset.erase_subset _ _) hz)
        choose! ym hymT hymMax using hsmaller
        obtain ⟨r₁, hr₁, hr₁sup⟩ :=
          Finset.exists_mem_eq_sup' hRne (fun r => psi q s (ym r))
        refine ⟨ym r₁, poly_mono_rays (Finset.erase_subset _ _) (hymT r₁ hr₁), ?_⟩
        intro z hz
        show psi q s z ≤ psi q s (ym r₁)
        obtain ⟨r₀, hr₀, y, hyT, hyle⟩ := hdrop z hz
        refine le_trans hyle (le_trans (hymMax r₀ hr₀ hyT) ?_)
        rw [← hr₁sup]
        exact Finset.le_sup' (fun r => psi q s (ym r)) hr₀
      · -- **Every direction of `Δ` is positively curved.**
        push Not at hA
        set C : Set Plane := convexHull ℝ (↑V : Set Plane) with hC
        set D : Set Plane := convexHull ℝ (↑R : Set Plane) with hD
        have hCc : IsCompact C := V.finite_toSet.isCompact_convexHull ℝ
        have hDc : IsCompact D := R.finite_toSet.isCompact_convexHull ℝ
        have hCne : C.Nonempty := by
          obtain ⟨v, hv⟩ := hVne; exact ⟨v, subset_convexHull ℝ _ hv⟩
        obtain ⟨c₀, hc₀⟩ := id hCne
        -- the minimum curvature on `Δ` is positive
        obtain ⟨δm, hδm, hδmin⟩ :=
          hDc.exists_isMinOn hDne (continuous_alongCurv q).continuousOn
        set m : ℝ := q.alongCurv δm with hm
        have hmpos : 0 < m := hA δm hδm
        -- a uniform bound on the slope, and on `psi`, over the compact parts
        obtain ⟨zB, hzB, hBmax⟩ :=
          (hCc.prod hDc).exists_isMaxOn (hCne.prod hDne) (continuous_absSlope q s).continuousOn
        set B : ℝ := |dot s zB.2 - dot (q.gradAt zB.1) zB.2| with hB
        have hBnn : 0 ≤ B := abs_nonneg _
        obtain ⟨cP, hcP, hPmax⟩ := hCc.exists_isMaxOn hCne (continuous_psi q s).continuousOn
        set P : ℝ := psi q s cP with hP
        set E : ℝ := max (P - psi q s c₀) 0 with hE
        have hEnn : 0 ≤ E := le_max_right _ _
        set l₀ : ℝ := max (2 * (B + 1) / m) (E + 1) with hl₀
        have hl₀pos : 0 < l₀ := lt_of_lt_of_le (by linarith) (le_max_right _ _)
        -- the compact truncation
        set K : Set Plane := (fun z : ℝ × Plane => z.1 • z.2) '' (Set.Icc 0 l₀ ×ˢ D) with hK
        have hKc : IsCompact K := (isCompact_Icc.prod hDc).image (by fun_prop)
        set A : Set Plane := (fun z : Plane × Plane => z.1 + z.2) '' (C ×ˢ K) with hA'
        have hAc : IsCompact A := (hCc.prod hKc).image (by fun_prop)
        have hzeroK : (0 : Plane) ∈ K := by
          obtain ⟨δ, hδ⟩ := hDne
          exact ⟨(0, δ), ⟨Set.left_mem_Icc.2 hl₀pos.le, hδ⟩, by simp⟩
        have hAne : A.Nonempty := ⟨c₀ + 0, ⟨(c₀, 0), ⟨hc₀, hzeroK⟩, rfl⟩⟩
        have hAsub : A ⊆ poly V R := by
          rintro _ ⟨⟨c, k⟩, ⟨hc, hk⟩, rfl⟩
          obtain ⟨⟨l, δ⟩, ⟨hl, hδ⟩, rfl⟩ := hk
          exact ⟨c, hc, l • δ, smul_mem_coneHull_of_mem_convexHull hδ hl.1, rfl⟩
        obtain ⟨xM, hxMA, hxMmax⟩ := hAc.exists_isMaxOn hAne (continuous_psi q s).continuousOn
        have hc₀A : c₀ ∈ A := ⟨(c₀, 0), ⟨hc₀, hzeroK⟩, by simp⟩
        refine ⟨xM, hAsub hxMA, ?_⟩
        rintro y ⟨c, hc, k, hk, rfl⟩
        show psi q s (c + k) ≤ psi q s xM
        by_cases hk0 : k = 0
        · subst hk0
          exact hxMmax ⟨(c, 0), ⟨hc, hzeroK⟩, by simp⟩
        obtain ⟨l, hlpos, δ, hδ, hkl⟩ := exists_scaling_of_mem_coneHull hk hk0
        by_cases hle : l ≤ l₀
        · refine hxMmax ⟨(c, k), ⟨hc, ?_⟩, rfl⟩
          exact ⟨(l, δ), ⟨⟨hlpos.le, hle⟩, hδ⟩, hkl.symm⟩
        · push Not at hle
          have hcP' : psi q s c ≤ P := hPmax hc
          have hslope : dot s δ - dot (q.gradAt c) δ ≤ B := by
            refine le_trans (le_abs_self _) ?_
            exact hBmax (Set.mk_mem_prod hc hδ)
          have hcurv : m ≤ q.alongCurv δ := hδmin hδ
          have hkey : psi q s (c + k) ≤ P + l * B - l ^ 2 * m / 2 := by
            rw [hkl, psi_along_dir]
            nlinarith [hlpos, hcurv, hslope, hcP']
          have hb1 : 2 * (B + 1) / m ≤ l₀ := le_max_left _ _
          have hb2 : E + 1 ≤ l₀ := le_max_right _ _
          have hml : 2 * (B + 1) ≤ m * l₀ := by
            rw [div_le_iff₀ hmpos] at hb1
            linarith
          have hfinal : P + l * B - l ^ 2 * m / 2 < psi q s c₀ := by
            have hEge : P - psi q s c₀ ≤ E := le_max_left _ _
            have h1 : m * l₀ < m * l := mul_lt_mul_of_pos_left hle hmpos
            have h2 : B + 1 ≤ m * l / 2 := by linarith
            have h3 : l * 1 ≤ l * (m * l / 2 - B) :=
              mul_le_mul_of_nonneg_left (by linarith) hlpos.le
            have h4 : l ^ 2 * m / 2 - l * B = l * (m * l / 2 - B) := by ring
            linarith
          exact le_trans hkey (le_of_lt (lt_of_lt_of_le hfinal (hxMmax hc₀A)))

/-- **Frank–Wolfe for a `QuaPiece`.** If `psi` is bounded above on the piece then
it attains its supremum there — with or without recession directions. -/
theorem exists_isMaxOn_of_bddAbove {p : QuaPiece} {s : Plane} {M : ℝ}
    (hM : ∀ x ∈ p.T, psi p.q s x ≤ M) :
    ∃ x ∈ p.T, IsMaxOn (psi p.q s) p.T x :=
  attain_aux p.q s p.rays.card p.verts p.rays rfl p.verts_nonempty M hM

/-- **Where the conjugate is finite, every piece attains its supremum.**

This is the hypothesis `selection` consumes, and it is the last thing Phase 7
owed: `le_conj` bounds each piece's objective by `f*(s)`, which is a real number
exactly when `f*(s) ≠ ⊤`, and Frank–Wolfe turns that bound into a maximiser. -/
theorem attained_of_conj_ne_top {f : QuaPol} {s : Plane} (h : f.conj s ≠ ⊤) :
    f.Attained s := by
  intro p hp
  refine exists_isMaxOn_of_bddAbove (M := (f.conj s).toReal) ?_
  intro x hx
  have hb : ((psi p.q s x : ℝ) : EReal) ≤ f.conj s := le_conj hp hx s
  rwa [← EReal.coe_toReal h (QuaPol.conj_ne_bot f s), EReal.coe_le_coe_iff] at hb

/-! ### Sanity: an unbounded piece whose conjugate really is `⊤`

`CLAUDE.md` → Verification, point 3, for the Phase 7 definitions. The piece is
the nonnegative `s₁`-axis carrying the zero quadratic; its conjugate is `+∞` at
`(1,0)`, because `⟨(1,0), (t,0)⟩ = t` is unbounded on it. A `coneHull` that had
collapsed to `{0}`, or a `T` that had dropped its rays, would make this false.

It lives here rather than beside `conj_isQuaCon` because two later files need
it: `QuaCon.lean` to show the `⊤` cell is inhabited, and `Biconj.lean` for the
counterexample of `TODO.md` C5. -/

namespace Sanity

/-- One piece: the nonnegative `s₁`-axis, carrying the zero quadratic. -/
noncomputable def rayPiece : QuaPiece :=
  ⟨{0}, Finset.singleton_nonempty _, {(1, 0)}, 0⟩

/-- The `QuaPol` with that single piece. -/
noncomputable def rayPol : QuaPol := ⟨{rayPiece}, Finset.singleton_nonempty _⟩

lemma mem_rayPiece {t : ℝ} (ht : 0 ≤ t) : ((t, 0) : Plane) ∈ rayPiece.T := by
  refine ⟨0, subset_convexHull ℝ _ (by simp [rayPiece]), t • ((1, 0) : Plane),
    smul_mem_coneHull (mem_coneHull_of_mem (show ((1, 0) : Plane) ∈ rayPiece.rays by
      simp [rayPiece])) ht, ?_⟩
  simp

lemma eval_rayPol {t : ℝ} (ht : 0 ≤ t) : rayPol.eval (t, 0) = 0 := by
  have hm := mem_rayPiece ht
  simp only [rayPol, QuaPol.eval, Finset.inf_singleton]
  rw [if_pos hm]
  simp [rayPiece]

/-- **The conjugate really is `⊤` there.** -/
theorem conj_rayPol_eq_top : rayPol.conj (1, 0) = ⊤ := by
  by_contra hne
  have hbot := QuaPol.conj_ne_bot rayPol (1, 0)
  set c : ℝ := (rayPol.conj (1, 0)).toReal with hc
  have hcoe : rayPol.conj (1, 0) = (c : EReal) := (EReal.coe_toReal hne hbot).symm
  have hle : ((|c| + 1 : ℝ) : EReal) ≤ rayPol.conj (1, 0) := by
    rw [QuaPol.conj_def]
    refine le_iSup_of_le ((|c| + 1, 0) : Plane) ?_
    rw [eval_rayPol (by positivity)]
    simp [dot]
  rw [hcoe, EReal.coe_le_coe_iff] at hle
  linarith [le_abs_self c]

end Sanity

end QuaConProof

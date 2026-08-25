/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.Realization
import QuaConProof.Bary

/-!
# A realised elliptical edge

`Shapes.lean` shows an elliptical tie conic needs an interior branch.
`Realization.lean` gives the certificate machinery. This file closes the loop
with an explicit `QuaPol` for which **infinitely many points of the plane have
two interior branches simultaneously active**, and all of them lie on one
non-degenerate ellipse.

That upgrades "an ellipse is in the candidate list" to "an ellipse is genuinely
met", which is the statement that bears on whether the conjugate can be stored
as a `QuaPar` (a structure requiring `b² - 4ac = 0` of every edge conic).

## Why this example and not the one in the documents

`../CONJ_FIELD_PROOF.md` §7.5's three-piece example has an elliptical edge, but
its tie point is `(-17/62 + √612030/186, 3/2)` — irrational. In fact that
ellipse, `93s₁² + 38s₁s₂ + 39s₂² - 6s₁ - 482s₂ - 1003 = 0`, has **no rational
point at all**: reduced mod `49` the homogenised form has no primitive zero, a
`7`-adic obstruction. So no rational witness exists there and a computational
certificate is impossible; realising *that* edge would need an
intermediate-value argument.

The example below is built instead so that the tie conic is a **circle through
the origin**, which is rationally parametrised, giving infinitely many rational
witnesses.

## The example

Two disjoint triangular pieces, both with positive definite Hessians:

* `q₁ = x₁² + x₂²` on the triangle `(-1,-1), (1,-1), (0,2)`;
* `q₂ = 2x₁² + 2x₂² + 8x₂ + 8` on `(-2,-4), (2,-4), (0,-3/2)`.

Their interior branches are `I₁ = (s₁² + s₂²)/4` and
`I₂ = (s₁² + s₂²)/8 - 2s₂`, and

`I₁ - I₂ = (s₁² + s₂² + 16s₂)/8`,

the circle `s₁² + (s₂+8)² = 64` — `disc = -1/16 < 0` and `det₃ = -1/2 ≠ 0`, so a
non-degenerate ellipse. Parametrising it from the origin by the slope `1/n`
gives the rational points

`sₙ = (-16n/(n²+1), -16/(n²+1))`,

and for `n ≥ 12` both maximisers lie genuinely inside their own pieces, so both
branches really are attained there — this is not an artefact of a branch
overshooting.

The pieces are disjoint, hence maximally **non-adjacent**, which is exactly the
comparison [JOGO] Theorem 6's proof does not cover.
-/

namespace QuaConProof

open scoped Classical

namespace Witness

/-! ### Membership in a triangle, from barycentric weights -/

/-- A convex combination of three points lies in their hull. -/
lemma mem_convexHull_triple {va vb vc x : Plane} {la lb lc : ℝ}
    (ha : 0 ≤ la) (hb : 0 ≤ lb) (hc : 0 ≤ lc) (hsum : la + lb + lc = 1)
    (hx : x = la • va + lb • vb + lc • vc) :
    x ∈ convexHull ℝ (↑({va, vb, vc} : Finset Plane) : Set Plane) := by
  refine mem_convexHull_of_weights (lam := fun z =>
    (if z = va then la else 0) + (if z = vb then lb else 0) + (if z = vc then lc else 0))
    ?_ ?_ ?_
  · intro z _
    positivity
  · have hva : va ∈ ({va, vb, vc} : Finset Plane) := by simp
    have hvb : vb ∈ ({va, vb, vc} : Finset Plane) := by simp
    have hvc : vc ∈ ({va, vb, vc} : Finset Plane) := by simp
    rw [Finset.sum_add_distrib, Finset.sum_add_distrib,
      Finset.sum_ite_eq' _ va, Finset.sum_ite_eq' _ vb, Finset.sum_ite_eq' _ vc,
      if_pos hva, if_pos hvb, if_pos hvc]
    exact hsum
  · have hexp : ∀ z : Plane,
        ((if z = va then la else 0) + (if z = vb then lb else 0)
          + (if z = vc then lc else 0)) • z
        = (if z = va then la else 0) • z + (if z = vb then lb else 0) • z
          + (if z = vc then lc else 0) • z := by
      intro z; module
    have hva : va ∈ ({va, vb, vc} : Finset Plane) := by simp
    have hvb : vb ∈ ({va, vb, vc} : Finset Plane) := by simp
    have hvc : vc ∈ ({va, vb, vc} : Finset Plane) := by simp
    simp only [hexp]
    rw [Finset.sum_add_distrib, Finset.sum_add_distrib]
    simp only [ite_smul, zero_smul]
    rw [Finset.sum_ite_eq' _ va, Finset.sum_ite_eq' _ vb, Finset.sum_ite_eq' _ vc,
      if_pos hva, if_pos hvb, if_pos hvc]
    exact hx.symm

/-! ### The two pieces -/

/-- `q₁ = x₁² + x₂²`, positive definite. -/
def q1 : Quad := ⟨1, 0, 1, 0, 0, 0⟩

/-- `q₂ = 2x₁² + 2x₂² + 8x₂ + 8`, positive definite. -/
def q2 : Quad := ⟨2, 0, 2, 0, 8, 8⟩

noncomputable def V1 : Finset Plane := {(-1, -1), (1, -1), (0, 2)}
noncomputable def V2 : Finset Plane := {(-2, -4), (2, -4), (0, -3/2)}

lemma V1_nonempty : V1.Nonempty := ⟨(-1, -1), by simp [V1]⟩
lemma V2_nonempty : V2.Nonempty := ⟨(-2, -4), by simp [V2]⟩

noncomputable def p1 : QuaPiece := ⟨V1, V1_nonempty, ∅, q1⟩
noncomputable def p2 : QuaPiece := ⟨V2, V2_nonempty, ∅, q2⟩

@[simp] lemma p1_rays : p1.rays = ∅ := rfl
@[simp] lemma p2_rays : p2.rays = ∅ := rfl

/-- The two-piece `QuaPol`. Its pieces are disjoint, hence non-adjacent. -/
noncomputable def f : QuaPol := ⟨{p1, p2}, ⟨p1, by simp⟩⟩

lemma mem_pieces_iff {p : QuaPiece} : p ∈ f.pieces ↔ p = p1 ∨ p = p2 := by
  simp [f]

/-! ### The tie conic is a non-degenerate ellipse -/

lemma interiorBranch_q1 : interiorBranch q1 = ⟨1/4, 0, 1/4, 0, 0, 0⟩ := by
  simp only [interiorBranch, Quad.hessDet, q1, Quad.mk.injEq]
  norm_num

lemma interiorBranch_q2 : interiorBranch q2 = ⟨1/8, 0, 1/8, 0, -2, 0⟩ := by
  simp only [interiorBranch, Quad.hessDet, q2, Quad.mk.injEq]
  norm_num

/-- **The tie conic of the two interior branches is a non-degenerate ellipse.**
It is the circle `s₁² + (s₂+8)² = 64`. -/
theorem tie_isEllipse : (interiorBranch q1 - interiorBranch q2).kind = ConicKind.ellipse := by
  rw [Quad.kind_eq_ellipse_iff]
  rw [interiorBranch_q1, interiorBranch_q2]
  refine ⟨?_, ?_, ?_⟩
  · rw [Ne, Quad.eq_zero_iff]; norm_num
  · simp only [Quad.det3]; norm_num
  · simp only [Quad.disc]; norm_num

/-! ### The rational points of the tie conic -/

/-- `sₙ = (-16n/(n²+1), -16/(n²+1))`, offset so that `n ≥ 12` always. -/
noncomputable def sw (k : ℕ) : Plane :=
  (-16 * ((k : ℝ) + 12) / (((k : ℝ) + 12) ^ 2 + 1), -16 / (((k : ℝ) + 12) ^ 2 + 1))

lemma sw_den_pos (k : ℕ) : (0 : ℝ) < ((k : ℝ) + 12) ^ 2 + 1 := by positivity

lemma sw_injective : Function.Injective sw := by
  intro m n h
  have h2 : (-16 : ℝ) / (((m : ℝ) + 12) ^ 2 + 1) = -16 / (((n : ℝ) + 12) ^ 2 + 1) :=
    congrArg Prod.snd h
  have hm := sw_den_pos m
  have hn := sw_den_pos n
  have : ((m : ℝ) + 12) ^ 2 = ((n : ℝ) + 12) ^ 2 := by
    field_simp at h2; linarith
  have hm0 : (0 : ℝ) ≤ (m : ℝ) := Nat.cast_nonneg m
  have hn0 : (0 : ℝ) ≤ (n : ℝ) := Nat.cast_nonneg n
  have hmn : (m : ℝ) = (n : ℝ) := by nlinarith
  exact_mod_cast hmn

/-- Each `sₙ` lies on the tie conic: the two interior branches agree there. -/
theorem tie_at_sw (k : ℕ) :
    (interiorBranch q1).eval (sw k) = (interiorBranch q2).eval (sw k) := by
  have hd := sw_den_pos k
  rw [interiorBranch_q1, interiorBranch_q2]
  simp only [Quad.eval, sw]
  field_simp
  ring

/-! ### Both maximisers lie inside their pieces -/

lemma interiorPoint_q1_sw (k : ℕ) :
    interiorPoint q1 (sw k)
      = (-8 * ((k : ℝ) + 12) / (((k : ℝ) + 12) ^ 2 + 1), -8 / (((k : ℝ) + 12) ^ 2 + 1)) := by
  have hd := sw_den_pos k
  simp only [interiorPoint, Quad.hessDet, q1, sw, Prod.mk.injEq]
  constructor <;> (field_simp; ring)

/-- **The maximiser of piece 1 lies inside piece 1**, for every `n ≥ 12`.

The barycentric weights over `(-1,-1), (1,-1), (0,2)` are
`(D+4+12n)/(3D)`, `(D+4-12n)/(3D)`, `(D-8)/(3D)` with `D = n²+1`; the middle one
is `(n² - 12n + 5)/(3D)`, which is where `n ≥ 12` is needed. -/
theorem interiorPoint_q1_mem (k : ℕ) : interiorPoint q1 (sw k) ∈ p1.T := by
  rw [QuaPiece.T_of_rays_empty p1_rays]
  have hd := sw_den_pos k
  set n : ℝ := (k : ℝ) + 12 with hn
  have hk0 : (0 : ℝ) ≤ (k : ℝ) := Nat.cast_nonneg k
  have hn12 : (12 : ℝ) ≤ n := by rw [hn]; linarith
  have hD : (0 : ℝ) < n ^ 2 + 1 := by positivity
  refine mem_convexHull_triple (la := (n ^ 2 + 1 + 4 + 12 * n) / (3 * (n ^ 2 + 1)))
    (lb := (n ^ 2 + 1 + 4 - 12 * n) / (3 * (n ^ 2 + 1)))
    (lc := (n ^ 2 + 1 - 8) / (3 * (n ^ 2 + 1))) ?_ ?_ ?_ ?_ ?_
  · apply div_nonneg (by nlinarith) (by linarith)
  · apply div_nonneg (by nlinarith) (by linarith)
  · apply div_nonneg (by nlinarith) (by linarith)
  · field_simp; ring
  · rw [interiorPoint_q1_sw]
    simp only [← hn, Prod.ext_iff, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
      smul_eq_mul]
    constructor <;> (field_simp; ring)

/-! ### The certificate, and the conclusion -/

/-- The witness is a Stage 1 `QuaPol`: neither piece has a recession direction. -/
lemma bounded : f.Bounded := by
  intro p hp
  rcases mem_pieces_iff.1 hp with rfl | rfl <;> rfl

lemma convex_pieces : ∀ p ∈ f.pieces, 0 < p.q.a ∧ 0 < p.q.hessDet := by
  intro p hp
  rcases mem_pieces_iff.1 hp with rfl | rfl
  · exact ⟨by norm_num [p1, q1], by norm_num [p1, q1, Quad.hessDet]⟩
  · exact ⟨by norm_num [p2, q2], by norm_num [p2, q2, Quad.hessDet]⟩

/-- **The conjugate at `sₙ` is the first interior branch.** -/
theorem conj_at_sw (k : ℕ) :
    f.conj (sw k) = (((interiorBranch q1).eval (sw k) : ℝ) : EReal) := by
  refine conj_eq_interiorBranch (p := p1) (by simp [f]) convex_pieces
    (interiorPoint_q1_mem k) ?_
  intro p' hp'
  rcases mem_pieces_iff.1 hp' with rfl | rfl
  · exact le_rfl
  · exact le_of_eq (tie_at_sw k).symm

lemma interiorBranch_q1_mem_cand : interiorBranch q1 ∈ cand f :=
  mem_cand (p := p1) (by simp [f]) (QuaPiece.interiorBranch_mem_branches
    (by norm_num [p1, q1, Quad.hessDet]))

lemma interiorBranch_q2_mem_cand : interiorBranch q2 ∈ cand f :=
  mem_cand (p := p2) (by simp [f]) (QuaPiece.interiorBranch_mem_branches
    (by norm_num [p2, q2, Quad.hessDet]))

/-- Both interior branches are active at every `sₙ`. -/
theorem both_active (k : ℕ) :
    interiorBranch q1 ∈ active f (sw k) ∧ interiorBranch q2 ∈ active f (sw k) := by
  constructor
  · exact mem_active_iff.2 ⟨interiorBranch_q1_mem_cand, (conj_at_sw k).symm⟩
  · refine mem_active_iff.2 ⟨interiorBranch_q2_mem_cand, ?_⟩
    rw [conj_at_sw k, tie_at_sw k]

/-- **An elliptical edge is realised.**

There is a `QuaPol` and two of its candidate quadratics such that

* their tie conic is a **non-degenerate ellipse**;
* **infinitely many** points of the plane have both of them active;
* and every such point lies on that ellipse.

At each of the exhibited points both maximisers lie genuinely inside their own
pieces, so both branches are attained: this is a real meeting of two faces of
`f*`, not a branch overshooting. -/
theorem ellipse_realised :
    (interiorBranch q1 - interiorBranch q2).kind = ConicKind.ellipse
  ∧ {s : Plane | interiorBranch q1 ∈ active f s ∧ interiorBranch q2 ∈ active f s}.Infinite
  ∧ {s : Plane | interiorBranch q1 ∈ active f s ∧ interiorBranch q2 ∈ active f s}
      ⊆ {s : Plane | (interiorBranch q1).eval s = (interiorBranch q2).eval s} := by
  refine ⟨tie_isEllipse, ?_, ?_⟩
  · exact Set.infinite_of_injective_forall_mem (f := sw) sw_injective
      (fun k => both_active k)
  · intro s hs
    exact mem_eqLocus_of_active hs.1 hs.2

end Witness

end QuaConProof

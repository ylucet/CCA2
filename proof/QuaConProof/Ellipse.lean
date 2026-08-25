/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.Conic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Linarith
import Mathlib.Algebra.Order.Field.Basic

/-!
# What `disc < 0` actually means: the elliptic normal form

`TODO.md` item A5. `Conic.lean` proves `disc` and `det₃` are invariants and
`Shapes.lean` reads a name off their signs, but so far nothing said the set so
named *has* that shape. This file closes the elliptic case, and with it the last
place where the development said "ellipse" without a theorem behind it.

## No rotation is needed

The plan called for `Quad.rotate` and `sin² + cos² = 1`, to reach the
axis-aligned normal form. Completing the square twice is shorter and gives the
same theorem: at the centre `c` (the solution of `∇q = 0`, which exists exactly
when `disc ≠ 0`),

    q(p) = a·u² + (-disc/4a)·v² + q(c),    u = (p₁ - c₁) + (b/2a)(p₂ - c₂),
                                           v = p₂ - c₂

so `u, v` is a **shear** of the coordinates rather than a rotation. That is
enough, because an ellipse is exactly the image of a circle under *any*
invertible affine map, not only under a rotation — which is how `IsEllipse` is
stated. Every positive definite binary form is `A(u + Dv)² + Bv²` with
`A, B > 0` (its Cholesky factorisation), so no generality is lost.

## The real/imaginary split

The sign of `det₃` alone cannot make it. With `disc < 0` the zero set is

| condition | the set |
|---|---|
| `a·det₃ < 0` | an ellipse |
| `det₃ = 0` | the single point `c` |
| `a·det₃ > 0` | **empty** — the imaginary ellipse |

`a·det₃` is invariant under scaling the whole quadratic (`a` by `t`, `det₃` by
`t³`, so the product by `t⁴ > 0`), so it is a legitimate test. Everything runs
through one identity, `eval_center_mul_disc`: `q(c)·disc = -det₃`.

## Main results

* `Quad.eval_eq_center_form` — the completed square.
* `isEllipse_zeroSet`, `zeroSet_eq_center_of_det3_eq_zero`,
  `zeroSet_eq_empty_of_disc_neg` — the three cases.
-/

namespace QuaConProof

namespace Quad

/-- The **centre** of a central conic: the solution of `∇q = 0`, which exists and
is unique exactly when `disc ≠ 0`. -/
noncomputable def center (q : Quad) : Plane :=
  ((2 * q.c * q.d - q.b * q.e) / q.disc, (2 * q.a * q.e - q.b * q.d) / q.disc)

/-- `disc < 0` forces `a ≠ 0`, since `4ac > b² ≥ 0`. -/
lemma a_ne_zero_of_disc_neg {q : Quad} (h : q.disc < 0) : q.a ≠ 0 := by
  intro ha
  rw [disc, ha] at h
  nlinarith [sq_nonneg q.b]

/-- The centre is stationary, first coordinate. -/
lemma center_stationary_fst {q : Quad} (hd : q.disc ≠ 0) :
    2 * q.a * q.center.1 + q.b * q.center.2 + q.d = 0 := by
  have hrw : 2 * q.a * q.center.1 + q.b * q.center.2 + q.d
      = (2 * q.a * (2 * q.c * q.d - q.b * q.e) + q.b * (2 * q.a * q.e - q.b * q.d)
          + q.d * q.disc) / q.disc := by
    simp only [center]
    field_simp
  rw [hrw]
  have hnum : 2 * q.a * (2 * q.c * q.d - q.b * q.e) + q.b * (2 * q.a * q.e - q.b * q.d)
      + q.d * q.disc = 0 := by
    simp only [disc]; ring
  rw [hnum, zero_div]

/-- The centre is stationary, second coordinate. -/
lemma center_stationary_snd {q : Quad} (hd : q.disc ≠ 0) :
    q.b * q.center.1 + 2 * q.c * q.center.2 + q.e = 0 := by
  have hrw : q.b * q.center.1 + 2 * q.c * q.center.2 + q.e
      = (q.b * (2 * q.c * q.d - q.b * q.e) + 2 * q.c * (2 * q.a * q.e - q.b * q.d)
          + q.e * q.disc) / q.disc := by
    simp only [center]
    field_simp
  rw [hrw]
  have hnum : q.b * (2 * q.c * q.d - q.b * q.e) + 2 * q.c * (2 * q.a * q.e - q.b * q.d)
      + q.e * q.disc = 0 := by
    simp only [disc]; ring
  rw [hnum, zero_div]

/-- **The completed square, at any stationary point.** -/
theorem eval_eq_stationary_form {q : Quad} (ha : q.a ≠ 0) {v : Plane}
    (hv1 : 2 * q.a * v.1 + q.b * v.2 + q.d = 0)
    (hv2 : q.b * v.1 + 2 * q.c * v.2 + q.e = 0) (p : Plane) :
    q.eval p = q.a * ((p.1 - v.1) + (q.b / (2 * q.a)) * (p.2 - v.2)) ^ 2
      + (-q.disc / (4 * q.a)) * (p.2 - v.2) ^ 2 + q.eval v := by
  have hdv : q.d = -(2 * q.a * v.1) - q.b * v.2 := by linarith
  have hev : q.e = -(q.b * v.1) - 2 * q.c * v.2 := by linarith
  simp only [eval, disc, hdv, hev]
  field_simp
  ring

/-- **The completed square, at the centre.** -/
theorem eval_eq_center_form {q : Quad} (hd : q.disc ≠ 0) (ha : q.a ≠ 0) (p : Plane) :
    q.eval p = q.a * ((p.1 - q.center.1) + (q.b / (2 * q.a)) * (p.2 - q.center.2)) ^ 2
      + (-q.disc / (4 * q.a)) * (p.2 - q.center.2) ^ 2 + q.eval q.center :=
  eval_eq_stationary_form ha (center_stationary_fst hd) (center_stationary_snd hd) p

/-- **The value at the centre, and the only place `det₃` is needed.**
`q(c)·disc = -det₃`. -/
theorem eval_center_mul_disc {q : Quad} (hd : q.disc ≠ 0) :
    q.eval q.center * q.disc = -q.det3 := by
  have hrw : q.eval q.center
      = (q.a * (2 * q.c * q.d - q.b * q.e) ^ 2
          + q.b * (2 * q.c * q.d - q.b * q.e) * (2 * q.a * q.e - q.b * q.d)
          + q.c * (2 * q.a * q.e - q.b * q.d) ^ 2
          + (q.d * (2 * q.c * q.d - q.b * q.e) + q.e * (2 * q.a * q.e - q.b * q.d)) * q.disc
          + q.f * q.disc ^ 2) / q.disc ^ 2 := by
    simp only [eval, center]
    field_simp
    ring
  rw [hrw, div_mul_eq_mul_div, div_eq_iff (pow_ne_zero 2 hd)]
  simp only [disc, det3]
  ring

end Quad

/-- Multiplying the completed square by `a` clears the `a` from the second
coefficient and makes the first a square. -/
private lemma mul_a_form {a Δ U V : ℝ} (ha : a ≠ 0) :
    a * (a * U ^ 2 + (-Δ / (4 * a)) * V ^ 2) = a ^ 2 * U ^ 2 + (-Δ / 4) * V ^ 2 := by
  field_simp
  try ring

/-- A set is an **ellipse**: the unit level set of a positive definite quadratic
form centred at a point, written in Cholesky form `A(u + D·v)² + B·v² = 1` with
`A, B > 0`.

Every positive definite binary form factors that way, so this is the general
ellipse in general position — equivalently, the image of the unit circle under an
invertible affine map. -/
def IsEllipse (C : Set Plane) : Prop :=
  ∃ (c : Plane) (A B D : ℝ), 0 < A ∧ 0 < B ∧
    C = {p : Plane | A * ((p.1 - c.1) + D * (p.2 - c.2)) ^ 2 + B * (p.2 - c.2) ^ 2 = 1}

/-- The unit circle is an ellipse, by the definition above. -/
example : IsEllipse {p : Plane | p.1 ^ 2 + p.2 ^ 2 = 1} :=
  ⟨(0, 0), 1, 1, 0, one_pos, one_pos, by ext p; norm_num⟩

/-- **`disc < 0` and `a·det₃ < 0` means ellipse.**

This is the theorem `Quad.kind` has been standing in for: the zero set really is
an ellipse, not merely a quadratic whose discriminant is negative. -/
theorem isEllipse_zeroSet {q : Quad} (hd : q.disc < 0) (h3 : q.a * q.det3 < 0) :
    IsEllipse {p : Plane | q.eval p = 0} := by
  have ha : q.a ≠ 0 := Quad.a_ne_zero_of_disc_neg hd
  have hd0 : q.disc ≠ 0 := ne_of_lt hd
  set K : ℝ := -q.eval q.center with hKdef
  -- `K · disc = det₃`, which is where every sign below comes from
  have hKd : K * q.disc = q.det3 := by
    have := Quad.eval_center_mul_disc hd0
    rw [hKdef]
    linarith [this]
  have hKa : 0 < q.a * K := by
    rcases lt_trichotomy (q.a * K) 0 with h | h | h
    · exfalso
      have hprod : 0 < (q.a * K) * q.disc := mul_pos_of_neg_of_neg h hd
      have heq : q.a * q.det3 = (q.a * K) * q.disc := by rw [← hKd]; ring
      linarith
    · exfalso
      have hK0 : K = 0 := by
        rcases mul_eq_zero.1 h with h' | h'
        · exact absurd h' ha
        · exact h'
      rw [hK0, zero_mul] at hKd
      rw [← hKd, mul_zero] at h3
      exact lt_irrefl 0 h3
    · exact h
  have hKne : K ≠ 0 := by
    intro h
    rw [h, mul_zero] at hKa
    exact lt_irrefl 0 hKa
  have hKK : 0 < K ^ 2 := by positivity
  refine ⟨q.center, q.a / K, (-q.disc / (4 * q.a)) / K, q.b / (2 * q.a), ?_, ?_, ?_⟩
  · have hrw : q.a / K = (q.a * K) / K ^ 2 := by field_simp; try ring
    rw [hrw]
    exact div_pos hKa hKK
  · have hrw : (-q.disc / (4 * q.a)) / K = ((-q.disc) * (q.a * K)) / (4 * q.a ^ 2 * K ^ 2) := by
      field_simp
      try ring
    rw [hrw]
    refine div_pos (mul_pos (by linarith) hKa) ?_
    have : 0 < q.a ^ 2 := by positivity
    positivity
  · ext p
    simp only [Set.mem_ofPred_eq]
    rw [Quad.eval_eq_center_form hd0 ha p]
    have hrw : (q.a / K) * ((p.1 - q.center.1) + (q.b / (2 * q.a)) * (p.2 - q.center.2)) ^ 2
          + ((-q.disc / (4 * q.a)) / K) * (p.2 - q.center.2) ^ 2
        = (q.a * ((p.1 - q.center.1) + (q.b / (2 * q.a)) * (p.2 - q.center.2)) ^ 2
          + (-q.disc / (4 * q.a)) * (p.2 - q.center.2) ^ 2) / K := by
      field_simp
    rw [hrw, div_eq_one_iff_eq hKne, hKdef]
    constructor <;> intro h <;> linarith

/-- **`disc < 0` and `det₃ = 0` means a single point** — the centre. -/
theorem zeroSet_eq_center_of_det3_eq_zero {q : Quad} (hd : q.disc < 0) (h3 : q.det3 = 0) :
    {p : Plane | q.eval p = 0} = {q.center} := by
  have ha : q.a ≠ 0 := Quad.a_ne_zero_of_disc_neg hd
  have hd0 : q.disc ≠ 0 := ne_of_lt hd
  have hc0 : q.eval q.center = 0 := by
    have := Quad.eval_center_mul_disc hd0
    rw [h3, neg_zero] at this
    rcases mul_eq_zero.1 this with h | h
    · exact h
    · exact absurd h hd0
  have hdd : 0 < -q.disc / 4 := by linarith
  have haa : 0 < q.a ^ 2 := by positivity
  ext p
  simp only [Set.mem_ofPred_eq, Set.mem_singleton_iff]
  rw [Quad.eval_eq_center_form hd0 ha p, hc0, add_zero]
  constructor
  · intro h
    set U : ℝ := (p.1 - q.center.1) + (q.b / (2 * q.a)) * (p.2 - q.center.2) with hU
    set V : ℝ := p.2 - q.center.2 with hV
    have key := mul_a_form (a := q.a) (Δ := q.disc) (U := U) (V := V) ha
    rw [h, mul_zero] at key
    have hVne : V = 0 := by
      by_contra hne
      have hVpos : 0 < V ^ 2 := lt_of_le_of_ne (sq_nonneg V) (Ne.symm (pow_ne_zero 2 hne))
      nlinarith [sq_nonneg U]
    have hUne : U = 0 := by
      by_contra hne
      have hUpos : 0 < U ^ 2 := lt_of_le_of_ne (sq_nonneg U) (Ne.symm (pow_ne_zero 2 hne))
      rw [hVne] at key
      nlinarith
    refine Prod.ext_iff.2 ⟨?_, ?_⟩
    · rw [hU, hVne, mul_zero, add_zero] at hUne
      linarith
    · rw [hV] at hVne
      linarith
  · intro h
    rw [h]
    ring

/-- **`disc < 0` and `a·det₃ > 0` means empty** — the imaginary ellipse. -/
theorem zeroSet_eq_empty_of_disc_neg {q : Quad} (hd : q.disc < 0) (h3 : 0 < q.a * q.det3) :
    {p : Plane | q.eval p = 0} = ∅ := by
  have ha : q.a ≠ 0 := Quad.a_ne_zero_of_disc_neg hd
  have hd0 : q.disc ≠ 0 := ne_of_lt hd
  have hdd : 0 < -q.disc / 4 := by linarith
  have haa : 0 < q.a ^ 2 := by positivity
  set K : ℝ := -q.eval q.center with hKdef
  have hKd : K * q.disc = q.det3 := by
    have := Quad.eval_center_mul_disc hd0
    rw [hKdef]
    linarith [this]
  have hKa : q.a * K < 0 := by
    rcases lt_trichotomy (q.a * K) 0 with h | h | h
    · exact h
    · exfalso
      have hK0 : K = 0 := by
        rcases mul_eq_zero.1 h with h' | h'
        · exact absurd h' ha
        · exact h'
      rw [hK0, zero_mul] at hKd
      rw [← hKd, mul_zero] at h3
      exact lt_irrefl 0 h3
    · exfalso
      have hprod : (q.a * K) * q.disc < 0 := mul_neg_of_pos_of_neg h hd
      have heq : q.a * q.det3 = (q.a * K) * q.disc := by rw [← hKd]; ring
      linarith
  ext p
  simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false]
  rw [Quad.eval_eq_center_form hd0 ha p]
  intro h
  set U : ℝ := (p.1 - q.center.1) + (q.b / (2 * q.a)) * (p.2 - q.center.2) with hU
  set V : ℝ := p.2 - q.center.2 with hV
  have hval : q.a * U ^ 2 + (-q.disc / (4 * q.a)) * V ^ 2 = K := by
    rw [hKdef]; linarith
  have key := mul_a_form (a := q.a) (Δ := q.disc) (U := U) (V := V) ha
  rw [hval] at key
  nlinarith [sq_nonneg U, sq_nonneg V]

/-! ### Sanity

`CLAUDE.md` → Verification, point 3. The documents' own elliptical edge, which
`Conic.lean` only knew was "elliptic and non-degenerate", is now an ellipse. -/

namespace Sanity

/-- `U₃|U₆` of `../doc/QuaConExample.md` §3.3 really is an ellipse. -/
theorem U3U6_isEllipse : IsEllipse {p : Plane | U3U6.eval p = 0} := by
  refine isEllipse_zeroSet ?_ ?_
  · norm_num [Quad.disc, U3U6]
  · norm_num [Quad.det3, U3U6]

end Sanity

end QuaConProof

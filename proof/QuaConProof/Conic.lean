/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.Quad

/-!
# Conics in the plane, and their discriminant classification

A **conic** is the zero set of a nonzero real quadratic in two variables. This is
the notion of "edge" the main theorem uses: `PROJECT_PLAN.md` §0.8 asks that
every multi-active cell of the conjugate's subdivision lie inside a conic, and
§0.7 asks that conics be classified into line / parabola / ellipse / hyperbola
and the degenerate cases.

**A line is a conic.** That is not a loophole, it is the point: the informal
Theorem 3 (`../CONJ_FIELD_PROOF.md`) shows that *adjacent* pieces of the input
always contribute degenerate edges — pairs of lines — while non-adjacent pieces
can contribute genuine ellipses and hyperbolas. Both are legal edges of a
`QuaCon`.

## Degeneracy is algebraic, never metric

The two invariants are

* `Quad.disc q = b² - 4ac`, which fixes the *type* (elliptic, parabolic,
  hyperbolic), and
* `Quad.det3 q`, four times the determinant of the conic's `3×3` matrix, which
  is zero exactly on the *degenerate* conics.

Only the **vanishing** of `det3` is independent of how the equation is scaled:
`det3` is cubic in the coefficients (`Quad.det3_smul`), so multiplying the whole
equation by `-1` flips its sign while `disc`, being quadratic, is unchanged
(`Quad.disc_smul`). The one place the sign is used is telling a real ellipse from
an empty one, where the scale-invariant test is `(a + c) * det3 < 0`.

`det3 = 0` is the whole definition of degenerate. It is emphatically **not** a
statement about curvature or about how straight the curve looks: `DECISIONS.md`,
2026-08-22, records two edges of a real example that are drawn as visually
straight segments and are nonetheless genuine hyperbolas, with `det3` of
absolute value greater than `2 × 10¹¹`. Their arcs merely sit far out on the
hyperbola where it hugs an asymptote. The `Sanity` section below pins exactly
those two, so that a classification predicate which confused the two notions
could not pass.

## Main results

* `IsConic` — the definition.
* `IsConic.ne_univ` — a conic is never the whole plane. This is what makes the
  notion non-vacuous, and it rests on `Quad.eq_zero_of_eval_eq_zero`.
* `isConic_eqLocus` — where two *distinct* coefficient vectors agree is a conic.
  This is the lemma the main theorem consumes.
* `Quad.disc_translate`, `Quad.det3_translate` — both invariants are unchanged
  by a translation of the plane, so the classification does not depend on where
  the origin is put.
-/

namespace QuaConProof

open Quad

/-- A subset of the plane is a **conic** when it is the zero set of a nonzero
real quadratic in two variables.

Degenerate conics are included on purpose: a line, a pair of lines, a single
point and the empty set all satisfy this, and all of them occur as edges of the
conjugate of a `QuaPol`. -/
def IsConic (C : Set Plane) : Prop :=
  ∃ q : Quad, q ≠ 0 ∧ C = {s : Plane | q.eval s = 0}

namespace Quad

/-- The discriminant `b² - 4ac`. Its **sign** fixes the type of the conic:
negative elliptic, zero parabolic, positive hyperbolic. -/
def disc (q : Quad) : ℝ := q.b ^ 2 - 4 * q.a * q.c

/-- Four times the determinant of the conic's `3×3` matrix

`⎡ a   b/2  d/2 ⎤`
`⎢ b/2  c   e/2 ⎥`
`⎣ d/2 e/2   f  ⎦`

The factor `4` clears the halves, so that a quadratic with integer coefficients
has an integer `det3`; only the vanishing and the sign of this quantity are ever
used, and the factor `4 > 0` disturbs neither. See `det3_eq_four_mul_det`. -/
def det3 (q : Quad) : ℝ :=
  4 * q.a * q.c * q.f + q.b * q.d * q.e - q.a * q.e ^ 2 - q.c * q.d ^ 2 - q.f * q.b ^ 2

/-- `det3` really is four times the `3×3` determinant, expanded along the first
row. This pins the definition against the textbook one. -/
lemma det3_eq_four_mul_det (q : Quad) :
    q.det3 = 4 * (q.a * (q.c * q.f - (q.e / 2) * (q.e / 2))
                - (q.b / 2) * ((q.b / 2) * q.f - (q.e / 2) * (q.d / 2))
                + (q.d / 2) * ((q.b / 2) * (q.e / 2) - q.c * (q.d / 2))) := by
  simp only [det3]; ring

/-- A conic is **degenerate** when its `3×3` determinant vanishes. -/
def IsDegenerate (q : Quad) : Prop := q.det3 = 0

/-- Neither invariant sees a translation of the plane. Since every conic can be
recentred, this is what makes `disc` and `det3` legitimate classifiers rather
than artefacts of the coordinate origin. -/
@[simp] lemma disc_translate (q : Quad) (v : Plane) : (q.translate v).disc = q.disc := by
  simp only [disc, translate_a, translate_b, translate_c]

@[simp] lemma det3_translate (q : Quad) (v : Plane) : (q.translate v).det3 = q.det3 := by
  simp only [det3, translate, eval]; ring

/-- Scaling a quadratic by `t` scales `disc` by `t²`, so it cannot change sign. -/
lemma disc_smul (q : Quad) (t : ℝ) :
    (⟨t * q.a, t * q.b, t * q.c, t * q.d, t * q.e, t * q.f⟩ : Quad).disc = t ^ 2 * q.disc := by
  simp only [disc]; ring

/-- Scaling a quadratic by `t` scales `det3` by `t³`, so it cannot change whether
it vanishes. -/
lemma det3_smul (q : Quad) (t : ℝ) :
    (⟨t * q.a, t * q.b, t * q.c, t * q.d, t * q.e, t * q.f⟩ : Quad).det3 = t ^ 3 * q.det3 := by
  simp only [det3]; ring

end Quad

/-! ### Basic properties of conics -/

/-- The zero set of a nonzero quadratic is a conic. -/
lemma isConic_zeroSet {q : Quad} (hq : q ≠ 0) : IsConic {s : Plane | q.eval s = 0} :=
  ⟨q, hq, rfl⟩

/-- **A conic is never the whole plane.**

This is the non-vacuity of the definition, and it is exactly where the choice to
store coefficient vectors rather than functions pays: `q ≠ 0` is a statement
about six real numbers, and `Quad.eq_zero_of_eval_eq_zero` converts it into a
point of the plane off the conic. -/
theorem IsConic.ne_univ {C : Set Plane} (h : IsConic C) : C ≠ Set.univ := by
  obtain ⟨q, hq, rfl⟩ := h
  obtain ⟨s, hs⟩ := Quad.exists_eval_ne_zero hq
  intro huniv
  have : s ∈ {s : Plane | q.eval s = 0} := huniv ▸ Set.mem_univ s
  exact hs this

/-- **Where two distinct quadratics agree is a conic.**

The lemma the main theorem consumes: a cell of the conjugate's subdivision on
which two candidate quadratics are simultaneously active is contained in the set
where they are equal, and that set is a conic.

Distinctness is of the *coefficient vectors*. That is why the candidate family is
deduplicated — `PROJECT_PLAN.md` §0.7. -/
theorem isConic_eqLocus {p q : Quad} (hpq : p ≠ q) :
    IsConic {s : Plane | p.eval s = q.eval s} := by
  refine ⟨p - q, ?_, ?_⟩
  · rw [Ne, Quad.sub_eq_zero_iff]; exact hpq
  · ext s
    simp only [Set.mem_ofPred_eq, Quad.eval_sub, sub_eq_zero]

/-- Consequently two distinct quadratics disagree somewhere: their equality locus
is a proper subset of the plane. -/
theorem exists_eval_ne_of_ne {p q : Quad} (hpq : p ≠ q) : ∃ s : Plane, p.eval s ≠ q.eval s := by
  have h := (isConic_eqLocus hpq).ne_univ
  by_contra hc
  push Not at hc
  exact h (Set.eq_univ_of_forall hc)

/-! ### Sanity checks

`CLAUDE.md` → Verification, point 3. The last two entries are the ones that
matter most: they are the real edges from `../doc/QuaCon.svg` row 3 that *look*
like straight lines and are not. A predicate that classified by appearance, or
that mixed up `disc` and `det3`, would fail here. -/

namespace Sanity

open Quad

/-- The unit circle, `s₁² + s₂² = 1`, is a conic. -/
example : IsConic {s : Plane | s.1 ^ 2 + s.2 ^ 2 = 1} := by
  refine ⟨⟨1, 0, 1, 0, 0, -1⟩, ?_, ?_⟩
  · rw [Ne, Quad.eq_zero_iff]; norm_num
  · ext s; simp only [Set.mem_ofPred_eq, eval]; constructor <;> intro h <;> linarith

/-- The unit circle is elliptic (`disc < 0`) and non-degenerate (`det3 ≠ 0`). -/
example : (⟨1, 0, 1, 0, 0, -1⟩ : Quad).disc = -4 := by norm_num [disc]
example : (⟨1, 0, 1, 0, 0, -1⟩ : Quad).det3 = -4 := by norm_num [det3]

/-- The parabola `s₂ = s₁²`, i.e. `s₁² - s₂ = 0`: `disc = 0`, non-degenerate. -/
example : IsConic {s : Plane | s.2 = s.1 ^ 2} := by
  refine ⟨⟨1, 0, 0, 0, -1, 0⟩, ?_, ?_⟩
  · rw [Ne, Quad.eq_zero_iff]; norm_num
  · ext s; simp only [Set.mem_ofPred_eq, eval]; constructor <;> intro h <;> linarith
example : (⟨1, 0, 0, 0, -1, 0⟩ : Quad).disc = 0 := by norm_num [disc]
example : (⟨1, 0, 0, 0, -1, 0⟩ : Quad).det3 = -1 := by norm_num [det3]

/-- The rectangular hyperbola `s₁s₂ = 1`: `disc > 0`, non-degenerate. -/
example : IsConic {s : Plane | s.1 * s.2 = 1} := by
  refine ⟨⟨0, 1, 0, 0, 0, -1⟩, ?_, ?_⟩
  · rw [Ne, Quad.eq_zero_iff]; norm_num
  · ext s; simp only [Set.mem_ofPred_eq, eval]; constructor <;> intro h <;> linarith
example : (⟨0, 1, 0, 0, 0, -1⟩ : Quad).disc = 1 := by norm_num [disc]
example : (⟨0, 1, 0, 0, 0, -1⟩ : Quad).det3 = 1 := by norm_num [det3]

/-- The crossing line pair `s₁s₂ = 0` — hyperbolic type, but **degenerate**.
Contrast with `s₁s₂ = 1` just above: the same `disc`, and `det3` tells them
apart. -/
example : (⟨0, 1, 0, 0, 0, 0⟩ : Quad).disc = 1 := by norm_num [disc]
example : (⟨0, 1, 0, 0, 0, 0⟩ : Quad).det3 = 0 := by norm_num [det3]

/-- A single line `s₂ = 0`, the all-zero-quadratic-part degenerate case. -/
example : IsConic {s : Plane | s.2 = 0} := by
  refine ⟨⟨0, 0, 0, 0, 1, 0⟩, ?_, ?_⟩
  · rw [Ne, Quad.eq_zero_iff]; norm_num
  · ext s; simp only [Set.mem_ofPred_eq, eval]; constructor <;> intro h <;> linarith
example : (⟨0, 0, 0, 0, 1, 0⟩ : Quad).disc = 0 := by norm_num [disc]
example : (⟨0, 0, 0, 0, 1, 0⟩ : Quad).det3 = 0 := by norm_num [det3]

/-! #### The two traps

Both edges below are drawn as visually straight segments in `../doc/QuaCon.svg`
row 3, and both are genuine, non-degenerate hyperbolas. Their active arcs sit
`14.8` and `2.3` transverse semi-axes out from the centre, where the hyperbola
hugs an asymptote, giving a sagitta of `0.13%` and `0.86%` of the chord — under
one line width as drawn. See `DECISIONS.md`, 2026-08-22.

The coefficients come from the difference of two face quadratics of `f*` for the
five-piece example, cleared to integers. -/

/-- `U₁ | U₆`, the boundary between the interior branches of pieces 5 and 3:
`2731 s₁² - 4598 s₁s₂ - 2189 s₂² - 107420 s₁ - 9988 s₂ + 421996 = 0`.
Hyperbolic **and non-degenerate**, despite looking straight. -/
def U1U6 : Quad := ⟨2731, -4598, -2189, -107420, -9988, 421996⟩

example : U1U6.disc = 45054240 := by norm_num [disc, U1U6]
example : 0 < U1U6.disc := by norm_num [disc, U1U6]
example : U1U6.det3 = 4 * 260148962304 := by norm_num [det3, U1U6]
example : ¬ U1U6.IsDegenerate := by norm_num [IsDegenerate, det3, U1U6]

/-- `U₁ | U₂`, piece 1's outer-edge branch against piece 5's interior:
`5969 s₁² + 5390 s₁s₂ - 847 s₂² + 67532 s₁ - 22748 s₂ - 353236 = 0`.
Hyperbolic **and non-degenerate**. -/
def U1U2 : Quad := ⟨5969, 5390, -847, 67532, -22748, -353236⟩

example : U1U2.disc = 49275072 := by norm_num [disc, U1U2]
example : 0 < U1U2.disc := by norm_num [disc, U1U2]
example : U1U2.det3 = 4 * 2474882726976 := by norm_num [det3, U1U2]
example : ¬ U1U2.IsDegenerate := by norm_num [IsDegenerate, det3, U1U2]

/-- `U₃ | U₆`, an edge of the same subdivision that really is an **ellipse**:
`93 s₁² + 38 s₁s₂ + 39 s₂² - 6 s₁ - 482 s₂ - 1003 = 0`. This is the arc that
`QuaPar` cannot store, since `QuaPar` demands `disc = 0`. -/
def U3U6 : Quad := ⟨93, 38, 39, -6, -482, -1003⟩

example : U3U6.disc = -13064 := by norm_num [disc, U3U6]
example : U3U6.disc < 0 := by norm_num [disc, U3U6]
example : U3U6.det3 = 4 * (-8650208) := by norm_num [det3, U3U6]
example : ¬ U3U6.IsDegenerate := by norm_num [IsDegenerate, det3, U3U6]

/-- The ellipse is a **real** one, not the empty "imaginary ellipse": the test is
`(a + c) * det3 < 0`, which is the one classification criterion whose *sign*
(not merely whose vanishing) is used. See the note on sign conventions above. -/
example : (U3U6.a + U3U6.c) * U3U6.det3 < 0 := by norm_num [det3, U3U6]

/-! #### The genuinely degenerate edges

`U₂|U₄` and `U₁|U₅` are adjacent-piece edges, the case of the informal Theorem 3.
Both have `det3 = 0` exactly, and both are `QuaPar`-legal. They are the contrast
that gives the two traps above their force: degeneracy is common in this very
subdivision, and it is still not what `U₁|U₆` and `U₁|U₂` are. -/

/-- `U₂ | U₄`: hyperbolic type, but **degenerate** — a crossing pair of lines. -/
def U2U4 : Quad := ⟨39, -88, -41, 88, 82, -41⟩

example : U2U4.disc = 14140 := by norm_num [disc, U2U4]
example : 0 < U2U4.disc := by norm_num [disc, U2U4]
example : U2U4.IsDegenerate := by norm_num [IsDegenerate, det3, U2U4]

/-- `U₁ | U₅`: `disc = 0` **and** `det3 = 0` — parallel or repeated lines, the
flattest degenerate case. Note it is not a parabola: a parabola needs
`det3 ≠ 0`. -/
def U1U5 : Quad := ⟨961, -3410, 3025, -19220, 34100, 96100⟩

example : U1U5.disc = 0 := by norm_num [disc, U1U5]
example : U1U5.IsDegenerate := by norm_num [IsDegenerate, det3, U1U5]

/-! #### The parabola, which needs an INDEFINITE piece

No pair of the seven faces `U₁..U₇` is a parabola — all five pieces of that
example are convex, and for a convex piece the argmax moves continuously, so the
single-piece conjugate is polyhedral. A parabolic edge needs the argmax to *jump*,
which needs an indefinite `Q`.

Panel (3b) of `../doc/QuaCon.svg` supplies one: `q = x₁x₂` on the triangle
`(0,0), (1,1), (2,0)`, whose Hessian `[[0,1],[1,0]]` is indefinite. Its conjugate
has the edge `¼(s₁+s₂)² = 2s₁`, where the branch along the segment
`[(0,0),(1,1)]` meets the vertex branch at `(2,0)`. Cleared to integers that is
`s₁² + 2s₁s₂ + s₂² - 8s₁ = 0`. -/
def Parab : Quad := ⟨1, 2, 1, -8, 0, 0⟩

example : Parab.disc = 0 := by norm_num [disc, Parab]
example : Parab.det3 = -64 := by norm_num [det3, Parab]

/-- It is a **genuine** parabola: parabolic type *and* non-degenerate. Without the
second half it could have been a repeated line, which is the other `disc = 0`
case — compare `U1U5` above, which has `disc = 0` and *is* degenerate. -/
theorem Parab_isGenuineParabola : Parab.disc = 0 ∧ ¬ Parab.IsDegenerate := by
  constructor
  · norm_num [disc, Parab]
  · norm_num [IsDegenerate, det3, Parab]

/-- All four non-degenerate types are realised by the examples in this file:
ellipse `U3U6`, parabola `Parab`, hyperbola `U1U6`, and the degenerate line pair
`U2U4`. -/
theorem all_four_types_occur :
    U3U6.disc < 0 ∧ Parab.disc = 0 ∧ 0 < U1U6.disc ∧ U2U4.IsDegenerate := by
  refine ⟨by norm_num [disc, U3U6], by norm_num [disc, Parab],
          by norm_num [disc, U1U6], by norm_num [IsDegenerate, det3, U2U4]⟩

/-- Not every conic is `QuaPar`-legal, stated as a theorem about this example:
`U₃|U₆` has nonzero discriminant, so it is neither a parabola nor a line, and it
is non-degenerate, so it is not a line pair either. It is a genuine ellipse, and
`QuaPar` — which demands `disc = 0` of every edge — cannot store it. -/
theorem U3U6_not_quaParLegal : U3U6.disc ≠ 0 ∧ ¬ U3U6.IsDegenerate := by
  constructor
  · norm_num [disc, U3U6]
  · norm_num [IsDegenerate, det3, U3U6]

end Sanity

end QuaConProof

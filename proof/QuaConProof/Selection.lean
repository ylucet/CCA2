/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.Candidates

/-!
# Towards the selection lemma: the second-order expansion, and step S8

`PROJECT_PLAN.md` §0.6. This file builds the selection lemma from the bottom up.
It starts with **S8**, the step flagged as most likely to fail, so that nothing
is built on top of it before it is known to work.

## The one identity everything rests on

`psi_along_dir`: for a quadratic `q`, moving from `x` in the direction `d`,

`ψ(x + t·d) = ψ(x) + t·⟨s - ∇q(x), d⟩ - t²·(dᵀHd)/2`

**exactly** — no remainder, because `q` is a quadratic. Every case of the
selection lemma is this identity plus a scalar fact:

| case | what is fed in | what comes out |
|---|---|---|
| S5, first order | `x` maximises, both signs of `t` allowed | `⟨s - ∇q(x), d⟩ = 0` |
| S7, second order | the same | `dᵀHd ≥ 0` |
| S7, `α > 0` | first order | `x` is the stationary point, so `edgeBranch` |
| S7, `α = 0` | first order | `ψ` is constant along `d`, so `vertexBranch` |
| **S8** | `∇q(x) = s` and `Hd = 0` | `ψ` is **constant** along `x + ℝd` |

## What S8 is for

At a three-vertex face the first-order condition forces `∇q(x) = s` outright. If
`H` is nonsingular that pins `x` and gives `interiorBranch`. If `H` is singular,
`ψ` is constant along the kernel direction, so sliding along it costs nothing;
the barycentric coordinates are affine in the parameter, so at the first moment
one of them vanishes there is *another* maximiser sitting on a proper face. The
case then reduces to a smaller face, and the three-vertex-singular branch never
needs a candidate of its own.

`psi_const_along_kernel` below is exactly that "sliding costs nothing" step, and
`exists_kernel_of_hessDet_eq_zero` supplies the direction.
-/

namespace QuaConProof

namespace Quad

/-- The Hessian applied to a vector: `H d`, where `H = [[2a, b], [b, 2c]]`. -/
def hessApply (q : Quad) (d : Plane) : Plane :=
  (2 * q.a * d.1 + q.b * d.2, q.b * d.1 + 2 * q.c * d.2)

@[simp] lemma hessApply_fst (q : Quad) (d : Plane) :
    (q.hessApply d).1 = 2 * q.a * d.1 + q.b * d.2 := rfl

@[simp] lemma hessApply_snd (q : Quad) (d : Plane) :
    (q.hessApply d).2 = q.b * d.1 + 2 * q.c * d.2 := rfl

/-- The curvature along `d` is `⟨Hd, d⟩`. -/
lemma alongCurv_eq_dot_hessApply (q : Quad) (d : Plane) :
    q.alongCurv d = dot (q.hessApply d) d := by
  simp only [alongCurv, dot, hessApply_fst, hessApply_snd]; ring

/-- The gradient of a quadratic is affine: `∇q(x + t·d) = ∇q(x) + t·Hd`. -/
lemma gradAt_add_smul (q : Quad) (x d : Plane) (t : ℝ) :
    q.gradAt (x + t • d) = q.gradAt x + t • q.hessApply d := by
  rw [Prod.ext_iff]
  simp only [gradAt, hessApply, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
    smul_eq_mul]
  constructor <;> ring

/-- A vector is in the kernel of the Hessian iff `H d = 0`; then the curvature
along it vanishes too. -/
lemma alongCurv_eq_zero_of_hessApply_eq_zero {q : Quad} {d : Plane}
    (h : q.hessApply d = 0) : q.alongCurv d = 0 := by
  rw [alongCurv_eq_dot_hessApply, h]
  simp [dot]

end Quad

/-! ### The exact second-order expansion -/

/-- **The identity everything rests on.** For a quadratic there is no remainder:
`ψ(x + t·d) = ψ(x) + t·⟨s - ∇q(x), d⟩ - t²(dᵀHd)/2`. -/
theorem psi_along_dir (q : Quad) (s x d : Plane) (t : ℝ) :
    psi q s (x + t • d)
      = psi q s x + t * (dot s d - dot (q.gradAt x) d) - t ^ 2 * q.alongCurv d / 2 := by
  simp only [psi, dot, Quad.gradAt, Quad.alongCurv, Quad.eval, Prod.fst_add, Prod.snd_add,
    Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
  ring

/-- `psi_along_line` is the special case `x = v`, `d = w - v`. -/
lemma psi_along_dir_eq_psi_along_line (q : Quad) (v w s : Plane) (t : ℝ) :
    psi q s (v + t • (w - v))
      = psi q s v + t * edgeSlope q v w s - t ^ 2 * edgeCurv q v w / 2 :=
  psi_along_line q v w s t

/-! ### S8, part one: a singular Hessian has a nonzero kernel vector -/

/-- **A singular Hessian kills some nonzero direction.**

In the plane this is explicit: `(b, -2a)` works whenever `(a, b) ≠ 0`, since
`H(b, -2a) = (0, b² - 4ac)` and `4ac - b² = 0`; and if `a = b = 0` the Hessian is
`[[0,0],[0,2c]]`, which kills `(1, 0)`. -/
theorem exists_kernel_of_hessDet_eq_zero {q : Quad} (h : q.hessDet = 0) :
    ∃ d : Plane, d ≠ 0 ∧ q.hessApply d = 0 := by
  by_cases hab : q.a = 0 ∧ q.b = 0
  · refine ⟨(1, 0), ?_, ?_⟩
    · simp [Prod.ext_iff]
    · simp only [Quad.hessApply, hab.1, hab.2, Prod.mk_eq_zero]
      constructor <;> ring
  · refine ⟨(q.b, -2 * q.a), ?_, ?_⟩
    · rw [Ne, Prod.ext_iff]
      simp only [Prod.fst_zero, Prod.snd_zero, not_and]
      intro hb
      have ha : q.a ≠ 0 := by
        intro ha; exact hab ⟨ha, hb⟩
      intro hcon
      exact ha (by linarith)
    · have hd : 4 * q.a * q.c - q.b ^ 2 = 0 := h
      simp only [Quad.hessApply, Prod.mk_eq_zero]
      constructor <;> nlinarith [hd]

/-! ### S8, part two: sliding along the kernel costs nothing -/

/-- **Step S8.** At a stationary point of a quadratic with singular Hessian, `ψ`
is *constant* along any kernel direction.

Both correction terms of `psi_along_dir` vanish: the first-order term because
`∇q(x) = s`, and the second-order term because `Hd = 0` forces `dᵀHd = 0`.

This is what lets the three-vertex singular case reduce to a smaller face: the
whole line `x + ℝd` consists of maximisers, so one may slide along it until a
barycentric coordinate hits zero without losing maximality. -/
theorem psi_const_along_kernel {q : Quad} {s x d : Plane}
    (hstat : q.gradAt x = s) (hker : q.hessApply d = 0) (t : ℝ) :
    psi q s (x + t • d) = psi q s x := by
  rw [psi_along_dir, hstat, Quad.alongCurv_eq_zero_of_hessApply_eq_zero hker]
  ring

/-- The same statement with the hypothesis in the form it will actually arrive
in: the Hessian is singular, so *some* nonzero direction has this property. -/
theorem exists_dir_psi_const {q : Quad} {s x : Plane} (hstat : q.gradAt x = s)
    (hsing : q.hessDet = 0) :
    ∃ d : Plane, d ≠ 0 ∧ ∀ t : ℝ, psi q s (x + t • d) = psi q s x := by
  obtain ⟨d, hd0, hdk⟩ := exists_kernel_of_hessDet_eq_zero hsing
  exact ⟨d, hd0, fun t => psi_const_along_kernel hstat hdk t⟩


/-! ### S2: the per-piece supremum is attained

`QuaPiece.T` is compact by construction (`QuaPiece.isCompact_T`), and `psi` is a
polynomial, hence continuous. The extreme value theorem does the rest. This is
the step the V-representation was chosen to make free. -/

/-- `psi q s` is continuous: it is a polynomial in the two coordinates. -/
lemma continuous_psi (q : Quad) (s : Plane) : Continuous (psi q s) := by
  unfold psi dot Quad.eval
  fun_prop

/-- **S2.** The supremum of `psi` over a piece is attained. -/
theorem exists_isMaxOn_piece (p : QuaPiece) (s : Plane) :
    ∃ x ∈ p.T, IsMaxOn (psi p.q s) p.T x :=
  p.isCompact_T.exists_isMaxOn p.T_nonempty (continuous_psi p.q s).continuousOn

/-! ### S1: the infimum over pieces is attained, and the conjugate is finite -/

/-- A finite infimum in `EReal` that is not `⊤` is attained. -/
private lemma inf_attained {α : Type*} [DecidableEq α] (s : Finset α) (g : α → EReal)
    (h : s.inf g ≠ ⊤) :
    ∃ a ∈ s, s.inf g = g a := by
  induction s using Finset.induction_on with
  | empty => simp at h
  | insert a s _ ih =>
      rw [Finset.inf_insert] at h ⊢
      rcases le_total (g a) (s.inf g) with hle | hle
      · exact ⟨a, Finset.mem_insert_self _ _, inf_eq_left.2 hle⟩
      · rw [inf_eq_right.2 hle] at h ⊢
        obtain ⟨b, hb, hbe⟩ := ih h
        exact ⟨b, Finset.mem_insert_of_mem hb, hbe⟩

/-- Where a `QuaPol` is finite, some piece realises its value. -/
theorem exists_piece_eq_eval {f : QuaPol} {x : Plane} (h : f.eval x ≠ ⊤) :
    ∃ p ∈ f.pieces, x ∈ p.T ∧ f.eval x = ((p.q.eval x : ℝ) : EReal) := by
  obtain ⟨p, hp, hpe⟩ := inf_attained f.pieces _ (by rwa [QuaPol.eval] at h)
  rw [← QuaPol.eval] at hpe
  by_cases hx : x ∈ p.T
  · exact ⟨p, hp, hx, by rw [hpe]; simp [hx]⟩
  · rw [hpe] at h; simp [hx] at h

/-- **Stage 1: the conjugate is finite everywhere.**

Every piece is compact, so its supremum is a real number, and there are finitely
many pieces. Hence `dom f* = ℝ²` and the `⊤` cell of `PROJECT_PLAN.md` §0.5 is
empty — it only becomes real once unbounded pieces are admitted (Phase 7).

This is stronger than the plan anticipated, and it makes `selection` hold
unconditionally. -/
theorem conj_ne_top (f : QuaPol) (s : Plane) : f.conj s ≠ ⊤ := by
  have hmax : ∀ p : QuaPiece, ∃ Mp : ℝ, ∀ x ∈ p.T, psi p.q s x ≤ Mp := by
    intro p
    obtain ⟨x, _, hxmax⟩ := exists_isMaxOn_piece p s
    exact ⟨psi p.q s x, fun y hy => hxmax hy⟩
  choose Mp hMp using hmax
  set M : ℝ := f.pieces.sup' f.pieces_nonempty Mp with hM
  have hbound : ∀ x : Plane, ((dot s x : ℝ) : EReal) - f.eval x ≤ ((M : ℝ) : EReal) := by
    intro x
    by_cases hx : f.eval x = ⊤
    · rw [hx, EReal.sub_top]; exact bot_le
    · obtain ⟨p, hp, hxT, hev⟩ := exists_piece_eq_eval hx
      rw [hev, ← EReal.coe_sub]
      exact EReal.coe_le_coe_iff.2
        (le_trans (hMp p x hxT) (Finset.le_sup' Mp hp))
  intro htop
  have hle := iSup_le hbound
  rw [← QuaPol.conj_def, htop, top_le_iff] at hle
  exact EReal.coe_ne_top M hle

/-- **S1 + S2 combined.** The conjugate is attained: some piece and some point of
that piece realise it, and that point maximises `psi` over the piece. -/
theorem exists_maximiser (f : QuaPol) (s : Plane) :
    ∃ p ∈ f.pieces, ∃ x ∈ p.T, IsMaxOn (psi p.q s) p.T x ∧
      ((psi p.q s x : ℝ) : EReal) = f.conj s := by
  -- pick a maximiser on each piece, then the best piece
  have hmax : ∀ p : QuaPiece, ∃ x, x ∈ p.T ∧ IsMaxOn (psi p.q s) p.T x := by
    intro p; obtain ⟨x, hx, hxm⟩ := exists_isMaxOn_piece p s; exact ⟨x, hx, hxm⟩
  choose xm hxmT hxmMax using hmax
  obtain ⟨p, hp, hpsup⟩ :=
    Finset.exists_mem_eq_sup' f.pieces_nonempty (fun p => psi p.q s (xm p))
  refine ⟨p, hp, xm p, hxmT p, hxmMax p, le_antisymm ?_ ?_⟩
  · -- `psi ≤ f*`: the value is attained, at `xm p` on the piece `p`
    rw [QuaPol.conj_def]
    refine le_iSup_of_le (xm p) ?_
    have hle : f.eval (xm p) ≤ ((p.q.eval (xm p) : ℝ) : EReal) :=
      QuaPol.eval_le_of_mem hp (hxmT p)
    calc ((psi p.q s (xm p) : ℝ) : EReal)
        = ((dot s (xm p) : ℝ) : EReal) - ((p.q.eval (xm p) : ℝ) : EReal) := by
          rw [← EReal.coe_sub]; rfl
      _ ≤ ((dot s (xm p) : ℝ) : EReal) - f.eval (xm p) :=
          EReal.sub_le_sub le_rfl hle
  · -- `f* ≤ psi`: no term of the supremum exceeds it
    rw [QuaPol.conj_def]
    refine iSup_le fun y => ?_
    by_cases hy : f.eval y = ⊤
    · rw [hy, EReal.sub_top]; exact bot_le
    · obtain ⟨p', hp', hyT, hev⟩ := exists_piece_eq_eval hy
      rw [hev, ← EReal.coe_sub]
      refine EReal.coe_le_coe_iff.2 (le_trans (hxmMax p' hyT) ?_)
      rw [← hpsup]
      exact Finset.le_sup' (fun p => psi p.q s (xm p)) hp'

/-! ### The scalar facts S5 and S7 will need -/

/-- **First-order condition.** If `t·A - t²·B/2 ≤ 0` for every `t` in a
neighbourhood of `0` — both signs — then `A = 0`.

This is the scalar content of S5: a maximum in the *relative interior* of a face
can be approached from both directions, so the linear term must vanish. -/
theorem eq_zero_of_forall_small {A B ε : ℝ} (hε : 0 < ε)
    (h : ∀ t : ℝ, |t| < ε → t * A - t ^ 2 * B / 2 ≤ 0) : A = 0 := by
  -- stepping forwards bounds `A` above, stepping backwards bounds it below
  have fwd : ∀ u : ℝ, 0 < u → u < ε → A ≤ u * B / 2 := by
    intro u hu huε
    have hb := h u (by rw [abs_of_pos hu]; exact huε)
    have : u * A ≤ u * (u * B / 2) := by nlinarith
    exact le_of_mul_le_mul_left this hu
  have bwd : ∀ u : ℝ, 0 < u → u < ε → -(u * B / 2) ≤ A := by
    intro u hu huε
    have hb := h (-u) (by rw [abs_of_neg (by linarith : -u < 0)]; simpa using huε)
    have : u * -(u * B / 2) ≤ u * A := by nlinarith
    exact le_of_mul_le_mul_left this hu
  -- hence `|A| ≤ u|B|/2` for every small `u > 0`, which forces `A = 0`
  have habs : ∀ u : ℝ, 0 < u → u < ε → |A| ≤ u * |B| / 2 := by
    intro u hu huε
    have h1 : A ≤ u * B / 2 := fwd u hu huε
    have h2 : -(u * B / 2) ≤ A := bwd u hu huε
    have h3 : u * B / 2 ≤ u * |B| / 2 := by
      have : B ≤ |B| := le_abs_self B
      nlinarith
    have h4 : -(u * |B| / 2) ≤ -(u * B / 2) := by
      have : -|B| ≤ B := neg_abs_le B
      nlinarith
    exact abs_le.2 ⟨by linarith, by linarith⟩
  by_contra hA
  have hApos : 0 < |A| := abs_pos.2 hA
  have hB1 : (0 : ℝ) < |B| + 1 := by positivity
  set u := min (ε / 2) (|A| / (|B| + 1)) with hu_def
  have hu0 : 0 < u := lt_min (by linarith) (by positivity)
  have huε : u < ε := lt_of_le_of_lt (min_le_left _ _) (by linarith)
  have hbound := habs u hu0 huε
  have hu2 : u ≤ |A| / (|B| + 1) := min_le_right _ _
  have : u * |B| ≤ |A| / (|B| + 1) * |B| :=
    mul_le_mul_of_nonneg_right hu2 (abs_nonneg B)
  have hlt : |A| / (|B| + 1) * |B| < |A| := by
    rw [div_mul_eq_mul_div, div_lt_iff₀ hB1]
    nlinarith [abs_nonneg B]
  linarith


/-! ### The remaining gap, and the selection lemma

Everything above is proved. What is left is exactly the case analysis S3–S9 for a
**single piece**, isolated here so that the `sorry` sits on a precise
mathematical statement rather than on the top-level theorem. -/

/-- **The core of the selection lemma (S3–S9), for one piece.**

If `x` maximises `ψ` over the piece, then one of the piece's candidate branches
takes the value `ψ(x)` at `s`.

Route, with the pieces already in hand marked:

* **S3** take `W ⊆ verts` of minimal cardinality with `x ∈ convexHull ↑W`
  (`Caratheodory.minCardFinsetOfMemConvexHull`, already in mathlib, together with
  `affineIndependent_minCardFinsetOfMemConvexHull`). Minimality forces every
  barycentric coordinate of `x` to be strictly positive, so `x` is in the
  *relative interior* of the simplex on `W`, and `W.card ≤ 3` in the plane.
* **S4** `x` maximises over `convexHull ↑W` too, since that is a subset.
* **S5** for every direction `d` of `affineSpan ℝ W`, `x ± t·d` stays in the
  simplex for small `t`; feeding `psi_along_dir` into `eq_zero_of_forall_small`
  ✓ gives `⟨s - ∇q(x), d⟩ = 0`.
* **S6** `W.card = 1`: `x` is a vertex, and `vertexBranch_eval` ✓ finishes.
* **S7** `W.card = 2`: S5 along `w - v` gives `L = t·α`. The second-order part of
  `psi_along_dir` ✓ gives `α ≥ 0`. If `α > 0` then `x` is the stationary point and
  `edgeBranch_eval` ✓ finishes, with `edgeBranch ∈ branches` because
  `edgePairs` keeps exactly the pairs with `0 < edgeCurv`. If `α = 0` then `L = 0`
  and `ψ` is constant along the segment, so the value is `vertexBranch` at `v`.
* **S8** `W.card = 3`: affine independence in the plane makes `affineSpan ℝ W = ⊤`,
  so S5 gives `∇q(x) = s` outright. If `hessDet ≠ 0` then `gradAt_interiorPoint` ✓
  and `interiorBranch_eval` ✓ finish. If `hessDet = 0`, `exists_dir_psi_const` ✓
  supplies a direction along which `ψ` is constant; the barycentric coordinates
  are affine in the parameter, so at the first parameter where one vanishes there
  is another maximiser on a proper face, and the induction below applies to it.
* **S9** strong induction on `W.card`, the value being the same at every
  maximiser.

The ticked steps are proved in this file or in `Candidates.lean`. What remains is
the barycentric bookkeeping of S3, S5 and the descent in S8, and the induction
that ties them together. -/
theorem exists_branch_eq_max (p : QuaPiece) (s : Plane) {x : Plane} (hx : x ∈ p.T)
    (hmax : IsMaxOn (psi p.q s) p.T x) :
    ∃ b ∈ p.branches, b.eval s = psi p.q s x := by
  sorry

/-- **Selection.** Some candidate quadratic attains the conjugate at every point.

At Stage 1 this is unconditional: `conj_ne_top` shows the conjugate is finite
everywhere, because every piece is compact. -/
theorem selection (f : QuaPol) (s : Plane) :
    ∃ q ∈ cand f, ((q.eval s : ℝ) : EReal) = f.conj s := by
  obtain ⟨p, hp, x, hxT, hxmax, hval⟩ := exists_maximiser f s
  obtain ⟨b, hb, hbe⟩ := exists_branch_eq_max p s hxT hxmax
  exact ⟨b, mem_cand hp hb, by rw [hbe]; exact hval⟩

/-! ### Sanity checks for S8 -/

namespace Sanity

/-- `q(x) = x₁²` has a singular Hessian. -/
example : (⟨1, 0, 0, 0, 0, 0⟩ : Quad).hessDet = 0 := by norm_num [Quad.hessDet]

/-- Its kernel direction is the `x₂` axis. -/
example : (⟨1, 0, 0, 0, 0, 0⟩ : Quad).hessApply (0, 1) = 0 := by
  simp only [Quad.hessApply, Prod.mk_eq_zero]; norm_num

/-- **S8 in a concrete instance.** For `q(x) = x₁²` and `s = (2, 0)`, the point
`(1, 0)` is stationary, and `ψ` is constant along the whole vertical line through
it — value `1` everywhere. A candidate list that expected a unique interior
maximiser here would be wrong; S8 is what handles it. -/
example (t : ℝ) :
    psi (⟨1, 0, 0, 0, 0, 0⟩ : Quad) (2, 0) ((1, 0) + t • (0, 1))
      = psi (⟨1, 0, 0, 0, 0, 0⟩ : Quad) (2, 0) (1, 0) := by
  refine psi_const_along_kernel ?_ ?_ t
  · simp only [Quad.gradAt, Prod.mk.injEq]; norm_num
  · simp only [Quad.hessApply, Prod.mk_eq_zero]; norm_num

/-- and that common value is `1`, by hand: `⟨(2,0),(1,t)⟩ - 1² = 2 - 1`. -/
example : psi (⟨1, 0, 0, 0, 0, 0⟩ : Quad) (2, 0) (1, 0) = 1 := by
  norm_num [psi, dot, Quad.eval]

end Sanity

end QuaConProof

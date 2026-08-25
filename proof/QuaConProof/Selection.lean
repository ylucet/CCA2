/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.Candidates
import QuaConProof.Bary

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

/-- **S2.** On a piece with no recession directions the supremum of `psi` is
attained, by compactness. With rays this needs a real theorem — see
`QuaPol.Attained` below. -/
theorem exists_isMaxOn_piece {p : QuaPiece} (hr : p.rays = ∅) (s : Plane) :
    ∃ x ∈ p.T, IsMaxOn (psi p.q s) p.T x :=
  (QuaPiece.isCompact_T hr).exists_isMaxOn p.T_nonempty (continuous_psi p.q s).continuousOn

/-- **Attainment at `s`**: every piece's objective attains its supremum. This is
the hypothesis every step below actually consumes; boundedness of the pieces is
only one way to get it (`attained_of_bounded`), and Phase 7 supplies the other. -/
def QuaPol.Attained (f : QuaPol) (s : Plane) : Prop :=
  ∀ p ∈ f.pieces, ∃ x ∈ p.T, IsMaxOn (psi p.q s) p.T x

/-- Bounded pieces are compact, so their suprema are attained. -/
theorem attained_of_bounded {f : QuaPol} (hb : f.Bounded) (s : Plane) : f.Attained s :=
  fun p hp => exists_isMaxOn_piece (hb p hp) s

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
theorem conj_ne_top {f : QuaPol} (hb : f.Bounded) (s : Plane) : f.conj s ≠ ⊤ := by
  have hmax : ∀ p : QuaPiece, ∃ Mp : ℝ, p ∈ f.pieces → ∀ x ∈ p.T, psi p.q s x ≤ Mp := by
    intro p
    by_cases hp : p ∈ f.pieces
    · obtain ⟨x, _, hxmax⟩ := exists_isMaxOn_piece (hb p hp) s
      exact ⟨psi p.q s x, fun _ y hy => hxmax hy⟩
    · exact ⟨0, fun h => absurd h hp⟩
  choose Mp hMp using hmax
  set M : ℝ := f.pieces.sup' f.pieces_nonempty Mp with hM
  have hbound : ∀ x : Plane, ((dot s x : ℝ) : EReal) - f.eval x ≤ ((M : ℝ) : EReal) := by
    intro x
    by_cases hx : f.eval x = ⊤
    · rw [hx, EReal.sub_top]; exact bot_le
    · obtain ⟨p, hp, hxT, hev⟩ := exists_piece_eq_eval hx
      rw [hev, ← EReal.coe_sub]
      exact EReal.coe_le_coe_iff.2
        (le_trans (hMp p hp x hxT) (Finset.le_sup' Mp hp))
  intro htop
  have hle := iSup_le hbound
  rw [← QuaPol.conj_def, htop, top_le_iff] at hle
  exact EReal.coe_ne_top M hle

/-- **S1 + S2 combined.** The conjugate is attained: some piece and some point of
that piece realise it, and that point maximises `psi` over the piece. -/
theorem exists_maximiser {f : QuaPol} {s : Plane} (hA : f.Attained s) :
    ∃ p ∈ f.pieces, ∃ x ∈ p.T, IsMaxOn (psi p.q s) p.T x ∧
      ((psi p.q s x : ℝ) : EReal) = f.conj s := by
  -- pick a maximiser on each piece, then the best piece
  have hmax : ∀ p : QuaPiece, ∃ x, p ∈ f.pieces → x ∈ p.T ∧ IsMaxOn (psi p.q s) p.T x := by
    intro p
    by_cases hp : p ∈ f.pieces
    · obtain ⟨x, hx, hxm⟩ := hA p hp; exact ⟨x, fun _ => ⟨hx, hxm⟩⟩
    · exact ⟨0, fun h => absurd h hp⟩
  choose xm hxm using hmax
  have hxmT : ∀ p ∈ f.pieces, xm p ∈ p.T := fun p hp => (hxm p hp).1
  have hxmMax : ∀ p ∈ f.pieces, IsMaxOn (psi p.q s) p.T (xm p) := fun p hp => (hxm p hp).2
  obtain ⟨p, hp, hpsup⟩ :=
    Finset.exists_mem_eq_sup' f.pieces_nonempty (fun p => psi p.q s (xm p))
  refine ⟨p, hp, xm p, hxmT p hp, hxmMax p hp, le_antisymm ?_ ?_⟩
  · -- `psi ≤ f*`: the value is attained, at `xm p` on the piece `p`
    rw [QuaPol.conj_def]
    refine le_iSup_of_le (xm p) ?_
    have hle : f.eval (xm p) ≤ ((p.q.eval (xm p) : ℝ) : EReal) :=
      QuaPol.eval_le_of_mem hp (hxmT p hp)
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
      refine EReal.coe_le_coe_iff.2 (le_trans (hxmMax p' hp' hyT) ?_)
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


/-! ### Small algebraic helpers -/

lemma dot_sub_left (a b c : Plane) : dot (a - b) c = dot a c - dot b c := by
  simp only [dot, Prod.fst_sub, Prod.snd_sub]; ring

/-- **The stationary point is unique when the Hessian is nonsingular.** -/
theorem interiorPoint_unique {q : Quad} {s y : Plane} (hD : q.hessDet ≠ 0)
    (hy : q.gradAt y = s) : y = interiorPoint q s := by
  have h1 : 2 * q.a * y.1 + q.b * y.2 + q.d = s.1 := congrArg Prod.fst hy
  have h2 : q.b * y.1 + 2 * q.c * y.2 + q.e = s.2 := congrArg Prod.snd hy
  have hD' : 4 * q.a * q.c - q.b ^ 2 ≠ 0 := hD
  refine Prod.ext_iff.2 ⟨?_, ?_⟩
  · simp only [interiorPoint, Quad.hessDet]
    rw [eq_div_iff hD']
    linear_combination (2 * q.c) * h1 - q.b * h2
  · simp only [interiorPoint, Quad.hessDet]
    rw [eq_div_iff hD']
    linear_combination (2 * q.a) * h2 - q.b * h1

/-! ### S5 and S7: the first- and second-order conditions -/

/-- **The first- and second-order conditions along a two-sided direction.**

If every weight of `x` is strictly positive then the step along `d` may be taken
in both directions, so the linear term of `psi_along_dir` must vanish and its
quadratic term must be nonnegative.

The direction is given by its weight representation, so this one statement
covers a step towards another vertex and a step along a recession direction. -/
theorem foc_soc {p : QuaPiece} {s x d : Plane} {W S : Finset Plane}
    {lam mu gam nu : Plane → ℝ}
    (hWV : W ⊆ p.verts) (hSR : S ⊆ p.rays)
    (h1 : (∑ v ∈ W, lam v) = 1)
    (hx : (∑ v ∈ W, lam v • v) + (∑ r ∈ S, mu r • r) = x)
    (hpos : ∀ v ∈ W, 0 < lam v) (hmupos : ∀ r ∈ S, 0 < mu r)
    (hrep : IsDirRep W S gam nu d)
    (hmax : IsMaxOn (psi p.q s) p.T x) :
    dot s d - dot (p.q.gradAt x) d = 0 ∧ 0 ≤ p.q.alongCurv d := by
  obtain ⟨ε, hε, hstep⟩ := exists_perturb_radius gam nu hpos hmupos
  set A : ℝ := dot s d - dot (p.q.gradAt x) d with hA
  set B : ℝ := p.q.alongCurv d with hB
  have key : ∀ t : ℝ, |t| < ε → t * A - t ^ 2 * B / 2 ≤ 0 := by
    intro t ht
    obtain ⟨hl, hm⟩ := hstep t ht
    have hmem : x + t • d ∈ p.T := mem_T_of_perturb hWV hSR h1 hx hrep hl hm
    have hle : psi p.q s (x + t • d) ≤ psi p.q s x := hmax hmem
    rw [psi_along_dir] at hle
    linarith
  refine ⟨eq_zero_of_forall_small hε key, ?_⟩
  have h2' := key (ε / 2) (by rw [abs_of_pos (by linarith)]; linarith)
  rw [eq_zero_of_forall_small hε key] at h2'
  by_contra hBneg
  push Not at hBneg
  have hpos2 : 0 < (ε / 2) ^ 2 := by positivity
  have : (ε / 2) ^ 2 * B < 0 := mul_neg_of_pos_of_neg hpos2 hBneg
  linarith

/-! ### The core of the selection lemma, and the selection lemma itself

The case analysis S3-S9 for a **single piece**, and the multi-piece assembly on
top of it. Both are proved; nothing in this development is `sorry`ed. -/

private theorem branch_aux (p : QuaPiece) (s : Plane) :
    ∀ n : ℕ, ∀ (W S : Finset Plane) (lam mu : Plane → ℝ) (x : Plane),
      W.card + S.card = n → W ⊆ p.verts → S ⊆ p.rays → W.Nonempty →
      (∀ v ∈ W, 0 ≤ lam v) → (∑ v ∈ W, lam v) = 1 → (∀ r ∈ S, 0 ≤ mu r) →
      (∑ v ∈ W, lam v • v) + (∑ r ∈ S, mu r • r) = x →
      IsMaxOn (psi p.q s) p.T x →
      ∃ b ∈ p.branches, b.eval s = psi p.q s x := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro W S lam mu x hcard hWV hSR hWne h0 h1 hm0 hx hmax
    -- **A vertex weight vanishes**: drop that vertex from the support.
    by_cases hz : ∃ v₀ ∈ W, lam v₀ = 0
    · obtain ⟨v₀, hv₀, hzero⟩ := hz
      obtain ⟨hne', h0', h1', hx'⟩ := erase_vert_data hv₀ h0 h1 hzero hx
      refine ih ((W.erase v₀).card + S.card) ?_ (W.erase v₀) S lam mu x rfl
        (fun z hz' => hWV (Finset.mem_of_mem_erase hz')) hSR hne' h0' h1' hm0 hx' hmax
      rw [← hcard]
      exact Nat.add_lt_add_right (Finset.card_erase_lt_of_mem hv₀) _
    -- **A ray weight vanishes**: drop that ray.
    by_cases hz2 : ∃ r₀ ∈ S, mu r₀ = 0
    · obtain ⟨r₀, hr₀, hzero⟩ := hz2
      refine ih (W.card + (S.erase r₀).card) ?_ W (S.erase r₀) lam mu x rfl hWV
        (fun z hz' => hSR (Finset.mem_of_mem_erase hz')) hWne h0 h1
        (fun z hz' => hm0 z (Finset.mem_of_mem_erase hz')) ?_ hmax
      · rw [← hcard]
        exact Nat.add_lt_add_left (Finset.card_erase_lt_of_mem hr₀) _
      · rw [Finset.sum_erase _ (by rw [hzero, zero_smul])]; exact hx
    -- **A zero recession direction**: it contributes nothing, so drop it.
    by_cases hz3 : (0 : Plane) ∈ S
    · refine ih (W.card + (S.erase 0).card) ?_ W (S.erase 0) lam mu x rfl hWV
        (fun z hz' => hSR (Finset.mem_of_mem_erase hz')) hWne h0 h1
        (fun z hz' => hm0 z (Finset.mem_of_mem_erase hz')) ?_ hmax
      · rw [← hcard]
        exact Nat.add_lt_add_left (Finset.card_erase_lt_of_mem hz3) _
      · rw [Finset.sum_erase _ (by simp)]; exact hx
    -- **Every weight is positive and every ray nonzero.**
    push Not at hz hz2
    have hpos : ∀ v ∈ W, 0 < lam v := fun v hv => lt_of_le_of_ne (h0 v hv) (Ne.symm (hz v hv))
    have hmupos : ∀ r ∈ S, 0 < mu r := fun r hr => lt_of_le_of_ne (hm0 r hr) (Ne.symm (hz2 r hr))
    -- Either the support offers a direction, or it is a single vertex.
    have hcase : (∃ (v₀ w₀ : Plane) (gam₁ nu₁ : Plane → ℝ),
          v₀ ∈ W ∧ (v₀, w₀) ∈ p.edgeCandidates ∧ w₀ - v₀ ≠ 0 ∧
            IsDirRep W S gam₁ nu₁ (w₀ - v₀))
        ∨ (∃ v, W = {v} ∧ S = ∅) := by
      rcases Finset.eq_empty_or_nonempty S with hSe | hSne
      · by_cases h2pt : ∃ v ∈ W, ∃ w ∈ W, v ≠ w
        · obtain ⟨v, hv, w, hw, hvw⟩ := h2pt
          exact Or.inl ⟨v, w, _, _, hv, QuaPiece.mem_edgeCandidates_vert (hWV hv) (hWV hw),
            sub_ne_zero.2 (Ne.symm hvw), isDirRep_vert_sub hv hw⟩
        · push Not at h2pt
          obtain ⟨v, hv⟩ := hWne
          exact Or.inr ⟨v, Finset.eq_singleton_iff_unique_mem.2
            ⟨hv, fun u hu => h2pt u hu v hv⟩, hSe⟩
      · obtain ⟨r₀, hr₀⟩ := hSne
        obtain ⟨v, hv⟩ := hWne
        have hr0ne : r₀ ≠ 0 := fun h => hz3 (h ▸ hr₀)
        refine Or.inl ⟨v, v + r₀, fun _ => 0, fun z => if z = r₀ then (1 : ℝ) else 0, hv,
          QuaPiece.mem_edgeCandidates_ray (hWV hv) (hSR hr₀), ?_, ?_⟩
        · rw [add_sub_cancel_left]; exact hr0ne
        · rw [add_sub_cancel_left]; exact isDirRep_ray hr₀
    rcases hcase with ⟨v₀, w₀, gam₁, nu₁, hv₀, hcand, hd₁ne, hrep₁⟩ | ⟨v, hWeq, hSeq⟩
    · -- **There is a direction.** Read off its first- and second-order conditions.
      obtain ⟨hfoc₁, hsoc₁⟩ := foc_soc hWV hSR h1 hx hpos hmupos hrep₁ hmax
      -- Is the whole support parallel to it?
      have hsplit2 : ((∀ u ∈ W, cross (w₀ - v₀) (u - v₀) = 0) ∧
            (∀ r ∈ S, cross (w₀ - v₀) r = 0))
          ∨ ∃ (d₂ : Plane) (gam₂ nu₂ : Plane → ℝ),
              cross (w₀ - v₀) d₂ ≠ 0 ∧ IsDirRep W S gam₂ nu₂ d₂ := by
        by_cases hW2 : ∀ u ∈ W, cross (w₀ - v₀) (u - v₀) = 0
        · by_cases hS2 : ∀ r ∈ S, cross (w₀ - v₀) r = 0
          · exact Or.inl ⟨hW2, hS2⟩
          · push Not at hS2
            obtain ⟨r₁, hr₁, hcr⟩ := hS2
            exact Or.inr ⟨r₁, _, _, hcr, isDirRep_ray hr₁⟩
        · push Not at hW2
          obtain ⟨u, hu, hcu⟩ := hW2
          exact Or.inr ⟨u - v₀, _, _, hcu, isDirRep_vert_sub hv₀ hu⟩
      rcases hsplit2 with ⟨hW2, hS2⟩ | ⟨d₂, gam₂, nu₂, hcross, hrep₂⟩
      · -- **S7: the face lies on a line.** Edge branch, or vertex branch if flat.
        have hxline : cross (w₀ - v₀) (x - v₀) = 0 := cross_eq_zero_of_support h1 hx hW2 hS2
        obtain ⟨t, htx⟩ := exists_smul_of_cross_eq_zero hd₁ne hxline
        have hxv : x = v₀ + t • (w₀ - v₀) := by rw [← htx]; abel
        have hgrad : p.q.gradAt x = p.q.gradAt v₀ + t • p.q.hessApply (w₀ - v₀) := by
          rw [hxv]; exact Quad.gradAt_add_smul _ _ _ _
        have hslope : edgeSlope p.q v₀ w₀ s = t * edgeCurv p.q v₀ w₀ := by
          have hd : dot (p.q.gradAt x) (w₀ - v₀)
              = dot (p.q.gradAt v₀) (w₀ - v₀) + t * p.q.alongCurv (w₀ - v₀) := by
            rw [hgrad, Quad.alongCurv_eq_dot_hessApply]
            simp only [dot, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
              smul_eq_mul]
            ring
          simp only [edgeSlope, edgeCurv]
          rw [hd] at hfoc₁
          linarith
        rcases lt_trichotomy (edgeCurv p.q v₀ w₀) 0 with hlt0 | heq0 | hgt0
        · exact absurd hsoc₁ (not_le.2 hlt0)
        · refine ⟨vertexBranch p.q v₀, QuaPiece.vertexBranch_mem_branches (hWV hv₀), ?_⟩
          rw [vertexBranch_eval, hxv, psi_along_line, hslope, heq0]
          ring
        · have hcoef : edgeSlope p.q v₀ w₀ s / edgeCurv p.q v₀ w₀ = t := by
            rw [hslope]; field_simp
          have hep : edgePoint p.q v₀ w₀ s = x := by rw [edgePoint, hcoef, hxv]
          refine ⟨edgeBranch p.q v₀ w₀,
            QuaPiece.edgeBranch_mem_branches_of_cand hcand hgt0, ?_⟩
          rw [edgeBranch_eval p.q v₀ w₀ s (ne_of_gt hgt0), hep]
      · -- **The face is two-dimensional**: the conditions pin `∇q(x) = s`.
        obtain ⟨hfoc₂, -⟩ := foc_soc hWV hSR h1 hx hpos hmupos hrep₂ hmax
        have hr0 : s - p.q.gradAt x = 0 :=
          eq_zero_of_dot_eq_zero_of_cross_ne_zero (by rw [dot_sub_left]; exact hfoc₁)
            (by rw [dot_sub_left]; exact hfoc₂) hcross
        have hstat : p.q.gradAt x = s := (sub_eq_zero.1 hr0).symm
        by_cases hdet : p.q.hessDet = 0
        · -- **S8: singular Hessian.** Slide along the kernel onto a proper face.
          obtain ⟨d, hd0, hdk⟩ := exists_kernel_of_hessDet_eq_zero hdet
          obtain ⟨b₁, b₂, hb⟩ := exists_combo_of_cross_ne_zero hcross d
          have hrepd : IsDirRep W S (fun z => b₁ * gam₁ z + b₂ * gam₂ z)
              (fun z => b₁ * nu₁ z + b₂ * nu₂ z) d := by
            have hcomb := (hrep₁.smul b₁).add (hrep₂.smul b₂)
            rwa [← hb] at hcomb
          have hne : (∃ v ∈ W, b₁ * gam₁ v + b₂ * gam₂ v ≠ 0)
              ∨ (∃ r ∈ S, b₁ * nu₁ r + b₂ * nu₂ r ≠ 0) := by
            by_contra hc
            push Not at hc
            have hzW : ∀ v ∈ W, (b₁ * gam₁ v + b₂ * gam₂ v) • v = (0 : Plane) := by
              intro v hv
              rw [hc.1 v hv, zero_smul]
            have hzS : ∀ r ∈ S, (b₁ * nu₁ r + b₂ * nu₂ r) • r = (0 : Plane) := by
              intro r hr
              rw [hc.2 r hr, zero_smul]
            refine hd0 ?_
            rw [← hrepd.2, Finset.sum_eq_zero hzW, Finset.sum_eq_zero hzS, add_zero]
          obtain ⟨t₀, hl, hm, hzero⟩ :=
            exists_descent_gen (lam := lam) (mu := mu)
              (gam := fun z => b₁ * gam₁ z + b₂ * gam₂ z)
              (nu := fun z => b₁ * nu₁ z + b₂ * nu₂ z) hpos hmupos hne
          have hymem : x + t₀ • d ∈ p.T := mem_T_of_perturb hWV hSR h1 hx hrepd hl hm
          have hpsiy : psi p.q s (x + t₀ • d) = psi p.q s x :=
            psi_const_along_kernel hstat hdk t₀
          have hmaxy : IsMaxOn (psi p.q s) p.T (x + t₀ • d) := by
            intro z hz'
            show psi p.q s z ≤ psi p.q s (x + t₀ • d)
            rw [hpsiy]
            exact hmax hz'
          have h1' : (∑ v ∈ W, (lam v + t₀ * (b₁ * gam₁ v + b₂ * gam₂ v))) = 1 :=
            perturb_sum_eq_one h1 hrepd t₀
          have hx' : (∑ v ∈ W, (lam v + t₀ * (b₁ * gam₁ v + b₂ * gam₂ v)) • v)
              + (∑ r ∈ S, (mu r + t₀ * (b₁ * nu₁ r + b₂ * nu₂ r)) • r) = x + t₀ • d :=
            perturb_combo hx hrepd t₀
          rcases hzero with ⟨v₁, hv₁, hz₁⟩ | ⟨r₁, hr₁, hz₁⟩
          · obtain ⟨hne', h0'', h1'', hx''⟩ :=
              erase_vert_data (lam := fun z => lam z + t₀ * (b₁ * gam₁ z + b₂ * gam₂ z))
                (mu := fun z => mu z + t₀ * (b₁ * nu₁ z + b₂ * nu₂ z))
                hv₁ hl h1' hz₁ hx'
            have hlt : (W.erase v₁).card + S.card < n := by
              rw [← hcard]
              exact Nat.add_lt_add_right (Finset.card_erase_lt_of_mem hv₁) _
            obtain ⟨b, hbmem, hbe⟩ :=
              ih ((W.erase v₁).card + S.card) hlt
                (W.erase v₁) S _ _ (x + t₀ • d) rfl
                (fun z hz' => hWV (Finset.mem_of_mem_erase hz')) hSR hne' h0'' h1'' hm hx''
                hmaxy
            exact ⟨b, hbmem, by rw [hbe, hpsiy]⟩
          · have hlt : W.card + (S.erase r₁).card < n := by
              rw [← hcard]
              exact Nat.add_lt_add_left (Finset.card_erase_lt_of_mem hr₁) _
            obtain ⟨b, hbmem, hbe⟩ :=
              ih (W.card + (S.erase r₁).card) hlt
                W (S.erase r₁) _ _ (x + t₀ • d) rfl hWV
                (fun z hz' => hSR (Finset.mem_of_mem_erase hz')) hWne hl h1'
                (fun z hz' => hm z (Finset.mem_of_mem_erase hz'))
                (by rw [Finset.sum_erase _ (by rw [hz₁, zero_smul])]; exact hx') hmaxy
            exact ⟨b, hbmem, by rw [hbe, hpsiy]⟩
        · -- **nonsingular Hessian**: the interior branch, by uniqueness of the
          -- stationary point
          refine ⟨interiorBranch p.q, QuaPiece.interiorBranch_mem_branches hdet, ?_⟩
          rw [interiorBranch_eval p.q s hdet, ← interiorPoint_unique hdet hstat]
    · -- **S6: a single vertex, no rays.**
      have h1' : lam v = 1 := by rw [hWeq, Finset.sum_singleton] at h1; exact h1
      have hxv : x = v := by
        rw [hWeq, hSeq, Finset.sum_singleton, Finset.sum_empty, add_zero, h1', one_smul] at hx
        exact hx.symm
      have hvV : v ∈ p.verts := hWV (by rw [hWeq]; exact Finset.mem_singleton_self v)
      exact ⟨vertexBranch p.q v, QuaPiece.vertexBranch_mem_branches hvV,
        by rw [vertexBranch_eval, hxv]⟩

/-- **The core of the selection lemma (S3-S9), for one piece.**

If `x` maximises `ψ` over the piece, then one of the piece's candidate branches
takes the value `ψ(x)` at `s`.

The proof is a strong induction on the size of the support of `x` — the number
of vertices carrying a positive weight, plus the number of recession directions
carrying a positive weight. A vanishing weight, or a zero recession direction,
shrinks the support. Otherwise `x` may be stepped in both directions along every
difference of two support vertices and along every support ray, so the first-
and second-order conditions hold there (`foc_soc`). If all of those directions
are parallel, `x` lies on a line through a vertex and the value is an
`edgeBranch` — of a vertex pair or of a vertex-and-ray pair, indifferently — or a
`vertexBranch` when the curvature vanishes. If they are not, the conditions force
`∇q(x) = s`, giving `interiorBranch` when the Hessian is nonsingular and, when it
is singular, a slide along `ker H` onto a proper face at no cost to the value
(`psi_const_along_kernel`, `exists_descent_gen`). A single vertex with no rays is
the base case. -/
theorem exists_branch_eq_max {p : QuaPiece} (s : Plane) {x : Plane} (hx : x ∈ p.T)
    (hmax : IsMaxOn (psi p.q s) p.T x) :
    ∃ b ∈ p.branches, b.eval s = psi p.q s x := by
  obtain ⟨c, hc, k, ⟨mu, hmu0, hmusum⟩, hck⟩ := hx
  obtain ⟨lam, h0, h1, h2⟩ := exists_weights hc
  exact branch_aux p s (p.verts.card + p.rays.card) p.verts p.rays lam mu x rfl
    (fun _ hz => hz) (fun _ hz => hz) p.verts_nonempty h0 h1 hmu0
    (by rw [h2, hmusum]; exact hck) hmax

/-- **Selection.** Some candidate quadratic attains the conjugate, wherever every
piece attains its own supremum.

For bounded pieces that is automatic (`attained_of_bounded`); with recession
directions it is the content of Phase 7. -/
theorem selection {f : QuaPol} {s : Plane} (hA : f.Attained s) :
    ∃ q ∈ cand f, ((q.eval s : ℝ) : EReal) = f.conj s := by
  obtain ⟨p, hp, x, hxT, hxmax, hval⟩ := exists_maximiser hA
  obtain ⟨b, hb, hbe⟩ := exists_branch_eq_max s hxT hxmax
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

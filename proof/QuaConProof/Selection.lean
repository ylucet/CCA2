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

/-- **The first- and second-order conditions at an interior weight.**

If all weights are positive and `v ≠ w` are vertices, then stepping along
`w - v` is possible in both directions, so the linear term of `psi_along_dir`
must vanish and its quadratic term must be nonnegative. -/
theorem foc_soc {p : QuaPiece} {s x : Plane} {W : Finset Plane} {lam : Plane → ℝ}
    (hWV : W ⊆ p.verts)
    (h0 : ∀ u ∈ W, 0 ≤ lam u) (h1 : (∑ u ∈ W, lam u) = 1) (h2 : (∑ u ∈ W, lam u • u) = x)
    (hpos : ∀ u ∈ W, 0 < lam u)
    (hmax : IsMaxOn (psi p.q s) p.T x)
    {v w : Plane} (hv : v ∈ W) (hw : w ∈ W) (hvw : v ≠ w) :
    dot s (w - v) - dot (p.q.gradAt x) (w - v) = 0 ∧ 0 ≤ p.q.alongCurv (w - v) := by
  set A : ℝ := dot s (w - v) - dot (p.q.gradAt x) (w - v) with hA
  set B : ℝ := p.q.alongCurv (w - v) with hB
  set ε : ℝ := min (lam v) (lam w) with hε
  have hε0 : 0 < ε := lt_min (hpos v hv) (hpos w hw)
  have hsub : convexHull ℝ (↑W : Set Plane) ⊆ p.T :=
    convexHull_mono (Finset.coe_subset.2 hWV)
  have key : ∀ t : ℝ, |t| < ε → t * A - t ^ 2 * B / 2 ≤ 0 := by
    intro t ht
    have htv : t ≤ lam v := le_trans (le_abs_self t) (le_of_lt (lt_of_lt_of_le ht (min_le_left _ _)))
    have htw : -t ≤ lam w := le_trans (neg_le_abs t) (le_of_lt (lt_of_lt_of_le ht (min_le_right _ _)))
    have hmem : x + t • (w - v) ∈ p.T :=
      hsub (mem_convexHull_perturb h0 h1 h2 hv hw hvw htv htw)
    have hle : psi p.q s (x + t • (w - v)) ≤ psi p.q s x := hmax hmem
    rw [psi_along_dir] at hle
    linarith
  refine ⟨eq_zero_of_forall_small hε0 key, ?_⟩
  have h2' := key (ε / 2) (by rw [abs_of_pos (by linarith)]; linarith)
  rw [eq_zero_of_forall_small hε0 key] at h2'
  by_contra hBneg
  push Not at hBneg
  have hpos2 : 0 < (ε / 2) ^ 2 := by positivity
  have : (ε / 2) ^ 2 * B < 0 := mul_neg_of_pos_of_neg hpos2 hBneg
  linarith

/-! ### The core of the selection lemma, and the selection lemma itself

The case analysis S3-S9 for a **single piece**, and the multi-piece assembly on
top of it. Both are proved; nothing in this development is `sorry`ed. -/

private theorem branch_aux (p : QuaPiece) (s : Plane) :
    ∀ n : ℕ, ∀ (W : Finset Plane) (x : Plane), W.card = n → W ⊆ p.verts → W.Nonempty →
      x ∈ convexHull ℝ (↑W : Set Plane) → IsMaxOn (psi p.q s) p.T x →
      ∃ b ∈ p.branches, b.eval s = psi p.q s x := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro W x hcard hWV hWne hxW hmax
    obtain ⟨lam, h0, h1, h2⟩ := exists_weights hxW
    by_cases hzero : ∃ v ∈ W, lam v = 0
    · -- **A weight vanishes**: the point already lies in the hull of a smaller set.
      obtain ⟨v₀, hv₀W, hv₀0⟩ := hzero
      have hxe : x ∈ convexHull ℝ (↑(W.erase v₀) : Set Plane) :=
        mem_convexHull_erase h0 h1 h2 hv₀0
      have hlt : (W.erase v₀).card < n := by
        rw [← hcard]; exact Finset.card_erase_lt_of_mem hv₀W
      have hnee : (W.erase v₀).Nonempty := by
        rcases Finset.eq_empty_or_nonempty (W.erase v₀) with he | he
        · rw [he] at hxe; simp at hxe
        · exact he
      exact ih _ hlt (W.erase v₀) x rfl
        (fun u hu => hWV (Finset.mem_of_mem_erase hu)) hnee hxe hmax
    · -- **Every weight is positive**: `x` is in the relative interior of the face.
      push Not at hzero
      have hpos : ∀ u ∈ W, 0 < lam u := fun u hu =>
        lt_of_le_of_ne (h0 u hu) (Ne.symm (hzero u hu))
      by_cases h2pt : ∃ v ∈ W, ∃ w ∈ W, v ≠ w
      · obtain ⟨v, hv, w, hw, hvw⟩ := h2pt
        have hDne : (w - v) ≠ 0 := sub_ne_zero.2 (Ne.symm hvw)
        obtain ⟨hfoc, hsoc⟩ := foc_soc hWV h0 h1 h2 hpos hmax hv hw hvw
        by_cases hcol : ∀ u ∈ W, cross (w - v) (u - v) = 0
        · -- **S7: the face lies on a line.** Edge branch, or vertex branch if flat.
          have hxline : cross (w - v) (x - v) = 0 := cross_eq_zero_of_forall h1 h2 hcol
          obtain ⟨t, htx⟩ := exists_smul_of_cross_eq_zero hDne hxline
          have hxv : x = v + t • (w - v) := by rw [← htx]; abel
          have hgrad : p.q.gradAt x = p.q.gradAt v + t • p.q.hessApply (w - v) := by
            rw [hxv]; exact Quad.gradAt_add_smul _ _ _ _
          have hslope : edgeSlope p.q v w s = t * edgeCurv p.q v w := by
            have hd : dot (p.q.gradAt x) (w - v)
                = dot (p.q.gradAt v) (w - v) + t * p.q.alongCurv (w - v) := by
              rw [hgrad, Quad.alongCurv_eq_dot_hessApply]
              simp only [dot, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
                smul_eq_mul]
              ring
            simp only [edgeSlope, edgeCurv]
            rw [hd] at hfoc
            linarith
          rcases lt_trichotomy (edgeCurv p.q v w) 0 with hlt0 | heq0 | hgt0
          · exact absurd hsoc (not_le.2 hlt0)
          · -- flat along the line: the value is already the vertex value at `v`
            refine ⟨vertexBranch p.q v, QuaPiece.vertexBranch_mem_branches (hWV hv), ?_⟩
            rw [vertexBranch_eval, hxv, psi_along_line, hslope, heq0]
            ring
          · -- a genuine edge branch, and the pair is kept by `edgePairs`
            have hcoef : edgeSlope p.q v w s / edgeCurv p.q v w = t := by
              rw [hslope]; field_simp
            have hep : edgePoint p.q v w s = x := by rw [edgePoint, hcoef, hxv]
            refine ⟨edgeBranch p.q v w,
              QuaPiece.edgeBranch_mem_branches (hWV hv) (hWV hw) hgt0, ?_⟩
            rw [edgeBranch_eval p.q v w s (ne_of_gt hgt0), hep]
        · -- **The face is two-dimensional**: the conditions pin `∇q(x) = s`.
          push Not at hcol
          obtain ⟨u, huW, hu0⟩ := hcol
          have huv : v ≠ u := by
            intro h
            rw [← h] at hu0
            exact hu0 (by simp [cross])
          obtain ⟨hfoc2, -⟩ := foc_soc hWV h0 h1 h2 hpos hmax hv huW huv
          have hr : s - p.q.gradAt x = 0 := by
            refine eq_zero_of_dot_eq_zero_of_cross_ne_zero ?_ ?_ hu0
            · rw [dot_sub_left]; exact hfoc
            · rw [dot_sub_left]; exact hfoc2
          have hstat : p.q.gradAt x = s := (sub_eq_zero.1 hr).symm
          by_cases hdet : p.q.hessDet = 0
          · -- **S8: singular Hessian.** Slide along the kernel onto a proper face.
            obtain ⟨d, hd0, hdk⟩ := exists_kernel_of_hessDet_eq_zero hdet
            have hwu : w ≠ u := by
              intro h
              rw [h] at hu0
              exact hu0 (cross_self _)
            obtain ⟨b₁, b₂, hb⟩ := exists_combo_of_cross_ne_zero hu0 d
            classical
            set gam : Plane → ℝ :=
              fun z => if z = v then -(b₁ + b₂) else if z = w then b₁ else if z = u then b₂ else 0
              with hgamdef
            have hgv : gam v = -(b₁ + b₂) := by simp [hgamdef]
            have hgw : gam w = b₁ := by simp [hgamdef, Ne.symm hvw]
            have hgu : gam u = b₂ := by simp [hgamdef, Ne.symm huv, Ne.symm hwu]
            have hgo : ∀ z, z ≠ v → z ≠ w → z ≠ u → gam z = 0 := by
              intro z hzv hzw hzu; simp [hgamdef, hzv, hzw, hzu]
            have hsub3 : ({v, w, u} : Finset Plane) ⊆ W := by
              intro z hz
              simp only [Finset.mem_insert, Finset.mem_singleton] at hz
              rcases hz with rfl | rfl | rfl <;> assumption
            have hsum3 : ∀ (M : Type) (_ : AddCommMonoid M) (g : Plane → M),
                (∀ z, z ≠ v → z ≠ w → z ≠ u → g z = 0) →
                (∑ z ∈ W, g z) = g v + g w + g u := by
              intro M _ g hg
              rw [← Finset.sum_subset hsub3 (fun z _ hz => by
                simp only [Finset.mem_insert, Finset.mem_singleton, not_or] at hz
                exact hg z hz.1 hz.2.1 hz.2.2)]
              rw [Finset.sum_insert (by simp [hvw, huv]), Finset.sum_insert (by simp [hwu]),
                Finset.sum_singleton, add_assoc]
            have hgsum : (∑ z ∈ W, gam z) = 0 := by
              rw [hsum3 ℝ inferInstance gam hgo, hgv, hgw, hgu]; ring
            have hgd : (∑ z ∈ W, gam z • z) = d := by
              rw [hsum3 Plane inferInstance (fun z => gam z • z)
                (fun z a b c => by rw [hgo z a b c, zero_smul]), hgv, hgw, hgu, hb]
              module
            have hb12 : gam w ≠ 0 ∨ gam u ≠ 0 := by
              rw [hgw, hgu]
              by_contra hc
              push Not at hc
              rw [hc.1, hc.2] at hb
              simp only [zero_smul, add_zero] at hb
              exact hd0 hb
            have hgne : ∃ z ∈ W, gam z ≠ 0 := by
              rcases hb12 with h | h
              · exact ⟨w, hw, h⟩
              · exact ⟨u, huW, h⟩
            obtain ⟨t₀, v₀, hv₀W, hymem⟩ := exists_descent hpos h1 h2 hgsum hgd hgne
            have hpsiy : psi p.q s (x + t₀ • d) = psi p.q s x :=
              psi_const_along_kernel hstat hdk t₀
            have hmaxy : IsMaxOn (psi p.q s) p.T (x + t₀ • d) := by
              intro z hz
              show psi p.q s z ≤ psi p.q s (x + t₀ • d)
              rw [hpsiy]
              exact hmax hz
            have hlt : (W.erase v₀).card < n := by
              rw [← hcard]; exact Finset.card_erase_lt_of_mem hv₀W
            have hnee : (W.erase v₀).Nonempty := by
              by_cases hvv : v₀ = v
              · exact ⟨w, Finset.mem_erase.2 ⟨by rw [hvv]; exact Ne.symm hvw, hw⟩⟩
              · exact ⟨v, Finset.mem_erase.2 ⟨fun h => hvv h.symm, hv⟩⟩
            obtain ⟨b, hbmem, hbe⟩ := ih _ hlt (W.erase v₀) (x + t₀ • d) rfl
              (fun z hz => hWV (Finset.mem_of_mem_erase hz)) hnee hymem hmaxy
            exact ⟨b, hbmem, by rw [hbe, hpsiy]⟩
          · -- **nonsingular Hessian**: the interior branch, by uniqueness of the
            -- stationary point
            refine ⟨interiorBranch p.q, QuaPiece.interiorBranch_mem_branches hdet, ?_⟩
            rw [interiorBranch_eval p.q s hdet, ← interiorPoint_unique hdet hstat]
      · -- **S6: a single vertex.**
        push Not at h2pt
        obtain ⟨v, hv⟩ := hWne
        have hWeq : W = {v} :=
          Finset.eq_singleton_iff_unique_mem.2 ⟨hv, fun u hu => h2pt u hu v hv⟩
        subst hWeq
        rw [Finset.sum_singleton] at h1 h2
        have hxv : x = v := by rw [← h2, h1, one_smul]
        refine ⟨vertexBranch p.q v, QuaPiece.vertexBranch_mem_branches (hWV hv), ?_⟩
        rw [vertexBranch_eval, hxv]

/-- **The core of the selection lemma (S3-S9), for one piece.**

If `x` maximises `ψ` over the piece, then one of the piece's candidate branches
takes the value `ψ(x)` at `s`.

The proof is a strong induction on the number of vertices of the face carrying
`x`. If a barycentric weight vanishes, the face shrinks. Otherwise `x` is in the
relative interior, so one may step both ways along any `w - v` and read off the
first- and second-order conditions (`foc_soc`). If the face lies on a line the
value is an `edgeBranch`, or a `vertexBranch` when the curvature vanishes; if it
does not, the conditions force `∇q(x) = s`, giving `interiorBranch` when the
Hessian is nonsingular and, when it is singular, a slide along `ker H` onto a
proper face at no cost to the value (`psi_const_along_kernel`, `exists_descent`).
A single vertex is the base case. -/
theorem exists_branch_eq_max (p : QuaPiece) (s : Plane) {x : Plane} (hx : x ∈ p.T)
    (hmax : IsMaxOn (psi p.q s) p.T x) :
    ∃ b ∈ p.branches, b.eval s = psi p.q s x :=
  branch_aux p s p.verts.card p.verts x rfl (fun _ hz => hz) p.verts_nonempty hx hmax

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

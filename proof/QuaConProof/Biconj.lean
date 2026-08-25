/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.Convexity

/-!
# The biconjugate

`TODO.md` item A4, and the base of Track C. `f**` is the conjugate of `f*`,

    f**(x) = sup_s ⟨s,x⟩ - f*(s),

taken over all of the plane and valued in `EReal`, exactly as `conj` is.

## Fenchel–Young, stated so that it is usable

The textbook form `⟨s,x⟩ ≤ f(x) + f*(s)` is an `EReal` inequality with an
addition on the right, and both summands can be `⊤`. Everything below is stated
instead through the **affine minorant**

    A_s(x) = ⟨s,x⟩ - f*(s)     (a genuine affine function of `x`, for `s ∈ dom f*`)

which is what every later argument actually uses: `affineMinorant_le_eval` says
`A_s ≤ f`, and `affineMinorant_le_biconj` says `A_s ≤ f**`. The real-valued
Fenchel–Young inequality is recorded as `fenchel_young` for the record.

Note `f**` really can be `⊥`: a single unbounded piece carrying a concave
quadratic has `f* ≡ ⊤`, hence `f** ≡ ⊥`. So there is no `biconj_ne_bot` to
match `conj_ne_bot`.

## Main results

* `biconj` — the definition.
* `affineMinorant_le_eval`, `affineMinorant_le_biconj` — the two halves of
  Fenchel–Young in the form Track C consumes.
* `biconj_le_eval` — `f** ≤ f`.
-/

namespace QuaConProof

namespace QuaPol

/-- The **biconjugate**, `f**(x) = sup_s ⟨s,x⟩ - f*(s)`, in `EReal`. -/
noncomputable def biconj (f : QuaPol) (x : Plane) : EReal :=
  ⨆ s : Plane, ((dot s x : ℝ) : EReal) - f.conj s

lemma biconj_def (f : QuaPol) (x : Plane) :
    f.biconj x = ⨆ s : Plane, ((dot s x : ℝ) : EReal) - f.conj s := rfl

/-- The **affine minorant** attached to a point `s` of `dom f*`:
`A_s(x) = ⟨s,x⟩ - f*(s)`. -/
noncomputable def affineMinorant (f : QuaPol) (s x : Plane) : ℝ :=
  dot s x - (f.conj s).toReal

/-- **`A_s ≤ f`.** This is Fenchel–Young: no point of the graph of `f` lies below
the affine function cut out by `s`. -/
theorem affineMinorant_le_eval {f : QuaPol} {s : Plane} (hs : f.conj s ≠ ⊤) (x : Plane) :
    ((f.affineMinorant s x : ℝ) : EReal) ≤ f.eval x := by
  have hC : f.conj s = (((f.conj s).toReal : ℝ) : EReal) :=
    (EReal.coe_toReal hs (conj_ne_bot f s)).symm
  by_cases hx : f.eval x = ⊤
  · rw [hx]
    exact le_top
  · have hE : f.eval x = (((f.eval x).toReal : ℝ) : EReal) :=
      (EReal.coe_toReal hx (eval_ne_bot f x)).symm
    have hterm : ((dot s x : ℝ) : EReal) - f.eval x ≤ f.conj s := term_le_conj f s x
    rw [hE, hC, ← EReal.coe_sub, EReal.coe_le_coe_iff] at hterm
    rw [hE, EReal.coe_le_coe_iff, affineMinorant]
    linarith

/-- **`A_s ≤ f**`.** One term of the supremum defining the biconjugate. -/
theorem affineMinorant_le_biconj {f : QuaPol} {s : Plane} (hs : f.conj s ≠ ⊤) (x : Plane) :
    ((f.affineMinorant s x : ℝ) : EReal) ≤ f.biconj x := by
  have hC : f.conj s = (((f.conj s).toReal : ℝ) : EReal) :=
    (EReal.coe_toReal hs (conj_ne_bot f s)).symm
  refine le_trans (le_of_eq ?_) (le_iSup (fun t : Plane =>
    ((dot t x : ℝ) : EReal) - f.conj t) s)
  rw [hC, ← EReal.coe_sub, affineMinorant]

/-- **`f** ≤ f`.** Every affine minorant of `f` is a minorant, so their supremum
is one too. -/
theorem biconj_le_eval (f : QuaPol) (x : Plane) : f.biconj x ≤ f.eval x := by
  rw [biconj_def]
  refine iSup_le fun s => ?_
  by_cases hs : f.conj s = ⊤
  · rw [hs, EReal.sub_top]
    exact bot_le
  · have hC : f.conj s = (((f.conj s).toReal : ℝ) : EReal) :=
      (EReal.coe_toReal hs (conj_ne_bot f s)).symm
    have := affineMinorant_le_eval hs x
    rwa [affineMinorant, EReal.coe_sub, ← hC] at this

/-- **Fenchel–Young**, in the real-valued form: where both `f` and `f*` are
finite, `⟨s,x⟩ ≤ f(x) + f*(s)`. -/
theorem fenchel_young {f : QuaPol} {s x : Plane} (hs : f.conj s ≠ ⊤) (hx : f.eval x ≠ ⊤) :
    dot s x ≤ (f.eval x).toReal + (f.conj s).toReal := by
  have h := affineMinorant_le_eval hs x
  have hE : f.eval x = (((f.eval x).toReal : ℝ) : EReal) :=
    (EReal.coe_toReal hx (eval_ne_bot f x)).symm
  rw [hE, EReal.coe_le_coe_iff, affineMinorant] at h
  linarith

end QuaPol

end QuaConProof

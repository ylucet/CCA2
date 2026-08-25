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
* `convex_epigraph_biconj`, `convex_dom_biconj`, `lowerSemicontinuous_biconj` —
  `f**` is convex and lsc, by the `supAffine` lemmas of `Convexity.lean`.
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

/-! ### `f**` is convex and lower semicontinuous

`TODO.md` item C1. Both statements come straight from `Convexity.lean`: `f**` is
`supAffine f.conj`, and `conj_ne_bot` supplies the one hypothesis those lemmas
need. Nothing is reproved. -/

/-- `f**` is a supremum of affine functions, by definition. -/
lemma biconj_eq_supAffine (f : QuaPol) : f.biconj = supAffine f.conj := rfl

/-- **The epigraph inequality for `f**`.** -/
theorem biconj_le_of_combo (f : QuaPol) {x₁ x₂ : Plane} {a b t₁ t₂ : ℝ}
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : a + b = 1)
    (h₁ : f.biconj x₁ ≤ (t₁ : EReal)) (h₂ : f.biconj x₂ ≤ (t₂ : EReal)) :
    f.biconj (a • x₁ + b • x₂) ≤ ((a * t₁ + b * t₂ : ℝ) : EReal) :=
  supAffine_le_of_combo (conj_ne_bot f) ha hb hab h₁ h₂

/-- **The epigraph of `f**` is convex.** -/
theorem convex_epigraph_biconj (f : QuaPol) :
    Convex ℝ {p : Plane × ℝ | f.biconj p.1 ≤ ((p.2 : ℝ) : EReal)} :=
  convex_epigraph_supAffine (conj_ne_bot f)

/-- **The domain of `f**` is convex.** -/
theorem convex_dom_biconj (f : QuaPol) : Convex ℝ {x : Plane | f.biconj x ≠ ⊤} :=
  convex_dom_supAffine (conj_ne_bot f)

/-- **`f**` is lower semicontinuous.** -/
theorem lowerSemicontinuous_biconj (f : QuaPol) : LowerSemicontinuous f.biconj :=
  lowerSemicontinuous_supAffine (conj_ne_bot f)

/-- **Fenchel–Young**, in the real-valued form: where both `f` and `f*` are
finite, `⟨s,x⟩ ≤ f(x) + f*(s)`. -/
theorem fenchel_young {f : QuaPol} {s x : Plane} (hs : f.conj s ≠ ⊤) (hx : f.eval x ≠ ⊤) :
    dot s x ≤ (f.eval x).toReal + (f.conj s).toReal := by
  have h := affineMinorant_le_eval hs x
  have hE : f.eval x = (((f.eval x).toReal : ℝ) : EReal) :=
    (EReal.coe_toReal hx (eval_ne_bot f x)).symm
  rw [hE, EReal.coe_le_coe_iff, affineMinorant] at h
  linarith

/-! ### The key lemma: `f**` is affine on the hull of the maximisers

`TODO.md` item C2, and `../CONJ_FIELD_PROOF.md` §5, where it is written for a
triangle. Everything in Theorem 4's table except the 2-cell rows is a corollary
of this one statement.

Fix `s` with `f*(s)` finite and let `maxSet f s` be the points of the pieces at
which the supremum defining `f*(s)` is attained. Then `f**` **equals** the affine
function `A_s` on the whole convex hull of that set.

The proof is short and worth reading. `A_s ≤ f` everywhere by Fenchel–Young, so
`A_s ≤ f**`. At a maximiser `x` in piece `p`,

    A_s(x) = q_p(x) ≥ f(x) ≥ f**(x) ≥ A_s(x)

so all four agree there. `f**` is convex and `A_s` is affine, so `f** ≤ A_s`
propagates from the maximisers to their hull; with the reverse inequality, they
are equal. -/

/-- The points at which the supremum defining `f*(s)` is attained. Empty when
`f*(s) = ⊤`, since `psi` is real-valued. -/
def maxSet (f : QuaPol) (s : Plane) : Set Plane :=
  {x | ∃ p ∈ f.pieces, x ∈ p.T ∧ ((psi p.q s x : ℝ) : EReal) = f.conj s}

/-- Where `f*` is finite the maximiser set is nonempty — that is Frank–Wolfe. -/
theorem maxSet_nonempty {f : QuaPol} {s : Plane} (hs : f.conj s ≠ ⊤) :
    (maxSet f s).Nonempty := by
  obtain ⟨p, hp, x, hxT, _, hval⟩ := exists_maximiser (attained_of_conj_ne_top hs)
  exact ⟨x, p, hp, hxT, hval⟩

/-- The affine minorant is affine: it respects convex combinations. -/
lemma affineMinorant_combo (f : QuaPol) (s : Plane) {a b : ℝ} (hab : a + b = 1)
    (x₁ x₂ : Plane) :
    f.affineMinorant s (a • x₁ + b • x₂)
      = a * f.affineMinorant s x₁ + b * f.affineMinorant s x₂ := by
  simp only [affineMinorant, dot_combo_right]
  linear_combination ((f.conj s).toReal) * hab

/-- **At a maximiser, `f**` already equals the affine minorant.** -/
theorem biconj_eq_affineMinorant_of_mem_maxSet {f : QuaPol} {s : Plane}
    (hs : f.conj s ≠ ⊤) {x : Plane} (hx : x ∈ maxSet f s) :
    f.biconj x = ((f.affineMinorant s x : ℝ) : EReal) := by
  obtain ⟨p, hp, hxT, hval⟩ := hx
  -- the piece's quadratic at `x` *is* the affine minorant there
  have hq : p.q.eval x = f.affineMinorant s x := by
    have hC : f.conj s = (((f.conj s).toReal : ℝ) : EReal) :=
      (EReal.coe_toReal hs (conj_ne_bot f s)).symm
    rw [hC, EReal.coe_eq_coe_iff] at hval
    simp only [affineMinorant, ← hval, psi]
    ring
  refine le_antisymm ?_ (affineMinorant_le_biconj hs x)
  calc f.biconj x ≤ f.eval x := biconj_le_eval f x
    _ ≤ ((p.q.eval x : ℝ) : EReal) := eval_le_of_mem hp hxT
    _ = ((f.affineMinorant s x : ℝ) : EReal) := by rw [hq]

/-- **The key lemma.** On the convex hull of the maximisers, `f**` is the affine
function `A_s`. -/
theorem biconj_eq_affineMinorant_on_hull {f : QuaPol} {s : Plane} (hs : f.conj s ≠ ⊤)
    {x : Plane} (hx : x ∈ convexHull ℝ (maxSet f s)) :
    f.biconj x = ((f.affineMinorant s x : ℝ) : EReal) := by
  refine le_antisymm ?_ (affineMinorant_le_biconj hs x)
  -- the set where `f**` lies below the affine function is convex and contains
  -- every maximiser, so it contains their hull
  refine convexHull_min (t := {z : Plane | f.biconj z ≤ ((f.affineMinorant s z : ℝ) : EReal)})
    (fun z hz => le_of_eq (biconj_eq_affineMinorant_of_mem_maxSet hs hz)) ?_ hx
  intro x₁ h₁ x₂ h₂ a b ha hb hab
  show f.biconj (a • x₁ + b • x₂) ≤ ((f.affineMinorant s (a • x₁ + b • x₂) : ℝ) : EReal)
  rw [affineMinorant_combo f s hab]
  exact biconj_le_of_combo f ha hb hab h₁ h₂

/-- **and `s` is a subgradient of `f**` there.** That is what "`f**` is affine
with gradient `s` on this cell" means without differentiating. -/
theorem biconj_subgradient {f : QuaPol} {s : Plane} (hs : f.conj s ≠ ⊤)
    {x : Plane} (hx : x ∈ convexHull ℝ (maxSet f s)) (z : Plane) :
    f.biconj x + ((dot s z - dot s x : ℝ) : EReal) ≤ f.biconj z := by
  rw [biconj_eq_affineMinorant_on_hull hs hx, ← EReal.coe_add]
  have hrw : f.affineMinorant s x + (dot s z - dot s x) = f.affineMinorant s z := by
    simp only [affineMinorant]
    ring
  rw [hrw]
  exact affineMinorant_le_biconj hs z

/-! ### Rows 4 and 5 of Theorem 4

`TODO.md` items C3 and C4. Both are corollaries of the key lemma: a 1-cell of
`f*` (two branches active) gives a **ruled** 2-cell of `f**` on which `f**` is
affine along each ruling, and a 0-cell (three active) gives a 2-cell on which
`f**` is **affine outright**. In both cases the value is `⟨s,x⟩ - f*(s)` and the
gradient is `s`.

## Getting from an active candidate to a maximiser

The corners of those cells are the *maximisers*, not the active candidates, and
the two are not the same thing: an edge branch can be active while its stationary
point lies outside its own segment, in which case it is not attained there
(`DECISIONS.md`, 2026-08-21 — this is exactly why `selection` is not
`f* = max over cand`). The three bridge lemmas below say when an active candidate
does deliver a maximiser. For a **vertex** branch there is no side condition, the
vertex always lies in its piece; for the other two the stationary point must. -/

/-- An active vertex branch always contributes its vertex to `maxSet`. -/
lemma mem_maxSet_of_vertexBranch_active {f : QuaPol} {p : QuaPiece} {v s : Plane}
    (hp : p ∈ f.pieces) (hv : v ∈ p.verts)
    (hact : (((vertexBranch p.q v).eval s : ℝ) : EReal) = f.conj s) :
    v ∈ maxSet f s :=
  ⟨p, hp, p.subset_T hv, by rwa [vertexBranch_eval] at hact⟩

/-- An active edge branch contributes its stationary point, when that point lies
in the piece. -/
lemma mem_maxSet_of_edgeBranch_active {f : QuaPol} {p : QuaPiece} {v w s : Plane}
    (hp : p ∈ f.pieces) (hcurv : edgeCurv p.q v w ≠ 0)
    (hin : edgePoint p.q v w s ∈ p.T)
    (hact : (((edgeBranch p.q v w).eval s : ℝ) : EReal) = f.conj s) :
    edgePoint p.q v w s ∈ maxSet f s :=
  ⟨p, hp, hin, by rwa [edgeBranch_eval p.q v w s hcurv] at hact⟩

/-- An active interior branch contributes its stationary point, when that point
lies in the piece. -/
lemma mem_maxSet_of_interiorBranch_active {f : QuaPol} {p : QuaPiece} {s : Plane}
    (hp : p ∈ f.pieces) (hdet : p.q.hessDet ≠ 0)
    (hin : interiorPoint p.q s ∈ p.T)
    (hact : (((interiorBranch p.q).eval s : ℝ) : EReal) = f.conj s) :
    interiorPoint p.q s ∈ maxSet f s :=
  ⟨p, hp, hin, by rwa [interiorBranch_eval p.q s hdet] at hact⟩

/-- **Row 4 of Theorem 4: the ruled cell.** Two maximisers at `s` give a segment
on which `f**` is affine — a ruling — with value `⟨s,x⟩ - f*(s)`. -/
theorem biconj_affine_on_segment {f : QuaPol} {s : Plane} (hs : f.conj s ≠ ⊤)
    {x₁ x₂ : Plane} (h₁ : x₁ ∈ maxSet f s) (h₂ : x₂ ∈ maxSet f s)
    {x : Plane} (hx : x ∈ segment ℝ x₁ x₂) :
    f.biconj x = ((f.affineMinorant s x : ℝ) : EReal) := by
  refine biconj_eq_affineMinorant_on_hull hs ?_
  refine convexHull_mono (s := ({x₁, x₂} : Set Plane)) ?_ ?_
  · rintro z (rfl | rfl)
    · exact h₁
    · exact h₂
  · rwa [convexHull_pair]

/-- **Row 5 of Theorem 4: the affine cell over a vertex of `f*`.** Three
maximisers at `s` give a triangle on which `f**` is affine. -/
theorem biconj_affine_on_triangle {f : QuaPol} {s : Plane} (hs : f.conj s ≠ ⊤)
    {x₁ x₂ x₃ : Plane} (h₁ : x₁ ∈ maxSet f s) (h₂ : x₂ ∈ maxSet f s)
    (h₃ : x₃ ∈ maxSet f s)
    {x : Plane} (hx : x ∈ convexHull ℝ ({x₁, x₂, x₃} : Set Plane)) :
    f.biconj x = ((f.affineMinorant s x : ℝ) : EReal) := by
  refine biconj_eq_affineMinorant_on_hull hs (convexHull_mono ?_ hx)
  rintro z (rfl | rfl | rfl)
  · exact h₁
  · exact h₂
  · exact h₃

end QuaPol

end QuaConProof

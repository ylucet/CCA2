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

/-! ### Rows 1 to 3 of Theorem 4, and what they can and cannot say

`TODO.md` item C5, which flagged a risk: *"exactly one active candidate implies a
unique maximiser" is not obvious and may be false.* **It is false**, and
`Sanity.oneActive_manyMaximisers` below is the witness — one piece, one
candidate, and a whole line of maximisers. So rows 1 to 3 cannot be indexed by
the activity pattern; they have to be indexed by the branch that is attained, and
by the rank of *its* quadratic part.

What is provable, and is proved here:

| row | branch attained | the maximiser |
|---|---|---|
| 1 | interior, `hessDet ≠ 0` | `interiorPoint q s`, an **affine** function of `s` |
| 2 | edge along `v..w` | lies on the **fixed line** `v + ℝ(w-v)`, for every `s` |
| 3 | vertex `v` | the **constant** `v` |

Row 1's uniqueness holds within a positive definite piece, where `ψ` is strictly
concave (`maxOn_unique_of_posDef`). Across pieces it can fail for the obvious
reason — two pieces can both attain — and that is not a defect, it is the ruled
cell of row 4.

Rows 1 to 3 as statements about *cells* of `f*` need those cells to be
two-dimensional, which is the face-to-face regularity this development
deliberately does not claim (`DECISIONS.md`, 2026-08-21). That part is
scope-reduced, not proved. -/

/-- The set on which `f**` and the affine minorant agree — the subdifferential of
`f*` at `s`, written without derivatives. -/
def subgradSet (f : QuaPol) (s : Plane) : Set Plane :=
  {x | f.biconj x = ((f.affineMinorant s x : ℝ) : EReal)}

/-- **The subdifferential correspondence.** The whole hull of the maximisers lies
in `∂f*(s)`. This is C2 read as a statement about sets, and it is what indexes
the cells of `f**` by the cells of `f*`. -/
theorem convexHull_maxSet_subset_subgradSet {f : QuaPol} {s : Plane} (hs : f.conj s ≠ ⊤) :
    convexHull ℝ (maxSet f s) ⊆ subgradSet f s :=
  fun _ hx => biconj_eq_affineMinorant_on_hull hs hx

/-! #### Row 1: a positive definite piece has a unique maximiser -/

/-- A positive definite quadratic has **strictly** positive curvature in every
nonzero direction. -/
theorem alongCurv_pos {q : Quad} (ha : 0 < q.a) (hD : 0 < q.hessDet) {d : Plane}
    (hd : d ≠ 0) : 0 < q.alongCurv d := by
  have hD' : 0 < 4 * q.a * q.c - q.b ^ 2 := hD
  have key : 2 * q.a * q.alongCurv d
      = (2 * q.a * d.1 + q.b * d.2) ^ 2 + (4 * q.a * q.c - q.b ^ 2) * d.2 ^ 2 := by
    simp only [Quad.alongCurv]
    ring
  have hpos : 0 < (2 * q.a * d.1 + q.b * d.2) ^ 2 + (4 * q.a * q.c - q.b ^ 2) * d.2 ^ 2 := by
    by_cases h2 : d.2 = 0
    · have h1 : d.1 ≠ 0 := fun h1 => hd (Prod.ext_iff.2 ⟨h1, h2⟩)
      rw [h2]
      have hsq : 0 < (2 * q.a * d.1 + q.b * 0) ^ 2 := by
        have : 2 * q.a * d.1 + q.b * 0 ≠ 0 := by
          simp only [mul_zero, add_zero]
          positivity
        exact lt_of_le_of_ne (sq_nonneg _) (Ne.symm (pow_ne_zero 2 this))
      nlinarith
    · have hd2 : 0 < d.2 ^ 2 := by positivity
      nlinarith [sq_nonneg (2 * q.a * d.1 + q.b * d.2)]
  by_contra hcon
  push Not at hcon
  nlinarith

/-- **At most one maximiser in a positive definite piece.** `ψ` is strictly
concave there, so two maximisers would be beaten by their midpoint. -/
theorem maxOn_unique_of_posDef {p : QuaPiece} {s : Plane} (ha : 0 < p.q.a)
    (hD : 0 < p.q.hessDet) {x y : Plane} (hx : x ∈ p.T) (hy : y ∈ p.T)
    (hxm : IsMaxOn (psi p.q s) p.T x) (hyv : psi p.q s y = psi p.q s x) : y = x := by
  by_contra hne
  set d : Plane := y - x with hd
  have hd0 : d ≠ 0 := sub_ne_zero.2 hne
  have hα : 0 < p.q.alongCurv d := alongCurv_pos ha hD hd0
  -- `ψ(y) = ψ(x)` forces the slope along `d` to be half the curvature
  have hxy : x + (1 : ℝ) • d = y := by rw [hd]; module
  have hone := psi_along_dir p.q s x d 1
  rw [hxy, hyv] at hone
  have hslope : dot s d - dot (p.q.gradAt x) d = p.q.alongCurv d / 2 := by linarith
  -- then the midpoint strictly beats `x`
  have hmid : x + (1 / 2 : ℝ) • d ∈ p.T := by
    have : x + (1 / 2 : ℝ) • d = (1 / 2 : ℝ) • x + (1 / 2 : ℝ) • y := by rw [hd]; module
    rw [this]
    exact p.convex_T hx hy (by norm_num) (by norm_num) (by norm_num)
  have hhalf := psi_along_dir p.q s x d (1 / 2)
  rw [hslope] at hhalf
  have hle : psi p.q s (x + (1 / 2 : ℝ) • d) ≤ psi p.q s x := hxm hmid
  rw [hhalf] at hle
  nlinarith

/-- **Row 1, the affine part.** The interior maximiser is an affine function of
`s`, so the cell it sweeps out is an affine image. -/
theorem interiorPoint_combo (q : Quad) (hD : q.hessDet ≠ 0) {a b : ℝ} (hab : a + b = 1)
    (s₁ s₂ : Plane) :
    interiorPoint q (a • s₁ + b • s₂) = a • interiorPoint q s₁ + b • interiorPoint q s₂ := by
  simp only [interiorPoint, Prod.ext_iff, Prod.fst_add, Prod.snd_add, Prod.smul_fst,
    Prod.smul_snd, smul_eq_mul]
  constructor
  · field_simp
    linear_combination (2 * q.c * q.d - q.b * q.e) * hab
  · field_simp
    linear_combination (2 * q.a * q.e - q.b * q.d) * hab

/-! #### Rows 2 and 3: the degenerate images -/

/-- **Row 2.** The edge maximiser lies on the fixed line through `v` in direction
`w - v`, whatever `s` is — so the cell it sweeps out is contained in a line. -/
theorem edgePoint_mem_line (q : Quad) (v w s : Plane) :
    ∃ t : ℝ, edgePoint q v w s = v + t • (w - v) :=
  ⟨edgeSlope q v w s / edgeCurv q v w, rfl⟩

/-- **Row 3.** The vertex maximiser does not move with `s` at all: the cell it
sweeps out is the single point `v`. -/
theorem vertexBranch_maximiser_const (q : Quad) (v : Plane) (s : Plane) :
    (vertexBranch q v).eval s = psi q s v := vertexBranch_eval q v s

end QuaPol

/-! #### The risk C5 flagged, and it is real

`TODO.md` C5 asked whether "exactly one candidate is active at `s`" forces the
maximiser to be unique, and said to test it before building on it. **It is
false.** The witness is the piece already used for the `⊤` cell: the nonnegative
`s₁`-axis carrying the zero quadratic. At `s = 0` it has exactly one candidate,
that candidate is active, and *every* point of the axis is a maximiser.

So rows 1 to 3 of Theorem 4 cannot be indexed by the activity pattern. They are
indexed above by the branch attained and the rank of its quadratic part, which is
what the table in `../CONJ_FIELD_PROOF.md` actually says. -/

namespace Sanity

/-- The only candidate of the ray piece is the zero quadratic: its single vertex
gives the zero vertex branch, the zero quadratic has no direction of positive
curvature so there are no edge branches, and its Hessian is singular so there is
no interior branch. -/
theorem cand_rayPol : cand rayPol ⊆ {(0 : Quad)} := by
  intro g hg
  obtain ⟨p, hp, hgp⟩ := Finset.mem_biUnion.1 hg
  have hpe : p = rayPiece := by simpa [rayPol] using hp
  subst hpe
  simp only [QuaPiece.branches, Finset.mem_union] at hgp
  rcases hgp with (hv | he) | hi
  · -- the vertex branch of the zero quadratic at the origin is `0`
    obtain ⟨v, hv1, rfl⟩ := Finset.mem_image.1 hv
    have hv0 : v = 0 := by simpa [rayPiece] using hv1
    subst hv0
    simp [Finset.mem_singleton, vertexBranch, rayPiece, Quad.eq_zero_iff, Quad.eval]
  · -- no edge pair survives the positive-curvature filter
    obtain ⟨vw, hvw, _⟩ := Finset.mem_image.1 he
    have hcurv : 0 < edgeCurv rayPiece.q vw.1 vw.2 := (Finset.mem_filter.1 hvw).2
    simp only [rayPiece, edgeCurv, Quad.alongCurv] at hcurv
    norm_num at hcurv
  · -- the Hessian is singular
    have hdet : rayPiece.q.hessDet = 0 := by simp [rayPiece, Quad.hessDet]
    simp only [QuaPiece.interiorPart, if_pos hdet, Finset.notMem_empty] at hi

/-- **At most one candidate is active** at the origin. -/
theorem active_rayPol_subsingleton {g h : Quad} (hg : g ∈ active rayPol 0)
    (hh : h ∈ active rayPol 0) : g = h := by
  have hg' := cand_rayPol (active_subset_cand _ _ hg)
  have hh' := cand_rayPol (active_subset_cand _ _ hh)
  rw [Finset.mem_singleton] at hg' hh'
  rw [hg', hh']

/-- and yet the origin and `(1,0)` are both maximisers. -/
theorem maxSet_rayPol_not_subsingleton :
    ((0, 0) : Plane) ∈ QuaPol.maxSet rayPol 0 ∧ ((1, 0) : Plane) ∈ QuaPol.maxSet rayPol 0 ∧
      ((0, 0) : Plane) ≠ ((1, 0) : Plane) := by
  have hconj : rayPol.conj 0 = 0 := by
    refine le_antisymm (iSup_le fun x => ?_) ?_
    · by_cases hx : rayPol.eval x = ⊤
      · rw [hx, EReal.sub_top]
        exact bot_le
      · obtain ⟨p, hp, hxT, hev⟩ := exists_piece_eq_eval hx
        have hpe : p = rayPiece := by simpa [rayPol] using hp
        subst hpe
        rw [hev]
        simp [rayPiece, dot, Quad.eval]
    · refine le_iSup_of_le ((0, 0) : Plane) ?_
      rw [eval_rayPol (le_refl (0 : ℝ))]
      simp [dot]
  have hmem : ∀ t : ℝ, 0 ≤ t → ((t, 0) : Plane) ∈ QuaPol.maxSet rayPol 0 := by
    intro t ht
    refine ⟨rayPiece, by simp [rayPol], mem_rayPiece ht, ?_⟩
    rw [hconj]
    simp [psi, dot, rayPiece, Quad.eval]
  refine ⟨hmem 0 le_rfl, hmem 1 zero_le_one, ?_⟩
  simp [Prod.ext_iff]

/-- **The refutation, packaged.** There is a `QuaPol` and a point at which at
most one candidate is active and yet the maximiser is not unique. -/
theorem exists_oneActive_manyMaximisers :
    ∃ (f : QuaPol) (s : Plane),
      (∀ g h : Quad, g ∈ active f s → h ∈ active f s → g = h) ∧
      ∃ x y : Plane, x ∈ QuaPol.maxSet f s ∧ y ∈ QuaPol.maxSet f s ∧ x ≠ y := by
  obtain ⟨h0, h1, hne⟩ := maxSet_rayPol_not_subsingleton
  exact ⟨rayPol, 0, fun g h hg hh => active_rayPol_subsingleton hg hh, _, _, h0, h1, hne⟩

end Sanity


end QuaConProof

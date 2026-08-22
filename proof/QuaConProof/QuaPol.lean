/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.Conic
import Mathlib.Analysis.Convex.Topology
import Mathlib.Data.EReal.Operations

/-!
# `QuaPol` and its Fenchel conjugate

A **`QuaPol`** is a piecewise linear-quadratic function on the plane: finitely
many polyhedral pieces, each carrying one quadratic, and `+∞` off their union.
This is Rockafellar–Wets Definition 10.20 restricted to two dimensions, with the
pieces presented by their vertices.

Its **conjugate** is `f*(s) = sup_x ⟨s,x⟩ - f(x)`, valued in `EReal`.

## Design, and why it differs from `PROJECT_PLAN.md` §0.2

* **Pieces are convex hulls of finite vertex sets** (V-representation), not
  intersections of half-spaces. mathlib has no Minkowski–Weyl, and this way the
  hull is compact for free (`QuaPiece.isCompact_T`), which is what makes the
  per-piece supremum attained.
* **A piece's quadratic is a `Quad`.** The plan wrote a separate self-adjoint
  operator `Q`, a vector `β` and a scalar `γ`; reusing `Quad` makes symmetry of
  the Hessian *structural* rather than a hypothesis to carry, and lets the whole
  file share `Quad`'s algebra. In `Quad` terms the Hessian is `[[2a, b], [b, 2c]]`
  and the gradient is `∇q(x) = (2a x₁ + b x₂ + d, b x₁ + 2c x₂ + e)`.
* **`f` is the pointwise `min` over pieces**, with **no** consistency hypothesis:
  pieces may overlap, and `f` need not be continuous. The conjugate depends only
  on that `min`, so nothing is lost, and the theorem is strictly stronger.
* **`Finset.inf` over the pieces**, so an empty piece list gives `f ≡ ⊤` and
  hence `f* ≡ ⊥` — the correct convex-analytic reading of "no pieces".

## Main definitions

* `dot` — the pairing `⟨s,x⟩`. A plain function; no `InnerProductSpace` is used
  anywhere in this development (`DECISIONS.md`, 2026-08-22).
* `QuaPiece`, `QuaPol`, `QuaPol.eval`, `QuaPol.conj`.

## Main results

* `QuaPiece.isCompact_T` — every piece is compact. This is the hypothesis the
  selection lemma's attainment step will consume.
* `QuaPol.conj_pt` — the conjugate of a one-point piece is exactly the
  affine "vertex branch" `s ↦ ⟨s,v⟩ - q(v)`. This is the first candidate branch
  of `PROJECT_PLAN.md` §0.4, computed from the definition rather than asserted,
  and it is the main guard against a mis-stated `conj`.
-/

namespace QuaConProof

open scoped Classical

/-- The pairing `⟨s, x⟩ = s₁x₁ + s₂x₂` on the plane. -/
def dot (s x : Plane) : ℝ := s.1 * x.1 + s.2 * x.2

lemma dot_comm (s x : Plane) : dot s x = dot x s := by simp only [dot]; ring

@[simp] lemma dot_zero (s : Plane) : dot s 0 = 0 := by simp [dot]

/-- One polyhedral piece of a `QuaPol`: a nonempty finite vertex set, and the
quadratic carried on its convex hull. -/
structure QuaPiece where
  /-- the vertices whose convex hull is the piece -/
  verts : Finset Plane
  /-- a piece is nonempty; without this its conjugate would be `⊥` -/
  verts_nonempty : verts.Nonempty
  /-- the quadratic carried on the piece -/
  q : Quad

namespace QuaPiece

noncomputable instance : DecidableEq QuaPiece := Classical.decEq _

/-- The piece itself: the convex hull of its vertices. -/
def T (p : QuaPiece) : Set Plane := convexHull ℝ (↑p.verts : Set Plane)

lemma subset_T (p : QuaPiece) : (↑p.verts : Set Plane) ⊆ p.T := subset_convexHull ℝ _

lemma T_nonempty (p : QuaPiece) : p.T.Nonempty :=
  let ⟨v, hv⟩ := p.verts_nonempty
  ⟨v, p.subset_T hv⟩

/-- **Every piece is compact.** The convex hull of a finite set is compact, so
this needs no hypothesis — which is exactly why the V-representation was chosen.
The per-piece supremum in `conj` is therefore attained. -/
lemma isCompact_T (p : QuaPiece) : IsCompact p.T :=
  (p.verts.finite_toSet).isCompact_convexHull ℝ

/-- A one-vertex piece is the single point. -/
@[simp] lemma T_singleton (v : Plane) (h : ({v} : Finset Plane).Nonempty) (q : Quad) :
    (⟨{v}, h, q⟩ : QuaPiece).T = {v} := by
  simp only [T, Finset.coe_singleton, convexHull_singleton]

end QuaPiece

/-- A **`QuaPol`**: finitely many pieces. The function is the pointwise minimum
of the pieces' quadratics over the pieces containing the point, and `⊤` where no
piece does. -/
structure QuaPol where
  /-- the pieces -/
  pieces : Finset QuaPiece
  /-- there is at least one piece -/
  pieces_nonempty : pieces.Nonempty

namespace QuaPol

/-- The value of a `QuaPol`, in `EReal`.

`⊤` off the union of the pieces. Where several pieces overlap the smallest value
wins; no agreement between overlapping pieces is required. -/
noncomputable def eval (f : QuaPol) (x : Plane) : EReal :=
  f.pieces.inf (fun p => if x ∈ p.T then ((p.q.eval x : ℝ) : EReal) else ⊤)

/-- The **Fenchel conjugate**, `f*(s) = sup_x ⟨s,x⟩ - f(x)`.

Points where `f = ⊤` contribute `⊥` to the supremum, since `x - ⊤ = ⊥` in
`EReal`, so they do not constrain it — which is the intended reading. -/
noncomputable def conj (f : QuaPol) (s : Plane) : EReal :=
  ⨆ x : Plane, ((dot s x : ℝ) : EReal) - f.eval x

lemma conj_def (f : QuaPol) (s : Plane) :
    f.conj s = ⨆ x : Plane, ((dot s x : ℝ) : EReal) - f.eval x := rfl

/-- Off every piece, `f` is `⊤`. -/
lemma eval_eq_top_of_notMem {f : QuaPol} {x : Plane} (h : ∀ p ∈ f.pieces, x ∉ p.T) :
    f.eval x = ⊤ := by
  rw [eval, Finset.inf_eq_top_iff]
  intro p hp
  simp [h p hp]

/-- On a piece, `f` is at most that piece's quadratic. -/
lemma eval_le_of_mem {f : QuaPol} {x : Plane} {p : QuaPiece} (hp : p ∈ f.pieces)
    (hx : x ∈ p.T) : f.eval x ≤ ((p.q.eval x : ℝ) : EReal) := by
  refine le_trans (Finset.inf_le hp) ?_
  simp [hx]

/-! ### A one-piece, one-point `QuaPol`

The smallest nontrivial instance, and the one that pins down `conj`. -/

section Singleton

variable (v : Plane) (q : Quad)

/-- The `QuaPol` with a single piece, that piece being the single point `v`
carrying `q`. -/
noncomputable def pt : QuaPol :=
  ⟨{⟨{v}, Finset.singleton_nonempty v, q⟩}, Finset.singleton_nonempty _⟩

@[simp] lemma eval_pt (x : Plane) :
    (pt v q).eval x = if x = v then ((q.eval x : ℝ) : EReal) else ⊤ := by
  simp only [pt, eval, Finset.inf_singleton, QuaPiece.T_singleton, Set.mem_singleton_iff]

/-- **The conjugate of a one-point piece is the vertex branch.**

`f*(s) = ⟨s,v⟩ - q(v)`, computed from the definition of the supremum. This is
branch one of `PROJECT_PLAN.md` §0.4, and the check that `conj` means what it is
supposed to mean: a sign error or a swapped argument in `dot`, `eval` or `conj`
would show up here. -/
theorem conj_pt (s : Plane) :
    (pt v q).conj s = ((dot s v - q.eval v : ℝ) : EReal) := by
  refine le_antisymm (iSup_le fun x => ?_) (le_iSup_of_le v ?_)
  · by_cases hx : x = v
    · subst hx; simp [eval_pt]
    · simp [eval_pt, hx, EReal.sub_top]
  · simp [eval_pt]

end Singleton

/-! ### Sanity checks

`CLAUDE.md` → Verification, point 3. `conj_pt` above is the structural check;
these are the numeric ones, hand-computed. -/

namespace Sanity

/-- `dot` is not symmetric-by-accident in a way that would hide a swap:
`⟨(3,5), (7,11)⟩ = 21 + 55 = 76`. -/
example : dot (3, 5) (7, 11) = 76 := by norm_num [dot]

/-- The single point `v = (1,2)` carrying the constant `3`. At `s = (5,7)`,
`f*(s) = 5·1 + 7·2 - 3 = 16`. -/
example : (pt (1, 2) ⟨0, 0, 0, 0, 0, 3⟩).conj (5, 7) = (16 : ℝ) := by
  rw [conj_pt]
  norm_num [dot, Quad.eval]

/-- The same point carrying `q(x) = x₁² + x₂²`, so `q(v) = 1 + 4 = 5`. At
`s = (5,7)`, `f*(s) = 19 - 5 = 14`. -/
example : (pt (1, 2) ⟨1, 0, 1, 0, 0, 0⟩).conj (5, 7) = (14 : ℝ) := by
  rw [conj_pt]
  norm_num [dot, Quad.eval]

/-- At the origin the conjugate of a one-point piece is `-q(v)`, so a nonzero
constant is visible: here `-3`. This catches a dropped constant term. -/
example : (pt (1, 2) ⟨0, 0, 0, 0, 0, 3⟩).conj (0, 0) = (-3 : ℝ) := by
  rw [conj_pt]
  norm_num [dot, Quad.eval]

/-- A `QuaPol` is `⊤` away from its pieces: the one-point `QuaPol` at `(1,2)`
takes the value `⊤` at `(0,0)`. -/
example : (pt (1, 2) ⟨0, 0, 0, 0, 0, 3⟩).eval (0, 0) = ⊤ := by
  rw [eval_pt]
  norm_num

end Sanity

end QuaPol

end QuaConProof

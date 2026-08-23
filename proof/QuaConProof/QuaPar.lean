/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.Witness

/-!
# The conjugate's subdivision need not be parabolic

This is the question the whole project exists to answer. CCA2's `QuaPar` stores a
dual subdivision by requiring **every edge conic to satisfy `b² - 4ac = 0`** — a
parabola, or a line as the all-zero degenerate case. [JOGO] Theorem 6 asserts
that the conjugate always admits such a subdivision.

`conj_isQuaCon` establishes that every multi-active cell lies inside *a* conic.
This file settles which conics are forced, and exhibits a `QuaPol` for which the
parabolic restriction fails.

## What is forced, exactly

From `Shapes.lean`, for two simultaneously active candidates:

| pair | `disc` | `QuaPar`-legal? |
|---|---|---|
| vertex, vertex | `0` | **yes**, always |
| vertex, edge | `0` | **yes**, always |
| vertex, interior | `-1/hessDet ≠ 0` | **no**, never |
| edge, edge | `cross(δ₁,δ₂)²/(α₁α₂)` | only when the directions are parallel |
| edge, interior · interior, interior | unconstrained | not forced |

So the sound core of [JOGO] Theorem 6's argument is precisely
`disc_eq_zero_of_vertexBranch_flat` below: **a vertex branch against any flat
branch is parabolic**. That is what "when we compare two functions we always get
one of them as linear" delivers, and it is correct as far as it goes.

Two things break it, and the second is the serious one:

* **vertex against interior is never parabolic** — `disc = -1/hessDet` is nonzero
  as soon as the Hessian is nonsingular. But by
  `det3_interiorBranch_sub_vertexBranch_self` a *same-piece* such pair is
  degenerate with `disc < 0`, i.e. a single point. A critic may fairly say a
  point is not an edge, so this alone does not refute the claim.
* **two interior branches can tie along a genuine ellipse**, met at infinitely
  many points. That is `Witness.ellipse_realised`, and it is what
  `not_hasParabolicEdges` below turns into a refutation.

## The statement

`IsEdgePair` demands the two candidates be active together at **infinitely many**
points, so the refutation cannot be dismissed as a tangency.
-/

namespace QuaConProof

open scoped Classical

/-! ### What [JOGO]'s argument does establish -/

/-- **The sound core.** A vertex branch against any flat branch — another vertex
branch, or an edge branch — ties along a conic of vanishing discriminant: a
parabola or a line. Both are `QuaPar`-legal.

This is exactly the case [JOGO] Theorem 6's proof covers, and within it the
conclusion holds. -/
theorem disc_eq_zero_of_vertexBranch_flat (q : Quad) (v : Plane) {b : Quad}
    (hb : IsFlatBranch b) : (vertexBranch q v - b).disc = 0 := by
  cases hb with
  | vertex q₂ w => exact disc_vertexBranch_sub_vertexBranch q q₂ v w
  | edge q₂ z u h => exact disc_vertexBranch_sub_edgeBranch q q₂ v z u (ne_of_gt h)

/-- **A vertex branch against an interior branch is never parabolic.** The
discriminant is `-1/hessDet`, which is nonzero whenever the interior branch
exists at all. -/
theorem disc_ne_zero_vertexBranch_interiorBranch (q₁ q₂ : Quad) (v : Plane)
    (hD : q₁.hessDet ≠ 0) : (interiorBranch q₁ - vertexBranch q₂ v).disc ≠ 0 := by
  rw [disc_interiorBranch_sub_vertexBranch q₁ q₂ v hD]
  simp only [ne_eq, neg_eq_zero, one_div, inv_eq_zero]
  exact hD

/-! ### Edges, and the parabolic restriction -/

/-- Two candidates form an **edge** of the conjugate's subdivision when they are
distinct and simultaneously active at **infinitely many** points.

The infinitude is deliberate. `../doc/QuaConExample.md` §3.3 records three pairs
that "touch at a single point ... they are not edges", and
`det3_interiorBranch_sub_vertexBranch_self` explains why: a same-piece
interior-versus-vertex pair is always a degenerate conic. Requiring infinitely
many points rules those out, so a refutation built on this notion cannot be
answered by "that is only a tangency". -/
def IsEdgePair (f : QuaPol) (q₁ q₂ : Quad) : Prop :=
  q₁ ≠ q₂ ∧ {s : Plane | q₁ ∈ active f s ∧ q₂ ∈ active f s}.Infinite

/-- The conjugate's subdivision is **parabolic** when every edge conic has
vanishing discriminant. This is precisely the condition `QuaPar` enforces on its
`Ec` rows, and precisely what [JOGO] Theorem 6 claims always holds. -/
def HasParabolicEdges (f : QuaPol) : Prop :=
  ∀ q₁ q₂ : Quad, IsEdgePair f q₁ q₂ → (q₁ - q₂).disc = 0

/-- Every edge conic really is a conic, and the cells that make it up really do
lie inside it — the part of the picture that *is* unconditional. -/
theorem isConic_of_isEdgePair {f : QuaPol} {q₁ q₂ : Quad} (h : IsEdgePair f q₁ q₂) :
    IsConic {s : Plane | q₁.eval s = q₂.eval s} := isConic_eqLocus h.1

/-! ### The refutation -/

namespace Witness

lemma interiorBranch_ne : interiorBranch q1 ≠ interiorBranch q2 := by
  rw [interiorBranch_q1, interiorBranch_q2, Ne, Quad.ext_iff]
  norm_num

/-- The two interior branches of the witness form a genuine edge pair: distinct,
and simultaneously active at infinitely many points. -/
theorem isEdgePair_interiorBranches :
    IsEdgePair f (interiorBranch q1) (interiorBranch q2) :=
  ⟨interiorBranch_ne, ellipse_realised.2.1⟩

end Witness

/-- **The conjugate's subdivision need not be parabolic.**

For the two-piece `QuaPol` of `Witness.lean`, two interior branches are
simultaneously active at infinitely many points and their tie conic is a
non-degenerate **ellipse**, whose discriminant is `-1/16 ≠ 0`.

So no `QuaPar` — which requires `b² - 4ac = 0` of every edge conic — can store
this conjugate, and [JOGO] Theorem 6's conclusion does not hold in the generality
claimed. What its proof does establish is
`disc_eq_zero_of_vertexBranch_flat`. -/
theorem not_hasParabolicEdges_witness : ¬ HasParabolicEdges Witness.f := by
  intro h
  have hd := h _ _ Witness.isEdgePair_interiorBranches
  have hell := (Quad.kind_eq_ellipse_iff _).1 Witness.tie_isEllipse
  exact absurd hd (ne_of_lt hell.2.2)

/-- **The headline, stated existentially.** There is a `QuaPol` whose conjugate
does not admit a parabolic subdivision. -/
theorem exists_not_hasParabolicEdges : ∃ f : QuaPol, ¬ HasParabolicEdges f :=
  ⟨Witness.f, not_hasParabolicEdges_witness⟩

/-- and the offending edge is an ellipse, not merely a non-parabola. -/
theorem exists_elliptical_edge :
    ∃ (f : QuaPol) (q₁ q₂ : Quad), IsEdgePair f q₁ q₂ ∧ (q₁ - q₂).kind = ConicKind.ellipse :=
  ⟨Witness.f, interiorBranch Witness.q1, interiorBranch Witness.q2,
    Witness.isEdgePair_interiorBranches, Witness.tie_isEllipse⟩

end QuaConProof

/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.Rational

/-!
# Rational input: `RatQuaPol`, and the embedding into `QuaPol`

`TODO.md` items B1–B4, and `../CONJ_FIELD_PROOF.md` Lemma 1 and Theorem 1.

`Rational.lean` already carries the *computational* rational layer — `RatQuad`,
the three branch formulas over `ℚ`, and the census that reproduces the worked
example. What it did not carry is an input class: a rational `QuaPol` that can be
embedded into the real one, so that a statement about `cand` of the embedding is
a statement about the objects the main theorem is about.

This file supplies it. `RatQuaPiece` and `RatQuaPol` mirror `QuaPiece` and
`QuaPol` field for field, over `ℚ`, and `toQuaPol` embeds them.

## Computability

`RatQuaPiece` carries a nonemptiness proof, so its `DecidableEq` is supplied
classically, exactly as `QuaPiece`'s is, and `RatQuaPol` is therefore not a
computable object. That is deliberate and costs nothing: the *census* runs on
`Rational.lean`'s list-based `RatPiece`, which stays computable, while this class
exists to carry theorems.

## Main definitions

* `toPlane` — the embedding `ℚ × ℚ → Plane`.
* `RatQuaPiece`, `RatQuaPol`, `toQuaPiece`, `toQuaPol`.
-/

namespace QuaConProof

open scoped Classical

/-- A rational point of the plane, as a real one. -/
def toPlane (v : RatPt) : Plane := ((v.1 : ℝ), (v.2 : ℝ))

@[simp] lemma toPlane_fst (v : RatPt) : (toPlane v).1 = (v.1 : ℝ) := rfl
@[simp] lemma toPlane_snd (v : RatPt) : (toPlane v).2 = (v.2 : ℝ) := rfl

lemma toPlane_injective : Function.Injective toPlane := by
  intro v w h
  have h1 : (v.1 : ℝ) = (w.1 : ℝ) := congrArg Prod.fst h
  have h2 : (v.2 : ℝ) = (w.2 : ℝ) := congrArg Prod.snd h
  exact Prod.ext_iff.2 ⟨by exact_mod_cast h1, by exact_mod_cast h2⟩

/-- One piece of a rational `QuaPol`: a nonempty finite rational vertex set, a
finite set of rational recession directions, and a rational quadratic. Mirrors
`QuaPiece` field for field. -/
structure RatQuaPiece where
  /-- the vertices, over `ℚ` -/
  verts : Finset RatPt
  /-- a piece is nonempty -/
  verts_nonempty : verts.Nonempty
  /-- the recession directions, over `ℚ` -/
  rays : Finset RatPt
  /-- the quadratic carried on the piece -/
  q : RatQuad

namespace RatQuaPiece

noncomputable instance : DecidableEq RatQuaPiece := Classical.decEq _

/-- The real piece with the same data. -/
noncomputable def toQuaPiece (p : RatQuaPiece) : QuaPiece where
  verts := p.verts.image toPlane
  verts_nonempty := p.verts_nonempty.image _
  rays := p.rays.image toPlane
  q := p.q.toQuad

@[simp] lemma toQuaPiece_q (p : RatQuaPiece) : p.toQuaPiece.q = p.q.toQuad := rfl

@[simp] lemma toQuaPiece_verts (p : RatQuaPiece) :
    p.toQuaPiece.verts = p.verts.image toPlane := rfl

@[simp] lemma toQuaPiece_rays (p : RatQuaPiece) :
    p.toQuaPiece.rays = p.rays.image toPlane := rfl

lemma mem_toQuaPiece_verts {p : RatQuaPiece} {v : RatPt} (hv : v ∈ p.verts) :
    toPlane v ∈ p.toQuaPiece.verts :=
  Finset.mem_image.2 ⟨v, hv, rfl⟩

lemma mem_toQuaPiece_rays {p : RatQuaPiece} {r : RatPt} (hr : r ∈ p.rays) :
    toPlane r ∈ p.toQuaPiece.rays :=
  Finset.mem_image.2 ⟨r, hr, rfl⟩

end RatQuaPiece

/-- A rational `QuaPol`: finitely many rational pieces. -/
structure RatQuaPol where
  /-- the pieces -/
  pieces : Finset RatQuaPiece
  /-- there is at least one piece -/
  pieces_nonempty : pieces.Nonempty

namespace RatQuaPol

/-- The real `QuaPol` with the same data. -/
noncomputable def toQuaPol (F : RatQuaPol) : QuaPol where
  pieces := F.pieces.image RatQuaPiece.toQuaPiece
  pieces_nonempty := F.pieces_nonempty.image _

@[simp] lemma toQuaPol_pieces (F : RatQuaPol) :
    F.toQuaPol.pieces = F.pieces.image RatQuaPiece.toQuaPiece := rfl

lemma mem_toQuaPol_pieces {F : RatQuaPol} {p : RatQuaPiece} (hp : p ∈ F.pieces) :
    p.toQuaPiece ∈ F.toQuaPol.pieces :=
  Finset.mem_image.2 ⟨p, hp, rfl⟩

/-- A rational `QuaPol` with no rays is bounded. -/
lemma bounded_toQuaPol {F : RatQuaPol} (h : ∀ p ∈ F.pieces, p.rays = ∅) :
    F.toQuaPol.Bounded := by
  intro p hp
  obtain ⟨p', hp', rfl⟩ := Finset.mem_image.1 hp
  simp [RatQuaPiece.toQuaPiece, h p' hp']

end RatQuaPol

/-! ### Lemma 1: the branch formulas never leave `ℚ`

`../CONJ_FIELD_PROOF.md` Lemma 1. Each of the three branches is built from the
piece data by **field operations only**, so applying it over `ℚ` and then
embedding gives the same real quadratic as embedding and then applying it over
`ℝ`.

Note that none of the three needs a nondegeneracy hypothesis. The edge branch
divides by the curvature and the interior branch by the Hessian determinant, but
`x / 0 = 0` in both `ℚ` and `ℝ`, so the two sides agree on the degenerate values
too. The hypotheses `α ≠ 0` and `hessDet ≠ 0` are needed to know the branch
*means* something, not to know the diagram commutes. -/

namespace RatQuad

@[simp] lemma toQuad_eval (q : RatQuad) (x : RatPt) :
    (toQuad q).eval (toPlane x) = ((q.eval x : ℚ) : ℝ) := by
  simp only [toQuad, Quad.eval, eval, toPlane]
  push_cast
  ring

@[simp] lemma toQuad_hessDet (q : RatQuad) : (toQuad q).hessDet = ((q.hessDet : ℚ) : ℝ) := by
  simp only [toQuad, Quad.hessDet, hessDet]
  push_cast
  ring

lemma toQuad_gradAt (q : RatQuad) (x : RatPt) :
    (toQuad q).gradAt (toPlane x) = toPlane (q.gradAt x) := by
  simp only [toQuad, Quad.gradAt, gradAt, toPlane, Prod.mk.injEq]
  constructor <;> push_cast <;> ring

end RatQuad

lemma toPlane_sub (v w : RatPt) :
    toPlane w - toPlane v = toPlane (w.1 - v.1, w.2 - v.2) := by
  simp only [toPlane, Prod.ext_iff, Prod.fst_sub, Prod.snd_sub]
  constructor <;> push_cast <;> ring

lemma toQuad_ratEdgeCurv (q : RatQuad) (v w : RatPt) :
    edgeCurv (RatQuad.toQuad q) (toPlane v) (toPlane w) = ((ratEdgeCurv q v w : ℚ) : ℝ) := by
  simp only [edgeCurv, ratEdgeCurv, RatQuad.toQuad, Quad.alongCurv, RatQuad.alongCurv,
    toPlane, Prod.fst_sub, Prod.snd_sub]
  push_cast
  ring

/-- **Lemma 1, vertex branch.** -/
theorem toQuad_ratVertexBranch (q : RatQuad) (v : RatPt) :
    RatQuad.toQuad (ratVertexBranch q v) = vertexBranch (RatQuad.toQuad q) (toPlane v) := by
  ext <;>
    simp only [RatQuad.toQuad, ratVertexBranch, vertexBranch, RatQuad.eval, Quad.eval,
      toPlane] <;>
    push_cast <;> ring

/-- **Lemma 1, edge branch.** -/
theorem toQuad_ratEdgeBranch (q : RatQuad) (v w : RatPt) :
    RatQuad.toQuad (ratEdgeBranch q v w)
      = edgeBranch (RatQuad.toQuad q) (toPlane v) (toPlane w) := by
  ext <;>
    simp only [RatQuad.toQuad, ratEdgeBranch, edgeBranch, edgeCurv, ratEdgeCurv,
      Quad.alongCurv, RatQuad.alongCurv, Quad.gradAt, RatQuad.gradAt, RatQuad.eval, Quad.eval,
      dot, toPlane, Prod.fst_sub, Prod.snd_sub] <;>
    push_cast <;> ring

/-- **Lemma 1, interior branch.** -/
theorem toQuad_ratInteriorBranch (q : RatQuad) :
    RatQuad.toQuad (ratInteriorBranch q) = interiorBranch (RatQuad.toQuad q) := by
  ext <;>
    simp only [RatQuad.toQuad, ratInteriorBranch, interiorBranch, Quad.hessDet,
      RatQuad.hessDet] <;>
    push_cast <;> ring

/-! ### Theorem 1: the whole candidate family stays rational

`../CONJ_FIELD_PROOF.md` Theorem 1. A maximum **selects** among the functions of
Lemma 1; it never manufactures a new one. So for a rational input every candidate
quadratic is rational, and — the half `Shapes.lean` and `QuaPar.lean` actually
use — every tie conic `q₁ - q₂` has rational coefficients. -/

lemma toPlane_add (v r : RatPt) :
    toPlane v + toPlane r = toPlane (v.1 + r.1, v.2 + r.2) := by
  simp only [toPlane, Prod.ext_iff, Prod.fst_add, Prod.snd_add]
  constructor <;> push_cast <;> ring

namespace Quad

/-- A real quadratic **is rational** when it is the embedding of a rational one. -/
def IsRat (g : Quad) : Prop := ∃ r : RatQuad, g = RatQuad.toQuad r

lemma IsRat.sub {g h : Quad} (hg : g.IsRat) (hh : h.IsRat) : (g - h).IsRat := by
  obtain ⟨r, rfl⟩ := hg
  obtain ⟨t, rfl⟩ := hh
  refine ⟨r - t, ?_⟩
  ext <;> simp only [RatQuad.toQuad, Quad.sub_a, Quad.sub_b, Quad.sub_c, Quad.sub_d,
    Quad.sub_e, Quad.sub_f, RatQuad.sub_a, RatQuad.sub_b, RatQuad.sub_c, RatQuad.sub_d,
    RatQuad.sub_e, RatQuad.sub_f] <;> push_cast <;> ring

end Quad

/-- Every edge-branch candidate of an embedded piece is rational. -/
lemma isRat_edgeBranch_of_mem_edgeCandidates {p : RatQuaPiece} {vw : Plane × Plane}
    (h : vw ∈ p.toQuaPiece.edgeCandidates) :
    (edgeBranch p.toQuaPiece.q vw.1 vw.2).IsRat := by
  classical
  rcases Finset.mem_union.1 h with h1 | h2
  · obtain ⟨hv, hw⟩ := Finset.mem_product.1 h1
    obtain ⟨v', _, hveq⟩ := Finset.mem_image.1 hv
    obtain ⟨w', _, hweq⟩ := Finset.mem_image.1 hw
    refine ⟨ratEdgeBranch p.q v' w', ?_⟩
    rw [toQuad_ratEdgeBranch, hveq, hweq]
    rfl
  · obtain ⟨vr, hvr, hvreq⟩ := Finset.mem_image.1 h2
    obtain ⟨hv, hr⟩ := Finset.mem_product.1 hvr
    obtain ⟨v', _, hveq⟩ := Finset.mem_image.1 hv
    obtain ⟨r', _, hreq⟩ := Finset.mem_image.1 hr
    refine ⟨ratEdgeBranch p.q v' (v'.1 + r'.1, v'.2 + r'.2), ?_⟩
    rw [toQuad_ratEdgeBranch, ← toPlane_add, hveq, hreq]
    rw [← hvreq]
    rfl

/-- **Theorem 1.** Every candidate quadratic of a rational `QuaPol` is rational. -/
theorem isRat_of_mem_cand (F : RatQuaPol) {g : Quad} (hg : g ∈ cand F.toQuaPol) :
    g.IsRat := by
  classical
  obtain ⟨p, hp, hgp⟩ := Finset.mem_biUnion.1 hg
  obtain ⟨p', _, rfl⟩ := Finset.mem_image.1 hp
  simp only [QuaPiece.branches, Finset.mem_union] at hgp
  rcases hgp with (hv | he) | hi
  · obtain ⟨v, hv1, rfl⟩ := Finset.mem_image.1 hv
    obtain ⟨v', _, rfl⟩ := Finset.mem_image.1 hv1
    exact ⟨ratVertexBranch p'.q v', (toQuad_ratVertexBranch _ _).symm⟩
  · obtain ⟨vw, hvw, rfl⟩ := Finset.mem_image.1 he
    exact isRat_edgeBranch_of_mem_edgeCandidates (Finset.mem_filter.1 hvw).1
  · by_cases hdet : (p'.toQuaPiece).q.hessDet = 0
    · simp only [QuaPiece.interiorPart, if_pos hdet, Finset.notMem_empty] at hi
    · simp only [QuaPiece.interiorPart, if_neg hdet, Finset.mem_singleton] at hi
      subst hi
      exact ⟨ratInteriorBranch p'.q, (toQuad_ratInteriorBranch _).symm⟩

/-- **Theorem 1, the corollary that gets used.** Every tie conic of a rational
`QuaPol` has rational coefficients — so the classification of `Shapes.lean` and
the refutation of `QuaPar.lean` are decidable questions about rational data, not
merely real ones. -/
theorem isRat_tie_conic (F : RatQuaPol) {g h : Quad}
    (hg : g ∈ cand F.toQuaPol) (hh : h ∈ cand F.toQuaPol) : (g - h).IsRat :=
  (isRat_of_mem_cand F hg).sub (isRat_of_mem_cand F hh)

/-- and every active quadratic is a candidate, so `f*` is a rational quadratic on
every cell of a rational input. -/
theorem isRat_of_mem_active (F : RatQuaPol) {s : Plane} {g : Quad}
    (hg : g ∈ active F.toQuaPol s) : g.IsRat :=
  isRat_of_mem_cand F (active_subset_cand _ _ hg)

/-! ### Sanity

`CLAUDE.md` → Verification, point 3. The embedding must not quietly lose a
vertex, a ray, or a coefficient. -/

namespace Sanity

/-- The embedding of a rational point is the point with the same coordinates. -/
example : toPlane (1 / 2, 3) = ((1 / 2 : ℝ), (3 : ℝ)) := by
  simp [toPlane]

/-- A rational triangle carrying `x₁² + x₂²`. -/
noncomputable def ratTri : RatQuaPiece :=
  ⟨{(0, 0), (1, 0), (0, 1)}, ⟨(0, 0), by simp⟩, ∅, ⟨1, 0, 1, 0, 0, 0⟩⟩

/-- Its quadratic embeds to the expected real one. -/
example : ratTri.toQuaPiece.q = (⟨1, 0, 1, 0, 0, 0⟩ : Quad) := by
  simp [ratTri, RatQuaPiece.toQuaPiece, RatQuad.toQuad]

/-- Its three vertices survive the embedding. -/
example : toPlane ((0 : ℚ), (0 : ℚ)) ∈ ratTri.toQuaPiece.verts :=
  RatQuaPiece.mem_toQuaPiece_verts (p := ratTri) (v := (0, 0)) (by simp [ratTri])

/-- and the embedded vertex is the real point one expects. -/
example : toPlane ((0 : ℚ), (0 : ℚ)) = ((0 : ℝ), (0 : ℝ)) := by simp [toPlane]

/-- It has no rays, so the piece it embeds to is a polytope. -/
example : ratTri.toQuaPiece.rays = ∅ := by
  simp [ratTri, RatQuaPiece.toQuaPiece]

end Sanity

end QuaConProof

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

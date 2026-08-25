/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.FrankWolfe

/-!
# The main theorem: the conjugate of a `QuaPol` is a `QuaCon`

`PROJECT_PLAN.md` §0.8.

The conjugate `f*` of a `QuaPol` on the plane is a **`QuaCon`**: the plane splits
into finitely many cells, `f*` is a single quadratic on each, and every cell on
which two or more quadratics are simultaneously active lies inside a conic.

## What is proved here, and what is not

Five of the six conjuncts are consequences of how the cells are *defined* — by
activity pattern — and are proved outright in `Candidates.lean`. The whole
mathematical content sits in one lemma:

> `selection` : some candidate quadratic attains `f*` at every point.

`selection` is **proved**, in `Selection.lean`, from the per-piece core lemma
`exists_branch_eq_max` — the case analysis S3-S9 for a single piece.

**The development is complete.** `lake build` is green, there are no `sorry`s,
and `#print axioms conj_isQuaCon` lists only `propext`, `Classical.choice` and
`Quot.sound`.

Note what `selection` is **not**. It is *not* `f* = max over cand f`, which is
false: an edge branch overshoots the true supremum whenever its stationary point
falls outside the segment (`DECISIONS.md`, 2026-08-21). Selection says only that
the value is *attained by* some member of the list, which is what the cell
construction needs and what the informal Theorem 1 actually asserts.

## Main results

* `cell_empty_eq` — the empty-activity cell is exactly where `f*` is `⊤`.
* `dom_conj_eq_univ` — at Stage 1 the conjugate is finite everywhere.
* `active_nonempty` — at every point some candidate quadratic is active. This is
  the substance: `f*` really is one of finitely many quadratics at each point.
* `conj_isQuaCon` — the theorem.
-/

namespace QuaConProof

open scoped Classical

/-! ### The cells, completed -/

/-- **The empty-activity cell is exactly where the conjugate is `⊤`.**

Both inclusions are now substantive. Where `f*` is finite, Frank-Wolfe gives a
maximiser on every piece (`attained_of_conj_ne_top`) and `selection` then names
an active candidate, so the cell is not the empty one; and where `f*` is `⊤` no
candidate can be active, since candidates take real values.

At Stage 1 both sides are empty, because `conj_ne_top` says the conjugate is
finite everywhere. With unbounded pieces the `⊤` region can be a genuine
half-plane, and this conjunct is what describes it. -/
theorem cell_empty_eq (f : QuaPol) : cell f ∅ = {s : Plane | f.conj s = ⊤} := by
  ext s
  simp only [mem_cell_iff, Set.mem_ofPred_eq]
  constructor
  · intro hact
    by_contra htop
    obtain ⟨q, hq, hqs⟩ := selection (attained_of_conj_ne_top htop)
    have hmem : q ∈ active f s := mem_active_iff.2 ⟨hq, hqs⟩
    rw [hact] at hmem
    exact Finset.notMem_empty q hmem
  · intro htop
    refine Finset.eq_empty_iff_forall_notMem.2 fun q hq => ?_
    have he := (mem_active_iff.1 hq).2
    rw [htop] at he
    exact EReal.coe_ne_top _ he

/-- Restated for the record: **at Stage 1 the conjugate is finite everywhere**, so
its domain is the whole plane and no cell carries `+∞`. -/
theorem dom_conj_eq_univ {f : QuaPol} (hb : f.Bounded) :
    {s : Plane | f.conj s ≠ ⊤} = Set.univ :=
  Set.eq_univ_of_forall fun s => conj_ne_top hb s

/-- **Wherever the conjugate is finite, some candidate quadratic is active.**

The bite of the theorem, stated separately because the six conjuncts of
`conj_isQuaCon` are each about *cells* and this one is about *points*: there is a
finite list of quadratics, computed from the input, and at every `s` of
`dom f*` the conjugate equals one of them. -/
theorem active_nonempty {f : QuaPol} {s : Plane} (h : f.conj s ≠ ⊤) :
    (active f s).Nonempty := by
  obtain ⟨q, hq, hqs⟩ := selection (attained_of_conj_ne_top h)
  exact ⟨q, mem_active_iff.2 ⟨hq, hqs⟩⟩

/-! ### The theorem -/

/-- **The conjugate of a `QuaPol` is a `QuaCon`.**

The plane is partitioned into cells; on each of them `f*` is one fixed
quadratic; only finitely many cells are nonempty; the cell of empty activity is
exactly the region where `f*` is `+∞`; and every cell carrying two or more active
quadratics lies inside a conic — a line, a parabola, an ellipse, a hyperbola, or
one of the degenerate cases, as `Conic.lean` classifies them.

What this deliberately does **not** claim, per the agreed regularity level: any
statement about dimension, connectedness, arcs, or a face-to-face CW structure.
See `DECISIONS.md`, 2026-08-21. -/
theorem conj_isQuaCon (f : QuaPol) :
    -- cells are pairwise disjoint
    (∀ S T : Finset Quad, S ≠ T → Disjoint (cell f S) (cell f T))
    -- and cover the plane
  ∧ (⋃ S : Finset Quad, cell f S) = Set.univ
    -- on each cell, `f*` is one fixed quadratic
  ∧ (∀ (S : Finset Quad) (q : Quad), q ∈ S → ∀ s ∈ cell f S,
        f.conj s = ((q.eval s : ℝ) : EReal))
    -- only finitely many cells are nonempty
  ∧ {S : Finset Quad | (cell f S).Nonempty}.Finite
    -- the `+∞` region is exactly the empty-activity cell
  ∧ cell f ∅ = {s : Plane | f.conj s = ⊤}
    -- every multi-active cell sits inside a conic
  ∧ (∀ (S : Finset Quad) (q₁ q₂ : Quad), q₁ ∈ S → q₂ ∈ S → q₁ ≠ q₂ →
        cell f S ⊆ {s : Plane | q₁.eval s = q₂.eval s}
        ∧ IsConic {s : Plane | q₁.eval s = q₂.eval s}) := by
  refine ⟨fun S T h => cell_disjoint f h, iUnion_cell_eq_univ f,
    fun S q hq s hs => conj_eq_of_mem_cell hq hs, finite_nonempty_cells f,
    cell_empty_eq f, fun S q₁ q₂ h₁ h₂ hne =>
      ⟨cell_subset_eqLocus h₁ h₂, isConic_eqLocus hne⟩⟩

/-! ### Sanity: the `⊤` cell is not vacuous

`CLAUDE.md` -> Verification, point 3, for the Phase 7 definitions. At Stage 1 the
fifth conjunct of `conj_isQuaCon` was the equality of two empty sets. With a
recession direction it describes a genuine region, and here is one.

The piece is the nonnegative `s₁`-axis carrying the zero quadratic. Its conjugate
is `+∞` at `(1,0)`, because `⟨(1,0), (t,0)⟩ = t` is unbounded on the piece. A
`coneHull` that had collapsed to `{0}`, or a `T` that had dropped its rays, would
make this false. -/

namespace Sanity

/-- One piece: the nonnegative `s₁`-axis, carrying the zero quadratic. -/
noncomputable def rayPiece : QuaPiece :=
  ⟨{0}, Finset.singleton_nonempty _, {(1, 0)}, 0⟩

/-- The `QuaPol` with that single piece. -/
noncomputable def rayPol : QuaPol := ⟨{rayPiece}, Finset.singleton_nonempty _⟩

lemma mem_rayPiece {t : ℝ} (ht : 0 ≤ t) : ((t, 0) : Plane) ∈ rayPiece.T := by
  refine ⟨0, subset_convexHull ℝ _ (by simp [rayPiece]), t • ((1, 0) : Plane),
    smul_mem_coneHull (mem_coneHull_of_mem (show ((1, 0) : Plane) ∈ rayPiece.rays by
      simp [rayPiece])) ht, ?_⟩
  simp

lemma eval_rayPol {t : ℝ} (ht : 0 ≤ t) : rayPol.eval (t, 0) = 0 := by
  have hm := mem_rayPiece ht
  simp only [rayPol, QuaPol.eval, Finset.inf_singleton]
  rw [if_pos hm]
  simp [rayPiece]

/-- **The conjugate really is `⊤` there.** -/
theorem conj_rayPol_eq_top : rayPol.conj (1, 0) = ⊤ := by
  by_contra hne
  have hbot := QuaPol.conj_ne_bot rayPol (1, 0)
  set c : ℝ := (rayPol.conj (1, 0)).toReal with hc
  have hcoe : rayPol.conj (1, 0) = (c : EReal) := (EReal.coe_toReal hne hbot).symm
  have hle : ((|c| + 1 : ℝ) : EReal) ≤ rayPol.conj (1, 0) := by
    rw [QuaPol.conj_def]
    refine le_iSup_of_le ((|c| + 1, 0) : Plane) ?_
    rw [eval_rayPol (by positivity)]
    simp [dot]
  rw [hcoe, EReal.coe_le_coe_iff] at hle
  linarith [le_abs_self c]

/-- **So the empty-activity cell is inhabited**, and the fifth conjunct of
`conj_isQuaCon` is carrying weight rather than comparing two empty sets. -/
theorem cell_empty_rayPol_nonempty : (cell rayPol ∅).Nonempty :=
  ⟨(1, 0), by rw [cell_empty_eq]; exact conj_rayPol_eq_top⟩

/-- and `rayPol` is genuinely outside Stage 1. -/
theorem rayPol_not_bounded : ¬ rayPol.Bounded := by
  intro hb
  have := hb rayPiece (by simp [rayPol])
  simp [rayPiece] at this

end Sanity

end QuaConProof

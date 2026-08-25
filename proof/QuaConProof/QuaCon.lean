/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.Selection

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

namespace QuaPol

/-- A `QuaPol` never takes the value `⊥`: every piece contributes either a real
number or `⊤`, and a finite infimum of such is one of them. -/
lemma eval_ne_bot (f : QuaPol) (x : Plane) : f.eval x ≠ ⊥ := by
  rw [eval]
  induction f.pieces using Finset.induction_on with
  | empty => simp
  | insert p s _ ih =>
      rw [Finset.inf_insert]
      have hp : (if x ∈ p.T then ((p.q.eval x : ℝ) : EReal) else ⊤) ≠ ⊥ := by
        split <;> simp [EReal.coe_ne_bot]
      rcases le_total (if x ∈ p.T then ((p.q.eval x : ℝ) : EReal) else ⊤)
          (s.inf fun p => if x ∈ p.T then ((p.q.eval x : ℝ) : EReal) else ⊤) with h | h
      · rw [inf_eq_left.2 h]; exact hp
      · rw [inf_eq_right.2 h]; exact ih

/-- On a piece, `f` is finite: neither `⊤` nor `⊥`. -/
lemma eval_ne_top_of_mem {f : QuaPol} {x : Plane} {p : QuaPiece} (hp : p ∈ f.pieces)
    (hx : x ∈ p.T) : f.eval x ≠ ⊤ := by
  intro h
  have hle := eval_le_of_mem hp hx
  rw [h] at hle
  exact EReal.coe_ne_top _ (top_le_iff.1 hle)

/-- **The conjugate is never `⊥`.**

A `QuaPol` has at least one piece, and that piece at least one vertex, so `f` is
finite there; the supremum defining `f*` therefore has at least one real term.
This is what lets `cell ∅` be characterised as `{f* = ⊤}` rather than as
`{f* ∉ ℝ}`. -/
theorem conj_ne_bot (f : QuaPol) (s : Plane) : f.conj s ≠ ⊥ := by
  obtain ⟨p, hp⟩ := f.pieces_nonempty
  obtain ⟨v, hv⟩ := p.verts_nonempty
  have hmem : v ∈ p.T := p.subset_T hv
  have hcoe : ((f.eval v).toReal : EReal) = f.eval v :=
    EReal.coe_toReal (eval_ne_top_of_mem hp hmem) (eval_ne_bot f v)
  have hle : ((dot s v - (f.eval v).toReal : ℝ) : EReal) ≤ f.conj s := by
    rw [conj_def]
    refine le_iSup_of_le v ?_
    rw [EReal.coe_sub, hcoe]
  intro hbot
  rw [hbot, le_bot_iff] at hle
  exact EReal.coe_ne_bot _ hle

end QuaPol

/-! ### The cells, completed -/

/-- **The empty-activity cell is exactly where the conjugate is `⊤`.**

At Stage 1 both sides are empty: `QuaPol.conj_ne_top` says the conjugate is
finite everywhere, because every piece is compact, and `selection` says some
candidate is always active. The conjunct becomes substantive only in Phase 7,
when unbounded pieces make `dom f*` a proper subset. -/
theorem cell_empty_eq {f : QuaPol} (hb : f.Bounded) : cell f ∅ = {s : Plane | f.conj s = ⊤} := by
  have hR : {s : Plane | f.conj s = ⊤} = (∅ : Set Plane) := by
    ext s
    simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false]
    exact conj_ne_top hb s
  have hL : cell f ∅ = (∅ : Set Plane) := by
    ext s
    simp only [mem_cell_iff, Set.mem_empty_iff_false, iff_false]
    intro hact
    obtain ⟨q, hq, hqs⟩ := selection hb s
    have hmem : q ∈ active f s := mem_active_iff.2 ⟨hq, hqs⟩
    rw [hact] at hmem
    exact Finset.notMem_empty q hmem
  rw [hL, hR]

/-- Restated for the record: **at Stage 1 the conjugate is finite everywhere**, so
its domain is the whole plane and no cell carries `+∞`. -/
theorem dom_conj_eq_univ {f : QuaPol} (hb : f.Bounded) :
    {s : Plane | f.conj s ≠ ⊤} = Set.univ :=
  Set.eq_univ_of_forall fun s => conj_ne_top hb s

/-- **At every point of the plane some candidate quadratic is active.**

The bite of the theorem, stated separately because the six conjuncts of
`conj_isQuaCon` are each about *cells* and this one is about *points*: there is a
finite list of quadratics, computed from the input, and at every `s` the
conjugate equals one of them. -/
theorem active_nonempty {f : QuaPol} (hb : f.Bounded) (s : Plane) :
    (active f s).Nonempty := by
  obtain ⟨q, hq, hqs⟩ := selection hb s
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
theorem conj_isQuaCon {f : QuaPol} (hb : f.Bounded) :
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
    cell_empty_eq hb, fun S q₁ q₂ h₁ h₂ hne =>
      ⟨cell_subset_eqLocus h₁ h₂, isConic_eqLocus hne⟩⟩

end QuaConProof

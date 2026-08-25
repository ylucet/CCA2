/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.Shapes

/-!
# Certificates: pinning the conjugate at a point

`conj_isQuaCon` says every multi-active cell lies inside a conic, and
`Shapes.lean` says which conic. Neither says a given cell is **nonempty**. That
is the realisation gap: "an ellipse is in the candidate list" is not yet "`f*`
has an elliptical edge".

Closing it needs a way to *pin* `f*` at a point, and that is what this file
provides. The conjugate is a supremum, so pinning it takes two bounds:

* `le_conj` — any admissible point of any piece gives a lower bound. Free.
* `conj_le` — a bound valid on **every** piece bounds `f*`. This is where
  `exists_maximiser` earns its keep: the supremum is attained, so a bound on the
  attaining point is a bound on the whole thing.
* `conj_eq_of_isMaxOn` — the two together. A point that beats every point of
  every piece pins `f*` exactly.

For a **convex** piece there is a bound with no side condition at all:

* `psi_le_interiorBranch` — if the Hessian is positive definite, `ψ` is concave
  and its unconstrained maximum is the interior branch. So the interior branch
  bounds that piece's contribution *everywhere*, whether or not the maximiser
  lies inside the piece.

Combining these turns "`f*(s) = V` and exactly these branches are active" into a
finite list of inequalities — which, over `ℚ`, is decidable. That is the machinery
a realisation witness needs.

## What is still missing after this file

A concrete witness. For the elliptical edge the tie point is irrational
(`../CONJ_FIELD_PROOF.md` §7.5 gives `s* = (-17/62 + √612030/186, 3/2)`), so a
*rational* witness on that edge may not exist and the last step would need an
intermediate-value argument rather than a computation. `TODO.md` carries this.
-/

namespace QuaConProof

open scoped Classical

/-! ### The two bounds -/

/-- **Lower bound.** Every admissible point of every piece bounds `f*` below.
This is the easy half: each such point is one term of the supremum. -/
theorem le_conj {f : QuaPol} {p : QuaPiece} (hp : p ∈ f.pieces) {x : Plane}
    (hx : x ∈ p.T) (s : Plane) : ((psi p.q s x : ℝ) : EReal) ≤ f.conj s := by
  rw [QuaPol.conj_def]
  refine le_iSup_of_le x ?_
  have hle : f.eval x ≤ ((p.q.eval x : ℝ) : EReal) := QuaPol.eval_le_of_mem hp hx
  calc ((psi p.q s x : ℝ) : EReal)
      = ((dot s x : ℝ) : EReal) - ((p.q.eval x : ℝ) : EReal) := by
        rw [← EReal.coe_sub]; rfl
    _ ≤ ((dot s x : ℝ) : EReal) - f.eval x := EReal.sub_le_sub le_rfl hle

/-- **Upper bound.** A bound valid on every point of every piece bounds `f*`.

The supremum defining `f*` ranges over all of the plane, not just over the
pieces, so this is not immediate. It needs no attainment, though: off the pieces
`f` is `⊤` and the term is `⊥`, and on them `exists_piece_eq_eval` names a piece
realising the value. That is what keeps this usable when pieces are unbounded. -/
theorem conj_le {f : QuaPol} {s : Plane} {M : ℝ}
    (h : ∀ p ∈ f.pieces, ∀ x ∈ p.T, psi p.q s x ≤ M) : f.conj s ≤ (M : EReal) := by
  rw [QuaPol.conj_def]
  refine iSup_le fun x => ?_
  by_cases hx : f.eval x = ⊤
  · rw [hx, EReal.sub_top]; exact bot_le
  · obtain ⟨p, hp, hxT, hev⟩ := exists_piece_eq_eval hx
    rw [hev, ← EReal.coe_sub]
    exact EReal.coe_le_coe_iff.2 (h p hp x hxT)

/-- **The certificate.** A point of a piece that beats every point of every piece
pins the conjugate exactly. -/
theorem conj_eq_of_isMaxOn {f : QuaPol} {s : Plane} {p : QuaPiece} (hp : p ∈ f.pieces)
    {x : Plane} (hx : x ∈ p.T)
    (h : ∀ p' ∈ f.pieces, ∀ y ∈ p'.T, psi p'.q s y ≤ psi p.q s x) :
    f.conj s = ((psi p.q s x : ℝ) : EReal) :=
  le_antisymm (conj_le h) (le_conj hp hx s)

/-! ### Convex pieces bound themselves, everywhere -/

/-- A positive definite quadratic has nonnegative curvature in every direction.

`2a d₁² + 2b d₁d₂ + 2c d₂² = ((2a d₁ + b d₂)² + (4ac - b²) d₂²) / (2a)`. -/
theorem alongCurv_nonneg {q : Quad} (ha : 0 < q.a) (hD : 0 < q.hessDet) (d : Plane) :
    0 ≤ q.alongCurv d := by
  have hD' : 0 < 4 * q.a * q.c - q.b ^ 2 := hD
  simp only [Quad.alongCurv]
  nlinarith [sq_nonneg (2 * q.a * d.1 + q.b * d.2), sq_nonneg d.2, sq_nonneg d.1, ha, hD']

/-- **A convex piece's interior branch bounds it everywhere.**

If the Hessian is positive definite then `ψ` is concave, so its unconstrained
maximum — the interior branch — dominates `ψ` at *every* point of the plane, and
in particular at every point of the piece. No side condition: the maximiser need
not lie inside the piece.

This is what makes the certificate finitely checkable for a `QuaPol` all of whose
pieces are convex, which is the case in `../CONJ_FIELD_PROOF.md` §4.1. -/
theorem psi_le_interiorBranch {q : Quad} (ha : 0 < q.a) (hD : 0 < q.hessDet)
    (s x : Plane) : psi q s x ≤ (interiorBranch q).eval s := by
  have hD0 : q.hessDet ≠ 0 := ne_of_gt hD
  set x0 : Plane := interiorPoint q s with hx0
  have hstat : q.gradAt x0 = s := gradAt_interiorPoint q s hD0
  have hexp : psi q s (x0 + (1 : ℝ) • (x - x0))
      = psi q s x0 + 1 * (dot s (x - x0) - dot (q.gradAt x0) (x - x0))
        - 1 ^ 2 * q.alongCurv (x - x0) / 2 := psi_along_dir q s x0 (x - x0) 1
  rw [hstat] at hexp
  have hxx : x0 + (1 : ℝ) • (x - x0) = x := by
    simp only [one_smul]; abel
  rw [hxx] at hexp
  have hzero : dot s (x - x0) - dot s (x - x0) = 0 := sub_self _
  rw [hzero] at hexp
  have hcurv : 0 ≤ q.alongCurv (x - x0) := alongCurv_nonneg ha hD _
  have hval : psi q s x0 = (interiorBranch q).eval s :=
    (interiorBranch_eval q s hD0).symm
  rw [← hval]
  linarith [hexp, hcurv]

/-- **Every piece convex ⇒ the interior branches bound `f*`.**

If a bound `M` dominates every piece's interior branch at `s`, it dominates
`f*(s)`. Each hypothesis is a single inequality between explicit quadratics, so
over `ℚ` this is decidable. -/
theorem conj_le_of_convex {f : QuaPol} {s : Plane} {M : ℝ}
    (hconv : ∀ p ∈ f.pieces, 0 < p.q.a ∧ 0 < p.q.hessDet)
    (hbound : ∀ p ∈ f.pieces, (interiorBranch p.q).eval s ≤ M) :
    f.conj s ≤ (M : EReal) := by
  refine conj_le fun p hp x _ => ?_
  obtain ⟨ha, hD⟩ := hconv p hp
  exact le_trans (psi_le_interiorBranch ha hD s x) (hbound p hp)

/-- **The convex certificate.** If some piece's unconstrained maximiser actually
lies inside that piece, and no other piece's interior branch exceeds its value,
then `f*` is that interior branch at `s`.

Both hypotheses are checkable: the first is membership of an explicit point in a
convex hull, the second a finite list of inequalities between quadratics. -/
theorem conj_eq_interiorBranch {f : QuaPol} {s : Plane} {p : QuaPiece} (hp : p ∈ f.pieces)
    (hconv : ∀ p' ∈ f.pieces, 0 < p'.q.a ∧ 0 < p'.q.hessDet)
    (hin : interiorPoint p.q s ∈ p.T)
    (hmax : ∀ p' ∈ f.pieces, (interiorBranch p'.q).eval s ≤ (interiorBranch p.q).eval s) :
    f.conj s = (((interiorBranch p.q).eval s : ℝ) : EReal) := by
  obtain ⟨ha, hD⟩ := hconv p hp
  have hD0 : p.q.hessDet ≠ 0 := ne_of_gt hD
  have hval : psi p.q s (interiorPoint p.q s) = (interiorBranch p.q).eval s :=
    (interiorBranch_eval p.q s hD0).symm
  refine le_antisymm (conj_le_of_convex hconv hmax) ?_
  have := le_conj hp hin s
  rwa [hval] at this

/-! ### Reading off the activity set

With `f*` pinned, membership in `active` is a numeric test: a candidate is active
exactly when it evaluates to the pinned value. -/

/-- Once `f*(s)` is known to be the real number `V`, a candidate is active at `s`
precisely when it takes the value `V` there. -/
theorem mem_active_iff_eq {f : QuaPol} {s : Plane} {V : ℝ} (hV : f.conj s = (V : EReal))
    {q : Quad} : q ∈ active f s ↔ q ∈ cand f ∧ q.eval s = V := by
  rw [mem_active_iff, hV]
  constructor
  · rintro ⟨h1, h2⟩; exact ⟨h1, by exact_mod_cast h2⟩
  · rintro ⟨h1, h2⟩; exact ⟨h1, by exact_mod_cast h2⟩

/-- **A nonempty multi-active cell is a realised conic.** If two *distinct*
candidates are both active at `s`, then `s` lies on their tie conic and that
conic is genuinely met — the containment of `conj_isQuaCon` is not vacuous
there. -/
theorem mem_eqLocus_of_active {f : QuaPol} {s : Plane} {q₁ q₂ : Quad}
    (h₁ : q₁ ∈ active f s) (h₂ : q₂ ∈ active f s) :
    q₁.eval s = q₂.eval s := by
  have e₁ := (mem_active_iff.1 h₁).2
  have e₂ := (mem_active_iff.1 h₂).2
  exact_mod_cast e₁.trans e₂.symm

/-- The cell of the activity pattern at a witnessed point is, by construction,
nonempty — so exhibiting one `s` with a given activity set realises that cell. -/
theorem cell_nonempty_of_witness (f : QuaPol) (s : Plane) :
    (cell f (active f s)).Nonempty := ⟨s, mem_cell_active f s⟩

end QuaConProof

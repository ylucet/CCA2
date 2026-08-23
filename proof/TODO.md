# TODO

Near-term action items only. Longer-range roadmap lives in `PROJECT_PLAN.md`;
narrative status lives in `CURRENT_STATE.md`.

**Dead ends do not belong here.** Anything tried and reverted, rejected on
method grounds, or deliberately not being done goes in `DECISIONS.md`. An item
under **Blocked** below is still an action item waiting on something; once it is
abandoned rather than blocked, move it to `DECISIONS.md` and delete it here.

**Format:** items use checkboxes -- `- [ ]` open, `- [x]` done. The boxes under
**Next up** are the machine-readable record of what is left in this project, so
tick an item in the same turn you finish it. Unattended runs (the `/overnight`
skill) read this section to decide whether there is more work to do.

**Write boundary:** everything below happens inside `CCA2/proof/`. Never edit
`../TODO.md` or any other file in `../`.

## Next up

**The objective is met.** `conj_isQuaCon` is proved, `sorry`-free, with clean
axioms, for the Stage 1 input class. Everything below is optional extension, and
none of it is required by `CLAUDE.md`'s statement of the objective.

**Decide with Yves before starting any of it.** These are three different
projects, and the right next one depends on what the result is for.

- [ ] **Phase 7 -- unbounded pieces.** The one that makes the theorem cover
      Rockafellar-Wets 10.20 in full. Needs a recession cone on `QuaPiece`, a
      Frank-Wolfe style attainment result, and `dom f*` becoming a proper subset
      so the fifth conjunct stops being the equality of two empty sets
      (`DECISIONS.md`, 2026-08-22). Largest of the three.
- [ ] **Realisation, the remaining half.** `Witness.lean` proves infinitely many
      points have two interior branches active, all on one non-degenerate ellipse.
      What is still open is that a **single cell** is infinite (the set proved
      infinite is a union of cells), which needs the activity pattern to be
      locally constant. Separately, realising the *documents'* elliptical edge
      cannot be done by computation at all -- that conic has no rational point,
      by a 7-adic obstruction -- so it needs an intermediate-value argument.
- [ ] **Phase 2 remainder -- conic normal forms.** `disc` and `det3` are proved
      invariants with computed values, but there is still no theorem saying
      `disc < 0` makes a set an ellipse. Needs `Quad.rotate` and
      `sin^2+cos^2 = 1` threaded through a coefficient identity. Self-contained,
      and the only thing that would let the write-up say "ellipse" as a *proved*
      classification rather than as a discriminant sign.
- [ ] **Phase 6 -- a `QuaPol` witness.** `Conic.lean` carries the five integer
      classification targets including both traps and the parabola, but as bare
      coefficient vectors. Connecting `U3U6` back to an actual `QuaPol` whose
      conjugate has that edge would make the "`QuaPar` is too narrow" claim a
      theorem about the pipeline rather than about a hand-copied conic.
- [ ] **Phase 8 -- write-up.** `PROOF_NOTES.md` mapping each Lean lemma to its
      counterpart in `../CONJ_FIELD_PROOF.md`; then decide with Yves whether this
      becomes a paper or a note inside CCA2.

## Blocked

- [ ] Nothing is blocked.

## Done recently

- [x] 2026-08-22 -- **`Realization.lean`**: the certificate machinery. `le_conj`,
      `conj_le` (which needs attainment, so `exists_maximiser` earns its keep),
      and for convex pieces `psi_le_interiorBranch` -- the interior branch bounds
      its piece *everywhere*, with no side condition.
- [x] 2026-08-22 -- **`Witness.lean`: `ellipse_realised`.** An explicit two-piece
      `QuaPol` with infinitely many points where two interior branches are both
      active, all on one non-degenerate ellipse, with both maximisers genuinely
      inside their pieces. Needed a purpose-built example: the documents' own
      elliptical edge has **no rational point** (7-adic obstruction).

- [x] 2026-08-22 -- **`Shapes.lean`: which conic arises.** Type is decided by the
      ranks of the two branches' quadratic parts, across pieces. Headline:
      `not_flat_of_disc_neg`, an elliptical tie conic requires an interior branch
      -- the [JOGO] Theorem 6 gap as a theorem. Plus Theorem 3, which needed no
      continuity hypothesis after all.
- [x] 2026-08-22 -- **`Rational.lean`: computable classifier over ℚ** with the
      `kind_toQuad` bridge, and a census reproducing `doc/QuaConExample.md` 3.1
      (7 per piece, 23 after dedup), all ten curved edges of 3.3, and the four
      adjacent-pair degeneracies. Mostly kernel-checked; two counts use
      `native_decide` and are audited in `SORRY_LEDGER.md`.

- [x] 2026-08-22 -- **`conj_isQuaCon` PROVED.** `lake build` green, zero `sorry`,
      `#print axioms` lists only `propext`, `Classical.choice`, `Quot.sound`.
- [x] 2026-08-22 -- **S3-S9 complete**: `exists_branch_eq_max`, by strong
      induction on the number of vertices of the face carrying the maximiser.
      Compiled on the first attempt once the infrastructure was in place.
- [x] 2026-08-22 -- `Bary.lean`: `mem_convexHull_erase`, `mem_convexHull_perturb`,
      `exists_descent`, plus plane geometry by scalar cross product
      (`exists_smul_of_cross_eq_zero`, `eq_zero_of_dot_eq_zero_of_cross_ne_zero`,
      `exists_combo_of_cross_ne_zero`, `cross_eq_zero_of_forall`). Using the cross
      product instead of `AffineIndependent` / `Collinear` avoided a large slice
      of mathlib API and made the case split concrete.
- [x] 2026-08-22 -- `foc_soc` (S5 and the second-order condition) and
      `interiorPoint_unique`.
- [x] 2026-08-22 -- **S8's core proved**: `psi_along_dir`,
      `exists_kernel_of_hessDet_eq_zero`, `psi_const_along_kernel`,
      `exists_dir_psi_const`.
- [x] 2026-08-22 -- **S1 and S2 proved**: `exists_isMaxOn_piece`,
      `exists_piece_eq_eval`, `exists_maximiser`; and `conj_ne_top`.
- [x] 2026-08-22 -- **Phase 3**: `QuaPol.lean`, `Candidates.lean`, `QuaCon.lean`;
      the three branch identities, each differential-tested before being written
      into Lean and then proved as a Lean identity.
- [x] 2026-08-22 -- **Phases 1 and 2 core**: Lake project on mathlib `v4.33.0`
      from local cache; `Quad.lean`, `Conic.lean`, eleven sanity witnesses.
- [x] 2026-08-21 -- target statement agreed and scaffolded.

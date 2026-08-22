# TODO

Near-term action items only. Longer-range roadmap lives in `PROJECT_PLAN.md`;
narrative status lives in `CURRENT_STATE.md`.

**Dead ends do not belong here.** Anything tried and reverted, rejected on
method grounds, or deliberately not being done goes in `DECISIONS.md`. An item
under **Blocked** below is still an action item waiting on something; once it is
abandoned rather than blocked, move it to `DECISIONS.md` and delete it here.

**Format:** items use checkboxes — `- [ ]` open, `- [x]` done. The boxes under
**Next up** are the machine-readable record of what is left in this project, so
tick an item in the same turn you finish it. Unattended runs (the `/overnight`
skill) read this section to decide whether there is more work to do.

**Write boundary:** everything below happens inside `CCA2/proof/`. Never edit
`../TODO.md` or any other file in `../`.

## Next up

**Phase 4 - one lemma left: `exists_branch_eq_max`.** S1, S2 and S8's core are
proved; `selection` and the whole outer structure are proved from it. What
remains is the barycentric bookkeeping inside a single piece. `SORRY_LEDGER.md`
has the full entry, including which supporting lemmas are already in hand.

- [ ] **S3a** minimal Caratheodory subset => strictly positive barycentric
      coordinates. `Caratheodory.minCardFinsetOfMemConvexHull` gives the subset;
      the content is turning card-minimality into positivity (if a coordinate is
      zero, the point lies in the hull of the smaller set)
- [ ] **S3b** affine independence in the plane => `W.card <= 3`
- [ ] **S5** positive coordinates => `x +- t*d` stays in `conv W` for small `t`,
      for `d` in the direction space; then `psi_along_dir` and
      `eq_zero_of_forall_small` (both proved) give the first-order condition
- [ ] **S6/S7** `W.card = 1` and `= 2`. The branch identities they need are
      proved; S7 also needs `alpha >= 0` from the second-order part of
      `psi_along_dir`
- [ ] **S8 descent** barycentric coordinates are affine in the parameter, so the
      first zero gives a maximiser on a proper face. The mathematics is proved
      (`psi_const_along_kernel`); this is the bookkeeping around it
- [ ] **S9** strong induction on `W.card`, then `exists_branch_eq_max` sorry-free
- [ ] then `#print axioms conj_isQuaCon` clean and `SORRY_LEDGER.md` at zero

**Phase 2 remainder — the conic normal forms.** Self-contained; nothing is
blocked on it. `disc` and `det3` are currently invariants with computed values,
not a proved geometric classification.

- [ ] `Quad.rotate` and `disc`/`det3` invariance under it (needs `sin²+cos²=1`
      threaded through a coefficient identity)
- [ ] the normal-form theorems, then `ConicKind` / `Quad.kind`

**Phase 6 — witnesses.** The five integer classification targets in
`PROJECT_PLAN.md` Phase 6 are already in `Conic.lean`, including the parabola.
What remains is the harder statement:

- [ ] exhibit a concrete `QuaPol` whose conjugate genuinely has a non-degenerate
      elliptical edge, i.e. connect `U3U6` back to an actual `QuaPol` rather than
      carrying it as a bare coefficient vector

## Blocked

- [ ] Nothing is blocked. (Phase 7, unbounded pieces, is sequenced after Phase 5,
      not blocked by anything external.)

## Done recently

- [x] 2026-08-22 - **S8's core proved**, the step flagged as most likely to sink
      the plan: `psi_along_dir` (the exact second-order expansion, no remainder),
      `exists_kernel_of_hessDet_eq_zero`, `psi_const_along_kernel`,
      `exists_dir_psi_const`. Sliding along the kernel of a singular Hessian at a
      stationary point costs nothing, so the three-vertex singular case reduces to
      a smaller face and needs no candidate of its own.
- [x] 2026-08-22 - **S1 and S2 proved**: `exists_isMaxOn_piece` (extreme value
      theorem on the compact piece), `exists_piece_eq_eval`, `exists_maximiser`.
- [x] 2026-08-22 - **`conj_ne_top`**: at Stage 1 the conjugate is finite
      everywhere, so `dom f* = R^2` and `selection` needs no hypothesis. See
      `DECISIONS.md` - it also means the fifth conjunct of the main theorem is
      currently the equality of two empty sets, and is not yet load-bearing.
- [x] 2026-08-22 - `eq_zero_of_forall_small`, the scalar first-order condition.
- [x] 2026-08-22 - the single `sorry` moved from the top-level `selection` down to
      the precise per-piece statement `exists_branch_eq_max`.

- [x] 2026-08-22 — **Phase 3 complete.** `QuaPol.lean`, `Candidates.lean`,
      `QuaCon.lean`. `conj_isQuaCon` is stated in full and proved modulo the
      single `selection` sorry; the other five conjuncts are proved outright.
- [x] 2026-08-22 — the three candidate branch formulas differential-tested
      against direct optimisation before being written into Lean, then proved in
      Lean as identities (`vertexBranch_eval`, `edgeBranch_eval`,
      `interiorBranch_eval`), plus `psi_le_edgeBranch` for maximality along the
      line and `gradAt_interiorPoint` for the stationarity of the interior point.
- [x] 2026-08-22 — `QuaPol.conj_pt`: the conjugate of a one-point piece computed
      from the definition, the main guard against a mis-stated `conj`.

- [x] 2026-08-22 — **Phase 1 complete.** Lake project on mathlib `v4.33.0`,
      manifest pinned to the same revisions as `AI/lean` so `lake exe cache get`
      hit the local cache (8689 files, nothing downloaded). `lake build` green.
- [x] 2026-08-22 — **Phase 2 core complete, 0 `sorry`.** `Quad.lean` and
      `Conic.lean`: coefficient vectors, `eval`, the everywhere-vanishing lemma,
      translation, `IsConic`, `IsConic.ne_univ`, `isConic_eqLocus`, `disc`,
      `det3` pinned against the `3×3` determinant, invariance under translation
      and covariance under scaling. `#print axioms` clean on all eight top-level
      results.
- [x] 2026-08-22 — sanity witnesses for all four conic types plus the degenerate
      cases, including the two edges that *look* straight and are hyperbolas, and
      a parabola witness verified independently against brute-force maximisation
      before being written into Lean.

- [x] 2026-08-21 — target statement agreed with Yves and written down as
      `PROJECT_PLAN.md` Phase 0: full geometric subdivision, activity-pattern
      cells, cover/disjoint/quadratic-on-cell/conic-containment regularity,
      proved discriminant trichotomy, `EReal`-valued with the `+∞` region as
      `cell ∅`, V-representation input, bounded pieces first
- [x] 2026-08-21 — project scaffolded in `CCA2/proof/` per the `/create-project`
      skill, adapted for a formalization (no `data/`, `src/`, `validation/`,
      `manuscript/`); six design decisions recorded in `DECISIONS.md`

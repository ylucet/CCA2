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

**Phase 4 — the selection lemma. This is the project.** Everything else in
`conj_isQuaCon` is proved; `selection` is the single `sorry` and the only thing
between here and a complete theorem. `SORRY_LEDGER.md` carries the full entry.

Order matters: S8 first, because it is the step most likely to fail, and nothing
should be built on top of it until it is known to work.

- [ ] **S8** `|W| = 3`, singular Hessian: `∇q(x*) = s` and `H` singular, so `psi`
      is constant along `ker H`; barycentric coordinates are affine in the
      parameter, so at the largest admissible step one coordinate vanishes and the
      maximum is also attained on a proper face. Contradicts minimality of `W`,
      descending to `|W| ≤ 2`. Write standalone, with its own example.
- [ ] **S2** attainment: `QuaPiece.isCompact_T` is already proved; combine with
      `IsCompact.exists_isMaxOn` and continuity of `psi` (polynomial, so `fun_prop`)
- [ ] **S3** minimal Carathéodory subset — `Caratheodory.minCardFinsetOfMemConvexHull`
      and `affineIndependent_minCardFinsetOfMemConvexHull` are already in mathlib
      and give exactly this; derive strict positivity of the barycentric coordinates
- [ ] **S1** the supremum over the finite union of pieces splits as a finite max
- [ ] **S5** first-order condition on the direction space of `affineSpan ℝ W`
- [ ] **S6** `|W| = 1` → `vertexBranch`; **S7** `|W| = 2` → `edgeBranch` when
      `0 < edgeCurv`, `vertexBranch` when it vanishes. The three `*_eval`
      identities these need are already proved.
- [ ] **S9** assemble by induction on `|W|`; `selection` `sorry`-free
- [ ] then `#print axioms conj_isQuaCon` must come back clean, and
      `SORRY_LEDGER.md` must go to zero

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

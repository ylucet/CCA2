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

Phase 2 — quadratics and conics. **Definitions and structure are done and green;
the geometric classification is not.** See `DECISIONS.md`, 2026-08-22, second
entry, for exactly where the line falls.

- [ ] `Quad.rotate` — precomposition with a rotation, and `disc`/`det3`
      invariance under it. Needs `sin²+cos²=1` threaded through a coefficient
      identity; `linear_combination` with that identity is the thing to try first
- [ ] the normal-form theorems: `disc < 0 ∧ det3 ≠ 0 ∧ (a+c)*det3 < 0` ⇒ the zero
      set is a translated, rotated `x²/A² + y²/B² = 1`; likewise parabola and
      hyperbola, and the degenerate branches
- [ ] `ConicKind` and `Quad.kind`, once the normal forms exist to give them
      meaning — defining the predicate before the theorem would just be a
      restatement of `disc`/`det3`

Phase 3 — the statement, sorried

- [ ] `QuaConProof/QuaPol.lean` — `QuaPiece`, `QuaPol`, `eval`, `conj`.
      Pairing is `dot s x = s.1*x.1 + s.2*x.2` on `ℝ × ℝ`; no `InnerProductSpace`
- [ ] two sanity `example`s pinning `conj` against hand computation (one piece
      that is a single point; one piece that is a segment with `Q = I`)
- [ ] `QuaConProof/Candidates.lean` — `cand`, `active`, `cell`
- [ ] `QuaConProof/QuaCon.lean` — `conj_isQuaCon` stated in full, proof `sorry`
- [ ] seed `SORRY_LEDGER.md` from that state

Phase 4 — the selection lemma, highest-risk step first

- [ ] S8 `|W| = 3` singular-`Q` descent, standalone, with its own example
- [ ] S1 supremum over a finite union splits into a finite max of suprema
- [ ] S2 per-piece attainment by compactness of `convexHull ↑verts`
- [ ] S3 Carathéodory + minimal `W` ⇒ strictly positive barycentric coordinates
- [ ] S5 first-order condition at a relative-interior maximum
- [ ] S6 `|W| = 1`; S7 `|W| = 2` both `α` branches
- [ ] S9 strong induction on `|W|`; selection lemma `sorry`-free

Phases 5-8 are in `PROJECT_PLAN.md` and are not expanded here until Phase 4 is
under way.

## Blocked

- [ ] Nothing is blocked. (Phase 7, unbounded pieces, is sequenced after Phase 5,
      not blocked by anything external.)

## Done recently

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

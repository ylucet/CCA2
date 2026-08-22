# Sorry Ledger

Every `sorry` in the Lean sources, with what it needs and what depends on it.

**Rules.**

- A commit may contain `sorry` **only if every one of them is listed here**, and
  the ledger is updated in the *same* commit as the code. A stale ledger is worse
  than no ledger.
- Never add a `sorry` that is not part of the plan in `PROJECT_PLAN.md`. In
  particular, never use `sorry` to make a build green after a refactor — fix the
  proof or revert the refactor (`CLAUDE.md` → Sorry policy).
- A `sorry` that turns out to be *unprovable as stated* is not a ledger entry.
  It is a `DECISIONS.md` entry and a change to the statement.
- Regenerate the count with:

      grep -rn "sorry" QuaConProof/ | grep -v "^.*--" | wc -l

- The project is done when this file has zero entries **and**
  `#print axioms conj_isQuaCon` lists only `propext`, `Classical.choice`,
  `Quot.sound`.

## Count

**1 sorry**, on the per-piece core lemma. Verified 2026-08-22:

    $ grep -rn "  sorry$" QuaConProof/*.lean
    QuaConProof/Selection.lean:342:  sorry

    $ #print axioms QuaConProof.conj_isQuaCon
    [propext, sorryAx, Classical.choice, Quot.sound]

Exactly three declarations depend on `sorryAx` — `exists_branch_eq_max`,
`selection`, and through them `cell_empty_eq` and `conj_isQuaCon`. Everything
else is clean, including all of S1, S2 and S8's core.

## Open

### `QuaConProof/Selection.lean` — `exists_branch_eq_max`

- **Statement:** if `x` maximises `ψ = ⟨s,·⟩ - q` over one piece `p`, then some
  branch `b ∈ p.branches` satisfies `b.eval s = ψ(x)`. This is S3–S9 of
  `PROJECT_PLAN.md` §0.6 for a **single piece**; the multi-piece assembly
  (`selection`) is proved from it.
- **Needs:** the barycentric bookkeeping, and only that. Specifically:
  - **S3** `Caratheodory.minCardFinsetOfMemConvexHull` (mathlib) gives a minimal
    affinely independent `W ⊆ verts` with `x ∈ convexHull ↑W`. Minimality must be
    converted into "every barycentric coordinate of `x` is strictly positive",
    and affine independence in the plane into `W.card ≤ 3`.
  - **S5** from positive coordinates, `x ± t·d` stays in the simplex for small
    `t`; then `psi_along_dir` and `eq_zero_of_forall_small` — **both proved** —
    give `⟨s - ∇q(x), d⟩ = 0`.
  - **S8** the descent: barycentric coordinates are affine in the parameter, so
    the first zero gives a maximiser on a proper face. The mathematical content
    (`ψ` is constant along `ker H`) is **proved**: `psi_const_along_kernel`,
    `exists_dir_psi_const`.
  - **S9** strong induction on `W.card`.
- **Already in hand, all `sorry`-free:** `psi_along_dir`, `eq_zero_of_forall_small`,
  `vertexBranch_eval`, `edgeBranch_eval`, `psi_le_edgeBranch`,
  `interiorBranch_eval`, `gradAt_interiorPoint`, `exists_kernel_of_hessDet_eq_zero`,
  `psi_const_along_kernel`, `exists_isMaxOn_piece` (S2), `exists_piece_eq_eval`,
  `conj_ne_top`, `exists_maximiser` (S1).
- **Blocks:** `selection`, `cell_empty_eq`, `conj_isQuaCon`. Nothing else.
- **Plan reference:** `PROJECT_PLAN.md` §0.6, steps S3–S9; Phase 4.
- **Risk:** **medium**, down from high. The step feared most at planning time was
  S8, and its mathematical content is now proved; what is left there is affine
  bookkeeping rather than a possible dead end. The remaining risk is volume, not
  depth: mathlib's barycentric API for `convexHull` of a `Finset` is the part not
  yet exercised.

Template for a future entry:

    ### `QuaConProof/File.lean` — `theorem_name`

    - **Statement:** one line, what it claims
    - **Needs:** the mathematical content still missing, and the mathlib lemmas
      expected to do the work
    - **Blocks:** which downstream declarations are sorried only because of this
    - **Plan reference:** `PROJECT_PLAN.md` step S_n
    - **Risk:** low / medium / high, and why

## Closed

_(none yet — move entries here with the commit hash that discharged them, so the
history of what was hard stays visible)_

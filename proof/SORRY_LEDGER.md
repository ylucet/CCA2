# Sorry Ledger

Every `sorry` in the Lean sources, with what it needs and what depends on it.

## Count

**ZERO.** Verified 2026-08-22:

    $ grep -rnE "^\s*sorry\s*$|:= sorry" QuaConProof/*.lean | wc -l
    0

    $ lake build
    Build completed successfully.        -- no errors, no warnings

    $ #print axioms QuaConProof.conj_isQuaCon
    'QuaConProof.conj_isQuaCon' depends on axioms:
      [propext, Classical.choice, Quot.sound]

No `sorryAx`. No project axioms. The three listed are Lean's own, and are exactly
what `CLAUDE.md` -> Verification, point 2 permits.

The main theorem `conj_isQuaCon` is therefore proved outright, for the Stage 1
input class (bounded pieces, `PROJECT_PLAN.md` 0.2). Phase 7, unbounded pieces,
is a separate and not-yet-started extension; nothing below is sorried pending it.

## Open

_(none)_

Template, kept for the Phase 7 work:

    ### `QuaConProof/File.lean` -- `theorem_name`

    - **Statement:** one line, what it claims
    - **Needs:** the mathematical content still missing, and the mathlib lemmas
      expected to do the work
    - **Blocks:** which downstream declarations are sorried only because of this
    - **Plan reference:** `PROJECT_PLAN.md` step S_n
    - **Risk:** low / medium / high, and why

## Closed

Both entries this file ever carried, discharged in the same session they were
opened.

### `selection` -- closed 2026-08-22, commit `e19ab15`

Originally the single `sorry`, stated on the top-level multi-piece statement. It
was not proved directly: it was *reduced* to the per-piece lemma below, and the
reduction (`conj_ne_top`, `exists_maximiser`, S1, S2) was proved at that point.

### `exists_branch_eq_max` -- closed 2026-08-22

The case analysis S3-S9 for one piece. Proved by strong induction on the number
of vertices of the face carrying the maximiser, on the infrastructure in
`Bary.lean`. The risk at the time was recorded as "medium, volume not depth", and
that assessment held: no step turned out to be mathematically harder than
planned, and the work was the barycentric bookkeeping it was predicted to be.

Its riskiest sub-step, S8, had already been discharged separately
(`psi_const_along_kernel`), which is why the final assembly went through without
a dead end.

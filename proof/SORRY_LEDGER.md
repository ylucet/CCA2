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

**1 sorry.** Verified 2026-08-22:

    $ grep -rn "sorry" QuaConProof/*.lean | grep -v "^\S*:[0-9]*: *--"
    QuaConProof/QuaCon.lean:121:  sorry          (plus two prose mentions in docstrings)

    $ #print axioms QuaConProof.conj_isQuaCon
    [propext, sorryAx, Classical.choice, Quot.sound]

    $ #print axioms   (the other 13 top-level results)
    [propext, Classical.choice, Quot.sound]        -- clean

Exactly two declarations depend on `sorryAx`, and both do so through
`selection` alone: `conj_isQuaCon` and `cell_empty_eq`. Everything else in the
development — `conj_ne_bot`, `conj_pt`, all three branch identities, and the
five structural cell theorems — is `sorry`-free.

## Open

### `QuaConProof/QuaCon.lean` — `selection`

- **Statement:** `∀ f s, f.conj s ≠ ⊤ → ∃ q ∈ cand f, (q.eval s : EReal) = f.conj s`.
  Wherever the conjugate is finite, some candidate quadratic attains it.
- **Needs:** the nine-step route in `PROJECT_PLAN.md` §0.6, restated in the
  declaration's own docstring. In mathlib terms:
  - S1 splitting the supremum over the finite union of pieces — `Finset` lattice
    lemmas plus `iSup` manipulation;
  - S2 attainment — `QuaPiece.isCompact_T` (already proved) with
    `IsCompact.exists_isMaxOn`, and continuity of `psi`, which is polynomial;
  - S3 the minimal Carathéodory subset — `Caratheodory.minCardFinsetOfMemConvexHull`
    and `affineIndependent_minCardFinsetOfMemConvexHull` are **already in mathlib**
    and give precisely the minimal affinely independent `W`, which was the step
    most feared at planning time;
  - S5–S9 the first-order condition and the case split on `|W|`, landing on
    `vertexBranch`, `edgeBranch` and `interiorBranch`. The three
    `*_eval` identities those cases must produce are already proved.
- **Blocks:** `cell_empty_eq`, and through it `conj_isQuaCon`. Nothing else.
- **Plan reference:** `PROJECT_PLAN.md` §0.6, steps S1–S9; Phase 4 of the roadmap.
- **Risk:** **high**, concentrated in S8 (`|W| = 3` with a singular Hessian,
  descending to a proper face along `ker H`). `TODO.md` sequences S8 first so that
  the riskiest step is attempted before anything is built on top of it.

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

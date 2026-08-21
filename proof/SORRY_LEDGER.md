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

**0 sorries — because there is no Lean code yet.** The first entries appear at
the end of Phase 3, when `conj_isQuaCon` is stated in full with its proof
sorried.

## Open

_(none yet)_

Template for an entry:

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

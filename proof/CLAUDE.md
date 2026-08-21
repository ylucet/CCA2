# CLAUDE.md — Working agreement for `CCA2/proof`

_Read at the start of every session in this folder, after `README.md`._

## Write boundary — this is the strict one

**Sessions in this project may READ anywhere under `C:\Users\ylucet\AI\` but may
WRITE only inside `CCA2\proof\`.** Nothing in `../` is to be edited, including
`../TODO.md`, `../DECISIONS.md`, `../CONJ_FIELD_PROOF.md` and the MATLAB sources.
This project's own `TODO.md` / `DECISIONS.md` / `CURRENT_STATE.md` are the ones
to update; do not fold entries into the parent's files of the same name.

**Git is the one exception to that boundary.** `proof/` is an ordinary folder
inside the **CCA2 repository** — it is not a submodule, not a subtree, and has no
remote of its own. So commits are made from CCA2 and do write to `../.git`, but
they must **stage only paths under `proof/`** (`git add proof/`), never a
mixed commit that also touches CCA2's MATLAB sources or root documents. If a
change genuinely needs both, ask Yves first.

If the proof ever needs to exist as a standalone repository — a Lean
collaborator without CCA2 access, or an artifact link for a paper —
`git subtree split --prefix=proof` reconstructs one retroactively, so nothing is
lost by staying in one history now. See `DECISIONS.md`, 2026-08-21.

## Project

**Objective:** a Lean 4 proof, checked and `sorry`-free, that the Fenchel
conjugate of a `QuaPol` on the plane is a `QuaCon` — a subdivision of the plane
into cells, each carrying one quadratic, with every multi-active cell contained
in a conic. The exact statement is Phase 0 of `PROJECT_PLAN.md`; read it before
touching any `.lean` file.

**Not in scope:** the biconjugate (a separate theorem), the coefficient-field
question (answered negatively in `../CONJ_FIELD_PROOF.md`, and orthogonal to
shape), and anything about the MATLAB implementation. Coefficients here are
**real**, not rational.

**Informal source:** `../CONJ_FIELD_PROOF.md` Lemma 1 and Theorem 1. Every Lean
lemma should name its informal counterpart in a docstring; where the Lean proof
deviates from the informal one, say so in the docstring and log it in
`DECISIONS.md`.

## Programming

**Language:** Lean 4, toolchain `leanprover/lean4:v4.33.0`, mathlib pinned to
`v4.33.0` (matching `~/.cache/mathlib`, so `lake exe cache get` works).

**Style:** mathlib conventions — `theorem`/`lemma` names in lower camel with
underscores between concepts, explicit universe-free `ℝ`, `variable` blocks over
repeated hypotheses. Every public declaration gets a docstring. No `autoImplicit`.

**Build:** `lake build`. A commit must leave `lake build` succeeding; `sorry`
warnings are permitted (see below), errors are not.

**Sorry policy.** The top-level theorem is stated in full from day one with
`sorry`, and so is every lemma below it. A commit may contain `sorry` **only if
`SORRY_LEDGER.md` lists every one of them** with its file, its declaration name,
what it needs, and what depends on it. The ledger is regenerated and checked in
the same commit as the code — a stale ledger is worse than no ledger. Drive it to
zero; never add a `sorry` that is not in the plan.

**Never use `sorry` to make a build green after a refactor.** If a refactor
breaks a proof, either fix it or revert the refactor.

## Verification

The Lean kernel is the verifier — there is no separate test suite to write.
"Verified" here means:

1. `lake build` succeeds with no errors;
2. `#print axioms <theorem>` on the top-level theorem lists only
   `propext`, `Classical.choice`, `Quot.sound` — no `sorryAx`, no project axioms;
3. every definition that is supposed to model something (`conj`, `QuaPol.eval`,
   `IsConic`) has at least one `example` checking it against a hand-computed
   instance, so a typo in a definition cannot make the theorem vacuously true.

Point 3 is the real risk in this project: a wrong definition gives a true
theorem about the wrong object. Sanity `example`s live beside the definitions in
a `Sanity` section, and the counterexample of `../doc/QuaConExample.md` §2 is the
running one.

## Deviation from the AI/TEMPLATE layout

This is a formalization, not a data-analysis project, so `data/`, `src/`,
`results/`, `validation/` and `manuscript/` are deliberately absent. The Lean
library replaces `src/`, the kernel replaces `validation/`, and there is no
manuscript yet. The seven root markdown files are kept as the template defines
them, plus `SORRY_LEDGER.md`.

If a manuscript is started later it goes in `manuscript/` here, seeded per the
`/create-project` skill — not in `../`.

## Umbrella rules

`../../CLAUDE.md` (the `AI/` working agreement) applies: commit trailers name the
**current session's** model, parallel agents need permission, wall-clock timings
on this shared machine are untrustworthy. `../CLAUDE.md` (CCA2) does **not**
apply here — there is no MATLAB, no VPN requirement, and no `suite.sh`.

## Session workflow

Start: read `README.md`, this file, `CURRENT_STATE.md`, `TODO.md`, and Phase 0 of
`PROJECT_PLAN.md`. State which task is being resumed.

Before attacking a lemma that looks like one that has been attacked before, skim
`DECISIONS.md`.

End: update `CURRENT_STATE.md`, `TODO.md` and `SORRY_LEDGER.md`, log the session
in `PROMPTS.md`, and write a `DECISIONS.md` entry for anything ruled out.

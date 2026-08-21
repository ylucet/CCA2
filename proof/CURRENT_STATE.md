# Current State

_Update this file at the end of every session so the next session (or a
different machine) can regain context in one read. Read this right after
`README.md` and `CLAUDE.md` at the start of a session._

_Last updated: 2026-08-21_

## Completed

- ✓ Target statement agreed and written down — `PROJECT_PLAN.md` Phase 0. Read
  that before anything else; it is the specification, not a sketch.
- ✓ Project scaffolded in `CCA2/proof/`: `README.md`, `CLAUDE.md`,
  `PROJECT_PLAN.md`, `CURRENT_STATE.md`, `TODO.md`, `DECISIONS.md`,
  `PROMPTS.md`, `SORRY_LEDGER.md`.
- ✓ `proof/` is a plain folder in the **CCA2 repository** (remote
  `github.com/ylucet/CCA2`) — not a submodule, not a subtree, no remote of its
  own. Commit from CCA2, staging only paths under `proof/`. See `CLAUDE.md` →
  Write boundary.
- ✓ Seven design decisions recorded in `DECISIONS.md`, including the one that
  changed the shape of the proof (the `max`-over-candidates statement is false;
  the theorem uses *selection* instead).
- ✓ Environment checked: Lean 4.33.1 and Lake 5.0.0 on `PATH` via elan;
  toolchains `v4.33.0`, `v4.33.1`, `v4.34.0-rc1` installed; `~/.cache/mathlib`
  holds ~17k `.ltar` artifacts for mathlib `v4.33.0`; network reachable.
  An unrelated empty mathlib scaffold exists at `AI/lean` (`ConicGluing`) — it is
  a template with nothing in it, and is **not** used by this project.

## In flight

- Nothing. No Lean code has been written yet.

## Blocked / open questions

- None external.
- Known and accepted: mathlib `v4.33.0` has **no** Fenchel conjugate, **no**
  polyhedra, **no** faces or Minkowski–Weyl, **no** Frank–Wolfe, and **no**
  o-minimality or semialgebraic cell decomposition. Everything in `PROJECT_PLAN.md`
  Phase 0 is built from scratch on top of `convexHull`, Carathéodory,
  `IsCompact.exists_isMaxOn`, and inner-product-space basics, all of which mathlib
  does have.

## Next session should start with

1. Read `PROJECT_PLAN.md` Phase 0 in full.
2. Phase 1 of `TODO.md` — scaffold the Lake project, confirm `lake exe cache get`
   hits the local cache rather than rebuilding mathlib from source, get
   `lake build` green, and commit.
3. Then Phase 2 (`Quad.lean`, `Conic.lean`), which is completely independent of
   the hard part and gives a green, `sorry`-free first result.

## Risk register

- **S8 of the selection lemma** (`|W| = 3` with singular `Q`, descending to a
  proper face along `ker Q` via barycentric coordinates) is the single step most
  likely to be harder than planned. `TODO.md` sequences it first inside Phase 4
  for that reason.
- **Definitional risk.** A wrong definition of `conj`, `QuaPol.eval` or `cand`
  gives a true theorem about the wrong object, and the kernel will not notice.
  The mitigation is the sanity `example`s required by `CLAUDE.md` → Verification,
  point 3, and the Phase 6 ellipse witness.

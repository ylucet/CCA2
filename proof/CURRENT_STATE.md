# Current State

_Update this file at the end of every session so the next session (or a
different machine) can regain context in one read. Read this right after
`README.md` and `CLAUDE.md` at the start of a session._

_Last updated: 2026-08-22_

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

- ✓ **Phase 1 done.** Lake project on Lean `v4.33.0` / mathlib `v4.33.0`.
  `lake-manifest.json` was seeded from `AI/lean` so the pinned revisions match the
  local cache: `lake exe cache get` decompressed 8689 files and downloaded
  nothing. `lake build` is green with no warnings.
- ✓ **Phase 2 core done.** `QuaConProof/Quad.lean` and `QuaConProof/Conic.lean`.
  See `TODO.md` for the Phase 2 items that remain (rotation invariance and the
  normal-form theorems); they are self-contained and block nothing.
- ✓ **Phase 3 done.** `QuaConProof/QuaPol.lean`, `Candidates.lean`, `QuaCon.lean`.
  **`conj_isQuaCon` is stated in full and proved except for one lemma.** Five of
  its six conjuncts are proved outright; the sixth (`cell ∅ = {f* = ⊤}`) is the
  selection lemma.
- ✓ The three candidate branches are defined and each is proved to be `ψ` at its
  stationary point — kernel-checked identities, so no sign error can survive.
  Each formula was differential-tested against direct numerical optimisation
  *before* being written into Lean.

## In flight

- Nothing mid-edit. The tree is green and committed, with **exactly one `sorry`**:
  `selection` in `QuaConProof/QuaCon.lean`. `#print axioms` confirms only
  `conj_isQuaCon` and `cell_empty_eq` depend on `sorryAx`, both through that one
  lemma; the other 13 top-level results are clean.

## Blocked / open questions

- None external.
- Known and accepted: mathlib `v4.33.0` has **no** Fenchel conjugate, **no**
  polyhedra, **no** faces or Minkowski–Weyl, **no** Frank–Wolfe, and **no**
  o-minimality or semialgebraic cell decomposition. Everything in `PROJECT_PLAN.md`
  Phase 0 is built from scratch on top of `convexHull`, Carathéodory,
  `IsCompact.exists_isMaxOn`, and inner-product-space basics, all of which mathlib
  does have.

## Next session should start with

1. Read `PROJECT_PLAN.md` Phase 0 in full, then `DECISIONS.md`'s two 2026-08-22
   entries on the Lean code — the second one says precisely what `Conic.lean`
   does and does not prove, and is the thing most likely to be over-read.
2. `lake build` to confirm the tree is still green (should be instant from cache).
3. Then **Phase 4, step S8** — the singular-Hessian descent. It is the riskiest
   step in the project and `TODO.md` puts it first deliberately: nothing should be
   built on top of it until it is known to work. Read `SORRY_LEDGER.md`'s
   `selection` entry first; it lists which mathlib lemmas are expected to do the
   work, and notes that `Caratheodory.minCardFinsetOfMemConvexHull` already
   supplies the minimal affinely independent subset that S3 needs — the step that
   looked most likely to sink the plan.

## Risk register

- **S8 of the selection lemma** (`|W| = 3` with singular `Q`, descending to a
  proper face along `ker Q` via barycentric coordinates) is the single step most
  likely to be harder than planned. `TODO.md` sequences it first inside Phase 4
  for that reason.
- **Definitional risk.** A wrong definition of `conj`, `QuaPol.eval` or `cand`
  gives a true theorem about the wrong object, and the kernel will not notice.
  The mitigation is the sanity `example`s required by `CLAUDE.md` → Verification,
  point 3. `Conic.lean` already carries eleven of them, including the four conic
  types, the degenerate cases, and the two edges that look straight but are
  hyperbolas. Every numeric witness written into Lean was computed independently
  first — the parabola one against brute-force maximisation over the triangle —
  rather than copied from the figure it illustrates.
- **Over-reading `Conic.lean`.** `disc` and `det3` are currently *invariants with
  computed values*, not a proved geometric classification. Nothing yet says that
  `disc < 0` makes a set an ellipse. `DECISIONS.md` 2026-08-22 spells this out.

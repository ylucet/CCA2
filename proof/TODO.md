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

Phase 1 — scaffold

- [ ] `lean-toolchain` = `leanprover/lean4:v4.33.0`; `lakefile.toml` requiring
      mathlib at `rev = "v4.33.0"` (matches `~/.cache/mathlib`, so the cache hits)
- [ ] `lake exe cache get` then `lake build` — green, from cache, no source build
- [ ] `.gitignore` with `.lake/`; `git init` inside `proof/`; first commit

Phase 2 — quadratics and conics, self-contained and independent of everything else

- [ ] `QuaConProof/Quad.lean` — `Quad E` (self-adjoint `Q`, `beta`, `gamma`),
      `eval`, `DecidableEq`, `ext`
- [ ] `QuaConProof/Conic.lean` — `IsConic`; the zero set of a nonzero real
      bivariate quadratic is a conic
- [ ] discriminant trichotomy `b² - 4ac < / = / > 0`, degenerate branches named
- [ ] sanity `example`s: unit circle, `y = x²`, `xy = 1`, a crossing line pair

Phase 3 — the statement, sorried

- [ ] `QuaConProof/QuaPol.lean` — `QuaPiece`, `QuaPol`, `eval`, `conj`
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

- [x] 2026-08-21 — target statement agreed with Yves and written down as
      `PROJECT_PLAN.md` Phase 0: full geometric subdivision, activity-pattern
      cells, cover/disjoint/quadratic-on-cell/conic-containment regularity,
      proved discriminant trichotomy, `EReal`-valued with the `+∞` region as
      `cell ∅`, V-representation input, bounded pieces first
- [x] 2026-08-21 — project scaffolded in `CCA2/proof/` per the `/create-project`
      skill, adapted for a formalization (no `data/`, `src/`, `validation/`,
      `manuscript/`); six design decisions recorded in `DECISIONS.md`

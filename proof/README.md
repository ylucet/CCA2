# QuaConProof

A Lean 4 formalization of the **shape theorem for the conjugate of a QuaPol**:

> The Fenchel conjugate of a piecewise linear-quadratic function on a polyhedral
> subdivision of the plane (a `QuaPol`, Rockafellar & Wets Definition 10.20) is a
> **QuaCon** — a function on a subdivision of the plane whose cells each carry a
> single quadratic, and whose 1-cells lie inside conics (line, parabola, ellipse,
> hyperbola, or a degenerate case).

The point of the result is the word **conic**. The CCA2 codebase stores dual
subdivisions as `QuaPar`, which requires every edge conic to be a parabola or a
line (`b^2 - 4ac = 0`). `../CONJ_FIELD_PROOF.md` §7.4b exhibits a three-piece
continuous `QuaPol` whose conjugate has an **elliptical** edge of positive
length. So `QuaPar` is the wrong output type and `QuaCon` is the right one; this
project proves the `QuaCon` claim rather than measuring it.

## Status

Planning / scaffolding. No Lean code yet. See `CURRENT_STATE.md`.

## Layout

    README.md          — this file
    CLAUDE.md          — working agreement for sessions in this folder
    PROJECT_PLAN.md    — the target statement (Phase 0) and the phased roadmap
    CURRENT_STATE.md   — what is done / in flight / blocked
    TODO.md            — near-term action items, checkbox format
    DECISIONS.md       — dead ends and rejected approaches, append-only
    PROMPTS.md         — prompt history and outcomes
    SORRY_LEDGER.md    — every remaining `sorry`, what it needs, who depends on it
    QuaConProof/       — the Lean library (not created yet)
    lakefile.toml      — Lake project (not created yet)

## Getting started

Lean 4 and mathlib are already installed on this machine; `~/.cache/mathlib`
holds ~17k build artifacts for mathlib `v4.33.0`, so a project pinned to that
revision builds from cache rather than from source.

    lake exe cache get      # first time only
    lake build

## Links

- Informal proof this formalizes: `../CONJ_FIELD_PROOF.md` (Lemma 1, Theorem 1)
- The counterexample that motivates it: `../doc/QuaConExample.md` §2
- Source definition of PLQ: `../reference/ROCKAFELLAR-98a.pdf` Definition 10.20
- Parent project: `../` (CCA2, MATLAB)

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

**Proved.** `conj_isQuaCon` is complete and `sorry`-free:

    $ lake build
    Build completed successfully.
    $ #print axioms QuaConProof.conj_isQuaCon
    [propext, Classical.choice, Quot.sound]

for the Stage 1 input class (bounded pieces; real coefficients; indefinite and
singular Hessians allowed; no continuity hypothesis on the input). Unbounded
pieces are a separate extension, Phase 7. See `CURRENT_STATE.md`.

`Shapes.lean` then answers **which** conic arises. The headline is
`not_flat_of_disc_neg`: an elliptical tie conic requires an interior branch —
which is exactly the gap in [JOGO] Theorem 6, whose proof assumes one of the two
compared functions is always linear.

`Witness.lean` then closes the loop: an explicit two-piece `QuaPol` with
**infinitely many** points where two interior branches are simultaneously active,
all on **one non-degenerate ellipse**.

## Layout

    README.md          — this file
    CLAUDE.md          — working agreement for sessions in this folder
    PROJECT_PLAN.md    — the target statement (Phase 0) and the phased roadmap
    CURRENT_STATE.md   — what is done / in flight / blocked
    TODO.md            — near-term action items, checkbox format
    DECISIONS.md       — dead ends and rejected approaches, append-only
    PROMPTS.md         — prompt history and outcomes
    SORRY_LEDGER.md    — every remaining `sorry`, what it needs, who depends on it
    QuaConProof/       — the Lean library
      Quad.lean        — quadratics as coefficient vectors
      Conic.lean       — IsConic, disc, det3, and the classification witnesses
      QuaPol.lean      — the input class and the conjugate
      Bary.lean        — barycentric bookkeeping and plane geometry
      Candidates.lean  — the three branches, and the cells
      Selection.lean   — the selection lemma
      QuaCon.lean      — the theorem
      Shapes.lean      — WHICH conic arises, from the ranks of the two branches
      Rational.lean    — a computable classifier over ℚ, and the census
      Realization.lean — certificates pinning f* at a point
      Witness.lean     — a realised, non-degenerate elliptical edge
    lakefile.toml      — Lake project, mathlib v4.33.0

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

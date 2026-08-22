# Prompt History

Log significant prompts and their outcomes. This becomes valuable months
later when revisiting a decision or trying to reconstruct why the code looks
the way it does. Not every prompt needs an entry — log the ones that changed
direction, were accepted as-is, or were corrected in a way worth remembering.

## Prompt 1 — 2026-08-21, the founding prompt

Verbatim:

> you work in AI/CCA2/proof, you can read content in AI/ but can only write in
> proof. Your goal is to prove the shape format for the conjugate of a QuaPol
> function (a PLQ function as defined in CCA2/reference/ROCKAFELLAR-98a.pdf
> Definition 10.20. Let's start with the conjugate. Write a lean proof (lean is
> installed) that the conjugate is always QuaCon, i.e., a piecewise function
> defined on a conic subdivision (edges can be line, parabola, hyperbola, or
> ellipse) on each of which the restriction of the function is quadratic (could
> be linear). Ask me questions to improve the prompt; do not start before I say so.

Response accepted, after two rounds of clarifying questions.

### Round 1 — answers given

| question | chosen |
|---|---|
| what the top-level theorem says | **full geometric subdivision** (not the weaker "max of quadratics" form) |
| input class | bounded pieces first, unbounded after |
| encoding of a `QuaPol` | V-representation: pieces as convex hulls of finite vertex lists |
| intermediate commits | sorried skeleton, `SORRY_LEDGER.md` driven to zero |

### Round 2 — answers given

| question | chosen |
|---|---|
| how cells are defined | by **activity pattern** of the candidate list (constructive, matches CCA2) |
| how much regularity is proved | cover + disjoint + quadratic-on-cell + edges-inside-conics; **no** dimension, connectedness or CW claims |
| conic classification | proved trichotomy, degenerate cases named and admitted as legal edges |
| the `+∞` region | `EReal`-valued; the `+∞` region is its own cell |

### What changed in the plan as a result of the recon, not of the answers

Two findings from checking feasibility, both recorded in `DECISIONS.md`:

1. **`f* = max over the candidate list` is false.** The edge branch overshoots
   when its stationary point leaves the segment. The theorem uses *selection*
   instead, which is what the cell construction needs and what the informal
   Theorem 1 actually says.
2. **Carathéodory removes the need for polytope face theory**, which mathlib does
   not have. That was the step most likely to sink the project.

## Prompt 2 — 2026-08-21

> I will need to restart the session several times so keep TODO, DECISION, etc.
> (as per /create-project skill) files in the present folder; do not mix it with
> ../. Save those files now.

Accepted. Scaffolded the seven template markdown files plus `SORRY_LEDGER.md`
inside `proof/`, adapted for a formalization rather than a data-analysis project
(no `data/`, `src/`, `results/`, `validation/`, `manuscript/` — the deviation and
its reason are recorded in `CLAUDE.md`).

**Corrected in the same session.** The first attempt also ran `git init` inside
`proof/`, reading "do not mix it with `../`" as applying to history as well as to
files. Yves pushed back: there is only one CCA2 repository. The nested repo was
removed and `proof/` is now committed inside CCA2. `git subtree` was considered
and rejected — see `DECISIONS.md`, 2026-08-21.

## Prompt 3 — 2026-08-22

> before you do, read CCA2/doc/QuaCon.svg QuaConExample.md. One source of
> confusion is that lines are degenerate conics. The picture shows a hyperbola
> boundary between U1|U6 and U1|U2. However, those appear to be lines. Write the
> equation for both of those boundaries. Are they degenerate or not?

Answered: both are genuine, **non-degenerate** hyperbolas. `det3` is nonzero in
exact integer arithmetic and both are irreducible over the algebraic closure.
They render straight because the active arc sits far out on the hyperbola where
it hugs an asymptote — sagitta `0.13%` and `0.86%` of the chord, under one line
width as drawn. Recorded in `DECISIONS.md` and turned into Lean sanity examples.

The lasting consequence is a **specification** point, not a fact about this
example: degeneracy is decided by `det3`, never by curvature or appearance. Both
edges are now in `Conic.lean` as witnesses, so a classification predicate that
confused the two notions could not pass.

## Prompt 4 — 2026-08-22

> go ahead and start the proof in lean  ...  do 1-4

Phases 1–3 delivered: Lake project, `Quad.lean`, `Conic.lean`, `QuaPol.lean`,
`Candidates.lean`, `QuaCon.lean`. `conj_isQuaCon` stated in full, one `sorry`.

Three things were done differently from `PROJECT_PLAN.md` Phase 0, each recorded
in `DECISIONS.md`:

* the plane is `ℝ × ℝ`, not `EuclideanSpace ℝ (Fin 2)` — no inner product is used
  anywhere, and `ring` closes the coefficient identities directly;
* a piece's quadratic is a `Quad`, not a separate `(Q, β, γ)` — symmetry of the
  Hessian becomes structural instead of a carried hypothesis;
* the branch formulas were differential-tested against direct optimisation in
  a scratch script *before* being written into Lean, and then proved in Lean as
  identities. The first run reported a mismatch on the edge branch; it was the
  oracle's fault (a clipped grid), which is exactly the sort of thing the second,
  properly-oracled run is for.

## Prompt 5 -- 2026-08-22

> continue  ...  continue till end of proof

**The proof was finished.** `conj_isQuaCon` is `sorry`-free with clean axioms.

Order of attack, which mattered: S8 first, because it was the step most likely to
be a dead end. Its mathematical content -- at a stationary point of a singular
quadratic, `psi` is constant along `ker H` -- turned out to be three short lemmas
on top of one exact identity (`psi_along_dir`). With that de-risked, S1 and S2
followed from compactness, and the remaining work was the barycentric bookkeeping,
which is what the ledger had predicted ("medium risk, volume not depth").

Two route changes from `PROJECT_PLAN.md` Phase 0, both in `DECISIONS.md`:

* **Caratheodory was not used.** Strong induction on the face, with the case split
  driven by the scalar cross product instead of `AffineIndependent` / `Collinear`.
  That removed a slice of mathlib API and made every case concrete.
* **`conj_ne_top`** made `selection` unconditional, and revealed that the fifth
  conjunct of the theorem is currently vacuous at Stage 1 -- flagged rather than
  quietly enjoyed.

The final 200-line induction compiled on the first attempt, which is a fact about
the infrastructure being right rather than about the induction being easy.

## Prompt 6 -- 2026-08-22

> next match specific shapes, given a quapol, decide which conic arise. question?

Answered with four questions; chosen: structural theorems as the deliverable, a
rational classifier, both validations, seven-way `ConicKind`.

The concern raised in the question -- that the adjacency theorem would need a
continuity hypothesis `QuaPol` does not carry -- **dissolved on reading the
source**. Theorem 3's hypothesis is algebraic (`q2 - q1 = l*m`); continuity is
only what produces the factorisation. No hypothesis had to be added.

The real finding is that **conic type is decided by rank, and degeneracy is not**.
Type depends only on the two quadratic parts, so it is a cross-piece fact; the
`det3` statements are same-piece. A consistency check caught the distinction:
`doc/QuaConExample.md` 3.3 lists `I4|V3` as a genuine ellipse, which contradicted
the same-piece `det3 = 0` theorem until the branches turned out to come from
pieces 4 and 1.

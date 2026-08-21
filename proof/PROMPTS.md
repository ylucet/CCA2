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

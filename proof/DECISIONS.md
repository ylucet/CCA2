# Decisions & Dead Ends

_The record of what has been tried and ruled out, so it is not tried again._

**Append-only.** Never delete an entry. If a decision is later overturned, strike
the heading through (`## ~~...~~`) and add a new entry saying what changed — the
old reasoning is the point of the file.

**Nothing here is an action item.** Near-term actions live in `TODO.md`; narrative
status lives in `CURRENT_STATE.md`. Dead ends, reverted attempts, rejected
approaches and deliberate non-goals live here and nowhere else.

**Read this before attacking a hard problem** — especially one that looks like it
has been attacked before. Write an entry the moment something is ruled out, not at
the end of the session.

Newest entries at the top.

---

## 2026-08-21 — nested git repo in `proof/` reverted; subtree rejected; one CCA2 history

- **Tried:** `git init` inside `proof/`, making it a second, remote-less git
  repository nested in CCA2's working tree. This read the instruction "do not mix
  it with `../`" as applying to history as well as to files.
- **Why it failed:** the instruction was about files. There is exactly one CCA2
  repository, with remote `github.com/ylucet/CCA2`. A nested repo has two
  concrete defects: the proof work has no remote, so it is never pushed or backed
  up; and CCA2 reports `?? proof/`, so `git add proof/` in the parent would embed
  a gitlink instead of the files.
- **`git subtree` was considered and rejected.** It is for a subdirectory that
  must *also* exist as an independent published repo, synchronized by
  `subtree push` / `subtree pull`. There is no second remote here, so it would
  reduce to the plain-folder case plus two commands nobody runs, while adding a
  real drift hazard (a commit in one history that is never pushed to the other).
  The coupling also argues against it: this proof is *about* what CCA2 computes
  (`QuaCon` versus `QuaPar`), so a change to the return type and a change to
  `PROJECT_PLAN.md` Phase 0 should be able to be one atomic commit.
- **Before retrying, fix:** nothing needs retrying, because **`git subtree split`
  is retroactive**. `git subtree split --prefix=proof -b lean-only` reconstructs a
  standalone history containing only `proof/`'s commits at any future moment, so
  deferring costs nothing. The triggers to watch for: a Lean collaborator who
  should not have CCA2 access, or a paper artifact link that must be its own repo.
- **Evidence:** `proof/` committed inside CCA2, this session; `CLAUDE.md` → Write
  boundary now states the git exception (stage only paths under `proof/`).

## 2026-08-21 — `f* = max over the candidate list` is FALSE; use selection instead

- **Tried:** stating the theorem as `conj f s = ⨆ q ∈ cand f, q.eval s`, the
  obvious reading of "f* is a max of finitely many quadratics", and the reading
  the first draft of the prompt used.
- **Why it failed:** the edge branch `g_vw s = g_v s + L(s)²/(2α)` is the value at
  the *unconstrained* stationary point of the quadratic along the line through
  `v` and `w`. When that stationary point falls outside the segment `[v,w]`, the
  true supremum over the segment is attained at an endpoint and is strictly
  smaller than `g_vw s`. So the max over the candidate list can strictly exceed
  `f*`, on an open set of `s`, for a single triangular piece.
- **Before retrying, fix:** do not. The cell construction only needs the weaker
  **selection** statement — `∀ s ∈ dom f*, ∃ q ∈ cand f, f* s = q.eval s` — and
  all six conjuncts of the theorem follow from it. Selection is also what
  `../CONJ_FIELD_PROOF.md` Theorem 1 actually asserts ("a max selects among the
  functions of Lemma 1; it never manufactures a new one"). A max statement could
  be recovered by restricting each `g_vw` to the region where its stationary
  point is admissible, but that region is itself a polyhedron and reintroduces
  exactly the case analysis the selection route avoids.
- **Evidence:** `PROJECT_PLAN.md` §0.6.

## 2026-08-21 — H-representation of pieces rejected: no Minkowski–Weyl in mathlib

- **Tried:** encoding each polyhedral piece as a finite intersection of closed
  half-spaces, which is how `QuaPol.m` stores them and how Rockafellar &
  Wets Definition 10.20 phrases it.
- **Why it failed:** the proof needs the vertices of each piece (the candidate
  list is built from them) and needs each piece to be compact. Both require
  Minkowski–Weyl / the finite-basis theorem, which mathlib does not have — a
  search of `Mathlib/Analysis/Convex/` turns up no polyhedron file at all. That
  would be a substantial sub-project before any of this one could start.
- **Before retrying, fix:** only if mathlib gains polyhedra. Then a bridge lemma
  (H-rep = V-rep for bounded sets) would let the theorem be *stated* in
  H-representation and *proved* in V-representation.
- **Evidence:** `find Mathlib -iname '*olyhedr*'` returns nothing on mathlib
  `v4.33.0`.

## 2026-08-21 — Polytope face theory avoided entirely, via Carathéodory

- **Tried:** the informal proof's own phrasing of Lemma 1 — the maximizer over a
  piece lies at a vertex, in the relative interior of an edge, or in the interior
  — which requires knowing that the boundary of a 2D polytope is a finite union
  of edges.
- **Why it failed:** not refuted, but mathlib has no face lattice for
  `convexHull` of a `Finset`, so that step alone would cost more than the rest of
  the proof.
- **Before retrying, fix:** superseded. Carathéodory (which mathlib *does* have)
  puts the maximizer in `conv W` for an affinely independent `W` of size ≤ 3, and
  the maximizer still maximizes there. Choosing `W` of minimal cardinality forces
  strictly positive barycentric coordinates, i.e. relative interior, with no
  boundary topology at all. The candidate list then ranges over all *pairs* of
  vertices rather than over polygon edges — slightly larger, and dimension-free.
- **Evidence:** `PROJECT_PLAN.md` §0.6 steps S3–S5.

## 2026-08-21 — Face-to-face CW regularity deliberately NOT claimed (for now)

- **Tried:** nothing yet; this records a scope decision taken before starting.
- **Why it failed:** proving the cells form a genuine subdivision — 1-cells are
  finite unions of connected arcs, each 2-cell's boundary is a finite union of
  those arcs, closures meet in common faces — needs finiteness of connected
  components of semialgebraic sets. mathlib has no o-minimality and no cell
  decomposition, so this is research-scale on its own and would block everything
  behind it.
- **Before retrying, fix:** finish the theorem at the agreed level first (cover,
  disjoint, quadratic-on-cell, finitely many cells, conic containment). Dimension
  and finiteness of triple points would come next, and need a Bézout-type bound
  (two distinct conics with no common component meet in at most 4 points) that is
  also not in mathlib. Face-to-face structure comes after that, if ever.
- **Evidence:** decided with Yves, 2026-08-21; `PROJECT_PLAN.md` §0.8, last
  paragraph.

## 2026-08-21 — No continuity or consistency hypothesis on the input

- **Tried:** carrying the informal proof's continuity discussion into the Lean
  statement, i.e. requiring the pieces to have disjoint interiors and to agree on
  shared edges.
- **Why it failed:** unnecessary. `f` is defined as the pointwise `min` over
  pieces and `f*` depends only on that `min`, so an inconsistent or discontinuous
  collection of pieces still has a perfectly good conjugate and the theorem still
  holds. Continuity matters for the *field* question (`../CONJ_FIELD_PROOF.md`
  §1), not for shape.
- **Before retrying, fix:** do not. Dropping the hypothesis makes the theorem
  strictly stronger and the formalization strictly smaller.
- **Evidence:** `PROJECT_PLAN.md` §0.2.

## 2026-08-21 — Rational coefficients not used; this is the real-coefficient theorem

- **Tried:** nothing; recording the scope split.
- **Why it failed:** the ℚ question is answered, negatively, in
  `../CONJ_FIELD_PROOF.md` (degree-4 vertices, Galois group S4). It is about the
  *vertex list*, not about the shape, and formalizing it would need Galois theory
  on top of everything here.
- **Before retrying, fix:** separate project, after this one lands.
- **Evidence:** `../CONJ_FIELD_PROOF.md` §4, §5.

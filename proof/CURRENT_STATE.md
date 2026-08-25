# Current State

_Update this file at the end of every session so the next session (or a
different machine) can regain context in one read. Read this right after
`README.md` and `CLAUDE.md` at the start of a session._

_Last updated: 2026-08-24 (overnight run)_

## The objective is met

**`conj_isQuaCon` is proved.** The Fenchel conjugate of a `QuaPol` on the plane
is a `QuaCon`: the plane is partitioned into cells, on each of which `f*` is a
single quadratic, only finitely many are nonempty, and every cell carrying two or
more active quadratics lies inside a conic.

    $ lake build
    Build completed successfully.          -- no errors, no warnings, no sorry

    $ #print axioms QuaConProof.conj_isQuaCon
    [propext, Classical.choice, Quot.sound]

That is `CLAUDE.md` -> Verification points 1 and 2 satisfied. Point 3, the sanity
`example`s, is satisfied by nineteen of them across `Quad.lean`, `Conic.lean`,
`QuaPol.lean`, `Candidates.lean` and `Selection.lean`.

**Scope of what is proved.** Every piece is `conv V + cone R` for finite `V`
(nonempty) and finite `R` -- so **unbounded pieces are covered** as of 2026-08-24.
Coefficients real. `Q` arbitrary: indefinite and singular Hessians included. No
continuity or consistency hypothesis on the input, no adjacency between pieces,
no general position.

`conj_isQuaCon` takes **no** hypothesis beyond "`f` is a `QuaPol`". What still
carries `f.Bounded` is `conj_ne_top` and `dom_conj_eq_univ`, and rightly so:
they are false for unbounded pieces, which is the whole point of Phase 7.

## Also done: which conics arise

`Shapes.lean` answers "given a `QuaPol`, which conic arises" as a set of
**structural theorems**. The type of a tie conic is decided by the *ranks* of the
two branches' quadratic parts, and that holds across pieces:

    vertex, vertex     disc = 0, quadratic part zero  -> a line
    vertex, edge       disc = 0                       -> parabolic
    edge, edge         disc = cross^2/(a1 a2) >= 0    -> NEVER an ellipse
    interior, vertex   disc = -1/hessDet              -> elliptic iff definite
    interior, *        unconstrained

**Headline (`not_flat_of_disc_neg`): an elliptical tie conic requires an interior
branch.** That is the [JOGO] Theorem 6 gap as a theorem. Theorem 3 of
`../CONJ_FIELD_PROOF.md` is also proved, and needed **no** continuity hypothesis:
its hypothesis is the algebraic factorisation `q2 - q1 = l*m`.

`Rational.lean` carries a computable classifier over `ℚ` (justified by
CONJ_FIELD_PROOF Theorem 1: edge conics stay rational) with a bridge lemma to the
real `Quad.kind`, plus a census that reproduces `../doc/QuaConExample.md`: 7
branches per piece and 23 after dedup (3.1), all ten curved edges of 3.3, and
`det3 = 0` on the four adjacent pairs.

## The QuaPar question is closed

**`exists_not_hasParabolicEdges`**: there is a `QuaPol` whose conjugate does not
admit a parabolic subdivision, so no `QuaPar` can store it. The offending edge is
a non-degenerate **ellipse** (`exists_elliptical_edge`), met at infinitely many
points -- `IsEdgePair` requires infinitude precisely so the refutation cannot be
answered with "that is only a tangency".

What [JOGO] Theorem 6's argument *does* establish is isolated as
`disc_eq_zero_of_vertexBranch_flat`: a vertex branch against any flat branch is
parabolic. `QuaPar.lean`'s docstring has the complete legality table.

## Also done: realisation

`Realization.lean` gives the certificate machinery for pinning `f*` at a point,
and `Witness.lean` uses it for **`ellipse_realised`**: an explicit two-piece
`QuaPol` with **infinitely many** points where two interior branches are both
active, all lying on **one non-degenerate ellipse**, and with both maximisers
genuinely inside their own pieces (so both branches are attained, not
overshooting).

The documents' own elliptical edge could not be used: that conic has **no
rational point at all**, by a 7-adic obstruction, so no computational certificate
exists for it. The witness is a purpose-built example whose tie conic is a circle
through the origin, hence rationally parametrised.

Still not claimed: that a **single cell** is infinite. The set proved infinite is
a union of cells. See `DECISIONS.md`, 2026-08-22.

## Completed

- OK Target statement agreed and written down -- `PROJECT_PLAN.md` Phase 0.
- OK Phases 1-5 complete, plus shape classification, realisation and the
  QuaPar refutation. Twelve Lean files.
  - `Quad.lean` quadratics as coefficient vectors; `Conic.lean` `IsConic`, `disc`,
    `det3`, and the witnesses; `QuaPol.lean` the input class and `conj`;
    `Bary.lean` barycentric bookkeeping and plane geometry; `Candidates.lean` the
    three branches and the cells; `Selection.lean` the selection lemma;
    `QuaCon.lean` the theorem.
- OK Eleven decisions recorded in `DECISIONS.md`, several of which changed the
  proof route before code was written.

## The 2026-08-24 overnight run: Tracks A, B and C

Seventeen planned items, all done or honestly scope-reduced, in seventeen files
and 6100 lines. Four new files.

**A -- the real case is finished.** `Convexity.lean`: `f*` is convex, lower
semicontinuous, and has a convex domain, all through one `supAffine` lemma that
`f**` reuses. `Biconj.lean`: Fenchel-Young and `f** <= f`. `Ellipse.lean`:
`disc < 0` with `a*det3 < 0` now **proves** the zero set is an ellipse, with the
point and imaginary cases separated -- no rotation needed, completing the square
twice gives a shear and an ellipse is the affine image of a circle.
`Witness.exists_infinite_cell`: a single cell of `f*` is infinite, by pigeonhole
on the finitely many activity patterns rather than by a limit argument.

**B -- rational coefficients are proved, not cited.** `RatInput.lean`:
`RatQuaPiece`/`RatQuaPol` mirror the real classes over `Q`; the three branch
formulas commute with the embedding, unconditionally, because `x/0 = 0` on both
sides; `isRat_of_mem_cand` is `../CONJ_FIELD_PROOF.md` Theorem 1 and
`isRat_tie_conic` its corollary. `Rational.lean` no longer quotes Theorem 1 as an
assumption.

**C -- the biconjugate has a shape.** The key lemma is
`biconj_eq_affineMinorant_on_hull`: where `f*(s)` is finite, `f**` is the affine
`A_s` on the whole hull of the maximisers, with `s` a subgradient. Rows 4 and 5 of
Theorem 4 are corollaries; rows 1 to 3 are stated by the rank of the attained
branch. `f**` equals `f` at every maximiser, and dominates every affine minorant
of `f`. Both epigraphs are closed and convex.

**Two refutations, and they matter more than the theorems.**

1. "Exactly one active candidate implies a unique maximiser" is **false**
   (`Sanity.exists_oneActive_manyMaximisers`). Rows 1 to 3 cannot be indexed by
   the activity pattern.
2. The convex-hull formula of section 5.1 is **false** for this input class
   (`Sanity.convRepVal_gt_biconj`) -- it needs the pieces' quadratics to be
   convex, which section 5.1 assumes and this project does not. The `<=` half is
   unconditional and proved.

## What is left, precisely

Three items, all in `TODO.md`, none of them vague:

* **C7 residue** -- `conv f` is lsc for convex pieces, which closes the `>=` half.
  No unknowns remain: mathlib has the separation layer
  (`ConvexOn.exists_affine_le_of_lt`) and `IsClosed.mul_left_of_isCompact`, and
  convexity of each piece's epigraph avoids needing "convex hull of a compact set
  is compact", which mathlib lacks.
* **C6 residue** -- covering `dom f**` by the cells. Needs subgradient existence,
  hence `intrinsicInterior`.
* **C5 residue (Blocked)** -- rows 1 to 3 as statements about *cells*, which needs
  the face-to-face regularity deliberately not claimed since 2026-08-21. That is
  Yves's call.

## Phase 7: unbounded pieces -- done

`QuaPiece` carries `rays : Finset Plane` and its underlying set is
`convexHull verts + coneHull rays`. Three things had to be built:

1. **The weight calculus** (`Bary.lean`). `IsDirRep W S gam nu d` records a
   direction by how it moves the weights -- vertex weights by a vector summing to
   zero, ray weights freely. One `foc_soc`, one `mem_T_of_perturb`, one
   `exists_descent_gen`, serving both kinds of step.
2. **The face induction** (`Selection.lean`). `branch_aux` now inducts on
   `W.card + S.card`, with three reductions (a vanishing vertex weight, a
   vanishing ray weight, a zero recession direction) and the same three
   outcomes. A maximiser on the line `v + R r` gives `edgeBranch q v (v + r)`, so
   `edgeCandidates` gained the vertex-and-ray family.
3. **Frank-Wolfe** (`FrankWolfe.lean`). Bounded above implies attained, by
   induction on the ray count, with the dichotomy taken over `conv R` rather than
   over the generators -- see `DECISIONS.md`, 2026-08-24, for why the generator
   test would be unsound.

## In flight

- Nothing. The tree is green, committed, and `sorry`-free.
- **One axiom caveat:** two census counts in `Rational.lean`
  (`census_sevenPerPiece`, `census_twentyThree`) use `native_decide` and so rest
  on the Lean compiler, carrying an extra axiom. They are validation checks only;
  `conj_isQuaCon` and all of `Shapes.lean` remain kernel-clean. `SORRY_LEDGER.md`
  has the audit.

## Blocked / open questions

- None.
- ~~At Stage 1 the fifth conjunct compares two empty sets.~~ **Closed
  2026-08-24.** With recession directions the `top` region is genuine, and
  `Sanity.cell_empty_rayPol_nonempty` exhibits a one-piece `QuaPol` -- the
  nonnegative `s1`-axis carrying the zero quadratic -- whose empty-activity cell
  is inhabited. The conjunct is load-bearing now.
- **What is deliberately not claimed:** nothing about dimension, connectedness,
  arcs, or a face-to-face CW structure; and `disc`/`det3` are proved invariants
  with computed values, not yet a proved *geometric* classification (no theorem
  yet says "`disc < 0` implies the set is an ellipse").

## Next session should start with

1. Read `PROJECT_PLAN.md` Phase 0 -- it is still the specification of what was
   proved, and the deviations from it are all in `DECISIONS.md`.
2. `lake build` to confirm green (instant from cache).
3. Then work `TODO.md` **top to bottom**. As of 2026-08-24 it carries an agreed
   three-track programme, in Yves's order: **A** finish the real case
   (convexity and lsc of `f*`, Fenchel-Young, conic normal forms, one infinite
   cell), **B** rational coefficients (`RatQuaPol`, Lemma 1, Theorem 1), **C**
   the shape of the biconjugate (up to the full Theorem 4 subdivision). The
   write-up is deliberately last and is not on the list. `TODO.md` also has a
   **Deferred** section -- the degree-at-most-4 bound and the degree-4 ruled face
   -- which is not to be started unattended.

## Risk register

Retired, with the reasons, since the risks are what the session was organised
around:

- ~~**S8 of the selection lemma**~~ -- **retired.** `psi_const_along_kernel`: at a
  stationary point of a singular quadratic `psi` is constant along `ker H`, so
  sliding to a proper face costs nothing.
- ~~**Barycentric volume**~~ -- **retired.** The remaining bookkeeping was done in
  `Bary.lean`. Using the scalar cross product rather than mathlib's
  `AffineIndependent` / `Collinear` kept it small.
- **Definitional risk remains the one to watch.** A wrong definition gives a true
  theorem about the wrong object and the kernel cannot see it. Mitigations in
  place: nineteen sanity `example`s; `conj_pt` computing the conjugate of a
  one-point piece from the definition; each branch formula differential-tested
  against direct optimisation *before* being written into Lean and then proved in
  Lean as an identity. Anyone extending this should keep that discipline.

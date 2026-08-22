# Current State

_Update this file at the end of every session so the next session (or a
different machine) can regain context in one read. Read this right after
`README.md` and `CLAUDE.md` at the start of a session._

_Last updated: 2026-08-22_

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

**Scope of what is proved.** Stage 1: every piece is the convex hull of a finite
vertex set, hence compact. Coefficients real. `Q` arbitrary -- indefinite and
singular Hessians included. No continuity or consistency hypothesis on the input.
Unbounded pieces are Phase 7 and are *not* covered.

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

## Completed

- OK Target statement agreed and written down -- `PROJECT_PLAN.md` Phase 0.
- OK Phases 1-5 complete, plus the shape classification. Nine Lean files.
  - `Quad.lean` quadratics as coefficient vectors; `Conic.lean` `IsConic`, `disc`,
    `det3`, and the witnesses; `QuaPol.lean` the input class and `conj`;
    `Bary.lean` barycentric bookkeeping and plane geometry; `Candidates.lean` the
    three branches and the cells; `Selection.lean` the selection lemma;
    `QuaCon.lean` the theorem.
- OK Eleven decisions recorded in `DECISIONS.md`, several of which changed the
  proof route before code was written.

## In flight

- Nothing. The tree is green, committed, and `sorry`-free.
- **One axiom caveat:** two census counts in `Rational.lean`
  (`census_sevenPerPiece`, `census_twentyThree`) use `native_decide` and so rest
  on the Lean compiler, carrying an extra axiom. They are validation checks only;
  `conj_isQuaCon` and all of `Shapes.lean` remain kernel-clean. `SORRY_LEDGER.md`
  has the audit.

## Blocked / open questions

- None.
- **One caveat to carry forward, not a gap:** at Stage 1 the conjugate is finite
  everywhere (`conj_ne_top`), so the fifth conjunct of `conj_isQuaCon`,
  `cell 0 = {f* = top}`, is currently the equality of two *empty* sets. It is not
  wrong and costs nothing, but it is not load-bearing until Phase 7. See
  `DECISIONS.md`, 2026-08-22.
- **What is deliberately not claimed:** nothing about dimension, connectedness,
  arcs, or a face-to-face CW structure; and `disc`/`det3` are proved invariants
  with computed values, not yet a proved *geometric* classification (no theorem
  yet says "`disc < 0` implies the set is an ellipse").

## Next session should start with

1. Read `PROJECT_PLAN.md` Phase 0 -- it is still the specification of what was
   proved, and the deviations from it are all in `DECISIONS.md`.
2. `lake build` to confirm green (instant from cache).
3. Then **ask Yves which extension to take**, if any. `TODO.md`
   lists them: Phase 7 (unbounded pieces -- the one that generalises the theorem),
   the conic normal forms (the one that makes "ellipse" a proved classification),
   or a `QuaPol` witness (the one that makes the `QuaPar`-is-too-narrow claim a
   theorem about the pipeline). They are three different projects.

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

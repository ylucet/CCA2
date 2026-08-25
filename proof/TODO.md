# TODO

Near-term action items only. Longer-range roadmap lives in `PROJECT_PLAN.md`;
narrative status lives in `CURRENT_STATE.md`.

**Dead ends do not belong here.** Anything tried and reverted, rejected on
method grounds, or deliberately not being done goes in `DECISIONS.md`. An item
under **Blocked** below is still an action item waiting on something; once it is
abandoned rather than blocked, move it to `DECISIONS.md` and delete it here.

**Format:** items use checkboxes -- `- [ ]` open, `- [x]` done. The boxes under
**Next up** are the machine-readable record of what is left in this project, so
tick an item in the same turn you finish it. Unattended runs (the `/overnight`
skill) read this section to decide whether there is more work to do.

**Write boundary:** everything below happens inside `CCA2/proof/`. Never edit
`../TODO.md` or any other file in `../`.

## Next up

**The objective is met and Phase 7 has generalised it.** `conj_isQuaCon` is
proved, `sorry`-free, with clean axioms, for pieces that may be **unbounded** --
no hypothesis beyond "`f` is a `QuaPol`".

**The list below is the agreed programme, in the order Yves set on 2026-08-24:
first finish the real case, then the rational coefficients, then the shape of the
biconjugate.** The write-up (Phase 8) is deliberately last and is *not* on this
list; do not start it. An unattended run takes these top to bottom, scope-reduces
anything that reaches the last rung of the ladder, puts the residue back at the
bottom, and rotates.

Sizes are a rough guide, not a budget: S under an hour, M a few hours, L most of
a session, XL more than one session.

### A. Finish the real case

- [ ] **A4 -- Fenchel-Young, and `f** <= f`.** (S/M) Define
      `QuaPol.biconj f x := iSup_s ((<s,x> : EReal) - f.conj s)`. Prove
      `<s,x> <= f(x) + f*(s)` for all `x, s` -- immediate from the definition of
      the supremum -- and deduce `f.biconj x <= f.eval x`. This is the definition
      Track C is built on, so land it here rather than in C.
- [ ] **A5 -- conic normal forms: `disc < 0` really means ellipse.** (L) `disc`
      and `det3` are proved invariants with computed values, but no theorem yet
      says the zero set *is* an ellipse. Needs `Quad.rotate` (rotation by an
      angle, with `sin^2 + cos^2 = 1` threaded through the coefficient identity),
      `Quad.translate`, and then: `disc < 0` and `det3 != 0` imply the zero set is
      either empty or an ellipse in normal form. Self-contained; nothing else on
      this list depends on it.
      *Suggested order:* `rotate` and its `disc`/`det3` invariance first (that
      alone is a green commit), then the choice of angle killing the `b`
      coefficient, then translation killing `d` and `e`, then the normal form.
- [ ] **A6 -- a single cell is infinite.** (L) `Witness.ellipse_realised` proves
      infinitely many points carry two active interior branches on one ellipse,
      but that set is a **union** of cells. To get one infinite *cell*, show the
      activity pattern is locally constant along the arc -- the two interior
      branches are active at every `s_n`, and no third candidate can join on a
      neighbourhood, because the other candidates are strictly below there.
      *Route worth trying first:* rather than a topological argument, show the
      activity set is literally the same finite set at every `s_n` -- compute
      `active f (sw k)` for the witness and prove it independent of `k`. That
      turns a limit argument into an algebraic one and would close the item
      outright.

### B. Rational coefficients

Yves chose the **parallel `RatQuaPol` structure** over a predicate on the real
development: everything stays computable and decidable, and it extends what
`Rational.lean` already has (`RatQuad`, `ratVertexBranch`, `ratEdgeBranch`,
`ratInteriorBranch`, `ratCand`, and the `kind_toQuad` bridge).

- [ ] **B1 -- `RatQuaPiece` and `RatQuaPol`.** (S) Vertices, rays and quadratic
      all over `Q`, mirroring `QuaPiece`/`QuaPol` field for field, plus
      `toQuaPiece` / `toQuaPol` and `toPlane : Q x Q -> Plane`. Sanity `example`s
      checking `toQuaPol` against a hand-computed instance, per `CLAUDE.md`
      Verification point 3.
- [ ] **B2 -- Lemma 1: the branches commute with the embedding.** (M) For each of
      the three families,
      `toQuad (ratXBranch q v ...) = XBranch (toQuad q) (toPlane v) ...`,
      under the nondegeneracy hypothesis that family needs (`ratEdgeCurv != 0`,
      `RatQuad.hessDet != 0`). This is `../CONJ_FIELD_PROOF.md` Lemma 1: the
      branch formulas are built from the piece data by field operations only, so
      they never leave `Q`. Name Lemma 1 in the docstring.
- [ ] **B3 -- Theorem 1: `cand` of a rational `QuaPol` is rational.** (S/M)
      `ratBranches` maps onto `branches` and `ratCand` onto `cand`, so every
      candidate quadratic of `toQuaPol F` is `toQuad` of a rational one. Corollary
      -- and this is the half `Shapes.lean` and `QuaPar.lean` actually use --
      **every tie conic `q1 - q2` has rational coefficients**, since `RatQuad`
      subtraction commutes with `toQuad`.
- [ ] **B4 -- retire the standing assumption in `Rational.lean`.** (S) Its
      module docstring currently says the computable classifier is "justified by
      CONJ_FIELD_PROOF Theorem 1". After B3 that is a proved lemma in this
      repository: rewrite the docstring to cite it, and make the census a check
      of a theorem rather than of an assumption.

### C. The shape of the biconjugate

Yves chose the **full Theorem 4 subdivision** (`../CONJ_FIELD_PROOF.md` Theorem 4
and section 5.1), not just the foundations. The key is C2: everything in the
table except the 2-cell rows is a corollary of it.

- [ ] **C1 -- `biconj` is convex and lsc.** (M) With `biconj` defined in A4:
      it is a supremum of functions affine in `x`, so the A1/A2 arguments
      transpose verbatim. Keep them as separate lemmas rather than copying the
      proofs -- if a shared lemma "a sup of affine functions is convex and lsc"
      falls out, extract it and use it for both `conj` and `biconj`.
- [ ] **C2 -- KEY LEMMA: `f**` is affine on the hull of the maximisers.** (M/L)
      Define `maxSet f s := {x | exists p in f.pieces, x in p.T and
      psi p.q s x = f*(s)}`, the points attaining the conjugate at `s`. Then for
      every `s` with `f*(s)` finite, `f**` equals the affine function
      `A_s(x) = <s,x> - f*(s)` on `convexHull (maxSet f s)`.
      *Proof, and it is short:* `A_s <= f` everywhere by Fenchel-Young, hence
      `A_s <= f**`. At a maximiser `x` in piece `p`,
      `A_s(x) = q_p(x) >= f(x) >= f**(x) >= A_s(x)`, so all four are equal. `f**`
      is convex (C1) and agrees with the affine `A_s` at every point of
      `maxSet`, so `f** <= A_s` on the hull; with the reverse inequality, equal.
      This is `../CONJ_FIELD_PROOF.md` section 5, written there for a triangle.
- [ ] **C3 -- row 5 of Theorem 4: a vertex of `f*` gives an affine 2-cell.** (S)
      Corollary of C2 with three active candidates: `f**` is affine on the
      triangle of the three maximisers, with gradient `s`. Include the sharpness
      half -- the cell is *exactly* that hull -- if it comes cheaply; if not,
      state the containment and note the gap.
- [ ] **C4 -- row 4 of Theorem 4: an arc of `f*` gives a ruled 2-cell.** (S/M)
      Corollary of C2 with two active candidates: `f**` is affine along each
      ruling `[x_i(s), x_j(s)]`, with value `<s,x> - f*(s)`.
- [ ] **C5 -- rows 1 to 3: a 2-cell of `f*`, by the rank of its Hessian.** (L)
      Where exactly one candidate `g` is active, the maximiser should be unique
      and equal to the point where `grad g(s)` sends `s`; then
      `f**(grad g(s)) = <s, grad g(s)> - g(s)`, the Legendre transform. The image
      of the cell is a 2-cell, a segment or a single point according as
      `Hess g` has rank 2, 1 or 0 -- that is, according as `g` is an interior,
      edge or vertex branch, which `Shapes.lean` already ties to the rank.
      *Risk to check first:* "exactly one active candidate implies a unique
      maximiser" is **not** obvious and may be false in degenerate cases. Test it
      before building on it; if it fails, weaken to "the maximiser set is a
      single point when the active branch is the interior branch of a
      nonsingular piece" and record the boundary.
- [ ] **C6 -- assembly: the subdivision of `f**`.** (XL) The sets
      `R(C) = union over s in relint C of convexHull (maxSet f s)` cover
      `dom f**` and meet only on their boundaries, so they *are* the subdivision.
      Expect this to be the item that gets scope-reduced; the honest partial
      result is "every point of `dom f**` lies in some `R(C)`", which is already
      most of the value.
- [ ] **C7 -- the convex-hull formula of section 5.1.** (L) `f**(x) =
      inf { sum_k lam_k q_k(z_k) : sum_k lam_k z_k = x, z_k in T_k, lam in the
      simplex }`. An independent route to the same object -- it does not mention
      `f*` at all -- and therefore the best available cross-check on C2 to C6.
      For bounded pieces the infimum is attained (compact simplex times compact
      pieces), which is what makes the formula closed and equal to `f**`; with
      recession directions that needs care, so **do the bounded case first** and
      quarantine the unbounded one.

## Blocked

- [ ] Nothing is blocked.

## Deferred -- do not start these unattended

These are on the roadmap but were explicitly held back on 2026-08-24. Starting
one costs a night and returns nothing certain.

- **Theorem 2, vertices are algebraic of degree at most 4 over `Q`.** Needs conic
  elimination -- a resultant, or a minimal-polynomial degree argument -- and
  mathlib `v4.33.0` may not have the machinery. The certain part of Track B is
  B1 to B4; this is not it. If a cheap special case appears while doing B3 (one
  tie conic degenerating to a **line**, where substitution leaves a rational
  quadratic and the degree is at most 2), land that alone and say so.
- **The degree-4, non-rational face on a ruled cell** (`DECISIONS.md`
  2026-08-23, parent repo). Out of reach with the current machinery.
- **Phase 8, the write-up.** Deliberately last. Yves: "there is a lot to do
  before thinking about a paper."

## Done recently



- [x] 2026-08-24 -- **A1, A2, A3: `f*` is convex, lsc, and has a convex domain.**
      `Convexity.lean`. Convexity is stated in epigraph form (`conj_le_of_combo`,
      `convex_epigraph_conj`) so that every arithmetic step stays in `R`;
      `convex_dom_conj` says the `top` region is the complement of a convex set,
      which sharpens the fifth conjunct of `conj_isQuaCon`.

- [x] 2026-08-24 -- **Phase 7: unbounded pieces.** `QuaPiece` gained
      `rays : Finset Plane`, with `T = convexHull verts + coneHull rays`. The
      face induction covers ray supports through `IsDirRep`, the candidate list
      gained the vertex-and-ray edge branches, and `FrankWolfe.lean` proves that
      a quadratic bounded above on `conv V + cone R` attains its supremum -- by
      induction on the ray count, with the curvature dichotomy taken over
      `conv R` rather than over the generators. `conj_isQuaCon` is now
      hypothesis-free and its `top`-cell conjunct is load-bearing, witnessed by
      `Sanity.cell_empty_rayPol_nonempty`.

- [x] 2026-08-22 -- **`QuaPar.lean`: the question the project exists for, closed.**
      `exists_not_hasParabolicEdges` -- the conjugate's subdivision need not be
      parabolic, so no `QuaPar` can store it. `IsEdgePair` requires infinitely
      many common points so the refutation cannot be dismissed as a tangency.
      `disc_eq_zero_of_vertexBranch_flat` isolates what [JOGO]'s argument does
      establish, and the module docstring carries the complete legality table.

- [x] 2026-08-22 -- **`Realization.lean`**: the certificate machinery. `le_conj`,
      `conj_le`, and for convex pieces `psi_le_interiorBranch` -- the interior
      branch bounds its piece *everywhere*, with no side condition. (2026-08-24:
      `le_conj` moved to `Selection.lean`, and `conj_le` was reproved without
      attainment so that it survives unbounded pieces.)
- [x] 2026-08-22 -- **`Witness.lean`: `ellipse_realised`.** An explicit two-piece
      `QuaPol` with infinitely many points where two interior branches are both
      active, all on one non-degenerate ellipse, with both maximisers genuinely
      inside their pieces. Needed a purpose-built example: the documents' own
      elliptical edge has **no rational point** (7-adic obstruction).

- [x] 2026-08-22 -- **`Shapes.lean`: which conic arises.** Type is decided by the
      ranks of the two branches' quadratic parts, across pieces. Headline:
      `not_flat_of_disc_neg`, an elliptical tie conic requires an interior branch
      -- the [JOGO] Theorem 6 gap as a theorem. Plus Theorem 3, which needed no
      continuity hypothesis after all.
- [x] 2026-08-22 -- **`Rational.lean`: computable classifier over ℚ** with the
      `kind_toQuad` bridge, and a census reproducing `doc/QuaConExample.md` 3.1
      (7 per piece, 23 after dedup), all ten curved edges of 3.3, and the four
      adjacent-pair degeneracies. Mostly kernel-checked; two counts use
      `native_decide` and are audited in `SORRY_LEDGER.md`.

- [x] 2026-08-22 -- **`conj_isQuaCon` PROVED.** `lake build` green, zero `sorry`,
      `#print axioms` lists only `propext`, `Classical.choice`, `Quot.sound`.
- [x] 2026-08-22 -- **S3-S9 complete**: `exists_branch_eq_max`, by strong
      induction on the number of vertices of the face carrying the maximiser.
      Compiled on the first attempt once the infrastructure was in place.
- [x] 2026-08-22 -- `Bary.lean`: `mem_convexHull_erase`, `mem_convexHull_perturb`,
      `exists_descent`, plus plane geometry by scalar cross product
      (`exists_smul_of_cross_eq_zero`, `eq_zero_of_dot_eq_zero_of_cross_ne_zero`,
      `exists_combo_of_cross_ne_zero`, `cross_eq_zero_of_forall`). Using the cross
      product instead of `AffineIndependent` / `Collinear` avoided a large slice
      of mathlib API and made the case split concrete.
- [x] 2026-08-22 -- `foc_soc` (S5 and the second-order condition) and
      `interiorPoint_unique`.
- [x] 2026-08-22 -- **S8's core proved**: `psi_along_dir`,
      `exists_kernel_of_hessDet_eq_zero`, `psi_const_along_kernel`,
      `exists_dir_psi_const`.
- [x] 2026-08-22 -- **S1 and S2 proved**: `exists_isMaxOn_piece`,
      `exists_piece_eq_eval`, `exists_maximiser`; and `conj_ne_top`.
- [x] 2026-08-22 -- **Phase 3**: `QuaPol.lean`, `Candidates.lean`, `QuaCon.lean`;
      the three branch identities, each differential-tested before being written
      into Lean and then proved as a Lean identity.
- [x] 2026-08-22 -- **Phases 1 and 2 core**: Lake project on mathlib `v4.33.0`
      from local cache; `Quad.lean`, `Conic.lean`, eleven sanity witnesses.
- [x] 2026-08-21 -- target statement agreed and scaffolded.

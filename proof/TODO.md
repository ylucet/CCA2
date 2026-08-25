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

### B. Rational coefficients

Yves chose the **parallel `RatQuaPol` structure** over a predicate on the real
development: everything stays computable and decidable, and it extends what
`Rational.lean` already has (`RatQuad`, `ratVertexBranch`, `ratEdgeBranch`,
`ratInteriorBranch`, `ratCand`, and the `kind_toQuad` bridge).

### C. The shape of the biconjugate

Yves chose the **full Theorem 4 subdivision** (`../CONJ_FIELD_PROOF.md` Theorem 4
and section 5.1), not just the foundations. The key is C2: everything in the
table except the 2-cell rows is a corollary of it.

- [ ] **C7 residue -- the `>=` half, for CONVEX pieces.** (L) Two findings from
      2026-08-24 reshape this item, both in `DECISIONS.md`:
      * the general one-point-per-piece formula is **false** -- witness
        `Sanity.convRepVal_gt_biconj`, a concave quadratic on a segment -- so the
        statement must assume the pieces' quadratics are convex, as
        `../CONJ_FIELD_PROOF.md` §5.1 itself does;
      * with that hypothesis there are **no unknowns left**. mathlib supplies the
        separation layer (`ConvexOn.exists_affine_le_of_lt`) and
        `IsClosed.mul_left_of_isCompact` for "compact + closed is closed"; and
        because each `epi_k = graph_k + ({0} x [0,inf))` is convex, the convex
        hull of the finite union is the continuous image of
        `simplex x prod graph_k`, hence compact -- so the general "convex hull of
        a compact set is compact", which mathlib does **not** have, is not needed.
      *Plan:* define the convex-piece class; prove `epi (conv f)` closed by that
      route; conclude `conv f` lsc; then `affine_le_biconj` plus
      `exists_affine_le_of_lt` gives `conv f <= f**` in three lines.
- [ ] **C6 residue -- covering.** (XL) "Every point of `dom f**` lies in one of
      the cells" needs `f**` to have a subgradient at that point, which for a
      closed proper convex function holds only on the **relative interior** of the
      domain. `exists_affine_le_of_lt` does not help: it gives an affine minorant
      strictly below at the point, never one that touches. So this needs
      `intrinsicInterior`, which this development has never used. Do it after the
      C7 residue, and with the C5 residue in view -- they want the same layer.
## Blocked

- [ ] **C5 residue -- rows 1 to 3 as statements about cells.** The pointwise
      content is proved (see Done recently). Turning it into "a 2-cell of `f*`
      maps to a 2-cell of `f**`" needs the cells to be two-dimensional, i.e. the
      face-to-face regularity that `DECISIONS.md` 2026-08-21 deliberately does
      **not** claim. Whether to take that on is a decision only Yves can make,
      so this is parked rather than attempted.

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

- [x] 2026-08-24 -- **C7: the convex-hull formula, `<=` half.**
      `IsConvRep`, `convRepVal`, and `biconj_le_convRepVal`: every convex
      combination of points of the pieces bounds `f**` from above. The finite
      Jensen step is `Convex.sum_mem` applied to the epigraph C1 proved convex,
      so it needed nothing new. The `>=` half is not proved -- it says the
      right-hand side is already closed, which needs Fenchel-Moreau or a
      separation argument; it is now the last item of **Next up**.

- [x] 2026-08-24 -- **C6: the contact set, and `f**` as the affine envelope.**
      Scope-reduced. `biconj_eq_eval_of_mem_maxSet` -- at a maximiser `f**`
      equals `f` itself, so every cell C2 produces has its corners on the graph
      of `f` while its interior lies below. `affine_le_biconj` -- every affine
      minorant of `f` is a minorant of `f**`, which with
      `affineMinorant_le_biconj` characterises `f**` as the supremum of the
      affine functions below `f`, without Fenchel-Moreau. The covering and
      disjointness halves are back at the bottom of **Next up**.

- [x] 2026-08-24 -- **C5: rows 1 to 3, scope-reduced, and one refutation.**
      The risk the item flagged is real: "exactly one active candidate implies a
      unique maximiser" is **false**, machine-checked
      (`Sanity.exists_oneActive_manyMaximisers`). Rows 1 to 3 are therefore
      indexed by the branch attained: `interiorPoint_combo` (affine in `s`),
      `maxOn_unique_of_posDef` (unique in a positive definite piece),
      `edgePoint_mem_line` (row 2), and the constant vertex maximiser (row 3).
      Plus `convexHull_maxSet_subset_subgradSet`, the subdifferential
      correspondence C6 needs. The cell-level form is in **Blocked**.



- [x] 2026-08-24 -- **C3, C4: rows 4 and 5 of Theorem 4.** `Biconj.lean`.
      `biconj_affine_on_segment` (the ruled cell) and
      `biconj_affine_on_triangle` (the affine cell over a vertex of `f*`), both
      corollaries of C2. Plus the three bridges from an *active candidate* to a
      *maximiser*: unconditional for a vertex branch, and conditional on the
      stationary point lying in the piece for the other two -- which is the edge
      branch overshoot of `DECISIONS.md` 2026-08-21 showing up again.

- [x] 2026-08-24 -- **C2: the key lemma.** `Biconj.lean`,
      `biconj_eq_affineMinorant_on_hull`: wherever `f*(s)` is finite, `f**`
      equals the affine `A_s` on the whole convex hull of the maximiser set, and
      `s` is a subgradient there. Rows 4 and 5 of Theorem 4's table are
      corollaries. Compiled first try.

- [x] 2026-08-24 -- **C1: `f**` is convex and lsc.** The shared lemma did fall
      out: `supAffine g z = sup_y (<y,z> - g y)` in `Convexity.lean` carries
      convexity, convex domain and lower semicontinuity once, and both `conj`
      and `biconj` are instances of it -- `dot` being symmetric is what lets the
      two orientations share a definition. Nothing is proved twice.





- [x] 2026-08-24 -- **B2, B3, B4: Lemma 1 and Theorem 1 proved.** `RatInput.lean`.
      The three branch formulas commute with the embedding (and need **no**
      nondegeneracy hypothesis, since `x/0 = 0` on both sides);
      `isRat_of_mem_cand` is Theorem 1 -- every candidate of a rational input is
      rational -- and `isRat_tie_conic` is the corollary the shape results use.
      `Rational.lean`'s docstring now cites the theorem instead of assuming it.

- [x] 2026-08-24 -- **B1: `RatQuaPiece`, `RatQuaPol`, and the embedding.**
      `RatInput.lean`. Mirrors `QuaPiece`/`QuaPol` field for field over `Q`,
      with `toPlane`, `toQuaPiece`, `toQuaPol`, and `bounded_toQuaPol`. The
      census keeps running on `Rational.lean`'s computable list-based
      `RatPiece`; this class exists to carry theorems.

- [x] 2026-08-24 -- **A6: a single cell is infinite.** `Witness.lean`,
      `exists_infinite_cell`. Closed by **counting**, not by a limit argument:
      every `s_n` has an activity pattern, every pattern is a subset of the
      finite `cand f`, and there are infinitely many `s_n`, so some pattern
      recurs infinitely often and its cell is infinite. No continuity, no local
      constancy, no topology.

- [x] 2026-08-24 -- **A5: `disc < 0` really means ellipse.** `Ellipse.lean`.
      `IsEllipse` (the unit level set of a positive definite form in Cholesky
      position), the completed square at the centre, and the three cases:
      `a*det3 < 0` ellipse, `det3 = 0` the single centre point, `a*det3 > 0`
      empty. No rotation was needed -- completing the square twice gives a
      *shear* to normal form, and an ellipse is the affine image of a circle,
      not only the rotated one. `U3U6_isEllipse` upgrades the documents' own
      elliptical edge from "elliptic and non-degenerate" to a proved ellipse.

- [x] 2026-08-24 -- **A4: Fenchel-Young and the biconjugate.** `Biconj.lean`.
      `biconj f x = sup_s (<s,x> - f*(s))`, the affine minorant
      `A_s(x) = <s,x> - f*(s)` with `A_s <= f` and `A_s <= f**`, and
      `biconj_le_eval : f** <= f`. Everything is routed through `A_s` rather than
      through `<s,x> <= f(x) + f*(s)`, so no `EReal` addition of two possibly
      infinite terms ever appears.



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

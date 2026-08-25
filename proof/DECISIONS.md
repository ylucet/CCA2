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

## 2026-08-24 — "one active candidate implies a unique maximiser" is FALSE

- **Tried:** indexing rows 1 to 3 of `../CONJ_FIELD_PROOF.md` Theorem 4 by the
  *activity pattern* — the reading that a cell of `f*` carrying a single face
  function has a single-point subdifferential, so that `f**`'s cell there is the
  affine image of it. `TODO.md` C5 flagged the step as unproved and said to test
  it first.
- **Why it failed, and it is a refutation, not a search failure:** the piece
  `Sanity.rayPol` — the nonnegative `s1`-axis carrying the **zero** quadratic —
  has exactly one candidate at `s = 0`. Its single vertex gives the zero vertex
  branch; the zero quadratic has no direction of positive curvature, so the
  edge-branch filter is empty; its Hessian is singular, so there is no interior
  branch. That one candidate is active. And **every** point of the axis is a
  maximiser, so `maxSet` is a ray, not a point. Machine-checked as
  `Sanity.exists_oneActive_manyMaximisers` in `Biconj.lean`.
- **What replaced it:** rows 1 to 3 are indexed by **which branch is attained**
  and the rank of its quadratic part, which is what Theorem 4's table actually
  says. Interior branch: the maximiser is `interiorPoint q s`, affine in `s`
  (`interiorPoint_combo`), and unique inside a positive definite piece
  (`maxOn_unique_of_posDef`, by strict concavity). Edge branch: the maximiser
  lies on the fixed line `v + R(w-v)` for every `s` (`edgePoint_mem_line`).
  Vertex branch: the maximiser is the constant `v`.
- **What is left, and it is a decision not a difficulty:** rows 1 to 3 as
  statements about *cells* need those cells to be two-dimensional, which is the
  face-to-face regularity this project deliberately does not claim
  (`DECISIONS.md`, 2026-08-21). Whether to take that on is Yves's call; until
  then the pointwise statements above are the honest form.
- **Before retrying, fix:** do not weaken the hypothesis to "the cell is a
  singleton activity pattern". The counterexample kills every version of that.
- **Evidence:** `QuaConProof/Biconj.lean`, `Sanity.cand_rayPol`,
  `Sanity.maxSet_rayPol_not_subsingleton`.

## 2026-08-24 — Frank–Wolfe proved in V-representation; H-representation still not needed

- **Tried first, on paper:** the textbook route to Frank–Wolfe — H-representation
  plus Farkas, induction on the active constraints. Rejected for the same reason
  as 2026-08-21: mathlib `v4.33.0` has no polyhedra, no Farkas, no recession
  cones (`find Mathlib -iname '*olyhedr*' -o -iname '*arkas*'` returns nothing).
- **Also tried, on paper:** induction on the rays with the split
  `T = T' + cone{r}`. It fails on the `α_r > 0` branch: the inner maximum over
  `t ≥ 0` is a *piecewise* quadratic in the remaining variables, split by the
  hyperplane `L(y,r) = 0`, and intersecting a V-represented polyhedron with a
  half-space needs a vertex-enumeration lemma — a sub-project of its own.
- **What worked:** do the dichotomy on `Δ = conv R`, not on the generators, and
  induct on `R.card`.
  - Some `δ ∈ Δ` has `α(δ) ≤ 0`. Then `α(δ) < 0` is impossible (`ψ` would be
    unbounded along `x + ℝ₊δ`), so `α(δ) = 0`; finiteness then forces the slope
    along `δ` to be `≤ 0` everywhere, so pushing *back* along `δ` never lowers
    `ψ`, and pushing back until a ray weight vanishes lands in
    `conv V + cone (R \ {r})`. Induction.
  - Otherwise `α > 0` on all of `Δ`, which is compact, so `α ≥ m > 0` there and
    `ψ(c + λδ) ≤ P + λB - λ²m/2` uniformly. Beyond an explicit `λ₀` nothing can
    be optimal, so the supremum is attained on the compact truncation
    `conv V + [0,λ₀]·Δ`.
- **Why the dichotomy must be over `Δ` and not over `R`:** curvature is not
  preserved by conic combination. With `H = diag(1,-1)`, the rays `(1, 0.1)` and
  `(-1, 0.1)` both have `α > 0`, yet their sum `(0, 0.2)` has `α < 0`, so `ψ` is
  unbounded along a direction no generator reveals. A test on the generators
  would have been unsound, and silently so.
- **Consequence:** `conj_isQuaCon` now has **no** boundedness hypothesis, and its
  fifth conjunct is load-bearing — `Sanity.cell_empty_rayPol_nonempty` exhibits a
  one-piece `QuaPol` whose `⊤` cell is inhabited. The 2026-08-22 entry recording
  that the conjunct compared two empty sets is superseded for the general input
  class; it still describes what happens when `f.Bounded` holds.
- **Before retrying, fix:** n/a. If mathlib ever gains polyhedra, the H-side
  statement could be *added*, but nothing here would need reproving.
- **Evidence:** `QuaConProof/FrankWolfe.lean`; `lake build` green, 0 sorry,
  `#print axioms conj_isQuaCon` = [propext, Classical.choice, Quot.sound].

## 2026-08-24 — one induction for vertices and rays, via `IsDirRep`

- **Tried:** duplicating the selection lemma's face induction — one version
  stepping along `w - v` for vertices, one stepping along a recession direction.
- **Why it failed:** not refuted, but the two bodies differ only in which weight
  family the step touches, and the middle of the argument (first- and
  second-order conditions, the independence test, the kernel slide) is identical.
  Two copies would have to be kept in step forever.
- **What replaced it:** `IsDirRep W S gam nu d` in `Bary.lean` — a direction
  recorded by how it changes the weights, with vertex changes summing to zero and
  ray changes free. `isDirRep_vert_sub` and `isDirRep_ray` are the two
  generators; `foc_soc`, `mem_T_of_perturb` and `exists_descent_gen` are stated
  once, for the record. `branch_aux` inducts on `W.card + S.card`.
- **A trap worth recording:** a ray may be the **zero vector**, which would make
  the "pick a nonzero direction" step false. It is a third reduction in the
  induction (drop it — it contributes nothing), alongside a vanishing vertex
  weight and a vanishing ray weight.
- **Evidence:** `QuaConProof/Selection.lean`, `QuaConProof/Bary.lean`.

## 2026-08-22 -- the QuaPar question, closed: the subdivision need NOT be parabolic

- **Question:** CCA2's `QuaPar` requires `b^2 - 4ac = 0` of every edge conic, and
  [JOGO] Theorem 6 asserts the conjugate always admits such a subdivision. Does it?
- **Answer: no.** `exists_not_hasParabolicEdges`. For the two-piece `QuaPol` of
  `Witness.lean`, two interior branches are simultaneously active at infinitely
  many points and their tie conic is a non-degenerate ellipse, `disc = -1/16`.
- **`IsEdgePair` demands INFINITELY many common points**, deliberately. A weaker
  definition would have been refuted more cheaply but unfairly: a *same-piece*
  interior-versus-vertex pair always has `disc = -1/hessDet != 0` yet is a
  degenerate conic, i.e. a single point, and a critic could rightly say a point is
  not an edge. `../doc/QuaConExample.md` 3.3 records exactly three such tangencies.
  Requiring infinitude rules them out, so the refutation cannot be answered that
  way.
- **What [JOGO]'s argument DOES establish**, isolated as
  `disc_eq_zero_of_vertexBranch_flat`: a vertex branch against any flat branch --
  vertex or edge -- has `disc = 0`. That is what "when we compare two functions we
  always get one of them as linear" delivers, and within that case the conclusion
  is correct. The complete legality table is in `QuaPar.lean`'s module docstring;
  the two failures are vertex-versus-interior (never parabolic) and interior-
  versus-interior (unconstrained, and the one that gives a genuine ellipse).
- **Before retrying, fix:** n/a. The remaining sharpening would be to show a
  single *cell* is infinite rather than the union of cells that `IsEdgePair` uses;
  that needs the activity pattern to be locally constant and is in `TODO.md`.
- **Evidence:** `QuaConProof/QuaPar.lean`, `#print axioms` clean.

## 2026-08-22 -- realisation: the documents' elliptical edge has NO rational point (7-adic obstruction)

- **Tried:** witnessing the elliptical edge of `../doc/QuaConExample.md` 7.5's
  three-piece example by exhibiting a rational point on it and checking a
  certificate. That would have turned "an ellipse is in the candidate list" into
  "an ellipse is met" by pure computation.
- **Why it failed, and this is a theorem not a search failure:** that ellipse,
  `93 s1^2 + 38 s1 s2 + 39 s2^2 - 6 s1 - 482 s2 - 1003 = 0`, has **no rational
  point at all**. Homogenised and reduced mod `49` it has no primitive zero, so
  there is no `7`-adic point and hence no rational one. (Consistent with the tie
  point being `(-17/62 + sqrt(612030)/186, 3/2)`, and with sympy's ternary
  quadratic solver returning nothing, and with a brute search over denominators
  up to 200 finding nothing.)
- **Consequence:** realising *that* edge cannot be done by a computation. It needs
  an intermediate-value argument -- continuity of the branches plus openness of
  the certificate conditions -- which is a genuinely analytic proof.
- **What replaced it:** a purpose-built example whose tie conic is a **circle
  through the origin**, hence rationally parametrised.
  `QuaConProof/Witness.lean`: `q1 = x1^2 + x2^2` and `q2 = 2x1^2 + 2x2^2 + 8x2 + 8`
  on two disjoint triangles. Their interior branches differ by
  `(s1^2 + s2^2 + 16 s2)/8`, the circle `s1^2 + (s2+8)^2 = 64`, with
  `disc = -1/16` and `det3 = -1/2`: a non-degenerate ellipse. The slope-`1/n`
  parametrisation gives `s_n = (-16n/(n^2+1), -16/(n^2+1))`, and for `n >= 12`
  both maximisers lie genuinely inside their own pieces.
- **Result:** `ellipse_realised` -- infinitely many points of the plane have two
  interior branches simultaneously active, and all of them lie on one
  non-degenerate ellipse. Both branches are *attained*, not overshooting.
- **What is still NOT claimed:** that a single `cell` is infinite. The set proved
  infinite is `{s | both active}`, which is a union of cells; the activity pattern
  could in principle differ between the witnesses. Proving one cell is an arc
  still needs the connectedness the main theorem does not address.
- **Before retrying, fix:** to realise the documents' own example, do the IVT
  argument. To strengthen this one to a single arc, prove the activity set is
  locally constant.
- **Evidence:** `QuaConProof/Witness.lean`, `#print axioms` clean; the obstruction
  and the construction were both checked in scratch scripts first.

## 2026-08-22 -- which conics arise: type is decided by RANK, degeneracy is not

- **Question:** given a `QuaPol`, decide which conics arise.
- **What was proved** (`QuaConProof/Shapes.lean`). The discriminant of a tie set
  depends only on the two quadratic *parts*, and those are fixed by which kind of
  branch each side is: vertex branches have rank 0, edge branches rank 1 and
  positive semidefinite, interior branches rank 2. Hence, **across pieces**:

  | pair | `disc` of the difference |
  |---|---|
  | vertex, vertex | `0`, quadratic part vanishes: a line |
  | vertex, edge | `0`: parabolic |
  | edge, edge | `cross(d1,d2)^2 / (a1*a2) >= 0`: **never an ellipse** |
  | interior, vertex | `-1 / hessDet`: elliptic iff the Hessian is definite |
  | interior, edge / interior, interior | unconstrained |

- **The headline, `not_flat_of_disc_neg`:** an elliptical tie conic requires an
  **interior branch**. That is the [JOGO] Theorem 6 gap stated as a theorem: its
  proof assumes "when we compare two functions we always get one of them as
  linear", which is exactly the flat-versus-flat case, and under that assumption
  `disc_nonneg_of_flat` makes the subdivision genuinely parabolic. The assumption
  fails at Step 3b, where non-adjacent pieces contribute two interior branches.
- **Theorem 3 needed NO continuity hypothesis.** Its hypothesis is algebraic --
  `q2 - q1` factors as a product of two affine forms -- and continuity across a
  shared edge is merely what *produces* that factorisation. So `QuaPol` did not
  have to grow a hypothesis, which had been the concern when the option was
  offered. `det3_interiorBranch_sub_of_factorisation`.
- **Degeneracy is NOT determined by the quadratic parts**, so those theorems are
  same-piece-specific. Interior against a vertex branch of the *same* piece has
  `det3 = 0`, hence a single point -- which is why `doc/QuaConExample.md` 3.3
  finds three pairs that "touch at a single point and are not edges". A check
  during the work caught this: 3.3 also lists `I4|V3` as a genuine ellipse, which
  looked like a contradiction until the branches turned out to come from pieces 4
  and 1. The doc is right and the theorem is correctly scoped.
- **Before retrying, fix:** these are statements about *pairs of branches*, not
  about *realised edges*. Whether a given tie conic actually appears as a cell of
  positive length is the regularity question the main theorem does not address.
  Do not quote these as "f* has an elliptical edge" without that step.
- **Evidence:** `QuaConProof/Shapes.lean`; every formula symbolically verified in
  a scratch script before being written into Lean.

## 2026-08-22 -- the census: kernel-checked where possible, `native_decide` twice

- **Tried:** validating the rational classifier entirely by `decide`, so that
  every census check is kernel-checked.
- **Why it failed:** `decide` gets stuck reducing `Rat` equality in the kernel --
  it stops at `ctorIdx.beq` after unfolding the `DecidableEq` instances. The
  numbers involved (`det3` around `10^12`) are not the problem; the instance
  chain is.
- **What replaced it:** `norm_num` with the definitions unfolded, which produces
  an ordinary kernel proof. That covers the ten curved edges of 3.3, the four
  adjacent-pair `det3 = 0` instances, the `I1` reproduction and the parabola
  witness -- all `#print axioms` clean.
- **The two exceptions**, `census_sevenPerPiece` and `census_twentyThree`, run a
  whole `List` computation and use `native_decide`, so they carry an extra
  `native_decide` axiom and rest on the Lean **compiler**. That is confined to
  those two declarations; `conj_isQuaCon` and everything in `Shapes.lean` remain
  clean, and the audit in `SORRY_LEDGER.md` lists exactly which declarations are
  affected.
- **Before retrying, fix:** if the two counts ever need to be kernel-checked,
  the route is a `Finset`-free reformulation with `Decidable` instances that
  reduce, or `decide +kernel`. Not worth it for a validation check.
- **Evidence:** `QuaConProof/Rational.lean`; the same three counts were computed
  independently in a scratch script in exact rational arithmetic first.

## 2026-08-22 -- Caratheodory NOT used; induction on the face plus the scalar cross product

- **Tried:** the plan's S3, which routes through
  `Caratheodory.minCardFinsetOfMemConvexHull` to get a minimal *affinely
  independent* subset carrying the maximiser, and then cases on its cardinality
  being 1, 2 or 3.
- **Why it was dropped:** not refuted, but two frictions showed up on contact.
  Getting `W.card <= 3` out of `AffineIndependent` in the plane needs a
  finrank argument, and turning card-minimality into strict positivity of the
  barycentric weights needs the same erase-lemma that the direct induction needs
  anyway. So Caratheodory was buying only the affine independence.
- **What replaced it:** strong induction on `W.card` directly, with the case split
  driven by the **scalar cross product** `cross a b = a.1*b.2 - a.2*b.1` rather
  than by `AffineIndependent` / `Collinear`:
  - a vanishing weight drops a vertex (`mem_convexHull_erase`) -- this is the
    induction step, and it also means "all weights positive" may be assumed;
  - with all weights positive, `cross (w-v) (u-v) = 0` for every `u` says the face
    lies on a line, and then the maximiser is on that line
    (`cross_eq_zero_of_forall`, `exists_smul_of_cross_eq_zero`) and the edge or
    vertex branch applies;
  - otherwise some `u` gives `cross (w-v) (u-v) != 0`, two independent first-order
    conditions, and hence `grad q(x) = s`
    (`eq_zero_of_dot_eq_zero_of_cross_ne_zero`).
  No affine-space API is used anywhere, and the whole of `Bary.lean` is 282 lines.
- **Before retrying, fix:** only if the development is generalised past two
  dimensions, where the cross product stops being a scalar and Caratheodory
  becomes the right tool again.
- **Evidence:** `QuaConProof/Bary.lean`, `QuaConProof/Selection.lean`
  (`branch_aux`). `PROJECT_PLAN.md` Phase 4's S3 line is ticked with the deviation
  noted.

## 2026-08-22 -- the proof is complete; what it does and does not say

- **Tried:** nothing ruled out. This entry exists so that a later reader does not
  over-read the result, which is the failure mode this project has been most at
  risk of.
- **What is proved:** `conj_isQuaCon`, `sorry`-free, `#print axioms` clean. For a
  `QuaPol` whose pieces are convex hulls of finite vertex sets: the cells defined
  by activity pattern are pairwise disjoint and cover the plane, `f*` is one fixed
  quadratic on each, only finitely many are nonempty, `cell 0` is exactly where
  `f*` is `+inf`, and every cell with two distinct active quadratics lies inside a
  conic. Plus `active_nonempty`: at every point of the plane some candidate is
  active, which is the substance -- `f*` really is one of finitely many explicit
  quadratics everywhere.
- **What it does NOT say**, all three deliberate:
  1. **Nothing about regularity.** No dimension, no connectedness, no arcs, no
     face-to-face CW structure. A "cell" here is an activity level set and could a
     priori be disconnected or lower-dimensional. This was the agreed level; see
     the 2026-08-21 entry.
  2. **`IsConic` is containment, and the classification is by invariants.**
     `disc` and `det3` are proved invariants with computed values; there is still
     no theorem saying `disc < 0` makes a set an ellipse. So the theorem's
     "lies inside a conic" is exact, but "line / parabola / ellipse / hyperbola" is
     currently a statement about two numbers, not a proved normal form.
  3. **Bounded pieces only.** Stage 1. And because every piece is compact, the
     conjugate is finite everywhere, so conjunct five is currently the equality of
     two empty sets.
- **Before retrying, fix:** n/a. Each of the three is a `TODO.md` item.
- **Evidence:** `SORRY_LEDGER.md`; `CURRENT_STATE.md` "Blocked / open questions".

## 2026-08-22 — Stage 1 has `dom f* = ℝ²`, so the `⊤` cell is EMPTY until Phase 7

- **Tried:** nothing ruled out; this records a finding that makes one conjunct of
  the main theorem weaker than it looked.
- **What was found:** `QuaPol.conj_ne_top` is provable outright at Stage 1. Every
  piece is a convex hull of a finite set, hence compact, so each per-piece
  supremum is a real number, and there are finitely many pieces. Therefore the
  conjugate is **finite everywhere** and `dom f* = ℝ²`.
- **Consequence:** `PROJECT_PLAN.md` §0.5 treats `cell ∅` as "the distinguished
  `+∞` cell". At Stage 1 it is simply **empty**, and so is `{s | f* s = ⊤}`; the
  fifth conjunct of `conj_isQuaCon` is currently the equality of two empty sets.
  It is not wrong, and it costs nothing to carry, but it is **not yet
  load-bearing** and should not be quoted as evidence that the `+∞` region is
  handled. It becomes substantive only in Phase 7, when unbounded pieces make
  `dom f*` a proper subset of the plane.
- **Also simplified:** `selection` no longer needs a `conj f s ≠ ⊤` hypothesis,
  since that is now automatic. The version in `Selection.lean` is unconditional.
- **Before retrying, fix:** n/a. The `EReal` machinery stays, because Phase 7
  needs it; only the claim about what the fifth conjunct currently *proves* is
  narrowed.
- **Evidence:** `QuaConProof/Selection.lean` (`conj_ne_top`),
  `QuaConProof/QuaCon.lean` (`cell_empty_eq`, `dom_conj_eq_univ`).

## 2026-08-22 — a piece's quadratic is a `Quad`, not a separate `(Q, β, γ)`

- **Tried:** `PROJECT_PLAN.md` §0.2 as written, where `QuaPiece` carries a
  self-adjoint operator `Q : E →L[ℝ] E` with a proof field `Q_symm`, a vector
  `beta` and a scalar `gamma`.
- **Why it failed:** not refuted, but it duplicates `Quad` and makes symmetry a
  *carried hypothesis* rather than a structural fact. Reusing `Quad` — the same
  six-coefficient record already used for the conjugate's face functions — makes
  the Hessian `[[2a, b], [b, 2c]]` symmetric by construction, with no proof field
  to thread through every lemma, and lets the piece side and the dual side share
  one algebra (`Quad.eval`, `ring`).
- **Consequence, and it is a good one:** `vertexBranch`, `edgeBranch` and
  `interiorBranch` all map `Quad → Quad`, so the statement "the conjugate is
  built from quadratics" is typed rather than narrated.
- **Before retrying, fix:** only if the development ever generalises past two
  dimensions, where a coefficient record stops being convenient. Nothing in the
  current plan needs that.
- **Evidence:** `QuaConProof/QuaPol.lean`, `QuaConProof/Candidates.lean`.

## 2026-08-22 — the branch formulas were oracled BEFORE being written into Lean

- **Tried:** writing the three candidate-branch coefficient formulas
  (`PROJECT_PLAN.md` §0.4) straight into Lean from the algebra.
- **Why that is not enough on its own:** a wrong formula still typechecks. The
  kernel would then certify a true theorem about the wrong candidate family — the
  exact "definitional risk" in `CURRENT_STATE.md`. So each formula was first
  differential-tested in a scratch script against direct numerical optimisation of
  `ψ(x) = ⟨s,x⟩ - q(x)` over the relevant set, on thousands of random instances.
- **What the oracle caught:** the first edge-branch run reported a mismatch of
  `1.1e-1`. The formula was right; the *oracle* was wrong, because its grid was
  clipped to `t ∈ [-60, 60]` and the stationary point `t* = L/α` escapes that
  range. Re-oracled against `ψ(t*)` directly, and separately against
  `ψ(t) ≤ branch` for random `t`, the agreement is `5.8e-11`. The same run
  confirmed that for `α < 0` the restriction of `ψ` to the line is unbounded above
  in 1008 of 1008 instances, which is why `edgePairs` keeps only `0 < edgeCurv`.
- **Before retrying, fix:** keep doing this for every formula that enters a
  definition rather than a theorem. Definitions are where the kernel cannot help.
  Note the failure mode observed here: a red from a hand-built oracle is at least
  as likely to be the oracle's fault as the formula's, so diagnose before editing.
- **Evidence:** the formulas are now Lean identities —
  `vertexBranch_eval`, `edgeBranch_eval`, `interiorBranch_eval`,
  `psi_le_edgeBranch`, `gradAt_interiorPoint`, all `sorry`-free.

## 2026-08-22 — the plane is `ℝ × ℝ`, not `EuclideanSpace ℝ (Fin 2)`

- **Tried:** `PROJECT_PLAN.md` §0.1 as written, which fixes
  `E := EuclideanSpace ℝ (Fin 2)` because mathlib carries an `InnerProductSpace`
  instance on it.
- **Why it failed:** not refuted, but the reason for it evaporated on contact.
  Nothing in `Quad.lean` or `Conic.lean` uses an inner product at all, and the
  conjugate needs only the pairing `⟪s, x⟫ = s.1*x.1 + s.2*x.2`, which is a
  two-line definition on `ℝ × ℝ`. Against that, `EuclideanSpace ℝ (Fin 2)` is
  `PiLp 2 fun _ : Fin 2 => ℝ`, so every coefficient manipulation goes through
  `Fin 2` indexing and `PiLp` coercions, and `ring` cannot see through them.
  `ℝ × ℝ` gives `s.1`/`s.2` directly and `ring` closes the coefficient identities
  outright — `det3_translate`, a genuinely messy polynomial identity, is one
  `simp only [...]; ring`.
- **What is not lost:** `ℝ × ℝ` is a `NormedAddCommGroup` and `NormedSpace ℝ`, so
  `convexHull`, compactness of the hull of a finite set, continuity, and
  Carathéodory all apply unchanged — those are the mathlib facts the selection
  lemma (§0.6) actually needs. The sup norm rather than the Euclidean norm gives
  the same topology, and no step of the plan depends on which norm it is.
- **Before retrying, fix:** switch only if a later phase genuinely needs
  orthogonality (a rotation normal form for the conic trichotomy might). Even
  then the cheaper move is a local `dot` definition plus the two-by-two rotation
  written out, not a change of the ambient type.
- **Evidence:** `QuaConProof/Quad.lean`, `QuaConProof/Conic.lean`; `lake build`
  green, `#print axioms` clean on all eight top-level results, 0 `sorry`.

## 2026-08-22 — first Lean code: what is proved, and what "classification" does NOT yet mean

- **Tried:** nothing ruled out; this entry exists to stop a later session from
  over-reading `Conic.lean`.
- **What is actually proved:** `IsConic` is defined, a conic is never the whole
  plane (`IsConic.ne_univ`, resting on `Quad.eq_zero_of_eval_eq_zero`), the
  equality locus of two distinct coefficient vectors is a conic
  (`isConic_eqLocus` — the lemma the main theorem consumes), and `disc`/`det3`
  are defined, pinned against the textbook `3×3` determinant, and shown invariant
  under translation and covariant under scaling.
- **What is NOT proved:** the *geometric* trichotomy of `PROJECT_PLAN.md` §0.7 —
  that `disc < 0` with `det3 ≠ 0` makes the zero set an ellipse *in the sense of a
  normal form*, and so on. Right now `disc` and `det3` are invariants with
  computed values, not yet a proved classification. Getting there needs a
  rotation to diagonalise the quadratic part, which needs `sin² + cos² = 1`
  threaded through a messy coefficient identity. That is the next real piece of
  work in Phase 2 and it is deliberately not started.
- **Before retrying, fix:** n/a. Sequenced, not blocked.
- **Evidence:** `TODO.md` Phase 2, remaining items.

## 2026-08-22 — "looks straight" is not "is a line": both are decided by det3, and the figure's two suspect edges are genuine hyperbolas

- **Tried:** reading the edge type off `doc/QuaCon.svg` row 3 by eye. The
  `U1|U6` and `U1|U2` boundaries are drawn essentially straight, which suggests
  they are lines — and since a line is a degenerate conic, that would have been a
  legal reading.
- **Why it failed:** it conflates a *metric* property (small curvature over a
  short arc) with an *algebraic* one (the conic factors into linear pieces).
  Recomputed from the primal data in exact rational arithmetic, both are
  **non-degenerate**:

  | edge | equation | `B^2-4AC` | `det3` | verdict |
  |---|---|---|---|---|
  | `U1\|U6` = `I3\|I5` | `2731 s1^2 - 4598 s1 s2 - 2189 s2^2 - 107420 s1 - 9988 s2 + 421996 = 0` | `+45054240` | `+260148962304` | genuine hyperbola |
  | `U1\|U2` = `E1-2\|I5` | `5969 s1^2 + 5390 s1 s2 - 847 s2^2 + 67532 s1 - 22748 s2 - 353236 = 0` | `+49275072` | `+2474882726976` | genuine hyperbola |

  (Sign convention: `det3` is cubic in the coefficients, so it flips when the
  whole equation is multiplied by `-1`, while `disc` is quadratic and does not.
  Only the **vanishing** of `det3` is scale-invariant, and that is all the
  degeneracy test uses. The values above are for the equations exactly as
  written. The one criterion that does read the sign is real-versus-empty ellipse,
  where the scale-invariant form is `(a+c)*det3 < 0`.)

  Both are irreducible over `Qbar` (`sympy.factor_list(..., extension=True)`
  returns one factor), so neither is a line pair. They render straight because the
  *active* arc sits far out on the hyperbola where it hugs an asymptote:
  `|arc - centre| / transverse semi-axis` is `14.8` and `2.3`, radius of curvature
  over chord length is `169` and `11`, and the maximum sagitta is `0.13%` and
  `0.86%` of the chord — about `0.3` and `0.9` points on a `425`-point-wide panel,
  i.e. under one line width.
- **Before retrying, fix:** nothing to retry; this fixes a **specification**
  point. `IsConic`'s classification in `PROJECT_PLAN.md` §0.7 must be decided by
  `det3` (and `B^2-4AC` for the type), never by any notion of how the curve looks
  or how small its curvature is. Degeneracy is not a limiting case of curvature:
  the same figure contains nine genuinely degenerate pairs with `det3 = 0`
  exactly — five line pairs (`U2|U4`, `U2|U5`, `U3|U7`, `U4|U5`, `U6|U7`) and four
  with `B^2-4AC = 0` as well (`U1|U5`, `U2|U3`, `U3|U4`, `U4|U7`), the
  adjacent-piece case of the informal Theorem 3.
- **Evidence:** all 21 pairs of the seven faces classified exactly — 7 hyperbolas,
  5 ellipses, 5 line pairs, 4 parallel/repeated lines, and **no parabola**. The
  face quadratics were taken from `../doc/QuaCon.svg` row 3 and independently
  re-derived from the primal `(Q_k, beta_k, T_k)` data; the per-piece constrained
  sup was checked against dense sampling of each triangle to `~1e-4` (grid
  effect), and the arc points lie on their conics to `<= 1e-10`.

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

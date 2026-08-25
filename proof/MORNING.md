# Morning report — 2026-08-24 overnight run (QuaConProof)

Branch: **`overnight/2026-08-24`** — *not* a branch of my own. See the deviation
note below; it matters for how you review this.

## Deviations from `/overnight`, and why

**No new branch.** The skill says to isolate the run on `overnight/YYYY-MM-DD`.
This repository is shared with a **concurrent MATLAB session in `CCA2/`**, which
had already created and checked out `overnight/2026-08-24`. One working tree
means one `HEAD`: creating or switching a branch would have moved that session's
checkout under it mid-run. So this run commits onto the branch that was already
checked out, and every commit stages **only paths under `proof/`**, as
`proof/CLAUDE.md` requires. The two projects' commits interleave in `git log`;
`git log --oneline -- proof/` separates them.

**The dirty files at start are not mine.** `MORNING.md`, `conjCPLQ.m`,
`checkConjSymFree.m`, `conjConvexPolygon.m` in `CCA2/` belong to the MATLAB run.
The skill says to checkpoint a dirty tree; I have not, because they are outside
my write boundary. They are untouched.

**`CCA2/MORNING.md` is the other run's report.** This file, `proof/MORNING.md`,
is mine.

**One cross-session incident, resolved.** At 22:14 the MATLAB run committed with
`git add -A` and swept two of my in-progress files (`QuaConProof.lean`,
`QuaConProof/Ellipse.lean`) into its commit; it noticed and un-tracked them one
commit later (`968b62a`), leaving the files on disk untouched. My next commit
re-added them properly. **Nothing was lost** — all 17 Lean files are tracked,
`git status` on `proof/` is clean, and the build is green — but the history has a
small wobble around `8ea857b..968b62a`, and it is the second symptom of two
unattended runs sharing one working tree. If both projects are to run overnight
again, they need separate worktrees (`git worktree add`), not just separate
folders.

> **CORRECTED 2026-08-25, by the other session, and the last sentence above is
> wrong.** Both reports reached "use separate worktrees" independently, and
> neither checked it. Measured: `git worktree add <path> main` fails with
> `fatal: 'main' is already used by worktree at ...`, so worktrees FORCE two
> different branches — which is what caused the tangle in the first place and
> what this repo's standing "commit on `main`, don't branch" preference exists to
> avoid. **The fix is simply not to branch:** with both sessions on `main` there
> is no branch for either to be dragged onto, and interleaved linear history in a
> solo repo is not a defect. The `git add -A` sweep described above has its own
> fix — stage explicit paths — and is unaffected by any branching scheme.
> Worktrees do earn their keep where different branches are genuinely wanted at
> once (`arch/co1d` runs six method-paths that way); that is a different problem.
> See `DECISIONS.md` 2026-08-25.

**Pushed before starting, on request.** `origin/main` moved `b9243d3..64043b5`,
publishing Phase 7. Note that `main` did *not* contain the later plan commit
`86243b1` — it had been left behind when the MATLAB session branched — so that
commit and everything tonight sit on `overnight/2026-08-24`, unpushed, as the
guardrails require.

> **SUPERSEDED 2026-08-25.** `overnight/2026-08-24` was folded into `main` by
> fast-forward, all 43 commits from both runs were pushed to `origin/main` on
> request, and every other branch — local and remote — was deleted after
> verifying each had zero commits not already on `main`. `main` is now the only
> branch anywhere.

## What changed

- **`convf` defined** and `biconj_le_convf` packages the `<=` half as a statement
  about functions, so the next session has something to state the converse
  against.

- **A second refutation, and it is the most useful thing the night found.** The
  convex-hull formula of `../CONJ_FIELD_PROOF.md` §5.1 is **false** for this
  development's input class: `Sanity.convRepVal_gt_biconj` is a concave quadratic
  on a segment where the infimum is `0` and `f**` is at most `-1`. §5.1 derives
  the formula from "each piece's epigraph is convex", which needs the pieces'
  quadratics to be convex, and this project allows indefinite ones. So the `<=`
  half I proved is the whole truth in general, and the `>=` half belongs to the
  convex-piece case.
  With that hypothesis the remaining work has **no unknowns**: mathlib has the
  separation layer and `IsClosed.mul_left_of_isCompact`, and convexity of each
  piece's epigraph makes the relevant convex hull a continuous image of a
  compact set, so the missing "convex hull of a compact set is compact" is not
  needed after all. `TODO.md` carries the plan.

- **Second pass on the two residues, and it changed what they are.** Landed
  `isClosed_epigraph_conj` and `isClosed_epigraph_biconj` -- the epigraphs are
  closed and convex, which is the object every separation argument needs. Then
  went looking for the missing layer and found that **mathlib already has it**:
  `ConvexOn.exists_affine_le_of_lt` and `sSup_affine_eq`. So C7's `>=` half is
  *not* blocked on Fenchel-Moreau; the only missing hypothesis is that `conv f`
  is lower semicontinuous, which for bounded pieces is a compactness argument.
  Both residues rewritten in `TODO.md` with the route and the one real unknown,
  and `DECISIONS.md` records the correction so nobody writes "needs
  Fenchel-Moreau" again.

- **C7 done, `<=` half.** `biconj_le_convRepVal`: every convex combination of
  points of the pieces bounds `f**` above, i.e. `f** <= conv f`. The finite
  Jensen step fell out of `Convex.sum_mem` applied to the epigraph C1 had
  already proved convex -- no new machinery at all. The `>=` half (that the
  right-hand side is *closed*, hence equal to `f**` rather than an upper bound)
  needs Fenchel-Moreau; it is now the last open item.

- **C6 done, scope-reduced.** Two real theorems: `f**` equals `f` **at every
  maximiser** (so the cells C2 produces have their corners on the graph of `f`),
  and `f**` dominates every affine minorant of `f` -- which characterises it as
  the supremum of those, with no Fenchel-Moreau needed. The covering and
  disjointness halves went back to the bottom of Next up with a note on what
  they need (`intrinsicInterior`, and the regularity that is Blocked).

- **C5 done, scope-reduced, and it produced a refutation.** The risk the item
  flagged is real and now machine-checked: exactly one active candidate does
  **not** force a unique maximiser. The zero quadratic on the nonnegative
  `s1`-axis has one candidate and a whole ray of maximisers. Rows 1 to 3 are
  therefore indexed by the branch attained, not by the activity pattern -- which
  is what Theorem 4's table says anyway. `DECISIONS.md` entry written.

- **C3, C4 done.** Rows 4 and 5 of Theorem 4's table, both immediate from C2.
  Recorded honestly: the corners of those cells are the *maximisers*, not the
  active candidates, and the bridge between them is unconditional only for
  vertex branches -- an active edge branch can have its stationary point outside
  its own segment. Also hit git index contention with the MATLAB session for the
  first time; a retry loop handled it.

- **C2 done** -- the key lemma of Track C.
  `biconj_eq_affineMinorant_on_hull`: where `f*(s)` is finite, `f**` equals the
  affine minorant `A_s` on the convex hull of `maxSet f s`, with `s` as
  subgradient. The convexity step needed no `AffineMap` machinery: the set where
  `f**` lies below `A_s` is convex directly from `biconj_le_of_combo` plus the
  fact that `A_s` respects convex combinations, so `convexHull_min` finishes it.
  Compiled first try.

- **C1 done**, and it paid for itself. `Convexity.lean` was refactored around
  `supAffine g z = sup_y (<y,z> - g y)`: convexity, convex domain and lower
  semicontinuity are proved once there, and `f*` and `f**` are both instances.
  The `f**` versions are four one-line corollaries. Also dropped a `!= bot`
  hypothesis from the domain statement -- `bot` is bounded by `0` like anything
  else -- which matters because `f**` genuinely can be `bot`.

- **B2, B3, B4 done** (`RatInput.lean`). `../CONJ_FIELD_PROOF.md` Lemma 1 and
  Theorem 1 are now theorems of this repository: every candidate quadratic of a
  rational `QuaPol` is the embedding of a rational one, and so is every tie
  conic. Worth recording: the branch formulas commute with the embedding
  **unconditionally** -- the edge branch divides by the curvature and the
  interior branch by the Hessian determinant, but `x/0 = 0` in both `Q` and `R`,
  so the diagram commutes on the degenerate values too. The nondegeneracy
  hypotheses are needed to know the branch *means* something, not to know it is
  rational. `Rational.lean`'s docstring no longer quotes Theorem 1 as an
  assumption.

- **B1 done** (`QuaConProof/RatInput.lean`, new file). `RatQuaPiece` and
  `RatQuaPol` over `Q`, mirroring the real classes field for field, with the
  embedding `toQuaPol` and four sanity checks that no vertex, ray or coefficient
  is lost. `toPlane` is proved injective, which B2 and B3 will need.

- **A6 done** (`Witness.exists_infinite_cell`). One activity pattern containing
  both interior branches has an **infinite cell**, sitting inside the
  non-degenerate ellipse. The plan expected a local-constancy argument; a
  pigeonhole on `(cand f).powerset` does it in a dozen lines instead, and needs
  no topology at all. Compiled first try.

- **A5 done** (`QuaConProof/Ellipse.lean`, new file). `disc < 0` plus
  `a*det3 < 0` now *proves* the zero set is an ellipse, with the imaginary and
  point cases separated. Deviation from the plan, and it made the item much
  shorter: no `Quad.rotate`, no `sin^2+cos^2=1`. Completing the square twice
  gives a shear rather than a rotation, and `IsEllipse` is stated as "unit level
  set of a positive definite form in Cholesky position" -- the affine image of a
  circle -- so a shear is enough. The sign work all funnels through one identity,
  `q(centre) * disc = -det3`.

- **A4 done** (`QuaConProof/Biconj.lean`, new file). `biconj`, the affine
  minorant `A_s`, `A_s <= f`, `A_s <= f**`, `f** <= f`, and real-valued
  Fenchel-Young. Compiled first try. Note for the record: `f**` really can be
  `bot` (one unbounded piece with a concave quadratic gives `f* = top`
  everywhere), so there is no `biconj_ne_bot` to match `conj_ne_bot`.

- **A1, A2, A3 done** (`QuaConProof/Convexity.lean`, new file). `f*` is convex,
  lower semicontinuous, and `dom f*` is convex. Convexity is stated as an
  epigraph inequality rather than `f*(as1+bs2) <= a f*(s1) + b f*(s2)`, which
  keeps `EReal` multiplication out of every proof; the direct form is recorded as
  `conj_combo_le`. Build green, axioms clean.


## What is broken

- Nothing is red. One item is parked: the **cell-level** form of Theorem 4's rows
  1 to 3. Its pointwise content is proved; the step from there to "a 2-cell maps
  to a 2-cell" needs face-to-face regularity, which the project has deliberately
  not claimed since 2026-08-21.
- Also scope-reduced, and back at the bottom of Next up rather than parked: the
  covering half of C6 and the `>=` half of C7. Neither is red; both are precisely
  specified now. C7 needs one lemma (`conv f` is lsc) and mathlib supplies
  everything else; C6 needs subgradient existence, hence `intrinsicInterior`.

## Needs a decision

- **Should the project claim face-to-face regularity of the subdivision?**
  `DECISIONS.md` 2026-08-21 declined it deliberately, and nothing since has
  needed it — until Theorem 4's rows 1 to 3, which cannot be stated about
  *cells* without it. Options: (a) leave it, and state rows 1 to 3 pointwise as
  they now are — my choice, and what is committed; (b) add regularity as a
  hypothesis on the input class, which is cheap but narrows the theorem;
  (c) prove it, which is a project of its own. I took (a) because it is the
  reversible one.
- **Git contention.** This repository is shared with the MATLAB overnight run in
  `CCA2/`, and both sessions commit to the same branch and index. I hit a
  `index.lock` collision once; a retry loop handles it. Worth knowing before the
  next time two runs are started together.

## Where I stopped

**All seventeen planned items are done or honestly scope-reduced**, in 16 green
commits. Sixteen loop iterations, well under the cap. `lake build` is green with
**no warnings** and **0 sorry**, and `#print axioms` is clean on every headline
result. 17 files, 6101 lines.

Three items remain, none of them vague, all in `TODO.md`:

1. **C7 residue** — `conv f` is lsc for convex pieces. Fully scoped after
   tonight's two findings; there are no unknowns left, only work.
2. **C6 residue** — covering `dom f**`, which needs `intrinsicInterior`.
3. **C5 residue (Blocked)** — the cell-level form of Theorem 4's rows 1 to 3,
   waiting on your decision about regularity.

I stopped rather than starting the C7 build because it is a third pass on the
same item in one night and a ~300-line build with a new hypothesis class; the
ladder says come back to it on a later pass, and an unfinished file is the one
outcome an unattended run should not leave. Everything on disk is green.

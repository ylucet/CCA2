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

**Pushed before starting, on request.** `origin/main` moved `b9243d3..64043b5`,
publishing Phase 7. Note that `main` did *not* contain the later plan commit
`86243b1` — it had been left behind when the MATLAB session branched — so that
commit and everything tonight sit on `overnight/2026-08-24`, unpushed, as the
guardrails require.

## What changed

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

_(nothing yet)_

## Needs a decision

_(nothing yet)_

## Where I stopped

Run started. Preflight done: `lake build` green at `9cbf484`, 0 sorry,
`#print axioms conj_isQuaCon` clean. Working the `TODO.md` programme top to
bottom — A (finish the real case), B (rational coefficients), C (biconjugate
shape), 17 items.

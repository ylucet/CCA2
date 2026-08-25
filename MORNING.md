# Morning report — 2026-08-24 overnight run

Branch: `overnight/2026-08-24` (merge to `main` with a fast-forward; a parallel
`proof/` session commits to `main`, which is the second reason this ran on a
branch despite the standing "commit on main" preference).

Task: steps 1–8 of the sym-free `conj` plan. `conj` must end sym-free with the
symbolic Case C kept only as a fallback. `biconj` untouched.

## What changed

- **Conic level of the lattice** (`Conic.m`, `RatCon.m`, thin `RatPar.m`). The
  subdivision axis is now `Pol < Par < Conic`; mesh data moved to `RatCon`, and
  every existing type is unchanged (still a `RatPar`, still a `Par`, still
  parabola-only). Fast bucket 249/0, identical to its recorded baseline.
  The marker is `Conic.m` because `CON` is a Windows device name -- `Con.m` is
  the console and MATLAB cannot open it.
- **`ratQ`** -- exact rational arithmetic on coefficient vectors, with two
  canonical forms (a face has a scale, a conic does not). 17 tests.
- **`conicMeet`** -- two conics intersected through their exact integer Sylvester
  resultant; only the root is floating point. 12 tests.

## What is broken

- Nothing red. One bounded limit found and pinned rather than fixed: the exact
  integer layer cannot express a near-tangency closer than about 1e-3 (the
  resultant wants ~1e17 and the budget is 2^53). It RAISES rather than rounding,
  which is the designed behaviour; `DECISIONS.md` 2026-08-24 has the numbers.

## Needs a decision

- (nothing yet)

## Where I stopped

- (run in progress)

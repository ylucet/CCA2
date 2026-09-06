# Morning report — 2026-09-05 overnight run

Branch: `overnight/2026-09-05`

Started from `main` @ `885d4fd`, tree clean. The slow bucket (`--slow -j 4`) was started in the
background at the top of the run, once, as umbrella rule 9's daily gate.

**A note on the branch.** Stored feedback says "commit on main, don't branch" for this repo. This
run branches anyway because `/overnight`'s safety model rests on it — abandoning the night is one
`git switch`, and the morning review is one diff. Say if you would rather these landed on `main`
and I will merge them straight across.

## The list for tonight

From the task list in session, items 2-6:

- [x] 2. `shearedRecession` — what is left of the unbounded envelope
- [x] 3. A CHAIN of edges — done for the non-collinear case, which is the whole reachable one
- [ ] 4. `conjQ.m` (88.4%) and `ratQ.m` (88.1%) coverage
- [ ] 5. Classify `isConvexQ`'s uncovered guards in `SUPPORT_MATRIX.md`
- [ ] 6. The IMPLICIT ALGEBRAIC FACE — the envelope's missing type

## What changed

**Item 2 — `shearedRecession` is gone; it was never a separate algorithm.** The lifted recession
ray's height must be the MINIMUM of `q`'s rate over the piece, not the constant `L'd` (which is
that minimum only when `H*d = 0`). Verified against two closed forms: `x*y` on the first quadrant
gives `co f = 0`, and `x*y` on `{x>=0, y>=2}` gives `co f = 2x`. The second discriminates --
`L'd` would give `co f = 0` there, a strict underestimate that is still a convex minorant, so only
a TIGHTNESS check catches it. Sweep 196/196, 0 wrong, after adding an `indefXY` class that the
old diagonal `indef` could not reach.

**The daily slow gate is GREEN: 98/0 over 5 suites (765 s, `-j 4`).** That clears this session's
edits to `QuaPol`/`QuaPar`/`RatPol` (createDom, isEdgeConvex, isConvex) and the `biconjQ` rewrite,
none of which had been through the slow bucket before tonight.

**Item 3 — the CHAIN of edges is enveloped, and the refusal it replaced rested on a false premise.**
`ambiguousChain` claimed several function rows with no `F` left it undecided which edge took which.
For a chain, row *j* IS edge *j*'s function -- `QuaPol.eval` already reads it that way, and a chain
must carry one row per edge for `eval` to work at all. A chain of two non-collinear segments has a
2-D hull, so it is an ordinary lower-hull answer; only the GATHER changes, iterating edges instead
of faces. Closed-form fixture: `-x^2` on (0,0)-(2,0) and `0` on (0,0)-(0,2) gives `co f = -2x` on
the triangle, finite at (0.5,0.5) where `f` itself is `+infinity`.

Two findings alongside it. `chainNotImplemented` was UNREACHABLE -- `ambiguousChain` fired first
for every well-formed chain -- and is now narrowed to the COLLINEAR chain, which is genuinely thin
and still refused. And every chain in `examples2` carries `x^2` or `y^2` along its edges, so all of
them curve upward and the corpus alone cannot exercise the new path; the fixture above had to be
written.

## What is broken

_(nothing yet)_

## Needs a decision

_(nothing yet)_

## Where I stopped

_(run in progress)_

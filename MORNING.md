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
- [x] 4. `conjQ.m` and `ratQ.m` coverage
- [x] 5. Classify `isConvexQ`'s uncovered guards in `SUPPORT_MATRIX.md`
- [~] 6. The IMPLICIT ALGEBRAIC FACE — spiked, not built; one decision and one routine left

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

**Item 5 — `SUPPORT_MATRIX.md` gains section 3.1 for the EXACT envelope path.** §3 documented only
the legacy `biconjCPLQ`. The new section classifies every guard in `biconjQ`, `isConvexQ` and
`domainIsConvexD`. Four of `isConvexQ`'s five uncovered lines are internal invariants and one is an
unreachable skip, so its 92.4% is not a coverage gap. One genuine branch is left uncovered and is
named as such: the recession half-plane test (line 248), which needs a domain whose non-convexity
shows only at infinity. `domainIsConvexD` covers the same branch at the mesh level.

**Item 4 — coverage: `ratQ.m` 88.1% -> 93.8%, `conjQ.m` 88.4% -> 89.6%.** Same finding as
`QuaCon` earlier in the session: the gap was ARGUMENT VALIDATORS with no test. `ratQ` is the layer
everything exact rests on, so an untested guard there is the worst case -- it passes forever
whatever it does. Ten new tests cover the reachable ones (`badDenominator`, `badDivisor`,
`sizeMismatch`, `badSystem`, `cubicFace`, `badMap`, `degreeOverflow`, `notExact`, and `chk`'s own
2^53 ceiling), plus `detExact`'s 1x1 branch and two `conjQ` input guards.

The rest is internal invariants -- `signAtNoConverge`, `isolateNoConverge`,
`pseudoRemNotIntegral`, and `conjQ`'s ear-clipping guards -- classified in `SUPPORT_MATRIX.md`
section 3.1 rather than chased.

## What is broken

Nothing. Every bucket is green at the end of the run:

| bucket | result |
|---|---|
| slow (the daily gate, run once at the start) | **98 / 0** over 5 suites, 765 s, `-j 4` |
| fast | **533 / 0** over 37 suites |
| `.claude/biconjAllHSweep.m` | 196 answered, 196 checked, **0 wrong**, worst 1.1e-15 |

Two things to know rather than to fix:

* **The run is on a branch**, `overnight/2026-09-05`, against your stored preference for landing
  on `main`. `/overnight` requires it. Merge with `git switch main && git merge overnight/2026-09-05`
  -- it is a fast-forward.
* **`PLQ:biconjQ:unbounded` no longer exists for a 2-D domain** (it still names a 1-D ray domain
  in `thinEnvelope`). Anything outside this repo matching on that identifier would need updating;
  nothing inside it does.

## Needs a decision

**Item 6 — the implicit algebraic face. Spiked, deliberately not built.**

The spike (`.claude/bridgeEliminationSpike.m`, DECISIONS.md 2026-09-05 overnight) answered the
question that had to come first: the bridge's data IS exact over Z. The tie conic comes out as the
primitive integer vector `[5 0 5 -26 14 -10]`, the second conic likewise, and solving the pair
reproduces the closed-form face value with difference 0.000e+00. The branch selector is `t` in
`[0,1]`.

So the type is viable. What I did NOT do is build it, and the reason is a decision rather than a
difficulty:

  1. **The lattice changes.** `RETURN_TYPE.md`'s grid is (`Rat` > `Qua`) x (`Par` > `Pol`). An
     implicit algebraic face is a THIRD value on the function axis, above `Rat`, taking the grid to
     six cells. That document explicitly defers such a move until an operator PRODUCES the new
     kind rather than merely tolerating it — which is now the case, so the trigger it names has
     been met and the call is yours.
  2. **`eval` stops being a formula.** An implicit face evaluates by solving `p(x,y,.) = 0` for `z`
     and picking a branch. Every downstream operator has to accept that.

**The one thing that is unblocked and worth doing either way: a bivariate RESULTANT in `ratQ`.**
It is what eliminates `s` between the two conics, it is self-contained and testable on its own, and
it is needed whichever way the type decision goes. I would start there.

**Options, if you want one named:** (a) build the resultant next session and revisit the type after
seeing its output — my recommendation; (b) commit to the type now and build both together;
(c) leave both, and keep `convexPieceCoupling` refused.

## Needs a decision

## Where I stopped

Items 2, 3, 4 and 5 are done, committed and green. Item 6 is spiked and parked on your decision --
the only thing left that I cannot take myself.

**Four commits** on `overnight/2026-09-05`, on top of `main` @ `885d4fd`.

**Where I would start next**, in order:

1. **A bivariate resultant in `ratQ`.** Needed for item 6 whichever way the type decision goes,
   self-contained, and testable on its own against the spike's worked instance.
2. **The collinear chain** -- the one half of item 3 still refused. It is a 1-D fold of several
   quadratics along one span, and `thinEnvelope` already has the H-form machinery for a single
   segment.
3. **`isConvexQ`'s recession half-plane branch** (line 248), the one genuinely uncovered line
   that is not an internal invariant. It needs a `QuaPol` whose domain is non-convex only at
   infinity; `domainIsConvexDTest` has such a fixture at the mesh level to copy from.

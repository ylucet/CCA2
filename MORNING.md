# Morning report — 2026-09-06/07 overnight run (second)

Branch: `overnight/2026-09-06-2`, from `main` @ `080f878`. **10 commits, tree clean, nothing
pushed.** The first run of the night is already merged and pushed, so this branch starts from a
published base and merges as a fast-forward.

**All three items are done.** Both daily gates were run once, chained, and both are green:
**slow 98/0** and **verylong 37/0** (the 8 "incomplete" are deliberate `legacyPin` assumption
filters, not reds). Both had last run before the BigInteger Sturm change, so the gate earned its
cost.

## What changed

- **Item 1 — the rank-1 EDGE cell face, derived and verified.** Along the Hessian's null direction
  the objective is LINEAR, so the maximiser is on a facet or on the stationary LINE and nowhere
  else; a 2-D degenerate problem collapses into one 1-D concave quadratic per facet. The
  stationary case is a LINE, not a point, and exists only when `y` is parallel to `v`. 1554 points
  to 6.0e-13. `7ab553d`

- **Item 2 — all three ranks are ONE routine, not three.** Rank 2 turned out not to be a separate
  build: for any PSD `H` the candidates are the stationary set, one 1-D concave quadratic per
  facet, and the vertices as those facets' endpoints, and **the rank only changes what the
  stationary set IS** (a point, a line, or everything). The facet coefficients are rank-free.
  Verified over **every cell of five fixtures — 53 cells, all ranks — worst error 3.5e-11**.
  `58c546f`

- **Item 3 — split, and the split is the finding.** The curved-RAY row is now confirmed **by the
  guard FIRING**, which is stronger than the condition count it rested on, and a test pins it. The
  `notDegenerate` row **cannot** be settled that way — its own message says the condition cannot
  arise for an in-contract pair — so it is moved to Blocked with the reason, and
  `SUPPORT_MATRIX.md` now says which row rests on which kind of evidence. `ef5ad50`

- **A guard that could not report itself.** Found while attempting the above: `QuaPar`'s size
  guard passed a 1×5 logical against five `%i`, and **`error` does not expand an array across
  format specifiers the way `sprintf` does**. It raised `MATLAB:error:nonScalarInput` instead of
  `QuaPar:size`, masking the real refusal so completely the actual defect could not be read. Fixed;
  red-then-green test. `ef5ad50`

- **Beyond the three items — the envelope line is now finished except for assembly.** The
  non-bridge evaluator is **EXACT** (every candidate decision an integer sign, no division, an
  infinite endpoint spelled as an infinite denominator so unbounded facets need no branch), it
  **names its winner** (all four candidate kinds measured to occur), and the regions those
  candidates own are **POLYHEDRAL** — so unlike the bridge faces, they need no curved subdivision
  at all. `7c1bce0` `6e0e6b0` `4f590e6`

## What is broken

**Nothing.** fast **575/0** · slow **98/0** · verylong **37/0**. No known red.

**Worth reading, though: five measurement defects tonight, and every single disagreement this
session was the CHECKER, not the code.** Four under-reported (an oracle that could not reach a
boundary; a growth ratio that could not grow; an agreement that did not depend on the value it
claimed to validate; unconstrained maximisers used on a constrained problem). **One
over-reported** — pooling several distinct boundaries into one line fit manufactured "22 of 22
curved" — and then a **threshold three orders below the measurement's own resolution** manufactured
"23 of 26 curved" on the same question, while two arcs sat at 1e-16 under both verdicts. Each is
recorded in `DECISIONS.md` with what replaced it. Two cheap tells came out of it: **a number
identical across unrelated fixtures is measuring the harness**, and **a threshold must be compared
against the measurement's own resolution before it is used**.

## Needs a decision

1. **The implicit-face TYPE** — a third value on `RETURN_TYPE.md`'s function axis. Still yours, and
   now fully independent: every non-bridge face is quadratic or piecewise affine in `x` with
   polyhedral regions, so nothing else waits on it.

2. **`cPLQ/` is gitignored AND shadows a tracked file** (raised in the first report, unchanged).
   Options: **(a)** delete if superseded, **(b)** track it, **(c)** note the shadowing in
   `CLAUDE.md`. I would take (a) if dead, else (c) — not an unattended action.

3. **An exercise harness for every named guard** — one constructed input each, asserting the guard
   fires *and* reports its own identifier. A session of work; it would settle reachability and
   self-reporting for every refusal at once, and it is what the refuted static lint was a cheap
   substitute for. Worth doing, but it is a scope call.

## Where I stopped

At the token budget, not on a blocker. Every commit leaves the tree green.

**`TODO.md`'s Next up, ordered:**

1. **The ASSEMBLY** — walk the cell pairs, build each face, and stitch bridge and non-bridge faces
   into one mesh. Every ingredient is now built and verified; this is the join.
2. **Tidy `.claude/cellFaceRegionProbe.m`** — set its straightness threshold from the trace step
   rather than a fixed 1e-6, and split the one remaining pooled label pair. Hygiene: it does not
   change the answer, but the probe as it stands would mislead the next reader.

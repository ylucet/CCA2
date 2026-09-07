# Morning report — 2026-09-06 overnight run

Branch: `overnight/2026-09-06`, started from `main` @ `a57a3e2`. **11 commits, all green, tree
clean.** Nothing pushed.

**On the branch.** Stored feedback says "commit on main, don't branch" for this repo. This run
branched anyway because `/overnight`'s safety model rests on it. The history is linear, so
`git switch main && git merge overnight/2026-09-06` is a fast-forward whenever you want it.

**All four items you listed are done.** Two left a residue; one of those was then closed inside the
same run. The daily gates were run once each and both are green.

## What changed

**The headline, if you read one line: the verylong gate caught a six-day-old regression, and
running the bridge stages IN SEQUENCE caught two defects that every unit probe had passed.**

1. **Item 1 — two `N/R` guards are REACHABLE, one survives the re-check.** `SUPPORT_MATRIX.md`
   classified guards as unreachable on a scope claim refuted yesterday. Re-tested each guard's
   CONDITION on real `conjQ` output over eight two-piece kink fixtures, with both single-piece
   fixtures as a control: **curved ray edges 52** (control 0), **irreducible ellipse/hyperbola
   differences 44** (control 0), **faces with >1 curved edge 0**. First two reclassified; the third
   survives. `d3e0cba`

2. **A SIX-DAY-OLD REGRESSION.** `--verylong` came back 20 failed / 28 incomplete of 45 jobs, all
   one error. `dec3d7b` (2026-08-31, *"delete 18 uncovered functions with no caller"*) deleted
   `functionNDomain.times`, which HAS a caller — written `eqFun.*eqFun`. **An operator overload is
   called by its SYMBOL, never its name**, so a caller scan by identifier reports it dead. The only
   caller sits behind verylong, so fast/normal/slow stayed green for six days. Restored with a
   header saying why it is not dead code. **verylong is now 37 passed, 0 failed**, the 8
   "incomplete" being deliberate `legacyPin` assumption filters. `db06226`

3. **Item 2 — REGION LOGIC, done.** The bridge region's two boundary conics are built exactly by
   one `ratQ.substAffine` call each (residuals 1e-17); there is **no envelope arc**, so the region
   really is two sign conditions and an H-form face can hold it; and the **SIDES are decided
   exactly by a Sturm root count** — floating point only ever proposes a witness, every decision is
   integer. It also reproduces the annulus from a completely different route than yesterday's
   tracing, agreeing to every printed digit. `140e036` `58380c9` `3e6fdc7`

4. **Item 4 — the CLOSED conic edge convention.** `E(:,3) = 2` means CLOSED with no endpoints,
   validated by `RatCon.assertEdgeFlags` and refused by name in six consumers via
   `RatCon.assertNoClosedEdges`. The column had **no validator at all** before, and every consumer
   reads the flag as `==0`/`~=0`, so a 2 would have been read as a segment and indexed vertex 0.
   `e59287d`

5. **Item 3 — the pipeline END TO END, and it found two real defects.** `c256ae6`
   - **`quadSignOnPoly ~= 0` was the wrong prune test.** `+1` means `q >= 0`, which *allows* `q =
     0`: a TANGENTIAL tie is still a tie and still carries a bridge. It dropped a real bridge —
     cell pair (5,5) of the diagonal-kink triangle pair. The routine now returns `[s, strict]`.
   - **Both filter probes ignored the weighted basis's halves**, so the difference conic named a
     different curve from the one where the faces actually agree.
   - Re-measured, **the design conclusion is unchanged**: the 16-gon pair still drops from 235
     candidate eliminations to 39. What moved is the free whole-plane column, which loses two
     thirds of its claim.
   - End to end: 49 cell pairs → 16 meet → 8 survive; of 150 sampled points, 5 had a tied
     maximiser and **all 5 were claimed by a kept pair, none dropped**; value agreement 3.8e-8.

6. **Item 2's residue, closed: a BigInteger Sturm chain.** `polyPseudoRem` now falls back on
   `ratQ:overflow` around the **whole** double path — a first version guarded only the site where
   the overflow had been seen and the next fixture overflowed one line further down. It then
   exposed a **scaling bug in the membership quartic that the isotropic instance had been hiding**,
   because the sides never read `Phi`. Both instances now decide: **81/81** exact membership and
   **4000/4000** sides on each. `7c66200`

7. **Item 3's residue, scoped and started.** Each non-bridge face is the conjugate of a quadratic
   over a POLYHEDRAL cell — i.e. `conjQ`'s Case B minus two hypotheses. The counts are exactly
   `1 + n + n` per n-gon, so Case B covers 1 cell of `2n+1` today. The **rank-0 (vertex) cells are
   now in closed form**: one affine function on one polyhedral cone, verified on 13 cells over
   three fixtures, 1498 points to 3.6e-15, 0 declined. `e892961` `4e17758`

## What is broken

**Nothing.** No known red anywhere.

- fast **573/0** · normal **13/0** · slow **98/0** · verylong **37/0** (+8 deliberate
  `legacyPin` filters)
- `ratQTest` 80/0 with eight new tests; `QuaConTest` 23/0 with four.

Two things worth knowing rather than fixing:

- **Three checks in this run passed while measuring nothing**, and each was replaced rather than
  tuned. The 4000/4000 sides agreement did not validate `Phi` (the sides never read it). A
  "growth factor" comparing sampled sups at two radii gave 4.3 where the answer must be 10, because
  the sample is drawn in a disc and rejected against a cone. And the first draft of the end-to-end
  probe disagreed with the definition by 1.94 — **a wrong answer from the checker, not from the
  thing checked**, which looked exactly like a real defect for an hour. All three are recorded in
  `DECISIONS.md` with what replaced them.

## Needs a decision

1. **The implicit-face TYPE** — a third value on `RETURN_TYPE.md`'s function axis, `p(x,y,z)=0`
   plus a branch selector. Still yours, and everything measured since supports the one-axis
   estimate. **New tonight:** it is also *independent* of the remaining non-bridge work — those
   faces are all quadratic or piecewise affine in `x`, so none is implicit, and the two can proceed
   in either order.

2. **`cPLQ/` is gitignored AND shadowed.** There are two `functionNDomain.m` — the tracked one in
   the project root and an untracked older copy in `cPLQ/` — and the root one wins on the path. The
   `cPLQ/` copy still contains `times`, so anyone checking there would have concluded the deletion
   was harmless. Options: **(a)** delete `cPLQ/` if it is truly superseded, **(b)** track it,
   **(c)** leave it and note the shadowing in `CLAUDE.md`. I would take **(a)** if it is dead and
   **(c)** otherwise — deleting a directory is not something an unattended run should do.

3. **Should any coverage-driven deletion be gated on verylong?** The `times` break survived six days
   precisely because its only caller lives behind a bucket that is a daily gate rather than a
   development loop. Not a code question, so I have not acted on it.

## Where I stopped

At `/overnight`'s **40-iteration cap**, not on a blocker. The tree is clean and every commit leaves
it green.

**`TODO.md`'s Next up is ordered for the next session:**

1. **rank 1, the EDGE cells** (`n` of `2n+1`) — quadratic in one direction, a support function in
   the other. The mixed case, and the real remaining work on the envelope.
2. **rank 2, the INTERIOR cell** (1 of `2n+1`) — bounded and already Case-B-ready.
3. **Confirm the two reclassified guards through the wired symbolic path** —
   `.claude/reachabilityProbe.m('slow')` does it and takes hours; the condition evidence is already
   in hand, this would show the guard itself firing.

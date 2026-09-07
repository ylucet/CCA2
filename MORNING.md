# Morning report — 2026-09-06 overnight run

Branch: `overnight/2026-09-06`

Started from `main` @ `a57a3e2`, tree clean.

**On the branch.** Stored feedback says "commit on main, don't branch" for this repo. This run
branches anyway because `/overnight`'s safety model rests on it. Say the word and I will merge
straight across to `main` — the history is linear.

**Daily gate.** `--slow -j 4` already ran earlier in this session: **98 passed, 0 failed**, which
clears the stale gate the handoff flagged. `--verylong -j 4` was started in the background at the
top of this run — and it earned its keep, see below.

The four items came from Yves directly rather than from `TODO.md`'s **Next up**, and were written
into `TODO.md` as a **Next up** block at the start of the run.

## What changed

- **Item 1 DONE — two `N/R` guards are reachable, one survives.** `SUPPORT_MATRIX.md` classified
  guards as unreachable on a scope claim that was refuted yesterday. Re-tested each guard's
  CONDITION on real `conjQ` output (`.claude/guardConditionProbe.m`), over eight two-piece kink
  fixtures with both single-piece fixtures as a control: **curved ray edges 52** (control 0),
  **irreducible ellipse/hyperbola differences 44** (control 0), **faces with more than one curved
  edge 0**. The first two rows are reclassified; the third survives the re-check. `d3e0cba`

- **A SIX-DAY-OLD REGRESSION, caught by the verylong gate.** It came back 20 failed / 28 incomplete
  of 45 jobs, all one error. `dec3d7b` (2026-08-31, *"delete 18 uncovered functions with no
  caller"*) deleted `functionNDomain.times`, which HAS a caller — written as `eqFun.*eqFun`. **An
  operator overload is called by its SYMBOL, never its name**, so a caller scan by identifier
  reports it dead. Restored from `dec3d7b^` with a header saying why it is not dead code.
  `testConvex` went red to green; fast stayed 569/0. `db06226`

- **Item 2 DONE bar one residue.** The bridge region's boundary conics are built exactly by one
  `ratQ.substAffine` call each (residuals 1e-17); there is **no envelope arc**, so the region is
  two sign conditions and an H-form face can hold it; and the **SIDES are decided exactly by a
  Sturm root count** on a quartic — floating point only ever proposes a witness, every decision is
  integer. 4000/4000 agreement with the independent classifier on the isotropic instance.
  `140e036` `58380c9`

## What is broken

- **`ratQ.sturmChain` overflows at 3.7e18** on the anisotropic region-membership quartic, so item
  2's sides are exact for the isotropic instance and refuse for the anisotropic one. Not a new
  problem and not external: `sturmChain`'s own header records the pseudo-remainder squaring its
  leading coefficient each step, and `ratQ` already has the BigInteger fallback pattern for exactly
  this growth (`resultantVarBig`, `bigDivExact`, `polyDetBig`) — `polyPseudoRem`/`sturmChain` do
  not yet use it. Back on the list as item 2 (residue).

## Needs a decision

- **The implicit-face TYPE** — a third value on `RETURN_TYPE.md`'s function axis, `p(x,y,z)=0` plus
  a branch selector. Unchanged from yesterday and still yours. Everything measured since supports
  the one-axis estimate: the region's boundaries are conics, but `QuaCon` already stores general
  conics, and this run showed the region needs only two sign conditions.

- **The `cPLQ/` directory is gitignored and SHADOWED.** There are two `functionNDomain.m` — the
  tracked one in the project root and an untracked older copy in `cPLQ/` — and the root one wins on
  the path. That is how the `times` deletion could look harmless to anyone who checked `cPLQ/`,
  where the method is still present. Options: (a) delete `cPLQ/` if it is truly superseded,
  (b) track it, (c) leave it and note the shadowing in `CLAUDE.md`. I would take (a) if it is dead
  and otherwise (c) — deleting a directory is not something an unattended run should do.

## Where I stopped

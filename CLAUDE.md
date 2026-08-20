# CCA2

Project facts. The rules are in the umbrella `AI/CLAUDE.md`.

## What it is

**Conjugates, convex envelopes and derived operators for bivariate PLQ
functions**, in exact/closed form, not on a grid. MATLAB + Symbolic Math
Toolbox. Alpha; API still changing.

**The design target is the conjugate and biconjugate of a `QuaPol`.** Everything
downstream is a special case, and the special-ness is load-bearing: cases that
can never arise from a `QuaPol` conjugate — hyperbolic edges, for instance — are
**out of scope, not gaps**. Several `notImplemented` guards protect unreachable
branches; they are assertions, not missing features. `SUPPORT_MATRIX.md`
classifies every one.

## Running it

**MATLAB needs the UBC VPN** — the licence server is on campus. Without it
nothing here runs; say so rather than reasoning about what a run would show.

    bash .claude/suite.sh --fast       # 18 suites, one MATLAB, ~90 s
    bash .claude/suite.sh --normal     # 3 suites, ~215 s
    bash .claude/suite.sh --slow       # 5 symbolic suites, one MATLAB each
    bash .claude/suite.sh --verylong   # testcPLQ + testMaxMultiRegion
    bash .claude/suite.sh --slow -j 4  # 4 at a time, sharded per TEST
    bash .claude/suite.sh regionTest   # named suites only

For umbrella §7/§9: **fast** is closed-form numerics, **slow** is the symbolic
pipeline — that is the whole distinction, so moving a case to closed form also
moves its test out of the slow bucket. **verylong** is the pipeline endurance
run. One process per slow suite is deliberate: `solve`/`simplify` cannot be
interrupted from M-code, so killing the process is the only working timeout.

## Tools that already exist

- **`plqCheck.m`** — the definition verifiers for §7: `co f <= f` with equality
  at vertices; `f*(s) = sup_D <s,x> - f(x)` against a numeric sup; `f** <= f`
  and convex on sampled segments. The sampled sup is a *lower* bound, so
  `f* < sup_sampled` is a definite defect, the other direction within tolerance
  is expected.
- **`plqStage.m`** — the per-stage cache for §7's test splitting, keyed on
  (fixture, stage), invalidated by any `.m` mtime. Untracked derived data;
  deleting it must change nothing but runtime.
- **`.claude/suite.sh`** — buckets, timeouts, `-j`. Its comments say why each
  choice is what it is; read them before changing it.

## Documents — grep, don't read

`DECISIONS.md` 129 KB · `DESIGN.md` 145 KB · `SUPPORT_MATRIX.md` 83 KB ·
`TODO.md` 54 KB · `ALGORITHM.md` 14 KB.

- `DECISIONS.md` — what has been tried and ruled out. Read the **headings**
  (`grep '^## ' DECISIONS.md`), then the one entry covering what you will touch.
  Struck-through headings are overturned.
- `SUPPORT_MATRIX.md` — every guard as OK / GAP / not reachable / internal
  invariant. The answer to "bug or by design".
- `TODO.md` actions · `ALGORITHM.md` order of operations.
- `RETURN_TYPE.md`, `SCIP_READINESS.md`, `FARFIELD_FIX_PLAN.md` — narrower,
  older, partly stale.

## Session Handoff
@.claude/SESSION_HANDOFF.md

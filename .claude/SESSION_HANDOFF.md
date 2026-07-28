# Session Handoff

_Last updated: 2026-07-28T00:00:00Z_

## What happened this session

Four distinct pieces of work, in order:

1. **Wired `clipArcByHalfPlane` into `maxQuaPar`** — it now accepts **one** curved (parabolic-edge)
   input, closing its long-standing TODO. Operands are swapped so the curved one is always `g1`
   (so only arc-vs-half-plane clipping is ever needed); the polyhedral clip path is left
   byte-identical and the curved case goes through a new `clipPolyHalfPlaneCurved`. Also enforced
   the arrangement invariant for arcs and added an independent output validity check.
2. **Status audit → `SUPPORT_MATRIX.md` and `RETURN_TYPE.md`**, generated from the actual error
   guards rather than from `DESIGN.md`, classifying every guard OK / GAP / not-reachable / internal.
3. **Built the `RatPar` type lattice** (three commits: constructors → hierarchy → rename), giving
   every operator a single static return type.
4. **Removed all third-party copyrighted PDFs from git history** (`git filter-repo`), then
   force-pushed.

## Where things stand

- Branch: `main` @ `4caf619` — "Ignore the two publisher PDFs stripped from history"
- Pushed: yes (force-pushed after the history rewrite; upstream tracking restored)
- Tests: **251 passed / 1 failed** across 21 suites. The single failure,
  `testRegion/testCreation`, is a longstanding toolbox-compatibility issue unrelated to the
  conjugate pipeline — it predates all of this session's work.

### The type lattice (new)

Two axes — function type (`Rat` ⊃ `Qua`, pinning `den ≡ 1`) × subdivision type (`Par` ⊃ `Pol`,
pinning `Ec ≡ 0`) — as property-less abstract markers, so either axis can be queried alone
(`isa(f,'Qua')`, `isa(f,'Pol')`). The four named types are the grid's cells:

```
     RatPar < Rat & Par        (abstract; nothing produces one yet)
     /    \
RatPol    QuaPar               < RatPar & Pol  /  < RatPar & Qua
     \    /
     QuaPol                    < RatPol & QuaPar      (renamed from QuaPoly)

QuaParCPLQ < RatPar & Qua      same maths as QuaPar, no mesh; use isMeshed()
```

Three MATLAB constraints shaped this and are **load-bearing — do not "simplify" them away**:

- A property defined in two superclasses is **fatal and unresolvable** in the child
  (`conflictingSuperClassProperty` / `RedefinedProperty`). So `den` and `Ec` live on `RatPar`
  **alone**, and each type's pinned value is enforced by `set.` validators reading the object's
  traits. Declaring `den` on `RatPol` would make `QuaPol` unconstructible.
- A shared base constructor **runs once per inheritance path and clobbers** the earlier path's
  writes. Hence the **leaf-only constructor protocol** documented in `RatPar.m`: every constructor
  has a no-arg path that writes nothing.
- Method conflicts, unlike properties, **are** resolvable by a child override. All 39 methods
  defined in both `RatPol` and `QuaPar` are already overridden by `QuaPol`.

Degree is deliberately **not** a type (nothing dispatches on it) — see `RETURN_TYPE.md`.

### Repository hygiene (new)

`reference/` (13 papers) and `doc/s10589-…` + `doc/s10898-…` (the published [COAP]/[JOGO] versions
of record) were stripped from **all** git history and are now gitignored; local copies remain on
disk. Repo went 20 MB → 2.2 MB. Verified from a fresh clone: zero copyrighted PDFs on either
branch. `cplq-engine` was force-pushed too — its old tip still held 12 of them.

Backups kept outside the repo, safe to delete once you're satisfied:
`../CCA2-backup-54d3049.bundle` (full pre-rewrite history) and `../CCA2-reference-backup/`.

## Next steps

1. **Ask GitHub Support to garbage-collect the repo.** Force-pushing made the old PDF blobs
   unreferenced, but GitHub can keep them reachable by direct SHA until GC. This is the only
   remaining exposure before going public.
2. **Fix `biconj` for a single bounded triangle** — now blocker #1 for 0.1, and half the project's
   stated goal. `biconj` is `conj∘conj`, and the conjugate of a bounded-domain function is finite
   everywhere (an unbounded multi-face domain), which `conjCPLQ` rejects (`conjCPLQ.m:103`).
   Coverage is pinned by `conjCPLQTest.biconjCoverageByInputCase`: Case A works, **Case B fails**,
   Case C works via the symbolic path. Note the shape — the *symbolic* path supports the
   biconjugate while the faster *numeric* one does not. Fixing the unbounded multi-face case closes
   both this and `SUPPORT_MATRIX.md` §1.2.
3. **Refresh `SUPPORT_MATRIX.md`** — it still lists the return type as the top blocker (now
   resolved) and predates the `QuaPoly`→`QuaPol` rename.
4. **`maxQuaPar`: split a cell that already carries an arc** — 30 of 115 sampled valid splits hit
   this guard. Shares a root cause with two-curved-inputs: a piece has only one `Ec` slot, so
   **multi-arc pieces** is the natural next unit of work.
5. Then 0.1 tagging (user asked to be consulted first — **do not tag without being asked**).
6. Longstanding, unaffected: `partialConj` (unimplemented for every engine); `'pqp'`/`'graph'`
   engines; `RatPol.conj`/`biconj`/`add`; the `mergeL`/`removeTangent` exact-tie-point bug;
   `QuaPar.eval` returning `Inf`/wrong exactly *at* a result vertex (~1.4% — pre-existing on the
   polyhedral path too); exact `[LOCATELLI]` citation.

## Relevant files

- `RatPar.m` — hierarchy rationale, the constructor protocol, `set.den`/`set.Ec` validators,
  `kind()`, `isMeshed()`. **Read this before touching any constructor or class hierarchy.**
- `Rat.m` / `Qua.m` / `Par.m` / `Pol.m` — axis markers; each documents what behaviour belongs on it
  if the per-axis code is ever factored out of the concrete classes.
- `RatParTest.m` — 12 tests pinning the lattice contract (full `isa` table asserted entry by entry).
- `SUPPORT_MATRIX.md` — what works / what doesn't, generated from the error guards. §0 states the
  scope rule that keeps being misread: never-arising `RatPol`/`QuaPar` cases (e.g. hyperbolic
  edges) are **out of scope, not gaps**.
- `RETURN_TYPE.md` — why `RatPar` exists; the MATLAB probe results table.
- `maxQuaPar.m` — curved-edge support; header VALIDATION block has the sweep numbers.
- `conjCPLQ.m` — Cases A/B/C; line 103 is the unbounded multi-face gap behind next step 2.
- `.gitignore` — records *why* `reference/`, `cPLQ/` and the two publisher PDFs stay local.

# Session Handoff

_Last updated: 2026-07-13T00:00:00-07:00_

## What happened this session

Pure design/roadmap update — no code was written. On request, redirected `DESIGN.md`'s roadmap:
the project's focus is the **nonconvex-PLQ operator pipeline** (`conj`→`infConv`/`moreau`→
`lasryLions`/`proxAverage`), not growing `RatPol`/`RatPar` into a general-purpose toolbox. Set the
build order to (1) extend `add` to `QuaPar` [now flagged as the load-bearing prerequisite, since
`infConv` and `proxAverage` both call `add` on the `QuaPar` type `conj` produces from `QuaPoly`
input — the `RatPol` extension of `add` is deprioritized since nothing in this pipeline needs it],
(2) `infConv`, (3) `moreau`, (4) `lasryLions`, (5) `proxAverage`. Also **rewrote the `proxAverage`
formula**: replaced the old `−moreau(−λ₁M_μf−λ₂M_μg)` composition with a direct-conjugate
derivation — `T_μP = conj(λ(T_μf)*+(1-λ)(T_μg)*)` where `T_μh:=μh+½‖·‖²` — a sandwich of two
conjugations around a weighted `add`, with a third `conj` recovering the result; no call to
`moreau` at all. Documented that, like `infConv`, this is only valid for convex `f,g` (its last
step is a biconjugation, unlike `moreau`'s convexity-free single-conjugate identity). Updated
`DESIGN.md`'s Implementation status, II.4 method table, II.6 code+derivation, Part III mapping,
Part IV item 11, and II.7 file layout accordingly. Committed and pushed.

## Where things stand

- Branch: `main` @ `809722d` — "Redirect roadmap to nonconvex-PLQ pipeline; rework proxAverage
  formula".
- Pushed: yes, `origin/main` @ `809722d`.

## Next steps

1. **Extend `add` to `QuaPar`** (curved/`Ec`-edge clipping in `clipByFace`/`clipPolyHalfPlane`,
   generalizing `addQuaPoly.m`) — **load-bearing prerequisite** for steps 2 and 5 below. `RatPol`'s
   `add` extension stays open but deprioritized.
2. **`infConv(f,g,engine) = conj(add(conj(f,engine),conj(g,engine)),engine)`** — convex `f,g`
   only; see `DESIGN.md` II.6.
3. **`moreau`**, single conjugate via expand-the-square [HIRIART-URRUTY-07] — needs
   `addQuadratic`/`addScaledEnergy` (per-face coefficient update, not implemented on any class
   yet) on `QuaPoly` and `QuaPar`. Not built on `infConv`; valid for nonconvex `f` too.
4. **`lasryLions`**, pure composition of `moreau` (`−e_μ(−e_λ f)`) once step 3 exists.
5. **`proxAverage`** — new formula (see "What happened"), needs only steps 1 and 3 above; no
   call to `moreau`.
6. **`'pqp'`/`'graph'` conjugate engines are unimplemented** (`conj(f,'pqp'|'graph')` errors
   explicitly, not silently wrong) — porting Jakee Khan's parametric-QP code or Tasnuva Haque's
   entity-graph algorithm is future work, not this session's focus.
7. Still-open, untouched: the standalone `RatPol.conj` gap; the `delta>0, Delta=0` degeneracy
   open question (see `/home/ylucet/CCA2/3-edge.tex`'s Conclusion, outside this repo).
8. Worth a quick audit: other unfused copies of buggy shared code across
   `QuaPoly.m`/`QuaPar.m`/`RatPol.m` beyond `orderEdges` (found twice already, in prior sessions).
9. Test command (Frances, prefix every batch call — Maple toolbox conflicts with Symbolic):
   ```
   matlab -batch "restoredefaultpath; rehash toolboxcache; cd('/home/ylucet/CCA2/CCA2'); \
     res=runtests({'RatPolTest','convEnvCPLQTest','QuaParTest','conjCPLQTest','conjPieceCPLQTest','PLQVCTest','maxQuaParTest','addQuaPolyTest'}); \
     fprintf('\n==== SUMMARY ====\nTOTAL %d PASSED %d FAILED %d\n',numel(res),sum([res.Passed]),sum([res.Failed])); \
     exit(sum([res.Failed])>0)"
   ```

## Relevant files

- `DESIGN.md` — the authoritative design doc; this session's changes are in "Implementation
  status" (new "Next planned, in priority order" note), II.4 (method table), II.6 (`proxAverage`
  code block + derivation, reordered before/after `lasryLions`), Part III (mapping table), Part IV
  (new item 11), and II.7 (file layout lists the six not-yet-started files:
  `addQuaPar.m`, `addQuadratic.m`, `infConv.m`, `moreau.m`, `lasryLions.m`, `proxAverage.m`).
- `addQuaPoly.m` — the existing `QuaPoly.add` implementation; the pattern to generalize when
  extending `add` to `QuaPar` (step 1 above).
- `QuaPoly.m`, `QuaPar.m`, `RatPol.m` — the three operated classes; `scalarMul`/`negate` already
  implemented on all three, `add` only on `QuaPoly`.

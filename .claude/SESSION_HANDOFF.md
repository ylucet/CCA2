# Session Handoff

_2026-08-21_

## Blocked

- **Phase C1 per-term cost target** — EXTERNAL, needs Yves. Box terms are 0.01 s,
  so the old 40–60 s figure is stale; C2/C3 wait on the number.

## State

- Branch `main` @ `a59dc4e` — "Write the re-planned port down"
- Pushed: pending
- Tests (2026-08-20): fast 249/0 · normal 12/0 · slow 88/0 (−j 3) · verylong NOT run
- Known reds: `testMaxMultiRegion/testPCE2` — `plq_1piece`'s envelope route, not
  `plq_1p`'s; T6 territory, see DECISIONS 2026-08-19 (night)

## Next

1. **T3** — `symbolicFunction`'s payload to a RATIONAL coefficient vector. Order,
   counts and reasons: `TODO.md`, top section.
2. **T6** — re-run the `plq_1piece` fixture swap first: `testPCE2`'s domain is the
   triangle A.4 now gets right, so one of the three regressions may be gone.
3. `bash .claude/suite.sh --verylong` before the next tag (~73 min); not run since
   the A.4 and `exactQ` changes.

## Files

- `TODO.md` top section — the ordered port plan, incl. Row 7 spelled out
- `DECISIONS.md` 2026-08-20 — five entries; read the T1 and degree-4 ones
- `.claude/t1_multiquadratic_example.md` — why one extension is not enough
- `.claude/evalbench.c` — SCIP evaluation cost, 24 ns indexed vs 1670 scanned
- `CONJ_FIELD_PROOF.md` — **UNTRACKED**, another session's degree-4 proof
- `exactQ.m` — multiquadratic now; `signExact` vs the screened `sign`

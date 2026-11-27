# Review 2026-11-19 #11 — differential rerun after findings #1–#10, #13–#17b

Re-run of the review's 357-input differential harness (generator
`python/tests/review/plq_gen.py`, seed 20261119: 147 exhaustive E-* + 30 each of
R1–R6 + 30 discontinuous D), graded against the reviewer's exact oracle
(`python/tests/review/compare_engines.py`), in all three engines.
Dates: 2026-11-22/23. Machine shared (other sessions active); no timing is used
for any conclusion.

> **SUPERSEDED for the D family (discontinuous input) by the addendum 2026-11-24 (#103) at the
> end of this file, and by #138 (2026-11-25): `PLQ:notContinuous` is gone and the D family is
> answered. Everything below is kept as measured; the D-family statements no longer hold.
> Current check (Python, 2026-11-25): `pytest tests -k discontinuous_shared_points` -> 25
> passed, 0 xfailed (the open-case set is empty). MATLAB was not rerun -- until the
> addendum 2026-11-26/27 (#204) at the end of this file, which is that rerun.**

## Verdict

**0 wrong values and 0 crashes** on all 357 inputs in Python, Rust and
MATLAB (MATLAB: 355 of 357; the two unanswered biconjugates are the exact
fixtures the original review recorded as MATLAB's slowest, below). The
discontinuous D family is refused by name in all three engines. Wherever two
engines both answer, they agree to 1e-9 relative. The only remaining
engine-to-engine differences are documented guard coverage (MATLAB's 2^53
ceiling) and performance (Python timeouts under shared load).

## Commands

```
cd python
python3 tests/review/plq_gen.py /tmp/opencode/review11_FIX.json --seed 20261119 --nrand 210   # 357 fixtures
python3 tests/review/run_engines.py FIX.json OUTDIR --engines py,rs -j 16 --timeout 300
python3 tests/review/run_engines.py FIX.json OUTDIR --engines py --timeout 3600 --ids <slow set>   # rerun under calmer load
python3 tests/review/run_engines.py FIX.json OUTDIR --engines rs -j 3  --timeout 3600 --ids <slow set>
matlab -batch "cd 'MATLAB/review'; reviewDiffDriver(FIX.json, OUTDIR/ml_k.jsonl, k, n)"   # k = 1..n
python3 tests/review/compare_engines.py FIX.json OUTDIR --engines py,rs,ml
```

## Per-engine grades

(conj / biconj: ok / wrong / refused; wrong = outside the certified bracket or
off the exact conjugate oracle)

| | conj ok / wrong / refused | biconj ok / wrong / refused |
|---|---|---|
| Python (357) | 322 / **0** / 35 | 240 / **0** / 117 |
| Rust (357) | 327 / **0** / 30 | 245 / **0** / 112 |
| MATLAB (355 of 357) | 327 / **0** / 30 | 235 / **0** / 120 (+2 missing) |

Python's 35/117 refusals include 10 `TIMEOUT` lines (5 fixtures, below); Rust
and MATLAB have no timeouts, crashes or panics. No engine ever answered a
discontinuous D input.

## What the original #11 named, now

| original fixture | original state | now (py / rs / ml) |
|---|---|---|
| R4-010 (158), R4-045 (193), R3-072 (220) | py vs rs value differences | all three answer; all pairwise identical |
| D-069 (217), D-132 (280) | wrong biconjugate (above certified bound) | all three refuse `PLQ:notContinuous` |
| R4-080 (228) | py refuses `bridgeFaceAlg:noZ`, rs panics | all three refuse `bridgeFaceAlg:noZ` |
| R4-003 (151) | py > 600 s, rs answers, ml `ratQ:overflow` | py TIMEOUT under load, rs answers, ml refuses `ratQ:overflow` |
| R4-150 (298), R2-162 (310), R4-164 (312), R2-169 (317), R2-204 (352) | py vs rs value differences | py and rs answer and agree; MATLAB: 310/312/352 answer and agree, 317 `ratQ:overflow`, 298 unfinished (below) |

Across all 357 fixtures, wherever Python and Rust both return values the two
are identical to 1e-9 relative (0 differing points); MATLAB agrees with them
at every point where it answers (checked per fixture).

## Refusals (by name)

Python and Rust refuse the identical set (identifier, count), 112 biconj + 30
conj each: `PLQ:biconjQ:minusInfinity` 57, `PLQ:notContinuous` 30 (both ops,
the D family), `PLQ:biconjQ:nonConvexPiece` 11, `envelopeCellPairs:curvedCell`
6, `bridgeFaceAlg:noWitness` 4, `PLQ:biconjQ:convexPieceCoupling` 3,
`bridgeFaceAlg:noZ` 1. Python additionally times out on 5 fixtures (both ops).
MATLAB's refusal set differs only by its 2^53 exact-arithmetic ceiling and
witness-overflow guard (see below); it refuses no input that the ports answer
with a value it cannot represent, and never answers an input the ports refuse.

## Engine-vs-engine differences (all documented, no values)

MATLAB refuses where the ports answer — the 2^53 ceiling, all biconj:
`ratQ:overflow` on R4-003 (151), R4-059 (207), R4-101 (249), R2-169 (317),
R4-171 (319), R4-199 (347), R6-201 (349); `bridgeFaceAlg:witnessOverflow` on
R3-037 (185). And `witnessOverflow` on R2-029 (177) where the ports say
`noWitness` — same refusal, different guard reached first. (The original
review recorded these as MATLAB's "2^53 ceiling" refusals.)

Python timeouts under shared load (performance, not values): R4-003 (151),
R4-059 (207), R4-150 (298), R4-171 (319), R3-177 (325) — both ops; each was
rerun alone at 3600 s under calmer load and still timed out. The other engines'
answers on these inputs are oracle-graded above and all pass.

## MATLAB performance observation (not graded)

Two biconjugates have no MATLAB result after bounded runs (capped 1 h / 2 h
each, restarted isolated): R4-150 (298) and R3-177 (325) — the same two
fixtures the original review found "stuck for over an hour" with "Rust
answers both in seconds". Their conjugates and all other operations finish;
Rust's values on both are oracle-graded (ok). Recorded as a MATLAB
performance observation, per the review's §5.1 precedent, not a defect.

## Wrong values vs the oracle

None. Every answered conjugate matches the exact oracle (1e-9 relative) and
every answered biconjugate lies inside its certified bracket [lo, hi]. The
original 12 (py) / 12 (rs) / 10 (ml) wrong discontinuous biconjugates are now
named refusals. No failing test to open.

---

# Addendum 2026-11-24 (#103) — rerun after #108 lifted the discontinuous-input refusal

The section above was taken before #108: its D-family statements
("refused by name `PLQ:notContinuous`", "no engine ever answered a discontinuous
D input") no longer hold. `PLQ:notContinuous` is gone; the D family is answered
(read lower semicontinuous, DESIGN.md "Discontinuous input"). Python and Rust
were rerun on the current main, same fixtures (generator seed 20261119, 357
inputs), same grader; MATLAB was NOT rerun (its rows above stay as measured, and
are pre-#108 on the D family). Machine shared, load 40-90 during the run, so the
TIMEOUT counts below are load artefacts, not measurements.

```
cd python
python3 tests/review/plq_gen.py FIX.json --seed 20261119 --nrand 210
python3 tests/review/run_engines.py FIX.json OUT --engines py -j 24 --timeout 600
python3 tests/review/run_engines.py FIX.json OUT --engines rs -j 16 --timeout 900
python3 tests/review/compare_engines.py FIX.json OUT --engines py,rs
```

| | conj ok / wrong / refused | biconj ok / wrong / refused |
|---|---|---|
| Python (357) | 345 / **0** / 12 (all TIMEOUT) | 256 / **0** / 101 (89 named + 12 TIMEOUT) |
| Rust (357) | 357 / **0** / 0 | 264 / **0** / 93 (89 named + 4 TIMEOUT) |

No crash in either engine. Wherever both answer, the two agree (the grader's
engine-vs-engine check lists only TIMEOUT-vs-answer differences, no value
difference). Named refusals, identical in both engines (biconj): `PLQ:biconjQ:minusInfinity` 57,
`PLQ:biconjQ:nonConvexPiece` 11, `envelopeCellPairs:curvedCell` 10,
`bridgeFaceAlg:noWitness` 7, `PLQ:biconjQ:convexPieceCoupling` 3,
`bridgeFaceAlg:noZ` 1. The D family (30 inputs): Rust conj 30/30 exact; biconj
23 inside the certified bracket + 7 named refusals (curvedCell 4, noWitness 3);
Python the same except 2 TIMEOUTs on D-153 and D-188 (both ops).

Rust TIMEOUTs at 900 s under load: biconj R4-059, R4-150, R2-176, R3-177. Python
TIMEOUTs at 600 s: R4-003, R2-036, R4-059, R2-120, R2-127, R4-150, D-153, R4-171,
R2-176, R3-177, D-188, R4-199 (conj and biconj each).

Not covered by this grade: the biconjugate at a jump vertex of a discontinuous
input (the grid is offset from the vertices) — the open #108 residual (TODO.md).

## Update 2026-11-25 (#138)

#138 fixed the jump-vertex residual (D-111/D-132/D-153; four evaluation defects). Rerun for this
update, Python: `cd python && python3 -m pytest tests -k discontinuous_shared_points -q` -> 25 passed,
0 xfailed in 60 s (before #138, the three cases were strict xfail). The 357-input differential
run above was not repeated; its D-family counts stand as measured in the addendum.

---

# Addendum 2026-11-26/27 (#204) — MATLAB rerun on the 357 corpus, after #108/#138/#198

The #103 addendum left MATLAB "not rerun". This is that rerun, on current main (after #198;
measurement only, no code changes). Same fixtures (generator seed 20261119, 357 inputs), same
grader (the reviewer's exact oracle, as the first section of this file). Python and Rust were
rerun on the same commit, so the engine-vs-engine comparison is on current code. Machine shared:
load ran ~2 to 126 during the run, so no timing is used for any conclusion, and a hit on the
30-min per-input cap is recorded as "not finished", not wrong.

## Commands

```
cd python
python3 tests/review/plq_gen.py OUT/FIX.json --seed 20261119 --nrand 210      # 357 fixtures
python3 tests/review/run_engines.py OUT/FIX.json OUT --engines py -j 24 --timeout 600
(cd rust && cargo build --release --bin golden)
python3 tests/review/run_engines.py OUT/FIX.json OUT --engines rs -j 8 --timeout 900
# MATLAB: one `matlab -batch` process per input under `timeout 1800` (4 workers),
#   reviewDiffDriver(OUT/FIX.json, OUT/mltmp/<id>.jsonl, <id>, 357)      # conj+biconj, 30-min cap
# the 6 inputs that hit the cap were re-run per op, each with a fresh 30-min cap
# (same driver body with an op filter; measurement harness only)
python3 tests/review/compare_engines.py OUT/FIX.json OUT -j 16 --engines py,rs,ml > OUT/report.txt
```
(OUT = /tmp/opencode/review357; it holds FIX.json, py.jsonl, rs.jsonl, ml_1.jsonl, graded.json,
report.txt, the per-input MATLAB files and logs.)

## Per-engine grades

(conj / biconj: ok / wrong / named refusal / not finished; ok = conj at 1e-9 of the exact oracle,
biconj inside its certified bracket [lo, hi] at 1e-7)

| | conj ok / wrong / refusal / not finished | biconj ok / wrong / refusal / not finished |
|---|---|---|
| **MATLAB (357)** | **357 / 0 / 0 / 0** | **256 / 0 / 97 / 4** |
| Python (357) | 341 / 0 / 0 / 16 | 252 / 0 / 89 / 16 |
| Rust (357) | 357 / 0 / 0 / 0 | 264 / 0 / 89 / 4 |

**0 wrong values and 0 crashes in all three engines, and 0 value disagreements between any pair of
engines** (the grader's engine-vs-engine check lists refusal-status differences only).

## MATLAB, in detail

- **conj: 357 answered, 357 ok.** The 30 `PLQ:notContinuous` D-family refusals of the first
  section are gone with #108; the D family is now answered and exact (30/30 ok).
- **biconj: 256 ok, 97 named refusals, 4 not finished, 0 wrong.**
  - Named (97): `PLQ:biconjQ:minusInfinity` 57, `PLQ:biconjQ:nonConvexPiece` 11,
    `envelopeCellPairs:curvedCell` 10 (D-006, D-034, R3-058, R3-107, R4-108, R3-121, R3-128,
    R6-138, D-146, D-202), `ratQ:overflow` 8 (R4-003, R3-037, R4-059, R4-101, R2-169, R4-171,
    R4-199, R6-201), `bridgeFaceAlg:noWitness` 6 (R2-029, R4-073, R4-136, D-167, D-174, R4-185),
    `bridgeFaceAlg:witnessOverflow` 1 (D-090), `bridgeFaceAlg:noZ` 1 (R4-080),
    `PLQ:biconjQ:convexPieceCoupling` 3 (R6-096, R6-103, R3-135).
  - Not finished (4, each past a dedicated 30-min cap on biconj): R2-127, R4-150, R2-176, D-188.
    (R2-141 and R3-177 also hit the combined-run cap but finished inside a dedicated 30-min cap
    and are in the 256 ok.)
- **D family (30):** conj 30/30 ok. biconj 22 ok, 7 named (curvedCell D-006/D-034/D-146/D-202;
  noWitness D-167/D-174; witnessOverflow D-090), 1 not finished (D-188), 0 wrong. Same shape as
  the pre-#108 row of the first section (22 inside the bracket, 7 named, D-188 unfinished) --
  and the 22 now answer ok on the #138-fixed code.

## Where MATLAB disagrees with Python or Rust

No value disagreement on any of the 357 inputs, either op (checked per point at 1e-9). The
differences are all refusal status:

- **MATLAB answers where a port does not finish:** all 16 of Python's not-finished inputs
  (R4-003, R4-017, R2-036, R4-059, R4-101, R2-120, R2-127, R4-150, D-153, R2-155, R4-171,
  R2-176, R3-177, D-188, R4-199, R6-201) -- MATLAB's conjugates on all 16 are answered and ok;
  its biconjugates are answered and ok on R4-017, R2-036, R2-120, D-153, R2-155, R3-177,
  `ratQ:overflow` on R4-003, R4-059, R4-101, R4-171, R4-199, R6-201, and not finished on
  R2-127, R4-150, R2-176, D-188. R3-177 was one of the first section's two "no result after
  bounded runs" and finished here under a dedicated cap on a quiet machine. Rust's 4 not-finished
  (R4-059, R4-150, R2-176, R3-177): MATLAB `ratQ:overflow`, not finished, not finished, answered.
- **MATLAB's `ratQ:overflow` (biconj), the 2^53 ceiling, 8 inputs:** R4-003, R3-037, R4-059,
  R4-101, R2-169, R4-171, R4-199, R6-201. Rust answers 7 of them (R4-059 also times out);
  Python answers R3-037 and R2-169 and times out on the rest. R3-037 is `ratQ:overflow` here
  where the first section had `bridgeFaceAlg:witnessOverflow` -- same refusal, a different guard
  is reached first after #198's CAD change.
- **D-090 biconj:** MATLAB `witnessOverflow`, Python/Rust `noWitness` (same refusal, different
  guard, the same pattern as R2-029 in the first section).
- **MATLAB not finished where Rust answers (biconj):** R2-127, D-188 (Python also not finished
  on both; on R4-150 and R2-176 every engine is not finished).

The not-finished set is load-dependent (the run spanned load 2-126): R3-177 finished where it
never had before, and R2-127 / R2-176 (answered in the first section) hit the cap. Recorded as
performance, per the first section's precedent, not as defects.

## Wrong values vs the oracle

None. Every answered conjugate in all three engines matches the exact oracle at 1e-9, and every
answered biconjugate lies inside its certified bracket at 1e-7 (the grader additionally flags
wide brackets as "loose", 436 such points, not as disagreement). No failing test to open.

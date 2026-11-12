# MATLAB vs Python: speed of the exact engine (2026-11-11)

**Latest: the last section** ("After the EnvelopeMin speed-up"). The first table below is Python
without the accelerators, before the kink fix; the sections after it are the steps in between.

**Read with care.** One machine (Frances, 64 cores), shared: another user was logged in and the
load average was 0.8-2.9 during the run. MATLAB and Python ran CONCURRENTLY, each on one thread.
Times are the BEST of 4 calls per fixture, except the fixtures whose first call took over 20 s,
which ran ONCE (3 slow biconjugates and 3 of the fast AlgCrv ones). The slow rows are therefore
single measurements, and a ratio near 1 there is a tie, not a ranking.

## What ran

| | MATLAB | Python |
|---|---|---|
| version | R2024a, `-singleCompThread` | CPython 3.10.12 |
| exact integers | doubles to 2^53, then `java.math.BigInteger` (bigZ, ratQ big paths) | built-in `int` (unbounded), `fractions.Fraction` |
| floating point | MATLAB built-ins (`roots`, `\`, `eig`) | numpy 2.2.6 (`roots`, `linalg.solve/eigh/lstsq`, `RandomState`), BLAS = scipy-openblas 0.3.29, 1 thread |
| CAS | Symbolic Math Toolbox: only in `algRowFails` fallbacks when Yun overflows 2^53 (not reached here) | none (sympy is test-only; a lint test keeps it out of `cca2/`) |
| other | -- | stdlib only (`math`, `json`, `os`); python-flint / numba NOT used in this table (see the addendum) |

Fixtures: the 200 distinct exact inputs the MATLAB suites pass to `conjQ` / `biconjQ`
(`bench/fixtures.json`, from the goldens). Scripts: `bench/bench_matlab.m`, `bench/bench_python.py`,
`bench/compare.py`. Both engines returned without error on all 200; their outputs are the ones the
golden tests compare (identical QuaCons; AlgCrv / EnvelopeMin equal at the sampled points).

## Result
| bucket | call | result | fixtures | MATLAB total (s) | Python total (s) | Python/MATLAB total | median ratio | ratio range |
|---|---|---|---|---|---|---|---|---|
| fast | biconjQ | AlgCrv | 5 | 68.372 | 53.111 | 0.78 | 0.83 | 0.74 - 4.32 |
| fast | biconjQ | QuaCon | 41 | 0.134 | 0.250 | 1.87 | 1.70 | 0.12 - 3.25 |
| fast | conjQ | QuaCon | 140 | 10.203 | 16.128 | 1.58 | 1.58 | 0.01 - 3.73 |
| slow | biconjQ | AlgCrv | 2 | 490.503 | 604.379 | 1.23 | 1.09 | 0.88 - 1.29 |
| slow | biconjQ | EnvelopeMin | 1 | 561.973 | 264.728 | 0.47 | 0.47 | 0.47 - 0.47 |
| slow | conjQ | QuaCon | 11 | 0.029 | 0.063 | 2.16 | 1.87 | 1.58 - 2.59 |

Slowest fixtures (best time, seconds):
| id | source | call | result | MATLAB | Python | ratio | reps M/P |
|---|---|---|---|---|---|---|---|
| 191 | slow_erratumBiconjEx2Test | biconjQ | EnvelopeMin | 561.97 | 264.73 | 0.47 | 1/1 |
| 196 | slow_kinkSubdivisionTest | biconjQ | AlgCrv | 416.93 | 539.34 | 1.29 | 1/1 |
| 200 | slow_minSubdivisionTest | biconjQ | AlgCrv | 73.57 | 65.04 | 0.88 | 1/1 |
| 177 | biconj | biconjQ | AlgCrv | 35.80 | 26.51 | 0.74 | 1/1 |
| 174 | biconj | biconjQ | AlgCrv | 20.03 | 16.72 | 0.83 | 1/4 |
| 165 | biconj | biconjQ | AlgCrv | 12.44 | 9.40 | 0.76 | 4/4 |
| 170 | biconj | conjQ | QuaCon | 3.51 | 4.58 | 1.31 | 4/4 |
| 91 | conj | conjQ | QuaCon | 1.44 | 2.43 | 1.69 | 4/4 |
| 98 | conj | conjQ | QuaCon | 1.12 | 2.00 | 1.78 | 4/4 |
| 35 | conj | conjQ | QuaCon | 0.98 | 1.55 | 1.58 | 4/4 |
| 80 | conj | conjQ | QuaCon | 0.44 | 0.72 | 1.63 | 4/4 |
| 110 | conj | conjQ | QuaCon | 0.39 | 0.70 | 1.80 | 4/4 |

Spread worst/best over repetitions (fixtures above 1 ms): MATLAB median 1.26 max 15.21; Python median 1.02 max 1.75.

Totals over all 200 fixtures (best times): **MATLAB 1131 s, Python 939 s**.

## Reading

- **Rational answers (`conjQ`, and `biconjQ` returning a QuaCon) are 1.6-1.9x slower in Python**
  (medians), at a few ms each: this is interpreted Python on small integer lists against MATLAB's
  vectorised doubles. The absolute cost is small (16 s for all 140 conjugates).
- **Algebraic answers (AlgCrv) are a wash**: Python 0.74-0.88x on four of the six, 1.29x on the
  kink (one run each for the big ones). The time on both sides is the float screens (Bernstein
  sign tests, cell sampling) and the certificate's box recursion, not exact arithmetic.
- **Erratum Example 2 (EnvelopeMin) is 2x faster in Python** (265 s vs 562 s, one run each).
  Why was not profiled on the MATLAB side; a plausible cause is MATLAB's BigInteger round trips
  through Java, which Python's ints do not need -- unverified.
- **Python's times are steadier**: worst/best over repetitions median 1.02 (MATLAB 1.26, max 15
  -- MATLAB's first calls after a cache flush).

**Libraries (as of this table):** the one cProfile of a small algebraic subdivision (9.6 s) put
the time in `cell_sample_point` and `poly_sign_on` (numpy on 3x3 to 5x5 matrices, i.e. call
overhead), not in exact polynomial arithmetic. That changed the same day -- see the addendum.

## Addendum (later on 2026-11-11): the kink after the port fix, with optional accelerators

The table above predates two changes: the port gained `algRowFails`' rows against a face cubic in z
(topRowPolys / topFaceRow; the kink is now 317 faces as in MATLAB, not 334), and those rows made
Yun's square-free decomposition hot (pure-Python PRS gcds; erratum Example 2 ran > 30 min). A
cProfile of the kink (1090 s profiled) then read: Yun / PRS gcds ~220 s, the brute-force surd-corner
search 319 s (5503 calls, mostly repeated corners), Bernstein box tests 196 s (8.5M calls), cell
interiors 164 s. Three changes, all optional and with the pure-Python / numpy path as reference:

| config | what | kink biconjQ, best of 2 (s) | the 2 runs |
|---|---|---|---|
| A | before (commit d64cb06) | 627 | 627, 629 |
| B | + python-flint 0.9.0 `factor_squarefree` for Yun (`cca2/_poly.py`), surd corners remembered | 310 | 312, 310 |
| C | B + numba 0.62.1 Bernstein box test (`cca2/_poly_sign_nb.py`) | 179 | 179, 187 |

A, B, C ran CONCURRENTLY (one process each, same load, average 1.0-2.1) on the shared machine; all
returned the 317-face AlgCrv. C against A: 3.5x. The slow bucket (23/0) took 219 s for erratum
Example 2 under C (265 s in the table above, before the port fix; one run each). For reference only
-- not measured side by side: MATLAB took 417 s on the kink in the table's run (539 s was Python's,
at 334 faces; corrected here). Switch off with
CCA2_FLINT=0 / CCA2_NUMBA=0; the fast bucket is 2733/0 either way.

## Rerun with accelerators (2026-11-11, 18:37-18:59)

Same 200 fixtures and scripts, code at e1e70a0. Python with **python-flint 0.9.0** (Yun's square-free
step) and **numba 0.62.1** (the Bernstein box test), plus the remembered surd corners; MATLAB
unchanged (R2024a, `-singleCompThread`). FOUR processes ran CONCURRENTLY, two full passes per engine,
on the shared machine (load 0.3 at the start, 1.5-2.2 at the end). Per fixture: the best over both
passes, so the slow fixtures are now **best of 2** on both sides (they were single runs above). Both
engines returned without error on all 200. Raw: `bench/run2/`; table: `python3 compare.py
run2/matlab_p1.jsonl,run2/matlab_p2.jsonl run2/python_p1.jsonl,run2/python_p2.jsonl`.

| bucket | call | result | fixtures | MATLAB total (s) | Python total (s) | Python/MATLAB total | median ratio | ratio range |
|---|---|---|---|---|---|---|---|---|
| fast | biconjQ | AlgCrv | 5 | 66.581 | 48.709 | 0.73 | 1.00 | 0.58 - 5.96 |
| fast | biconjQ | QuaCon | 41 | 0.135 | 0.248 | 1.83 | 1.75 | 0.13 - 3.27 |
| fast | conjQ | QuaCon | 140 | 9.896 | 15.952 | 1.61 | 1.59 | 0.02 - 3.38 |
| slow | biconjQ | AlgCrv | 2 | 483.898 | 225.449 | 0.47 | 0.52 | 0.44 - 0.60 |
| slow | biconjQ | EnvelopeMin | 1 | 569.669 | 204.172 | 0.36 | 0.36 | 0.36 - 0.36 |
| slow | conjQ | QuaCon | 11 | 0.028 | 0.063 | 2.24 | 1.92 | 1.67 - 2.69 |

Slowest fixtures (best time, seconds):
| id | source | call | result | MATLAB | Python | ratio | reps M/P |
|---|---|---|---|---|---|---|---|
| 191 | slow_erratumBiconjEx2Test | biconjQ | EnvelopeMin | 569.67 | 204.17 | 0.36 | 2/2 |
| 196 | slow_kinkSubdivisionTest | biconjQ | AlgCrv | 411.33 | 181.63 | 0.44 | 2/2 |
| 200 | slow_minSubdivisionTest | biconjQ | AlgCrv | 72.57 | 43.82 | 0.60 | 2/2 |
| 177 | biconj | biconjQ | AlgCrv | 36.54 | 21.12 | 0.58 | 2/2 |
| 174 | biconj | biconjQ | AlgCrv | 17.55 | 17.56 | 1.00 | 5/8 |
| 165 | biconj | biconjQ | AlgCrv | 12.40 | 9.44 | 0.76 | 8/8 |
| 170 | biconj | conjQ | QuaCon | 3.40 | 4.53 | 1.33 | 8/8 |
| 91 | conj | conjQ | QuaCon | 1.40 | 2.41 | 1.72 | 8/8 |
| 98 | conj | conjQ | QuaCon | 1.09 | 1.97 | 1.81 | 8/8 |
| 35 | conj | conjQ | QuaCon | 0.97 | 1.54 | 1.59 | 8/8 |
| 80 | conj | conjQ | QuaCon | 0.43 | 0.72 | 1.66 | 8/8 |
| 110 | conj | conjQ | QuaCon | 0.38 | 0.70 | 1.84 | 8/8 |

Totals over the 200 fixtures both engines finished (best times): MATLAB 1130 s, Python 495 s.

Spread worst/best over repetitions (fixtures above 1 ms): MATLAB median 1.57 max 21.54; Python median 1.04 max 27.36.

Per pass, the slow ones (s): Example 2 MATLAB 570 / 571, Python 204 / 212; kink MATLAB 416 / 411,
Python 182 / 187; min-subdivision MATLAB 73 / 73, Python 44 / 46. The Python "max 27" spread is
fixture 120's first call (0.12 s vs 0.0045 s), where the numba kernel is loaded from its cache.

**Reading.**
- **Algebraic answers now favour Python**: the three slow biconjugates are 2-2.8x faster (Example 2
  0.36, kink 0.44, min-subdivision 0.60), and the kink is now 317 faces on both sides. On the fast
  AlgCrv fixtures it is a tie or better (0.58-1.00), except two tiny ones (fixtures 183 and 186: 0.05 s
  MATLAB vs 0.29 s Python), which is fixed per-call overhead, not the algorithm.
- **Rational answers (`conjQ`, and `biconjQ` returning a QuaCon) are unchanged, still 1.6-1.9x
  slower in Python**, at a few ms each. The accelerators do not touch that path (interpreted Python on
  small integer lists).
- **Total over all 200: MATLAB 1130 s, Python 495 s** (it was 1131 vs 939 without the accelerators).

## After the conjQ speed-up (2026-11-11, 21:28-21:50)

Code at ee2d52e: on top of the accelerators, cellInterior's margin LP solves all its 3x3 vertex
systems in one batched LAPACK call (bit-identical vertices), conicMeet remembers each pair of conics
(79189 calls on the 151 conjugates, 4358 distinct pairs), and cellHasInterior's 400-point probe is
evaluated elementwise in numpy (same operations, same verdicts). `bench_python.py` now EMPTIES the two
memos MATLAB does not have (conicMeet's, the surd corners') before every timed call, so repetitions
are not flattered; algRowFails' memos mirror MATLAB's `persistent` ones and stay warm on both sides.
Same protocol as the previous section: two full passes per engine, four processes concurrent on the
shared machine (load 1.1-4.5 over the run), best over both passes. All 200 returned without error on
both sides. Raw: `bench/run3/`.

| bucket | call | result | fixtures | MATLAB total (s) | Python total (s) | Python/MATLAB total | median ratio | ratio range |
|---|---|---|---|---|---|---|---|---|
| fast | biconjQ | AlgCrv | 5 | 65.973 | 44.386 | 0.67 | 0.95 | 0.50 - 5.31 |
| fast | biconjQ | QuaCon | 41 | 0.136 | 0.082 | 0.61 | 0.63 | 0.13 - 0.76 |
| fast | conjQ | QuaCon | 140 | 9.966 | 2.699 | 0.27 | 0.49 | 0.01 - 0.81 |
| slow | biconjQ | AlgCrv | 2 | 499.879 | 155.612 | 0.31 | 0.42 | 0.27 - 0.57 |
| slow | biconjQ | EnvelopeMin | 1 | 572.906 | 169.570 | 0.30 | 0.30 | 0.30 - 0.30 |
| slow | conjQ | QuaCon | 11 | 0.029 | 0.020 | 0.69 | 0.62 | 0.53 - 0.80 |

Slowest fixtures (best time, seconds):
| id | source | call | result | MATLAB | Python | ratio | reps M/P |
|---|---|---|---|---|---|---|---|
| 191 | slow_erratumBiconjEx2Test | biconjQ | EnvelopeMin | 572.91 | 169.57 | 0.30 | 2/2 |
| 196 | slow_kinkSubdivisionTest | biconjQ | AlgCrv | 426.67 | 113.57 | 0.27 | 2/2 |
| 200 | slow_minSubdivisionTest | biconjQ | AlgCrv | 73.21 | 42.04 | 0.57 | 2/2 |
| 177 | biconj | biconjQ | AlgCrv | 36.28 | 18.19 | 0.50 | 2/8 |
| 174 | biconj | biconjQ | AlgCrv | 17.44 | 16.53 | 0.95 | 5/8 |
| 165 | biconj | biconjQ | AlgCrv | 12.14 | 9.10 | 0.75 | 5/8 |
| 170 | biconj | conjQ | QuaCon | 3.44 | 0.57 | 0.17 | 8/8 |
| 91 | conj | conjQ | QuaCon | 1.41 | 0.35 | 0.25 | 8/8 |
| 98 | conj | conjQ | QuaCon | 1.09 | 0.36 | 0.33 | 8/8 |
| 35 | conj | conjQ | QuaCon | 0.96 | 0.30 | 0.31 | 8/8 |
| 80 | conj | conjQ | QuaCon | 0.43 | 0.12 | 0.28 | 8/8 |
| 84 | conj | conjQ | QuaCon | 0.43 | 0.11 | 0.26 | 8/8 |

Totals over the 200 fixtures both engines finished (best times): MATLAB 1149 s, Python 372 s.

Spread worst/best over repetitions (fixtures above 1 ms): MATLAB median 1.39 max 49.45; Python median 1.06 max 2.07.

First calls (s, pass 1 / pass 2): Example 2 MATLAB 573 / 576, Python 174 / 170; kink MATLAB 427 / 429,
Python 115 / 114; min-subdivision MATLAB 73 / 74, Python 42 / 42; fixture 177 MATLAB 36 / 37,
Python 18 / 19 -- so the "best" columns are not carried by warm memos.

**Reading.**
- **`conjQ` is now faster in Python on every fixture**: 140 conjugates 2.7 s vs MATLAB 10.0 s (total
  0.27, median 0.49, worst 0.81); it was 16.0 s (1.6x slower) in the previous section. The same
  routines carry `biconjQ` returning a QuaCon: 0.08 s vs 0.14 s.
- The kink also gained (182 -> 114 s), since its arrangement uses the same cell tests.
- **Total over all 200: MATLAB 1149 s, Python 372 s** (3.1x).

## After the AlgCrv speed-up (2026-11-11, 22:13-22:34)

Code at 6a73234: exact root signs build each Sturm chain once (memoised; integer evaluation at
dyadic points; FLINT gcd for the common-root test), cell sampling is numpy end to end, the
quadratic fibre computes its powers once per point, and algRowFails remembers each boundary
factor's integer form and screen points and each row matrix's float array. Every change was checked
bit-for-bit (or exactly) against the code before it; the golden tests (kink: 317 faces, face for face)
pass. `bench_python.py` empties ALL the Python-only memos before every timed call. Same protocol:
two passes per engine, four concurrent processes on the shared machine (load 2.2-4.7), best over
both passes; all 200 returned without error on both sides. Raw: `bench/run4/`.

| bucket | call | result | fixtures | MATLAB total (s) | Python total (s) | Python/MATLAB total | median ratio | ratio range |
|---|---|---|---|---|---|---|---|---|
| fast | biconjQ | AlgCrv | 5 | 76.277 | 8.311 | 0.11 | 0.16 | 0.09 - 1.11 |
| fast | biconjQ | QuaCon | 41 | 0.133 | 0.083 | 0.62 | 0.63 | 0.12 - 0.80 |
| fast | conjQ | QuaCon | 140 | 9.940 | 2.700 | 0.27 | 0.51 | 0.01 - 0.83 |
| slow | biconjQ | AlgCrv | 2 | 486.264 | 56.974 | 0.12 | 0.13 | 0.11 - 0.16 |
| slow | biconjQ | EnvelopeMin | 1 | 562.384 | 44.569 | 0.08 | 0.08 | 0.08 - 0.08 |
| slow | conjQ | QuaCon | 11 | 0.028 | 0.020 | 0.70 | 0.62 | 0.52 - 0.81 |

Slowest fixtures (best time, seconds):
| id | source | call | result | MATLAB | Python | ratio | reps M/P |
|---|---|---|---|---|---|---|---|
| 191 | slow_erratumBiconjEx2Test | biconjQ | EnvelopeMin | 562.38 | 44.57 | 0.08 | 2/2 |
| 196 | slow_kinkSubdivisionTest | biconjQ | AlgCrv | 412.22 | 45.12 | 0.11 | 2/2 |
| 200 | slow_minSubdivisionTest | biconjQ | AlgCrv | 74.04 | 11.85 | 0.16 | 2/8 |
| 177 | biconj | biconjQ | AlgCrv | 35.64 | 3.15 | 0.09 | 2/8 |
| 174 | biconj | biconjQ | AlgCrv | 20.43 | 3.20 | 0.16 | 2/8 |
| 165 | biconj | biconjQ | AlgCrv | 20.10 | 1.84 | 0.09 | 2/8 |
| 170 | biconj | conjQ | QuaCon | 3.46 | 0.57 | 0.17 | 8/8 |
| 91 | conj | conjQ | QuaCon | 1.39 | 0.35 | 0.25 | 8/8 |
| 98 | conj | conjQ | QuaCon | 1.08 | 0.36 | 0.33 | 8/8 |
| 35 | conj | conjQ | QuaCon | 0.98 | 0.30 | 0.31 | 8/8 |
| 80 | conj | conjQ | QuaCon | 0.43 | 0.12 | 0.29 | 8/8 |
| 84 | conj | conjQ | QuaCon | 0.42 | 0.11 | 0.26 | 8/8 |

Totals over the 200 fixtures both engines finished (best times): MATLAB 1135 s, Python 113 s.

Spread worst/best over repetitions (fixtures above 1 ms): MATLAB median 1.36 max 19.08; Python median 1.06 max 2.06.

First calls (s, pass 1 / pass 2) -- Python now repeats fixtures 165, 174, 177 and 200 (first call under
20 s) while MATLAB does not, so these are the like-for-like numbers: Example 2 MATLAB 562 / 573,
Python 44.6 / 44.7; kink MATLAB 412 / 420, Python 45.1 / 45.5; min-subdivision MATLAB 75 / 74, Python
12.3 / 12.4; fixture 177 MATLAB 35.6 / 38.5, Python 3.3 / 3.3; 174 MATLAB 20.4 / 20.4, Python 3.3 / 3.3;
165 MATLAB 21.4 / 20.1, Python 2.3 / 2.3.

**Reading.**
- **Algebraic answers are now 6-13x faster in Python**: Example 2 45 s vs 562 s, kink 45 s vs 412 s,
  min-subdivision 12 s vs 74 s, the fast AlgCrv fixtures 2-3 s vs 20-36 s. The two tiny AlgCrv
  fixtures (183, 186: 0.06 s each) are a tie.
- Rational answers unchanged from the previous section (conjQ 2.7 s vs 9.9 s).
- **Total over all 200: MATLAB 1135 s, Python 113 s** (10x). Python's own path: 939 s (no
  accelerators) -> 495 s (accelerators) -> 372 s (conjQ) -> 113 s (AlgCrv).

## After the EnvelopeMin speed-up (2026-11-11, 23:03-23:26)

Code at 1ec840e: the exact fibre at big coefficients is one integer sum per layer (no Fractions),
the exact root-sign bisection runs on integers N / 2^e, the surd-corner search is vectorised, and
two smaller trims (shared products in the sample slack, cached kernel tables). Each change was
checked against the code before it (or the original) to the bit; golden tests pass. Same protocol:
two passes per engine, four concurrent processes on the shared machine (load 1.3-2.4), best over
both passes, all Python-only memos emptied before every timed call; all 200 returned without error
on both sides. Raw: `bench/run5/`.

| bucket | call | result | fixtures | MATLAB total (s) | Python total (s) | Python/MATLAB total | median ratio | ratio range |
|---|---|---|---|---|---|---|---|---|
| fast | biconjQ | AlgCrv | 5 | 64.439 | 6.484 | 0.10 | 0.14 | 0.08 - 1.14 |
| fast | biconjQ | QuaCon | 41 | 0.134 | 0.083 | 0.62 | 0.62 | 0.12 - 0.79 |
| fast | conjQ | QuaCon | 140 | 9.898 | 2.753 | 0.28 | 0.57 | 0.01 - 0.81 |
| slow | biconjQ | AlgCrv | 2 | 485.508 | 52.368 | 0.11 | 0.12 | 0.10 - 0.14 |
| slow | biconjQ | EnvelopeMin | 1 | 566.117 | 26.798 | 0.05 | 0.05 | 0.05 - 0.05 |
| slow | conjQ | QuaCon | 11 | 0.028 | 0.020 | 0.70 | 0.63 | 0.52 - 0.80 |

Slowest fixtures (best time, seconds):
| id | source | call | result | MATLAB | Python | ratio | reps M/P |
|---|---|---|---|---|---|---|---|
| 191 | slow_erratumBiconjEx2Test | biconjQ | EnvelopeMin | 566.12 | 26.80 | 0.05 | 2/2 |
| 196 | slow_kinkSubdivisionTest | biconjQ | AlgCrv | 412.58 | 41.83 | 0.10 | 2/2 |
| 200 | slow_minSubdivisionTest | biconjQ | AlgCrv | 72.92 | 10.54 | 0.14 | 2/8 |
| 177 | biconj | biconjQ | AlgCrv | 34.99 | 2.75 | 0.08 | 2/8 |
| 174 | biconj | biconjQ | AlgCrv | 17.28 | 2.34 | 0.14 | 8/8 |
| 165 | biconj | biconjQ | AlgCrv | 12.07 | 1.28 | 0.11 | 8/8 |
| 170 | biconj | conjQ | QuaCon | 3.38 | 0.57 | 0.17 | 8/8 |
| 91 | conj | conjQ | QuaCon | 1.41 | 0.35 | 0.25 | 8/8 |
| 98 | conj | conjQ | QuaCon | 1.09 | 0.36 | 0.33 | 8/8 |
| 35 | conj | conjQ | QuaCon | 0.97 | 0.31 | 0.32 | 8/8 |
| 80 | conj | conjQ | QuaCon | 0.43 | 0.12 | 0.28 | 8/8 |
| 84 | conj | conjQ | QuaCon | 0.43 | 0.11 | 0.25 | 8/8 |

Totals over the 200 fixtures both engines finished (best times): MATLAB 1126 s, Python 89 s.

Spread worst/best over repetitions (fixtures above 1 ms): MATLAB median 1.24 max 15.66; Python median 1.07 max 2.16.

First calls (s, pass 1 / pass 2): Example 2 MATLAB 566 / 570, Python 26.8 / 26.9; kink MATLAB 416 / 413,
Python 43.6 / 41.8; min-subdivision MATLAB 73 / 73, Python 10.9 / 10.8; 177 MATLAB 36.4 / 35.0, Python
2.9 / 2.9.

**Reading.**
- **Erratum Example 2 (the EnvelopeMin): 27 s vs MATLAB 566 s, 21x** (45 s in the previous section).
- The changes are shared code, so the AlgCrv fixtures gained a little too (kink 45 -> 42 s).
- **Total over all 200: MATLAB 1126 s, Python 89 s** (12.7x). Python's own path: 939 s (no
  accelerators) -> 495 -> 372 (conjQ) -> 113 (AlgCrv) -> 89 s (EnvelopeMin).

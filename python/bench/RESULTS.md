# MATLAB vs Python: speed of the exact engine (2026-11-11)

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
| other | -- | stdlib only (`math`, `json`, `os`); python-flint / cypari2 NOT used |

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

**Libraries that are unlikely to help, from one cProfile of one algebraic subdivision (9.6 s):**
python-flint (FLINT) / cypari2 (PARI) accelerate exact polynomial arithmetic, which is not where
that profile spends its time; the top entries are `cell_sample_point` and `poly_sign_on` (numpy on
3x3 to 5x5 matrices, i.e. call overhead). What would
help: batching the Bernstein tests into one numpy call per box level, or a compiled kernel
(the Rust/C core of DECISION.md).

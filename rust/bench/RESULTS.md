# Rust benchmark -- 2026-11-13 (after the code-review fixes)

Protocol of python/bench (bench_python.py): 200 fixtures, first call + 3 repetitions unless the first
took > 20 s, best time per fixture, two passes, Rust-only memos emptied before every timed call (the
same set bench_python.py empties). Rust: single thread, release build, FLINT (from python-flint 0.9.0) for
the square-free parts exactly where the Python port uses it. MATLAB and Python columns are
CCA2/python/bench/run5 (2026-11-11), NOT re-measured. Shared machine: Rust ran with load ~0.5-2, run5
with four concurrent processes (load 1.3-2.4) -- read the ratios as indicative, not exact.

fixtures timed by all three: 200; result kind differing from the fixture: []

| bucket | call | result | n | MATLAB (s) | Python (s) | Rust (s) | MATLAB/Python | Python/Rust | MATLAB/Rust |
|---|---|---|---|---|---|---|---|---|---|
| fast | biconjQ | AlgCrv | 5 | 64.439 | 6.484 | 1.007 | 9.9 | 6.4 | 64 |
| fast | biconjQ | QuaCon | 41 | 0.134 | 0.083 | 0.005 | 1.6 | 17.8 | 29 |
| fast | conjQ | QuaCon | 140 | 9.898 | 2.753 | 0.285 | 3.6 | 9.6 | 35 |
| slow | biconjQ | AlgCrv | 2 | 485.508 | 52.368 | 8.997 | 9.3 | 5.8 | 54 |
| slow | biconjQ | EnvelopeMin | 1 | 566.117 | 26.798 | 5.357 | 21.1 | 5.0 | 106 |
| slow | conjQ | QuaCon | 11 | 0.028 | 0.020 | 0.001 | 1.4 | 17.2 | 25 |
| **all** | | | 200 | **1126.1** | **88.5** | **15.7** | **12.7** | **5.7** | **72** |

Slowest fixtures (best time, seconds):
| id | source | call | result | MATLAB | Python | Rust | Python/Rust |
|---|---|---|---|---|---|---|---|
| 191 | slow_erratumBiconjEx2Test | biconjQ | EnvelopeMin | 566.12 | 26.80 | 5.357 | 5.0 |
| 196 | slow_kinkSubdivisionTest | biconjQ | AlgCrv | 412.58 | 41.83 | 7.559 | 5.5 |
| 200 | slow_minSubdivisionTest | biconjQ | AlgCrv | 72.92 | 10.54 | 1.439 | 7.3 |
| 177 | biconj | biconjQ | AlgCrv | 34.99 | 2.75 | 0.517 | 5.3 |
| 174 | biconj | biconjQ | AlgCrv | 17.28 | 2.34 | 0.321 | 7.3 |
| 165 | biconj | biconjQ | AlgCrv | 12.07 | 1.28 | 0.149 | 8.6 |
| 170 | biconj | conjQ | QuaCon | 3.38 | 0.57 | 0.050 | 11.4 |
| 91 | conj | conjQ | QuaCon | 1.41 | 0.35 | 0.049 | 7.1 |
| 98 | conj | conjQ | QuaCon | 1.09 | 0.36 | 0.027 | 13.2 |
| 35 | conj | conjQ | QuaCon | 0.97 | 0.31 | 0.039 | 7.8 |

Per-fixture Python/Rust ratio: median 14.1, range 5.0 - 38.7.

# Rust benchmark -- 2026-11-13

Protocol of python/bench (bench_python.py): 200 fixtures, first call + 3 repetitions unless the first
took > 20 s, best time per fixture, two passes, Rust-only memos emptied before every timed call (the
same set bench_python.py empties). Rust: single thread, release build, FLINT (from python-flint 0.9.0) for
the square-free parts exactly where the Python port uses it. MATLAB and Python columns are
CCA2/python/bench/run5 (2026-11-11), NOT re-measured. Shared machine: Rust ran with load ~0.5-2, run5
with four concurrent processes (load 1.3-2.4) -- read the ratios as indicative, not exact.

fixtures timed by all three: 200; result kind differing from the fixture: []

| bucket | call | result | n | MATLAB (s) | Python (s) | Rust (s) | MATLAB/Python | Python/Rust | MATLAB/Rust |
|---|---|---|---|---|---|---|---|---|---|
| fast | biconjQ | AlgCrv | 5 | 64.439 | 6.484 | 1.009 | 9.9 | 6.4 | 64 |
| fast | biconjQ | QuaCon | 41 | 0.134 | 0.083 | 0.005 | 1.6 | 18.3 | 30 |
| fast | conjQ | QuaCon | 140 | 9.898 | 2.753 | 0.264 | 3.6 | 10.4 | 38 |
| slow | biconjQ | AlgCrv | 2 | 485.508 | 52.368 | 9.439 | 9.3 | 5.5 | 51 |
| slow | biconjQ | EnvelopeMin | 1 | 566.117 | 26.798 | 5.432 | 21.1 | 4.9 | 104 |
| slow | conjQ | QuaCon | 11 | 0.028 | 0.020 | 0.001 | 1.4 | 18.0 | 26 |
| **all** | | | 200 | **1126.1** | **88.5** | **16.1** | **12.7** | **5.5** | **70** |

Slowest fixtures (best time, seconds):
| id | source | call | result | MATLAB | Python | Rust | Python/Rust |
|---|---|---|---|---|---|---|---|
| 191 | slow_erratumBiconjEx2Test | biconjQ | EnvelopeMin | 566.12 | 26.80 | 5.432 | 4.9 |
| 196 | slow_kinkSubdivisionTest | biconjQ | AlgCrv | 412.58 | 41.83 | 8.001 | 5.2 |
| 200 | slow_minSubdivisionTest | biconjQ | AlgCrv | 72.92 | 10.54 | 1.438 | 7.3 |
| 177 | biconj | biconjQ | AlgCrv | 34.99 | 2.75 | 0.520 | 5.3 |
| 174 | biconj | biconjQ | AlgCrv | 17.28 | 2.34 | 0.320 | 7.3 |
| 165 | biconj | biconjQ | AlgCrv | 12.07 | 1.28 | 0.149 | 8.6 |
| 170 | biconj | conjQ | QuaCon | 3.38 | 0.57 | 0.050 | 11.4 |
| 91 | conj | conjQ | QuaCon | 1.41 | 0.35 | 0.030 | 11.7 |
| 98 | conj | conjQ | QuaCon | 1.09 | 0.36 | 0.027 | 13.3 |
| 35 | conj | conjQ | QuaCon | 0.97 | 0.31 | 0.039 | 7.9 |

Per-fixture Python/Rust ratio: median 15.4, range 4.9 - 44.1.

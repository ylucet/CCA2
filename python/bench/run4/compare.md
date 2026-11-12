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

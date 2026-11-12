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

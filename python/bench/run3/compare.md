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

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

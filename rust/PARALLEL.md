# Parallel computing in CCA2: what is used, what could be, and how much it would gain

_Study of 2026-11-14 on Frances. Nothing in the engine was changed. The numbers are measurements of
this repository at `0b2f384`; the speed-ups of section 6 are ESTIMATES built from them._

## 0. Answer

1. **No parallel computing is used in the computation, in any of the three implementations.** The
   Rust engine is single-threaded, and so are the Python port and the MATLAB reference. The 72x of
   Rust over MATLAB (5.7x over Python) is a single-core figure.
2. Parallelism exists only around the engine: the MATLAB slow test bucket runs several MATLAB
   processes (`.claude/suite.sh -j`), the Python tests can use `pytest -n auto`, and `cargo test`
   runs the two fast parity tests on two threads.
3. **Several independent calls can already run in parallel, with no engine change** (each call owns
   its data and its memos). Measured on the 200 bench fixtures: 16.5 s on one core, 9.5 s on 2 to 32
   threads. That is only 1.75x (2.1x with processes), because one call, the kink fixture, is half
   of the total.
4. **Inside one call the potential is large.** A serial run was instrumented to measure how much of
   the work sits in loops whose iterations do not depend on each other: on the six heaviest
   fixtures the work is 9.5 to 27 times the longest dependent chain.
5. **Estimate for Frances (32 cores / 64 hardware threads):** the heaviest call about 5-9x faster on
   16 to 32 threads; the whole 200-fixture benchmark from 16.5 s to about 1-2 s, i.e. **8-20x, plan
   on 10x**, which would put Rust near 700x over MATLAB. The price is a mechanical but wide refactor
   (section 7).

## 1. Was parallel computing used?

| Part | Parallel? | Evidence |
|---|---|---|
| Rust engine (`rust/src`, 15 800 lines) | no | `Cargo.toml` has no threading crate; no `std::thread`, `rayon`, `Arc`, `Mutex`, atomics or SIMD in `src/`. Shared data are `Rc` (45 type uses in 12 modules) and the 15 memo tables are `thread_local!` + `RefCell`, i.e. written for one thread. `bench/RESULTS.md` says "single thread". |
| FLINT (the one C kernel, `src/flint.rs`) | no | Only `fmpz_mpoly_factor_squarefree` is called; `flint_set_num_threads` is never called, so FLINT stays on its default of one thread. |
| Python port (`python/cca2`) | no | No `multiprocessing`, `concurrent.futures`, `threading`, `joblib`; the two numba kernels are `@njit(cache=True)` without `parallel=True` or `prange`. |
| MATLAB reference (`MATLAB/`) | no | No `parfor`, `parpool`, `parfeval`, `spmd`, `gpuArray`. |
| MATLAB test driver | yes, processes | `.claude/suite.sh --slow -j N` runs N MATLAB processes (`xargs -P`), default 4. It parallelises the TESTS, not a computation. |
| Python tests | optional | `pytest -n auto` (pytest-xdist), again per test. |
| Rust tests | incidental | `cargo test` runs the test functions of `tests/parity.rs` on separate threads. This works because a call touches nothing outside its own thread. |
| Benchmarks (`bench.rs`, `bench_python.py`, `bench_matlab.m`) | no | One process, one call at a time. The 1126 s / 88.5 s / 15.7 s figures are all serial. |

## 2. The machine

Frances: 2 x Intel Xeon Gold 6426Y (Sapphire Rapids), 16 cores per socket, 32 cores, 64 hardware
threads, 2 NUMA nodes, 251 GB RAM, AVX-512. Linux 6.8, glibc 2.35. Two GPUs, used by the local
language models and of no use here (section 8). The machine is shared: other projects usually hold a
few cores, and single timings scatter by up to 20%.

## 3. Where the time goes (serial)

The benchmark is top-heavy: of the 16.5 s that the 200 fixtures take (one call each), fixture 196
(kink, 317 faces) is 7.9 s, 191 (erratum Example 2) 5.6 s, 200 (min-subdivision) 1.6 s, and 177, 174
and 165 together 1.2 s. The other 194 fixtures share 0.3 s. So "parallel CCA2" means "parallel
`biconjQ` on the algebraic-curve route"; `conjQ` and the QuaCon biconjugates are already in the
millisecond range and are not worth a thread.

Sampling profile of the three heaviest calls (`prof`, section 9), inclusive time:

| Stage (function) | 196 | 191 | 200 |
|---|---|---|---|
| envelope subdivision (`esa`) | 91% | 84% | 82% |
| - final loop over the cells: drop dominated faces (`drop_below` -> `alg_face_below` -> `alg_row_fails_box`) | 59% | 42% | 16% |
| - final loop over the cells: winner faces (`winner_faces`) | 16% | 8% | 16% |
| - building the arrangement, curve by curve (`live_cell`) | 15% | 19% | 16% |
| - overlay of several envelopes (`min_overlay`) | - | 11% | 29% |
| envelope assembly (`envelope_assemble` -> `bridge_witness`) | 8% | 16% | 18% |

## 4. What could run in parallel

A parallel region is a place where the code does the same thing to several items and no item needs
the result of another. They were found by reading the call tree of section 3, and each was checked
for what it writes: all of them only append results in order, plus the memo tables (pure caches).

**S1: loops over cells, pairs and tuples** (independent; results collected in order)

| Loop | Where |
|---|---|
| one subdivision per envelope when there are several | `envelope_subdivision_alg.rs:60` |
| arrangement: every cell is tested against the next curve (the curves stay sequential) | `:103` |
| final loop: every cell gets its faces | `:139` |
| winner faces of one cell | `:205` |
| overlay: merge of cell pairs, then one job per tuple | `:1199`, `:1230` |
| envelope assembly, one job per cell pair; the atoms' conjugates and pairs | `envelope.rs:489`, `biconj_atoms.rs:43`, `:162`, `:166` |

**S2 = S1 + two inner loops** (independent as well, one level further down)

| Loop | Where |
|---|---|
| the corners of a region, each certified by its own blow-up | `alg_row_fails.rs:175` |
| the ~1800 tie sample points tried as bridge witnesses (stops at 12 hits; nearly always visits all) | `bridge.rs:478` |

**S3 = S2 + the two quadtree recursions** (`a && b && c && d` over the four sub-boxes: parallel only
speculatively, with cancellation) -- `alg_face_below.rs:112`, `alg_row_fails.rs:349`.

Not parallel: the curve-by-curve order of the arrangement, and the second loop of `drop_below`
(`envelope_subdivision_alg.rs:620`), where each test depends on what the previous ones removed.

### Work and span, measured

`T1` is the serial time of the call. `T_inf` is the time the call would take with unlimited workers:
the longest chain of dependent work. It was measured, not guessed: `span_instrumentation.patch`
times every branch of every region in a normal serial run and replaces each region's time by that of
its longest branch. `T1 / T_inf` is the largest speed-up any number of cores could give.

| Fixture | T1 (s) | S1: T_inf | T1/T_inf | S2: T_inf | T1/T_inf | S3: T_inf | T1/T_inf |
|---|---|---|---|---|---|---|---|
| 196 kink | 8.0 | 3.28 | 2.4 | 0.71 | 11 | 0.50 | 16 |
| 191 Example 2 | 5.6 | 0.70 | 8.2 | 0.37 | 15 | 0.21 | 27 |
| 200 min-subdivision | 1.6-1.8 | 0.24 | 6.5 | 0.10 | 16 | 0.09 | 18 |
| 177 | 0.56 | 0.19 | 3.0 | 0.06 | 9.5 | 0.06 | 10 |
| 174 | 0.39-0.48 | 0.23 | 2.8 | 0.03 | 12 | 0.03 | 12 |
| 165 | 0.27 | 0.11 | 2.4 | 0.01 | 27 | 0.01 | 27 |

Two things to read from it:

- **The cell loops alone (S1) are not enough for the kink.** One of its 161 cells takes 2.9 s and
  another 1.1 s, so parallel cells stop at 2.4x. The time inside those cells is the corner loop:
  adding it (S2) brings the chain down to 0.7 s.
- **S3 adds little** until more than about 20 cores are free, and it is the only part that needs
  speculation.

### The memo tables

The serial run profits from 15 memo tables: a later cell finds what an earlier one computed. In
parallel that sharing is either lost (tables kept per thread) or kept (tables made shared). The
worst case was measured by emptying every table at the start of every branch (S2):

| Fixture | work without any sharing | T_inf without any sharing |
|---|---|---|
| 196 | 8.0 -> 9.9 s (+24%) | 0.71 -> 0.85 s |
| 191 | 5.6 -> 7.7 s (+38%) | 0.37 -> 0.54 s |
| 200 | 1.6 -> 1.9 s (+21%) | 0.10 -> 0.11 s |
| 165 | 0.27 -> 1.4 s (5x) | 0.010 -> 0.013 s |

So per-thread tables would cost a quarter to a third of the gain on the big fixtures and nearly all
of it on fixture 165. Shared tables are the right design.

## 5. What Frances delivers

Measured with `parbench` (N identical calls of fixture 196 at once; `bench/parallel/scaling.txt`):

| N at once | threads, one process (s) | one process each (s) | capacity, threads | capacity, processes |
|---|---|---|---|---|
| 1 | 9.4 | 7.9 | 0.84 | 1 |
| 8 | 9.95 | 8.1 | 6.4 | 7.8 |
| 16 | 10.9 | - | 11.6 | - |
| 32 | 11.7 | 10.3 | 21.6 | 24.6 |
| 64 | 19.7 | 17.4 | 25.7 | 29.0 |

"Capacity" is the throughput in units of one undisturbed core.

- **Thread tax, about 20%.** The same call takes 7.9 s in a process that never had a second thread
  and 9.4 s as soon as one exists, even an idle one. glibc's allocator drops its single-thread fast
  path, and this code allocates constantly (big integers, short vectors). Processes do not pay it.
  A thread-friendly allocator (mimalloc, jemalloc) is the usual cure; it could not be tried here,
  because downloading a new crate was not permitted in this session.
- **32 threads are worth about 22 cores, 64 threads about 26.** The loss is the tax plus the
  frequency drop of a loaded Xeon; hyper-threading adds about 20% on top of the 32 cores.

## 6. Estimated speed-ups

For P threads the time is bracketed by

- optimistic: `max(T1 / capacity(P), T_inf / 0.84)` -- perfect scheduling, memo tables shared;
- conservative: `T1_cold / capacity(P) + T_inf_cold / 0.84` -- the textbook bound for a greedy
  scheduler, with no memo sharing at all.

Both use the measured capacity of section 5, so the thread tax is inside.

**One call, loops S2:**

| Fixture | serial (s) | 8 threads | 16 threads | 32 threads | 64 threads |
|---|---|---|---|---|---|
| 196 kink | 7.9 | 3.1-6.4x | 4.3-9.4x | 5.4-9.4x | 5.7-9.4x |
| 191 Example 2 | 5.6 | 3.0-6.4x | 4.3-12x | 5.6-13x | 5.9-13x |
| 200 min-subdivision | 1.6 | 3.7-6.4x | 5.4-12x | 7.2-13x | 7.7-13x |
| 177 | 0.54 | 2.8-6.4x | 3.9-8.0x | 5.0-8.0x | 5.3-8.0x |

With S1 only, the kink gets 1.4-2.0x whatever the thread count; Example 2 gets 4-7x. With S3 on top
of S2 the 32-thread figures become 6.9-13x (kink) and 8.5-22x (Example 2).

**The whole benchmark (200 fixtures, 16.5 s on one core):**

| Approach | 8 threads | 16 threads | 32 threads | 64 threads | Status |
|---|---|---|---|---|---|
| A. one call per thread, engine unchanged | 1.75x | 1.75x | 1.75x | - | **measured** (9.5 s) |
| A'. one call per process | 2.1x | 2.1x | 2.1x | 2.1x | measured for the longest call (7.9 s) |
| B. A + parallel loops S2 inside a call | 3.7-6.4x | 5.6-12x | 8.1-20x | 8.8-20x | estimate (0.8-2.0 s at 32) |
| C. B + speculative recursions S3 | 3.9-6.4x | 6.3-12x | 9.5-22x | 10-26x | estimate |

**What to plan on:** 10x for the benchmark and 5-6x for its heaviest single call, on 16 to 32
threads. Past 32 threads the gain is small, because the remaining time is the dependent chain.

What the estimate does not contain: the cost of atomic reference counts and of locking the shared
memo tables (a few percent each, usually), a better allocator (which would push the figures up), and
the other jobs on Frances (which push them down: with 8 cores free, read the 8-thread column).

For a different workload the rule is simple. Many small calls: approach A scales with the core
count, about 22x on 32 threads. One large `biconjQ`: only B helps. The benchmark is an awkward mix,
half of it being a single call.

## 7. What approach B would take

In order of payoff; each step can be checked with the existing parity tests, since every region
appends its results in the serial order and the memo tables are pure caches, so the output must stay
identical.

1. **Make the data shareable**: `Rc` -> `Arc` (12 modules), `OnceCell` -> `OnceLock` in
   `envelope.rs`. Mechanical. Costs a few percent on one thread.
2. **Share the 15 memo tables** (a locked or sharded map in place of `thread_local!` + `RefCell`).
   Without it step 4 loses a quarter to a third of its gain (section 4).
3. **A work-stealing pool** (the `rayon` crate): the regions are nested (cells contain corner loops),
   so plain threads per loop would oversubscribe. New dependency.
4. **Parallel loops**, S1 first, then the corner loop and the witness loop. An error must be the one
   the serial order would have met first, so that refusals keep their identifier.
5. **Replace the allocator**, to remove the 20% thread tax.
6. Only then, and only if more than about 20 cores are really free: the speculative recursions (S3).

Risks: FLINT is called from several threads (it is thread-safe for distinct objects, which is how it
is used; to be confirmed under load); deep recursions need a larger stack on worker threads
(`parbench` gives them 64 MB); memory grows with the number of threads, but the calls are small.

Approach A needs nothing: `parbench` is the proof, and a library user can already call `conj_q` /
`biconj_q` from several threads today.

## 8. Not worth doing

- **GPU.** The work is exact big-integer arithmetic and data-dependent recursion on tiny data. There
  is no large uniform array to hand to a GPU.
- **SIMD.** The only dense floating-point kernel is the Bernstein sign test (`poly_sign_on`), at
  most about 15% of a heavy call, on matrices of a few entries.
- **Threads for `conjQ` and the QuaCon biconjugates.** They take milliseconds; the 20% thread tax
  would cost more than a thread gains. Parallelise only the algebraic-curve route.
- **Parallelising the Python or MATLAB engine.** Batch runs (approach A) would work with
  `multiprocessing` or separate MATLAB processes and hit the same 2x ceiling; anything finer is
  better spent in Rust.

## 9. How to reproduce

    cargo build --release
    # profile of one call (perf is not permitted on Frances: perf_event_paranoid = 4)
    ./target/release/prof ../python/bench/fixtures.json 196 /tmp/196.stacks
    python3 tools/prof_report.py /tmp/196.stacks --tree 9 --min 2.5
    # all fixtures, one call per thread; 0 = main thread only
    ./target/release/parbench ../python/bench/fixtures.json /tmp/out.jsonl 8
    ./target/release/parbench ../python/bench/fixtures.json /tmp/out.jsonl 32 copies 196
    # work / span (instruments the engine: apply, measure, revert)
    git apply bench/parallel/span_instrumentation.patch && cargo build --release
    CCA2_PT_CLASSES=49 ./target/release/prof ../python/bench/fixtures.json 196 /tmp/x   # "SPAN" line
    CCA2_PT_CLASSES=49 CCA2_PT_COLD=49 ./target/release/prof ../python/bench/fixtures.json 196 /tmp/x
    git apply -R bench/parallel/span_instrumentation.patch && cargo build --release

Classes for `CCA2_PT_CLASSES` (a bit mask): 1 = S1, 16 = corner loop, 32 = witness loop (49 = S2),
2 and 4 = the two recursions (55 = S3), 8 and 64 = the two loops of `drop_below` (not counted above).
Raw data: `bench/parallel/span.tsv`, `bench/parallel/scaling.txt`.

Caveats on the numbers: shared machine, wall-clock times, minima of 2 or 3 repetitions. `T1 / T_inf`
is a ratio within one run and is stable to two digits across repetitions; the absolute seconds are
not.

## 10. Addendum (2026-11-14): the max / min reduction

Question raised: conj = (triangulate, envelope, conjugate each piece: independent) then a MAX
(a reduction); biconj likewise. Sections 4 and 6 did not treat the reduction as its own item.

**Where the reductions are in the code** (all serial left folds today):
`conj_q.rs:49-53` (over the pieces of f), `conj_q.rs:122` (over the triangles of one piece),
`conj_q.rs:877` (unbounded non-convex), and `min_overlay` (`envelope_subdivision_alg.rs:1181`, the
min over the envelopes of the atoms route). `max_q` itself is a double loop over (cell of g1, cell of
g2) -- independent -- followed by `assemble_qua_con_cells` (merge / dedupe, mostly serial).

**What the benchmark contains** (measured, instrumented serial run, FLINT linked):
- the number of envelopes in the `min_overlay` of the heavy fixtures is tiny: 191 has np = 3
  (cells 28/42/33, tuples 114 -> 271, ~0.65 s of the 5.6 s), 200 has np = 2 (14/14 cells, 40 tuples,
  0.01 s of the arrangement); 196 (the kink) has ONE envelope, no overlay at all;
- a balanced tree over 2-3 operands has the same depth as the left fold (2 levels, or 1), so a
  tree reduction gains NOTHING on these fixtures. The reduction's real parallelism is inside each
  step -- the (a, k) merge loop and the per-tuple loop -- and that is already in S1/S2 of section 4
  (and so in the measured `T_inf`). The estimates of section 6 therefore stand for the benchmark.
- `conj_q`: 140 fixtures in 0.29 s (about 2 ms each, 1-3 pieces). The fold is negligible and thread
  spawn plus the ~20% thread tax cost more than any gain.

**Estimate for a workload where it matters** (n pieces or triangles, n >> 3; NOT in the benchmark).
Let c be the cost of one piece (envelope + conjugate), m(s) the cost of a max of two operands of
total size s, P the threads, capacity(P) as in section 5 (22 at 32 threads, 26 at 64):
- serial: T1 = n c + sum_k m(size of the running max) -- the running operand grows, so the fold is
  the quadratic-looking part;
- map phase: n c / min(P, n) -- near-ideal, about 20x at 32 threads when n >= 64;
- reduction as a tree: log2 n levels, level l has n / 2^l maxes. Level 1 is wide (n/2 maxes, parallel
  across them), the last level is ONE max of two big operands, whose cost is only reduced by the
  inner (i, j) cell loop of `max_q`. Span of the reduction = sum over levels of m(2^l * s0) / (inner
  parallelism), dominated by the last 1-2 levels;
- a tree is also less work than the left fold (balanced operand sizes), so part of the gain is
  algorithmic and shows even on one thread -- measure it with a serial tree first.
Rule of thumb: the map phase scales ~ capacity(P) (20x at 32 threads); the reduction scales only to
~ n/2 workers at its first level and to the inner-loop width at its last, so overall 8-15x at 32
threads for n of 64-128 if the max operands stay smallish, and less if the final max dominates
(then Amdahl: 1 / (f_last + (1 - f_last) / 20), e.g. 6x for f_last = 0.15).
This is an ESTIMATE, not a measurement: no fixture of the suite has more than ~3 operands, and the
needed fixtures (a triangulated non-convex piece with dozens of triangles) would have to be built.

**Order of the work if it were done**: (1) swap the left fold for a tree serially and check the
parity tests (is the result identical? `assemble_qua_con_cells` order may change face order, which
parity compares); (2) only then parallelise the map and the levels.

## 11. Implementation and measurements (2026-11-15)

The parallelism of sections 4-7 was implemented (rayon work stealing, `src/par.rs`; shared memo tables,
`src/memo.rs`). **Results are bit-identical to the serial engine**: `tests/parallel.rs` runs the
digest tool (`src/bin/digest.rs`, a hash of the COMPLETE result of each of the 200 fixtures, every face
and coefficient in order) at 1, 2, 8 and 32 threads and compares with the digests recorded from the
serial engine before any change (`bench/parallel/digest_serial.tsv`); the 114 + 176 parity traces and
the 20 slow golden records also pass (`cargo test --release -- --include-ignored`).

### What runs in parallel

| Region | Where |
|---|---|
| conjQ: one conjugate per piece / per triangle | `conj_q.rs` (`par::map`) |
| conjQ: the max over them (reduction) | `conj_q.rs` -- left fold by default, see below |
| max / assembly: liveness of every cell, shedding of implied rows | `cells.rs` (`assemble_qua_con_cells`, `compact_cells`) -- the inside of every `max_q` |
| biconjQ atoms route: conjugates of the atoms, envelope of each pair, pruning, triangle faces | `biconj_atoms.rs` |
| envelope assembly: one job per cell pair | `envelope.rs` |
| `esa`: one subdivision per envelope; every cell against the next curve; the final per-cell loop (`drop_below`, `winner_faces`) | `envelope_subdivision_alg.rs` |
| the min over several envelopes (`min_overlay`): the merge of tuples and the per-tuple resolution | `envelope_subdivision_alg.rs` |
| corners of a region (each its own blow-up); the bridge witness samples (in chunks, stopping rule applied in order) | `alg_row_fails.rs`, `bridge.rs` |
| quadtree recursions (S3) | implemented, **off** (`CCA2_SPECULATE`), measured slower |

Every helper returns results in input order and, when several items fail, the error the serial loop
would have met first. The 14 memo tables are sharded `RwLock` maps (pure caches, so sharing changes
the time, never the result); `Rc` became `Arc`.

Switches: `CCA2_THREADS=n` (default all 64 hardware threads), `CCA2_THREADS=1` or `CCA2_PAR=0` for the
serial reference, `CCA2_REDUCE=tree`, `CCA2_SPECULATE=1`. `--features mimalloc` for the binaries.

### The max / min reduction: the honest result

1. **Tree vs fold.** A balanced tree gives a different face partition than the left fold (conj
   fixture 35: 67 faces instead of 71), so it would break "identical to serial". It is therefore
   opt-in (`CCA2_REDUCE=tree`), and on the stress test below it was never faster than the fold (3.7 s
   vs 1.85 s at n = 8 serially equal, 1.2 s vs 0.64 s at 32 threads): the cost of `max_q` grows with the
   face counts of both operands, so the last levels of the tree dominate and the depth gain is lost.
2. **Where the time really is.** `conjstress` (star polygons, `f = (x^2+y^2)/2`, 2n vertices, 2n-2
   triangles; `src/bin/conjstress.rs`): at n = 8 the 14 triangle conjugates take 6 ms and the max
   reduction 1.85 s (99.7%). The reduction is not a small tail, it is the whole cost.
3. **What worked: parallelism INSIDE each `max_q`** (liveness and implied-row loops of
   `assemble_qua_con_cells` / `compact_cells`), which keeps the fold and the output bit-identical:

| star n | triangles | faces | 1 thread (s) | 32 threads (s) | speed-up |
|---|---|---|---|---|---|
| 5 | 8 | 325 | 0.336 | 0.125 | 2.7x |
| 6 | 10 | 731 | 0.994 | 0.358 | 2.8x |
| 8 | 14 | 933 | 1.829 | 0.640 | 2.9x |

   (min of 3, shared machine.) The remaining serial part is `merge_adjacent_cells` and the fold's
   own dependency chain.
4. **A pre-existing defect found by this test, NOT caused by the parallel work** (it reproduces with
   `CCA2_PAR=0`, whose output is digest-identical to the pre-parallel engine): `conjQ` of the star
   polygon with n = 6 returns wrong values at 111 of 625 sample points (fold) -- the face picked at
   some points is far below the true maximum -- and for n = 8 the fold is wrong while the tree is
   right. n = 4 and 5 agree with the definition to 1e-16. The check is in `conjstress` (the definition:
   interior point `x = s` if inside the polygon, else the best boundary point). It needs its own
   investigation (suspect: `max_q` / cell assembly with many faces).

### Benchmark on Frances (200 fixtures, same protocol as bench/RESULTS.md)

Method: `tools/par_sweep.sh` -- 3 sweeps; in each, every thread count runs the whole benchmark (first
call + 3 repetitions, best per fixture, memos emptied before every call); the table takes the
**minimum per fixture over the sweeps** and sums. SHARED MACHINE: other jobs kept the load average at
3-8 during the runs, so read the figures as conservative; raw data `bench/parallel/runs/`, summary
`bench/parallel/sweep_summary.txt`. FLINT linked. Kind, face count and refusal identifier of all 200
fixtures are equal in all 21 runs.

| threads | all 200 fixtures (s) | speed-up | #196 kink | #191 Example 2 | #200 | #177 | #174 | #165 |
|---|---|---|---|---|---|---|---|---|
| 1 (serial, `CCA2_PAR=0`) | 14.81 | 1.0x | 6.79 | 5.21 | 1.51 | 0.54 | 0.33 | 0.16 |
| 2 | 8.42 | 1.8x | 3.96 | 2.74 | 0.79 | 0.30 | 0.24 | 0.11 |
| 4 | 4.80 | 3.1x | 2.06 | 1.59 | 0.45 | 0.25 | 0.17 | 0.09 |
| 8 | 3.44 | 4.3x | 1.47 | 1.22 | 0.32 | 0.15 | 0.12 | 0.05 |
| 16 | 3.02 | **4.9x** | 1.38 | 0.94 | 0.31 | 0.14 | 0.10 | 0.03 |
| 32 | 3.08 | 4.8x | 1.34 | 0.96 | 0.37 | 0.15 | 0.10 | 0.03 |
| 64 | 3.32 | 4.5x | 1.37 | 1.01 | 0.38 | 0.17 | 0.12 | 0.05 |

Per heavy call: kink 5.1x, Example 2 5.4x, #200 4.9x (best thread counts). Relative to MATLAB
(1126 s) the whole benchmark is now about 370x (15.7 s -> ~3 s at 16 threads).

### Against the estimate of section 6

Section 6 predicted 5-9x for the heaviest call and 8-20x for the benchmark (plan on 10x). Measured:
5x and 4.9x. The heaviest calls match the lower end; the benchmark misses because (a) the kink and
Example 2 still have a ~1 s critical path inside a few cells (profiled: `drop_below` of cell 3
0.9 s, `winner_faces` of cell 22 1.0 s, cell 8 0.6 s) whose loops are dependent (`keep[]` in
`drop_below`) or have few items, (b) the quadtree speculation that would shorten them wastes more than
it gains (measured both in `box_rec` and in `alg_face_below`), (c) 64 threads are 32 cores plus
hyper-threading, and the machine is shared. The 20% allocator "thread tax" is only partly removed by
mimalloc (a few percent). What would still help: restructure `drop_below` so its dependent loops can
be speculated cheaply, and cache `row_polys` across the repeated `alg_row_fails_box` calls on one cell.

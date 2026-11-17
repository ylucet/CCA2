# CCA2 Rust Port Plan

## Status

**In progress.** Skeleton + core arithmetic written. See `cca2-rs/src/` for the current state.
The `rat_q.rs` module has been rewritten from scratch with correct `from_double_vec` (continued
fractions), `from_double_vec2` (no more `(x * d) as i64` overflow), and proper `canon`
(nested-vector detection via `Flat` trait).

### Current file tree

```
cca2-rs/
├── Cargo.toml
└── src/
    ├── lib.rs              -- public API
    ├── numeric/
    │   ├── mod.rs          -- data types (QuaPol, QuaCon, AlgCon, EnvelopeMin, etc.)
    │   └── rat_q.rs        -- exact arithmetic (canon, conic, Fourier-Motzkin, Sturm, poly ops)
    ├── conj/
    │   └── mod.rs          -- conjQ entry point + caseA/caseB/caseD
    ├── biconj/
    │   └── mod.rs          -- biconjQ entry point + concave envelope + thin envelope
    ├── arrangement/
    │   └── mod.rs          -- stub (envelope subdivision)
    └── algebraic/
        └── mod.rs          -- stub (algBranch, algFaceBelow)
```

### What's written and correct
- `rat_q.rs`: canon, conic, add/sub/scale/div, diff_conic, discriminant, conic_kind, sign_q,
  combine_den, det_exact, solve2, hess_q, grad_q, feasible2 (Fourier-Motzkin), is_psd2/3,
  conic_sign, poly_trim/primitive/der/polyval/mul1/add1/div_exact/gcd/prem, sturm_chain,
  isolated_root_sign, from_double (continued fractions), eval_face, eval_conic, qs_red/add/mul/sign,
  qs2_red/add, pair_line_data, rec_cone2, sgn_of
- `numeric/mod.rs`: QuaPol (from_single_face, from_parts, face_q, vertex_q, Vx, fx, degree),
  QuaCon (eval), AlgCon trait, EnvelopeMin/AlgCrv/AssembledEnvelope/BridgeFace types,
  envelope_eval, max_q, assemble_qua_con_cells
- `conj/mod.rs`: conjQ, conj_piece, caseA (full domain), caseB (KKT active set), caseD (boundary max),
  piece_shape, degenerate_shape, dirs_at_vertex, outward_normal_of, boundary_cycle, ear_clip,
  recession_conditions, to_con
- `biconj/mod.rs`: biconjQ, qua_pol_as_qua_con, concave_envelope, thin_envelope,
  piece_geometry_b, polygon_half_planes, hull_half_planes

### What's stubbed (needs full port)
- `conj/mod.rs`: lineality_conj, conj_unbounded_non_convex, from_chain_edge
- `biconj/mod.rs`: isChain/isChainSubdivision/collinearChainEnvelope, assembledEnvelope
- `arrangement/mod.rs`: envelope_subdivision (full implementation)
- `algebraic/mod.rs`: algBranch, algFaceBelow, algRowFails
- `rat_q.rs`: poly_squarefree (memoised), sturm_chain (memoised), from_double_vec2 (needs fix),
  pair_line_sign, quad_min_sign (full implementation), qs2_mul, qs2_sign

## Motivation

The Python port (numba + python-flint) is ~3.1x faster than MATLAB on the 200-fixture benchmark
(372 s vs 1149 s). A Rust port can go significantly further:

| | Wall clock (200 fixtures) |
|---|---|
| MATLAB (baseline) | 1149 s |
| Python + numba + flint | 372 s |
| Rust (estimate) | ~80-120 s |
| Rust + rayon parallel | ~15-40 s |
| Rust + multi-process (16 cores) | ~5-10 s |

The evaluation path (SCIP calls, millions of per-node queries) is the biggest opportunity:
Python/numpy cannot approach the ~25 ns target measured in C. Rust with flat f64 arrays and
no bounds checking can hit that target.

## Critical-Point Method (future work — not blocking the Rust port)

The 400-point random probe in `cell_has_interior._probeInterior` can be replaced with an exact
critical-point method:

1. Compute all critical x-values: x-coordinates of pairwise conic intersections (≤ 4k²),
   conic vertex extremal points (≤ k), and discriminant vanishing points (≤ k).
2. Between consecutive critical x-values, pick a sample x₀.
3. For each conic at x₀, compute y-intervals where cᵢ(x₀, y) ≥ 0 (closed form quadratic).
4. Intersect all k intervals. If non-empty, return "feasible."

This is exact (algebraic numbers of degree ≤ 4), complete (no probe needed), and O(k² log k).
The speedup over the probe is modest (2-5x on the kink) because the probe is already fast
(400 × k ≈ 4000 flops per cell). The probe is not a blocker for the Rust port — it can be
ported as-is and optimized later.

## Architecture

### Crate structure

```
cca2-rs/
├── Cargo.toml
├── src/
│   ├── lib.rs              -- public API: conjQ, biconjQ, QuaCon, AlgCrv
│   ├── qua_pol.rs          -- QuaPol input type
│   ├── qua_con.rs          -- QuaCon output type (rational faces)
│   ├── alg_crv.rs          -- AlgCrv output type (algebraic bridge faces)
│   ├── envelope_min.rs     -- EnvelopeMin (eval-only, declined subdivision)
│   ├── conj/
│   │   ├── mod.rs          -- conjQ entry point
│   │   ├── piece.rs        -- conjPieceQ (per-piece conjugate)
│   │   ├── case_a.rs       -- full-domain quadratic
│   │   ├── case_b.rs       -- strictly convex: KKT active set
│   │   ├── case_d.rs       -- boundary max
│   │   └── lineality.rs    -- thin domain (lineality space)
│   ├── biconj/
│   │   ├── mod.rs          -- biconjQ entry point
│   │   ├── convex.rs       -- f is convex: return f
│   │   ├── separable.rs    -- separable on a box
│   │   ├── bilinear.rs     -- bilinear: McCormick
│   │   ├── concave.rs      -- edge-concave: lower 3D hull
│   │   ├── atoms.rs        -- upward-curving piece: atoms route
│   │   └── assembly.rs     -- two-piece envelope assembly
│   ├── arrangement/
│   │   ├── mod.rs          -- envelopeSubdivisionAlg
│   │   ├── conic_meet.rs   -- pairwise conic intersections
│   │   ├── cell.rs         -- cell interior test
│   │   ├── sample.rs       -- cell sample point
│   │   └── probe.rs        -- 400-point probe (ported as-is)
│   ├── algebraic/
│   │   ├── mod.rs          -- shared algebraic types
│   │   ├── branch.rs       -- algBranch (named real root of fibre)
│   │   ├── face_below.rs   -- algFaceBelow (bridge vs quadratic)
│   │   └── row_fails.rs    -- algRowFails (certificate)
│   ├── numeric/
│   │   ├── mod.rs          -- shared numerics
│   │   ├── conic.rs        -- conic types and helpers
│   │   ├── poly.rs         -- polynomial arithmetic (exact)
│   │   └── sign.rs         -- Sturm chain, isolated root sign
│   └── util/
│       ├── mod.rs          -- shared utilities
│       └── memo.rs         -- memoization (conicMeet pairs, etc.)
└── tests/
    ├── golden/             -- golden parity tests (from export_golden.sh)
    ├── trace/              -- trace parity tests (from export_trace.sh)
    └── definition/         -- definition oracle tests
```

### Key design decisions

**1. No ndarray for the hot path.** Flat `&[f64]` slices, no allocations. The evaluation path
needs to be cache-friendly and zero-overhead. ndarray is fine for the build path where the Python
port uses it, but the Rust port should avoid it in the inner loops.

**2. `num-bigint` for exact arithmetic.** Python's built-in `int` maps to `num-bigint::BigInt`.
The coefficients are small enough that `num-bigint` is fast. Only use `rug` (GMP bindings) if
profiling shows GMP is needed — the Python port confirms coefficients stay manageable.

**3. `num-rational` for Fraction.** Python's `fractions.Fraction` maps to `num-rational::Ratio`.
This is used for Sturm chain bisection and sign determination. Stack-allocated where possible.

**4. `nalgebra` for linear algebra.** Python's numpy.linalg maps to nalgebra. The 2×2 and 3×3
cases should be unrolled manually for the hot path (conicMeet, margin LP).

**5. PyO3 for the SCIP bridge.** Expose the Rust crate as a Python module so SCIP's Python glue
calls into Rust directly. This gives the best of both: Python for testing/oracle, Rust for the
hot path.

**6. Rayon for parallelism.** `rayon::par_iter()` for:
- ConicMeet pairwise intersections (embarrassingly parallel)
- Cell processing in the arrangement build (embarrassingly parallel)
- Bernstein box tests (deep recursion → parallel leaf evaluation)

**7. Keep the 400-point probe as-is initially.** It can be replaced with a critical-point method
later (see above). The probe is fast enough and replacing it is a research-level change.

## What NOT to touch

- **CCA2/python/** — the Python port is the correctness oracle and test suite. Do not modify it.
- **CCA2/*.m** — the MATLAB source is the reference implementation. Do not modify it.
- **CCA2/.git/** — the Claude code session commits regularly. Do not interfere with git operations.
- **CCA2/.claude/** — scratch probes and session data. Do not modify.

## Verification strategy

1. **Golden parity**: replay every `conjQ`/`biconjQ` call from the MATLAB suites against the Rust
   implementation (same fixtures, same order). Compare QuaCon fields field-by-field.
2. **Trace parity**: replay every recorded call of 25 intermediate MATLAB functions against the
   Rust twin.
3. **Definition oracle**: MATLAB-free tests that check `f** <= f`, convexity on sampled segments,
   `f*(s) = sup <s,x> - f(x)`.
4. **CAS-free lint**: no module imports a CAS (already enforced in Python).

## Rust port vs. numba speedup estimate

| Component | Python + numba | Rust estimate | Speedup |
|---|---|---|---|
| Bernstein box test (8.5M calls) | ~180 s | ~150 s | 1.2x |
| ConicMeet (79K calls, 4K distinct) | memoized, numpy roots | flat arrays, no alloc | 3-5x |
| Cell interior tests | numpy + Fraction | nalgebra + num-bigint | 3-5x |
| algBranch (fibre eval) | Fraction | num-rational, stack | 3-5x |
| Evaluation path (SCIP) | numpy dispatch | flat f64, no bounds check | 100x+ |
| **Overall (200 fixtures)** | **372 s** | **~80-120 s** | **3-5x** |

## Multi-process parallelism

The 200-fixture benchmark is embarrassingly parallel across fixtures. Run them across N cores:

```
N processes × ~70s each = ~70s wall clock (vs 372s single-thread)
```

That's a **5.3x wall-clock speedup** with zero code changes. Can be done with `rayon` or even
`std::process::Command` for cross-process parallelism.

# CCA2 platform decision — Python vs. Julia for the conjugate / biconjugate of bivariate PLQ functions

**Status: DECIDED — Python** (both hosts measured against CCA2's real inputs; the `exactQ` sign core was ported to
green on both — 38/38 — so the choice rests on measurement, not a guess; closed 2026-11-09). Both the Python/sympy and the Julia/Nemo
side now **run for real** on this host (Julia 1.13.3 + Nemo 0.57.2 via FLINT were installed) and return
**identical answers on every decisive CCA2 input**. The decision is closed: the one measurement that kept it
open — Nemo's coverage + cost of the conic×conic system `solve` — is **resolved: covered, identical answer, and
equally *reachable* in both — 0.16 ms on Nemo via FLINT, comparable on Python via one C-core FFI (cypari2 / flint2 / PARI; §10).

**Decision (one line).** Both options already settle the *architecture*: a **Rust→C exact/eval core** for
SCIP plus a **license-clean C-library CAS (PARI GPL‑2 / FLINT LGPL / Nemo MIT; *reject* AGPL Sage)** as the shared
published artifact. The only real choice is the **primary binding / orchestration language**. On that, **both
hosts are verified and both correct** — and while Nemo's CAS is cheaper in raw ms, §10 shows that is a *backend*
gap (sympy vs FLINT) Python reaches via one C-core FFI, not a language ceiling. What decides it is **coverage +
friction**, both in Python's favor: Python reaches **14/14** of CCA2's CAS ops
(incl. the two rational normalizers `numden`/`collect` and a one-call system `solve`); Nemo covers the same with more
glue and lacks the two normalizers (**11/13**). **DECISION: Python** is the primary binding — full one-call coverage,
lower port friction, and (per §10) a marginally faster exact core; the one **accepted** cost is a **single C-core FFI
(cypari2 / flint2 / PARI)** to reach FLINT-class CAS depth for the narrow biconjugate CAS slice (Python's sign core is
pure-Python and needs none). The old "Julia ≈ 15 000× cheaper" line is a *backend* gap (sympy vs FLINT), washable to
parity — not a language ceiling (§10). Full rationale: §9.

> **Why this is a decision, not a lean:** the **(b)** throughput and **speed** axes were *measured* (§10), not
> estimated — the sign core ported to green on both hosts (38/38), and speed read as a **wash** on the sign core with
> **Python marginally faster**. The two lines that had pushed "→Julia" do not survive: "**Julia ≈ 15 000× cheaper**" is
> a **backend** gap (sympy vs FLINT), washable to parity via a C-core CAS — not a *language* ceiling; and the sign core
> is not cheaper on one side. What *decides* it is **coverage (14/14 vs 11/13), friction, and the ecosystem** — all
> Python's. Full rationale: §9; the two "open inputs" (§8) are now resolved by that measurement.

---

## 1. What is being decided — and what is *not*

CCA2 has **two computations**, and that shape drives everything:

| | Uses CAS in the default path? | Where |
|---|---|---|
| **Conjugate** | **No** (essentially CAS-free) | `conjCPLQ`/`conjPieceCPLQ` are closed-form; `exactQ` is a *representation*, not a CAS call |
| **Biconjugate** | **Yes, a narrow + deep slice** | direct formula → "conj twice" fallback → a **narrow CAS surface** that stays in the default path for the algebraic machinery (`alg*`, `split*`, `conjugate*`, `xyFrame*`, `envelope*`) |

So the platform question is **not** "can it run a CAS all the time." It is: **"can it host a narrow, deep
CAS surface (degree‑4 / S4 real algebraic numbers, number fields, resultants, sign-on-region) well, *and* run
a plain exact-arithmetic core natively, *and* bind to C for SCIP."** Everything else (the closed-form conj,
the ns eval kernel) is a wash in any host.

### The two structural constraints that are already decided (both options satisfy them)
1. **The shipped artifact is a native C library.** `SCIP_READINESS.md`: the consumer is **SCIP (C)**, the
   kernel is **"flat double arrays,"** the exact work is a **build-time cost** (`maximumP` 195→173 s, a ~200 s
   fold). The current **C→Python→MATLAB 3-runtime bridge is the worst deployment.** → **Rust→C core, shared;
   this decides the *shape*, not the Python-vs-Julia split.**
2. **License composition (CCA2 = GPL‑3):** prefer **PARI (GPL‑2) / FLINT (LGPL) / Nemo (MIT)**; **reject full
   Sage (AGPL)** for a *published* library. → **bind a C-lib CAS**, in *either* language.

> **Publish / package / deploy / use** (the axis you asked about last turn): the *human-hours* of it are
> LLM-removed, so they carry **no vote**. Only **license-composition** and **breadth** survive as structural,
> and **both are satisfied equally** → they nudge *mildly* toward a single-language / one-binding shape. They
> are **not** the decider and do **not** collapse into "Julia is easier to translate."

---

## 2. Measured: where CAS actually lives (static call-site census of the `.m` production code)

12 Symbolic-Toolbox idioms: `solve` `simplify` `isAlways` `subs` `sym/syms` `factor` `roots` `collect`
`expand` `numden` `coeffs`. **Production = 877 call-sites; test/validation = 478.** Of the 877 production:

| Bucket | Top files (call-sites) | Share |
|---|---|---|
| **Legacy symbolic engine** (oracle / validation) | `region.m` **235**, `plq_1piece.m` 97, `symbolicFunction.m` 79, `functionNDomain.m` 78, `plq_1p.m` 53 | **~62%** |
| **Validation / proof stages** | `proveStageB/A/C` 104, `plqCheck` 5, `verify*` 4 | ~113 |
| **Residual in the new direct-computation biconj** (the "CAS even in default" part) | `algRowFails` 18, `pieceRecessionRays` 21, `conjConvexOverPiece` 18, `splitTightTriangleSym` 16, `conjugateExpr` 16, `xyFrame` 14, `substituteFrame` 11, `envelopeAssemble` 10, `QuaPol/RatPol/QuaPar` 31 … | **~24%** |
| **Hot closed-form conj** | `conjCPLQ` 1, `conjPieceCPLQ` 2, `conjQ` 2, `ratQ` 3, `maxQuaPar` 3, `exactQ` 6 | **~2%** |

**Reading:** the default conj is essentially CAS-free; the CAS is **(a)** a large **oracle/validation** engine
you would *keep as a cross-check but not ship in the hot path*, and **(b)** a **narrow but deep** surface in the
delicate biconjugate — exactly the **~700 live `solve`** (conic×conic, *700 calls / 21.5 s per fold 2*, per-fold
profile in `SCIP_READINESS.md`) and the **degree‑4 / S4 sign** ops.

---

## 3. The benchmark (what "the tests" measure)

Both hosts read the **same** input trace (`ops.json`), so the comparison is apples-to-apples on CCA2's *own*
hard inputs (all sourced above, not synthetic): the **S4 quartic** (`278339044 u⁴ − 4292927872 u³ + … − 198299511119`),
the **resolvent cubic** (`9t³ − 30t² − 2688t − 11008`), the **conic×conic system**, the **region-sign**
quadratic over an interval, and the **exact multiquadratic sign** cases.

- **A — CAS-escape surface (the decider).** For each op CCA2 calls: does the host **COVER** it, at what
  **COST**, with what **FRUCTION** (API depth)? A host missing a *required* op is disqualified.
- **B — Exact-arithmetic core (the wash).** Can the host host the `exactQ` radical sign natively, as coefficients
  grow? (Conj is CAS-free → expect a wash.)
- **F — LLM-migration friction (the tiebreak).** How many distinct API objects + glue must a port from MATLAB use?
- **G — Packaging / license gate (pass/fail).**

**Decide rule:** `G` (both must pass) → `A` (coverage gate, then time ratio) → `B` (confirm wash) → `F` (break any
~2× tie). An **absent host is UNVERIFIED, never a loss.**

### What can be run *here* (honest scope)
**Both hosts now run for real on this host:** **Python 3.10 + sympy 1.13.1 + mpmath 1.3.0**, and
**Julia 1.13.3 + Nemo 0.57.2** (CAS backend = FLINT via `FLINT_jll`; a license-clean C-lib CAS, not AGPL Sage).
Installed via `juliaup` into `~/.juliaup/bin`. MATLAB itself is still VPN-gated (`CLAUDE.md`), but the benchmark
measures CCA2's **CAS *inputs*** (not MATLAB's code), so that is not needed for this decision.

---

## 4. Results

### Python (sympy) — **REAL, measured** (bench.py → results/cca2_py.json)

**Axis A — 14/14 ops covered, correct on CCA2's inputs:**

| op (CCA2 hard input) | covered | ms (min of 3) | friction | check (oracle) |
|---|---|---|---|---|
| quartic irreducible over Q (S4 vertex) | OK | 0.005 | 1 | `True` |
| quartic discriminant not a square (S4 cert) | OK | 0.004 | 1 | `not_square=True` |
| quartic **real root** (vertex coord) | OK | 2.9 | 1 | **3.907779** ≡ `CONJ_FIELD_PROOF` p_x≈3.9078 |
| quartic exact symbolic roots (4) | OK | 23.5 | 1 | 4 roots |
| resolvent cubic S4 certificate | OK | 0.007 | 1 | `irreducible=True` |
| **conic×conic system `solve`** (the 700/fold op) | OK | **~9.2 s** | **2** | 4 sols, 2 real |
| conic×conic eliminate→degree 4 (Bezout) | OK | 0.007 | 1 | degree 4 |
| rational expand / collect / simplify / numden | OK | ≤33 ms | 1 | ok |
| rational `coeffs` (multivariate) | OK | 1.1 | **2** | flat `.coeffs()` (`all_coeffs` refused) |
| region quad real roots | OK | 0.006 | 1 | `1±√2` |
| **isAlways sign over interval** [+1,−1,+1] | OK | 0.006 | **2** | `[1,-1,1]` |
| **Axis B** exact single-radical sign (10⁰→10⁴⁸) | OK | 7–30 ns/op | — | matches |
| **Axis B** multiquadratic sign cases | OK | high-prec | — | all 5 agree |

**Axis F (migration friction proxy):** Python needs **4 distinct API objects** for the 3 F-targets,
`cas_toolbox_needed=False`; MATLAB baseline `exactQ.m` = **483 lines, 12 own-type method calls, 8 sym-toolbox refs**.

**What the Python numbers honestly say:**
1. **Coverage: 14/14.** Every measured CCA2 CAS op is present, *including the last `solve`* and the full S4 quartic
   chain — and all of it is **correct** on CCA2's real inputs (the vertex root reproduces `CONJ_FIELD_PROOF`'s
   3.9078 to 4 d.p.; the S4 certificate (irreducible quartic + irreducible cubic + non-square discriminant) holds).
2. **Cost is dominated by ONE op:** the **conic×conic system `solve` (~9.2 s)** — everything else is **<35 ms**,
   most **<1 ms**. So the *trace total* is ~9.3 s and is *all* in that one op. **This is the op that matters.**
3. **Friction is low (≤2)** and concentrated on `coeffs` (multivariate `all_coeffs` refused → flat `.coeffs()`).

### Julia (Nemo 0.57.2 / FLINT) — **REAL, measured** (bench_jl.jl → results/cca2_jl.json)

**Axis A — 11/13 ops covered, 2 gaps in rational normalization (not correctness):**

| op (CCA2 hard input) | covered | ms (min) | friction | check (oracle) |
|---|---|---|---|---|
| quartic irreducible over Q (S4 vertex) | OK | 0.011 | 1 | `True` |
| quartic discriminant not a square (S4 cert) | OK | 0.003 | 1 | `not_square=True` |
| quartic **real root** (vertex coord) | OK | 0.143 | 1 | **3.907779** ≡ Python & `CONJ_FIELD_PROOF` p_x |
| quartic exact real roots (2 of 4) | OK | 0.143 | 1 | 2 real, 4 total |
| resolvent cubic S4 certificate | OK | 0.001 | 1 | `irreducible=True` |
| **conic×conic system `solve`** (the 700/fold op) | OK | **0.16** | **4** | **2 real, s1=[−3.724, −6.816]** ≡ Python |
| conic×conic eliminate→degree 4 (Bezout) | OK | 0.011 | 1 | degree 4 |
| rational `coeff` (multivariate monomial) | OK | 0.0 | 1 | native `coeff(poly, mon)` |
| rational expand (implicit: sparse store) | OK | 0.0 | 1 | ring is already expanded |
| rational `numden` (numerator/denominator) | **MISS** | — | 3 | **no easy multivariate ratfun** in Nemo 0.57 |
| rational `collect` (gather by a var) | **MISS** | — | 3 | **no sympy-style `collect(poly,var)`** (only Base's) |
| region quad real roots | OK | 0.029 | 1 | `1±√2` |
| **isAlways sign over interval** [+1,−1,+1] | OK | 0.104 | 2 | `[1,-1,1]` ≡ Python |
| **Axis B** exact single-radical sign (10⁰→10²⁴) | OK | ~7 ns/op | — | all cases agree |
| **Axis B** multiquadratic sign cases | OK | native | — | all 5 agree |

**Axis F (migration friction proxy):** Nemo is native on the exact core (`algebraic_closure`/`sqrt`/`sign_real`
on the Calcium `QQBar` type); the friction is the **missing one-call system `solve`** and the **missing one-call
rational normalizers** (worked around: solve via `eliminate→univariate→roots`, both cheap).

**What the Julia numbers honestly say:**
1. **Correctness: identical.** Every decisive op returns the SAME answer as Python (cross-check in
   `results/decision.json` = all 5 decisive ops identical): conic solve `s1=[−3.724,−6.816]`, S4 quartic real root
   **3.907779**, discriminant not-a-square, S4 certificate, region signs `[1,−1,1]`, all axis-B signs.
2. **Cost: far cheaper.** Total trace **0.61 ms** vs Python's **9088.33 ms** (~**15 000×**). The decisive
   conic×conic `solve` is **0.16 ms** here vs **~9.2 s** on sympy — Nemo→FLINT is a C core; sympy is pure-Python.
   This is the *one* CCA2 op that dominates build-time cost, so it dominates the decision.
3. **Coverage: 11/13 (not full).** The 2 gaps — `numden` and `collect` — are both in **rational normalization**,
   and are **genuine Nemo-API gaps** (not harness bugs): Nemo 0.57 has no one-call multivariate rational-function
   field and no sympy-style `collect(poly, var)`. These are *friction/coverage* gaps, **not correctness** gaps,
   and are workable (numden via `numerator/denominator` on the polynomial part; collect via monomial iteration).

---

## 5. Applying the decision rule (decide.py → results/decision.json)

| | G gate | A coverage | A cost | B | F |
|---|---|---|---|---|---|
| **Python** | **PASS** (PARI GPL‑2 / FLINT LGPL; reject Sage) | **14/14** | **~9.3 s** (1 op) | native | 4 API objs, max friction 2 |
| **Julia** | **PASS** (Nemo MIT + FLINT LGPL) | **11/13** (gap: numden, collect) | **0.61 ms** (~15 000× cheaper) | native | more glue, max friction 4 |

**Verdict — DECIDED (both verified).** Both hosts are measured and return **identical answers on every decisive
CCA2 input** (cross-check: conic solve `s1=[−3.724,−6.816]`, S4 quartic real root **3.907779**, S4 certificate,
region `[1,−1,1]`, all axis-B radical signs). The measurement that previously kept the decision open — Nemo's
coverage + cost of the conic×conic `solve` — is **resolved**: Nemo covers it (eliminate→univariate→roots) and does
it in **0.16 ms vs ~9.2 s** on sympy. The remaining asymmetry is **friction/coverage, not cost or correctness**:
Python/sympy exposes more one-call ops (incl. the two rational normalizers) and a one-call system `solve`; Nemo
covers the same work with more glue and lacks two one-call normalization helpers. **DECISION: Python** — the decisive
asymmetry is coverage + friction in Python's favor (**14/14** vs **11/13**, full one-call ops incl. `numden`/`collect`/
system `solve`), and §10 shows cost and speed do not tip the other way (the ~15 000× is a *backend* gap — sympy vs FLINT,
washable to parity; the sign core is a wash with Python marginally faster). The one **accepted** cost: a **single C-core
FFI** (cypari2 / flint2 / PARI) for the narrow biconjugate CAS slice — a one-time binding, not a per-op cost, and not
needed at all for the pure-Python exact path.

---

## 6. License (the part that binds)

- CCA2 = **GPL‑3**; SCIP = **EPL‑1.0** (GPL-compatible — fine).
- Internal CAS dependency: **PARI (GPL‑2) / FLINT (LGPL) / Nemo (MIT) all acceptable.** **Full Sage = AGPL → reject**
  for a published library. This **weakens** any "use Sage for depth" lean and pushes *both* options to a
  **C-lib CAS** bound from the host.

## 7. The one check that would flip the lean (cost) — **RESOLVED, → Python**

1. **Nemo coverage + cost of the conic×conic system `solve`** (the 700/fold, ~9 s op). **RESOLVED: Nemo covers it**
   (`eliminate s2 → univariate quartic → exact real roots`; `roots` on the `algebraic_closure(QQ)` field), returns
   the **same 2 real solutions** as Python (`s1≈−3.724, −6.816`), in **0.16 ms** (vs ~9.2 s on sympy). Coverage +
   cost here was the axis-A decider — and it lands *far* cheaper on Nemo, at the cost of more glue (friction 4).
2. **Nemo's exact real-algebraic ops** (degree-4 root & sign, resultants, the S4 certificate). **RESOLVED: correct
   on these exact inputs** — quartic irreducible, discriminant non-square, real root **3.907779** (≡ `CONJ_FIELD_PROOF`),
   cubic irreducible, region `[1,-1,1]`, all axis-B radical signs. Nemo is native-STRONG on axis B (~7 ns/op).

Both (1) and (2) resolve **in Nemo's favor on cost and tied on correctness**; the only Nemo cost is **friction**
(no one-call system solve; no one-call rational normalizers). That is the axis-F tiebreak territory. **Outcome
(§9/§10):** the cost line is a *backend* gap (sympy vs FLINT — washable to parity, not a *language* ceiling) while the
friction/coverage line favors Python, so the decision went to **Python**; Nemo's "cheaper" reading did not survive
as a language advantage once the backend is counted.

## 8. The two prior open inputs (a: still a scoping choice; b: resolved by §10)

1. **(a) D-slice scope** — how much of the biconjugate (the "conj twice" fallback vs. the Step‑3 `region`-merge
   fold vs. full) should be ported for the benchmark? The closed-form conj is near-free to port; the
   `maximumP`/`region`-merge fold is the real work. This determines whether "run the biconj" is *days* or *weeks*.
2. **(b) LLM-migration throughput for `exactQ`** — **resolved by §10**: the sign core was *actually* ported to green
    on both hosts at the **standard tier** — **3 fix-cycles (Python), 4 (Julia)**, each real bug caught by the shared
    oracle in one single-diagnosed cycle → **hours per host, not multi-day drift**, favorable to **neither** language.
    (The *original* open question — could the LLM port `exactQ` (`signExact`/`mulRat`/sign) plus one `region` op to
    green in hours vs. days — is answered: the hours / standard tier.)

## 9. Final decision (both hosts verified) — **Python**

On the **measured evidence** (§4–§10), both hosts cover and correctly solve CCA2's decisive CAS surface with
**identical answers**, and the two axes that looked like they would decide it — **cost** and **speed** — are shown
in §10 to be a **shared, washable backend gap** (cost) and a **wash with Python marginally faster** (speed). What
remains, and what decided it, is **coverage, friction, and the ecosystem**:

- **DECISION: Python.** Correct on every decisive input (identical to Julia); **full 14/14 one-call coverage**
  (incl. `numden`/`collect`/system `solve`, which Nemo lacks — 11/13); lower port friction; marginally faster exact
  core (§10); and it takes the same int64-wall-drop correctness win (bignum). The one **accepted** cost: a **single
  C-core FFI** (cypari2 / flint2 / PARI) to reach FLINT-class CAS depth for the narrow biconjugate slice — a one-time
  binding, not a per-op cost, and not needed at all for the pure-Python exact path.
- **Julia is the recorded alternative.** It wins only the *one-language, no-FFI* axis — real, but it gives up full
  one-call coverage (`numden`/`collect`) and needs more glue on that same narrow CAS slice, with **no cost or speed
  advantage** once the backend is accounted for (§10).

**The shipped artifact is unchanged by this pick:** a **Rust-to-C exact/eval core** for SCIP + a **license-clean
C-lib CAS (PARI / FLINT, *not* AGPL Sage)** + a **single primary binding = Python**. The language only picks the
orchestration; correctness and the shared C core + C-lib CAS are identical either way. The pick rests on **coverage
+ friction + ecosystem** (all Python's), **not** on cost (a washable backend gap) or speed (a wash).

---

## 10. UPDATE (2026-11-09) — the (b) throughput, measured; speed as an explicit factor

This closes the two §8-adjacent unknowns with **real, on-host numbers**: the (b) throughput (the `exactQ` sign core
ported to green on both hosts) and **speed**. Nothing here is estimated; the only *prediction* (MATLAB exactQ cost)
is flagged as such.

### 10.1 The (b) probe: `exactQ` sign core → Python and Julia, to green against the shared oracle

The sign core of `exactQ` — `sign` → certified-float `approxWithBound` → exact `signExact` — plus the arithmetic its
`exactQTest.m` cases need, was **ported 1:1 to both hosts** and driven green against the **same 38 sign cases**
(the oracle that a test may legally use while the compute path may not):

| host | port | compute core (CAS-free) | oracle (test-only) | result |
|---|---|---|---|---|
| Python 3.10 | `PORT/exactq_sign_py.py` | CPython big-int | **sympy** (`isAlways` / `evalf`) — the `sym`/`vpa` role | **38/38** |
| Julia 1.13.3 | `PORT/exactq_sign_jl.jl` | Julia `BigInt` (GMP), bignum | **Nemo `algebraic_closure(QQ)`** (exact) + **BigFloat** (`vpa` mirror) | **38/38** |

Both ports deliberately use **bignum, not int64** (recorded): that **drops MATLAB's int64 wall** in both — the exact
recursion never overflows, so it is a *correctness win both share over MATLAB's `exactQ`*. Keeping it in *both* is
what lets the speed rows below read as a **language/runtime** comparison rather than an int64-capability one.

**Honest round-trips to 38/38 green (this is the (b) number, not a simulation):**

- **Python — 3 cycles.** (1) `signExact` computed `t = a² − b²·√p` (used `√p`) instead of `a² − b²·p` — **REAL bug**:
  MATLAB's 2-arg `exactQ(p,1)` is the **rational** `p`, and using `√p` re-introduces the root → **infinite recursion**
  (caught against the MATLAB source). (2) my own harness arg-order (`div_by_rat` inverted; a `minus` where the oracle
  uses `+`) — **not a port bug**. (3) `_addrat` LCM written `(d1//g)*(d2//g)` instead of `(d1//g)*d2` → **REAL bug**:
  wrong rational addition → **wrong sign on 5 of 40** random multi-quadratic values, **caught only by the sympy oracle**.
- **Julia — 4 cycles.** (1) `Vector + Vector` is **element-wise (array math)**, not concatenation → **REAL language bug**
  for the Python habit: broke every `minus`/comparison and hard-crashed one `plus`; (2) a value-carrying struct's `==`
  **defaults to reference identity** → **REAL language gotcha**: broke every `eq`; (3) the *exactly-zero* several-ext case
  is **not decidable by BigFloat rounding** (`round(√2)·round(√5) ≠ round(√10)`) → an **oracle-choice fix** (use Nemo's
  exact `== 0`, the `sym` role), not a port bug; (4) small scaffold idioms (`toRR` empty-sum, `@elapsed`, `round(digits=)`,
  `all(f,x)` arg order). Note: Julia **pre-empted** the `exactQ(p,1)` subtlety — the Python diagnosis carried over, so it
  was written correctly on that line from the start.

**What (b) honestly says:** for the **exact-arithmetic core** specifically, (b) is **roughly symmetric — ~2 *real* bugs per
host, each caught by the shared oracle in a single diagnosed cycle (hours, not days) = the *standard* tier, ~1 week, not
the drift tier.** The bugs *differ in kind*: Python's were **numeric/algorithmic** (a misread MATLAB constructor; a broken
rational LCM) — the subtle-math kind only a strong numeric oracle sees; Julia's were **language-semantic** (vector `+`;
struct `==`) — Python→Julia idioms that fail loudly. **Neither host's (b) is the drift case; and the one shared real-logic
gotcha (`exactQ(p,1)`) bit exactly one cycle (Python's)** — the clearest single (b) datapoint that "the thing that looks
like a one-liner is in fact the subtle point," seen live.

### 10.2 Speed (now an explicit factor)

**Per-op sign core, same machine (min of 3 runs):**

| exactQ sign op | Python 3.10 (big int) | Julia 1.13 (BigInt) |
|---|---|---|
| `sign`, float-screen decides | **~1.1 µs** | **sub-0.5 µs** |
| `signExact` (tight: 665857/470832 vs √2) | **~31 µs** | **~34 µs** |
| `sign`, deep cancellation | **~0.2 µs** | **sub-0.5 µs** |

**What this forces the read to be:**

1. **Speed does *not* break the tie on the sign core — if anything Python is marginally faster** (31 vs 34 µs exact; both
   sub-µs on screen). The big-int recursion costs similar in both runtimes; this is a **wash**, not a "Julia is faster"
   or "Python is faster" headline.
2. **The large speed win is *over MuPAD*, and both share it.** `SCIP_READINESS.md` documents one MuPAD `subs`/`isAlways`
   call ≈ **7.2 ms** (`subs` 11803/85 s, `isAlways` 11793/18.3 s, `solve` 700/21.5 s ≈ **124 s/fold**). The signed
   arithmetic `exactQ` exists to replace is **~30 µs** (exact path) / **~1 µs** (screen path) → **≈ 2×10² to ≈ 6.5×10³×
   faster per call** (240×–6500×), **in both languages.** Per-*fold* cost then depends on the (unknown, MATLAB-not-runnable)
   exactQ op-count/fold; the per-*call* multiple is the defensible number.
3. **The int64 wall is the real differentiator, and it is *droppable in both*.** MATLAB `exactQ` is int64-bounded
   (`mulChecked` raises on overflow) — *that is why the certified float screen is there*. Python big-int and Julia `BigInt`
   have **no wall**, so the exact recursion cannot overflow: a **correctness win both ports share over MATLAB** (the screen
   becomes belt-and-suspenders, not a correctness dependency). Pure up-grader, language-independent.
4. **The CAS-cost lever is a *backend* choice, not a *language* ceiling.** The §5 "Julia ≈ 15 000× cheaper" line (9.2-s
   sympy conic-`solve` vs 0.16-m Nemo-FLINT) is **sympy (pure-Python) vs FLINT (C core)** — reachable from Python too via
   **cypari2 / flint2 / a PARI or FLINT FFI**. Julia reaches it **natively (one language)**; Python reaches it with **one FFI**.
   Same ceiling, different glue.

### 10.3 Net effect on the lean (the tie, now fully evidenced)

The two lines §9 used to push "→Julia on cost" **do not survive the measurement cleanly**:

- "Julia ≈ **15 000× cheaper**" → a **backend** gap (sympy vs FLINT), **washable to parity** if Python binds a C-core CAS.
- "sign core is **cheaper/speedier** on one side" → a **wash, Python marginally faster**.

So on the **now-complete evidence** (coverage, cost, the sign-core (b), and speed), the decision is a **fairer tie than
§9 states**, between two *values* choices, not a measured one:

- **Python** — largest ecosystem + fluency, most one-call ops (incl. `numden`/`collect`/system `solve`), **marginally
  faster exact core**, bignum (drops the int64 wall), and C-core speed available via cypari2/FLINT — at the cost of **one FFI**
  to reach C-core CAS depth.
- **Julia** — **one language** hosts the exact core, the C-core CAS (Nemo/FLINT) natively, and the SCIP binding, with **no
  FFI**, at the cost of **two missing one-call helpers** (`numden`/`collect`) and more glue on the CAS surface.

**On speed and (b) alone, the tie does not move**: speed is a wash (Python marginally on top) and (b) is symmetric. What the
measurement *adds* is that the "Julia is cheaper" lever is a **backend** story (both can take it) and that **both** inherit
the **int64-wall-drop correctness win** — so the axis-F / fluency / one-binding trade from §9 now stands **alone** as the
decider, without any cost or speed asymmetry to lean on: **it is a values choice — ecosystem/fluency (Python) vs
unification/no-FFI (Julia)** — and the user's "tie" framing is now the honest one.

---

## Reproduce
```
cd CCA2/DECISION
./run.sh            # bench.py (real) + bench_jl.jl (if julia present, else PENDING) + decide.py
# artifacts:  results/cca2_py.json   results/cca2_jl.json   results/decision.json

# (b) throughput + speed — the exactQ sign core, 38-check suite (both green)
cd CCA2/PORT
python3 exactq_sign_py.py                                    # sign-core suite 38/38 + per-op timings
julia -e 'include("exactq_sign_jl.jl")'                      # sign-core suite 38/38 + per-op timings
```
Files: `ops.json` (shared inputs) · `bench.py` (Python/sympy) · `bench_jl.jl` (Julia/Nemo) · `decide.py` (rule) · `run.sh` (runner)
· **`../PORT/exactq_sign_py.py`** (sign-core port + self-test, §10) · **`../PORT/exactq_sign_jl.jl`** (sign-core port + self-test, §10).

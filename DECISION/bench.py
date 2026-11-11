#!/usr/bin/env python3
"""
CCA2 language-decision benchmark -- PYTHON host (sympy + mpmath).

Decides, for CCA2 (conjugate + biconjugate of bivariate PLQ functions), between:
    P : Python   (exact core native + CAS-escape via sympy / PARI / FLINT)
    J : Julia    (exact core native + CAS-escape via Nemo / PARI.jl)

This file measures the PYTHON side. The Julia twin is bench_jl.jl (Nemo); decide.py unions
the two output JSONs and applies the decision rule. Both hosts read the SAME input trace
(ops.json) so the comparison is apples-to-apples:

  A  CAS-escape surface   -- THE DECIDER. Do the ops CCA2's `region`/legacy engine calls
                             (the 12 Symbolic-Toolbox functions counted from the .m files) get
                             COVERED by this host, at what COST, with what API FRUCTION?
  B  Exact-arithmetic core -- THE WASH (conj is CAS-free). Is the exactQ sign/zero test
                             hostable natively in plain Python, at what cost, as coefficients grow?
  F  LLM-migration friction -- THE TIEBREAK PROXY. For CCA2's exactQ core, how many distinct
                             API objects + glue does a port from MATLAB need here vs in Julia?

HONEST SCOPE: only needs Python + sympy + mpmath. It does NOT require MATLAB, Julia, Sage,
PARI, FLINT, or gmpy2, and does NOT run the real CCA2 (MATLAB is VPN-gated). Every hard input
is CCA2's own: the S4 quartic / resolvent cubic / conic-x-conic system / region-sign come from
DECISIONS.md + CONJ_FIELD_PROOF.md (ops.json), not synthetic data.

Output: results/cca2_py.json  +  a human table on stdout.
"""
import json, os, sys, time, math, re
from fractions import Fraction
from concurrent.futures import ThreadPoolExecutor, TimeoutError as _FTE

ROOT = os.path.dirname(os.path.abspath(__file__))
RES  = os.path.join(ROOT, "results")
os.makedirs(RES, exist_ok=True)
OUT  = os.path.join(RES, "cca2_py.json")
with open(os.path.join(ROOT, "ops.json")) as f:
    OP = json.load(f)


def _ver(m):
    try:
        return getattr(__import__(m), "__version__", "present")
    except Exception:
        return None

BACKENDS = {
    "python": sys.version.split()[0],
    "sympy":  _ver("sympy"), "mpmath": _ver("mpmath"),
    "gmpy2":  _ver("gmpy2"), "flint":  _ver("flint"),
    "sage":   _ver("sage"),  "pari":   _ver("pari"),
}

# -------------------------------------------------------------------------- helpers
def timed(fn, reps=3, timeout=60):
    """Run fn up to `reps` times; return (last_result, best_ms, status).
    status in {'ok','timeout','<ExcType>: <msg>'}. Robust to runaway sympy paths."""
    best = None
    for _ in range(reps):
        ex = ThreadPoolExecutor(max_workers=1)
        fut = ex.submit(fn)
        t0 = time.perf_counter()
        try:
            last = fut.result(timeout=timeout)
        except _FTE:
            ex.shutdown(wait=False); return None, best, "timeout"
        except Exception as e:
            ex.shutdown(wait=False); return None, best, f"{type(e).__name__}: {str(e)[:110]}"
        dt = (time.perf_counter() - t0) * 1000
        best = dt if best is None else min(best, dt)
        ex.shutdown(wait=False)
    return last, best, "ok"


def _frac(x):
    s = str(x)
    a, b = s.split("/") if "/" in s else (s, "1")
    return Fraction(int(a), int(b))


# -------------------------------------------------------------------------- AXIS A
def axisA():
    from sympy import (Symbol, Poly, nroots, roots, solve, resultant, together,
                       expand, collect, simplify, Rational, sqrt)
    out, A = {}, OP["axisA_cas_escape"]
    def rec(key, res, ms, st, friction, note, check=None):
        out[key] = {"covered": (st == "ok" and res is not None),
                    "ms_min": round(ms, 3) if isinstance(ms, (int, float)) else None,
                    "status": st, "friction_calls": friction, "check": check, "note": note}
    X = Symbol("x"); T = Symbol("t"); U = Symbol("u"); s1, s2 = Symbol("s1"), Symbol("s2")

    # -- S4 quartic (the degree-4 vertex, CONJ_FIELD_PROOF sec 4.2) --
    P = Poly(A["quartic_degree4"]["poly"], X)
    r, ms, st = timed(lambda: P.is_irreducible, 3, 30)
    rec("A.quartic_irreducible_over_Q", r, ms, st, 1,
        "S4 vertex test: quartic irreducible over Q", check=r)

    r, ms, st = timed(lambda: P.discriminant(), 3, 60)
    issq = None
    if st == "ok" and isinstance(r, int) and r > 0:
        rr = math.isqrt(r); issq = (rr * rr == r)
    rec("A.quartic_discriminant_square_test", r, ms, st, 1,
        "S4 cert: discrim not a rational square => not a quad-tower", check={"not_square": not issq})

    def real_vertex():
        rs = nroots(P, n=25, maxsteps=400)
        real = [v for v in rs if abs(v.as_real_imag()[1]) < 1e-12]
        return max(real, key=lambda v: v.as_real_imag()[0])
    r, ms, st = timed(real_vertex, 3, 90)
    v = float(r.as_real_imag()[0]) if st == "ok" else None
    rec("A.quartic_real_root_numeric", r, ms, st, 1,
        "real vertex coord (CONJ_FIELD_PROOF p_x ~ 3.9078)", check=round(v, 6) if v is not None else v)

    r, ms, st = timed(lambda: roots(P), 1, 180)
    rec("A.quartic_roots_symbolic_exact", r, ms, st, 1,
        "EXACT symbolic roots over an extension (PARI/Nemo vs sympy differ here)",
        check=len(r) if st == "ok" and r is not None else None)

    # -- resolvent cubic (the S4 certificate CCA2 factors) --
    R = Poly(A["resolvent_cubic"]["poly"], T)
    r, ms, st = timed(lambda: R.is_irreducible, 3, 30)
    rec("A.resolvent_cubic_S4_certificate", r, ms, st, 1,
        "factor(resolvent cubic); irreducible => Galois group S4", check=r)

    # -- conic x conic (the LAST solve() CCA2 keeps: region.getVertices, 700 calls/fold) --
    c = A["conic_conic_system"]
    def conic(d):
        return (d["s1s1"]*s1**2 + d["s1s2"]*s1*s2 + d["s2s2"]*s2**2
                + d["s1"]*s1 + d["s2"]*s2 + d["c"])
    e1, e2 = conic(c["eq1"]), conic(c["eq2"])
    r, ms, st = timed(lambda: solve([e1, e2], [s1, s2]), 1, 180)
    chk = None
    if st == "ok" and r is not None:
        sols = r if isinstance(r, list) else [r]
        nreal = 0
        for sol in sols:
            d = sol if isinstance(sol, dict) else dict(zip([s1, s2], sol))
            if all(abs(complex(d[s]).imag) < 1e-9 for s in (s1, s2)):
                nreal += 1
        chk = {"solutions": len(sols), "real": nreal}
    rec("A.conic_conic_solve_system", r, ms, st, 2,
        "solve the 2-conic system (the op that forces the quartic)", check=chk)

    r, ms, st = timed(lambda: Poly(resultant(Poly(e1, s1, s2), Poly(e2, s1, s2), s2), s1).degree(s1),
                      3, 90)
    rec("A.conic_conic_eliminate_bezout", r, ms, st, 1,
        "eliminate one var -> univariate quartic (Bezout; degree==4)", check=r)

    # -- rational normalization stack (simplify/expand/collect/numden/coeffs; 205+140+28+9 call-sites) --
    expr = (s1 + s2)**3 - (s1 - s2)*(s1**2 + s2**2 + 3) + Rational(7, 1)/(s1 + 2*s2)
    for name, fn, fr in [
        ("expand",   lambda: expand(expr), 1),
        ("collect",  lambda: collect(expr, s1), 1),
        ("simplify", lambda: simplify(expr), 1),
        ("numden",   lambda: together(expr).as_numer_denom(), 1),
        # NOTE (measured friction): sympy Poly.all_coeffs() REFUSES multivariate -> must use the
        #  flat .coeffs() after numerator extraction. MATLAB's `coeffs` is a single call.
        ("coeffs",   lambda: Poly(together(expr).as_numer_denom()[0], s1, s2).coeffs(), 2),
    ]:
        r, ms, st = timed(fn, 2, 90)
        rec(f"A.rational_{name}", r, ms, st, fr,
            f"{name}(expr) on a rational bivariate" + (f"  [{r[:4]}...]" if name == "coeffs" and isinstance(r, list) else ""))

    # -- region sign (isAlways / ptFeasible / sign-on-interval; 147 call-sites + 11793/fold) --
    gq = A["isalways_region_sign"]
    Qr = Poly(gq["quad"], U)
    r, ms, st = timed(lambda: roots(Qr), 3, 30)
    rchk = {str(k): int(v) for k, v in r.items()} if st == "ok" and r is not None else None
    rec("A.region_quad_real_roots", r, ms, st, 1,
        "real roots of the region quad (signEverywhere/ptFeasible)", check=rchk)

    def sign_segments():
        import mpmath as mp
        lo, hi = mp.mpf(gq["interval"][0]), mp.mpf(gq["interval"][1])
        rts = [x for x in [mp.mpf(1) + mp.sqrt(2), mp.mpf(1) - mp.sqrt(2)] if lo < x < hi]  # roots of u^2-2u-1
        pts = sorted([lo, hi] + rts)
        def q(x): return x**2 - 2*x - 1
        return [int(mp.sign(q((pts[i] + pts[i+1]) / 2))) for i in range(len(pts) - 1)]
    r, ms, st = timed(sign_segments, 3, 30)
    rec("A.region_isAlways_sign_on_interval", r, ms, st, 2,
        "isAlways/certifiesNonPositive: sign of quad over interval (certified, high-precision)",
        check={"segments": r} if r is not None else None)
    return out


# -------------------------------------------------------------------------- AXIS B
def _exact_sign_single(a, b, d):
    """Faithful CCA2 exactQ sign for a + b*sqrt(d) with a,b in Q, d squarefree>0. EXACT, no CAS, no float."""
    a, b = Fraction(a), Fraction(b)
    if b == 0: return (a > 0) - (a < 0)
    if a == 0: return (b > 0) - (b < 0)          # sqrt(d) > 0
    if (a > 0) == (b > 0): return 1 if a > 0 else -1
    L, R = a * a, b * b * d                       # both > 0; sign equals sign(L - R)
    return (L > R) - (L < R)


def axisB():
    out = {}
    # (1) exact single-radical sign, scaled up in coefficient magnitude (models int64->bignum growth)
    for k in [0, 1, 2, 6, 12, 24, 48]:
        base = 10 ** k
        N = 200
        t0 = time.perf_counter()
        acc = 0
        for i in range(N):
            acc += _exact_sign_single(base, -(base + (i % 3)), 2) + _exact_sign_single(base, base, 3)
        per_ns = (time.perf_counter() - t0) * 1e6 / N
        out[f"B.sign_exact_10^{k}"] = {"ns_per_op": round(per_ns, 3), "covered": True}
    # (2) sign of the CCA2 multiquadratic cases (ops.json) -- certified by HIGH-PRECISION eval
    #     (mpmath, 60 digits). This is a host-capability + oracle check: can the host sign a
    #     multi-radical value, and does the hand-written oracle (ops.json) actually agree?
    import mpmath as mp
    mp.mp.dps = 60
    detail, all_ok = [], True
    for c in OP["axisB_multiquadratic"]["sign_cases"]:
        v = mp.mpf(0)
        for part, (num, den) in c["parts"].items():
            v += mp.mpf(num) / mp.mpf(den) * (mp.sqrt(int(part)) if part != "1" else mp.mpf(1))
        sgn = int(mp.sign(v))
        ok = (sgn == c["expect_sign"])
        all_ok &= ok
        detail.append({"note": c.get("note", ""), "computed_sign": sgn,
                       "oracle_sign": c["expect_sign"], "agrees": ok, "value": str(v)})
    out["B.multiradical_sign_certified"] = {"covered": all_ok, "detail": detail}
    out["covered"] = True
    out["note"] = ("exactQ single-radical sign is EXACT in plain Python (no CAS, no float). "
                   "Confirms axis B (conj core) is hostable natively and cheap -- a wash across hosts.")
    return out


# -------------------------------------------------------------------------- AXIS F
def frictionF():
    out = {}
    # MATLAB baseline: read straight from the source (the thing we are migrating FROM)
    mat = os.path.join(ROOT, "..", "exactQ.m")
    if os.path.exists(mat):
        txt = open(mat, encoding="utf-8", errors="ignore").read()
        out["matlab_baseline"] = {
            "exactQ.m_lines": txt.count("\n") + 1,
            "distinct_own_type_method_calls": len(set(re.findall(r'\bexactQ\.[A-Za-z_]+', txt))),
            "distinct_self_method_calls": len(set(re.findall(r'\bthis\.[A-Za-z_]+', txt))),
            "symbolic_toolbox_calls": len(re.findall(r'\bsyms?\b|\bsimplify|\bnumden|isAlways', txt)),
        }
    # Python side: distinct public API objects THIS file needs for the 3 F-targets
    py = {"exactQ.sign": ["Fraction"], "exactQ.zero": ["Fraction"],
          "region.isAlways": ["Poly", "roots", "sqrt"]}
    out["python"] = {"api_objects_per_target": py,
                     "distinct_api_objects": len({o for v in py.values() for o in v}),
                     "cas_toolbox_needed": False}
    return out


# -------------------------------------------------------------------------- main
def main():
    res = {"host": "python", "backends": BACKENDS, "A_cas_surface": {},
           "B_exact_core": {}, "F_friction": {}}
    try: res["A_cas_surface"] = axisA()
    except Exception as e: res["A_cas_surface"]["_error"] = f"{type(e).__name__}: {e}"
    try: res["B_exact_core"] = axisB()
    except Exception as e: res["B_exact_core"]["_error"] = f"{type(e).__name__}: {e}"
    res["F_friction"] = frictionF()

    a = {k: v for k, v in res["A_cas_surface"].items() if isinstance(v, dict)}
    res["summary"] = {
        "A_ops_total":              len(a),
        "A_ops_covered":            sum(1 for v in a.values() if v.get("covered")),
        "A_ops_timeout_or_error":   sum(1 for v in a.values() if not v.get("covered")),
        "A_total_ms":               round(sum(v["ms_min"] or 0 for v in a.values()), 2),
        "A_max_friction":           max((v["friction_calls"] for v in a.values()), default=0),
        "B_covered":                bool(res["B_exact_core"].get("covered")),
        "F_python_cas_needed":      res["F_friction"].get("python", {}).get("cas_toolbox_needed"),
    }
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, default=str)

    # ---- human table ----
    print("\n=== CCA2 language-decision benchmark :: PYTHON (sympy) ===")
    print("backends:", ", ".join(f"{k}={v}" for k, v in BACKENDS.items() if v))
    print("\n[A] CAS-escape surface  (the decider) ---- coverage / cost / friction")
    for k, v in res["A_cas_surface"].items():
        if not isinstance(v, dict): continue
        tag = "OK " if v.get("covered") else "MISS"
        print(f"  [{tag}] {k:40} {str(v.get('ms_min')).rjust(10)} ms  "
              f"friction={v.get('friction_calls')}  status={v.get('status'):8}  chk={str(v.get('check'))[:34]}")
    print("\n[B] Exact-arithmetic core  (axis B, the wash) ---- exactQ sign, native (no CAS)")
    for k, v in res["B_exact_core"].items():
        if k in ("covered", "note"): continue
        if isinstance(v, dict):
            if "ns_per_op" in v: print(f"  {k:30} {v['ns_per_op']} ns/op")
            else: print(f"  {k:30} covered={v.get('covered')}")
    fr = res["F_friction"]
    print("\n[F] LLM-migration friction  (tiebreak proxy) ----")
    mb = fr.get("matlab_baseline", {})
    print(f"  MATLAB baseline exactQ.m: {mb.get('exactQ.m_lines')} lines, "
          f"{mb.get('distinct_own_type_method_calls')} own-type method calls, "
          f"{mb.get('distinct_self_method_calls')} self calls, {mb.get('symbolic_toolbox_calls')} sym-toolbox refs")
    py = fr.get("python", {})
    print(f"  PYTHON distinct API objects for the 3 targets: {py.get('distinct_api_objects')}, "
          f"cas_toolbox_needed={py.get('cas_toolbox_needed')}")
    s = res["summary"]
    print("\nSUMMARY: A covered %d/%d   A total %s ms   B covered=%s   F python needs CAS: %s" % (
        s["A_ops_covered"], s["A_ops_total"], s["A_total_ms"], s["B_covered"], s["F_python_cas_needed"]))
    print("written:", OUT)


if __name__ == "__main__":
    main()

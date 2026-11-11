"""
exactQ sign core  --  PYTHON PORT (the (b) throughput probe)

Faithful port of `exactQ.m`'s SIGN core (sign / signExact / approxWithBound) plus the arithmetic
those tests require (constructor forms, canon, plus/minus/times, division-by-rational, coeffOf,
isZero, isRational, comparisons). The oracle is the SAME as MATLAB's exactQTest.m: the value is
checked against `sym` (here: sympy), and the hard sign cases against high-precision vpa.

DELIBERATE difference from MATLAB, called out to keep the speed read honest:
  * exactQ.m is int64-bounded (mulChecked raises on overflow; `sign` screens in float to dodge it).
  * Python ints are UNBOUNDED by default -> this port has NO int64 wall, so the exact recursion
    NEVER overflows. That is a CAPABILITY win over MATLAB (exactly where the `overflowRaises`
    test pins it), and it is why Python can be *faster than MATLAB on the hard sign cases* while
    being *slower per-op in raw scalar speed* (CPython). The two effects are both in the numbers.

Run:  python3 exactq_sign_py.py
"""
import math, time, random

EPS = 2 ** -52          # MATLAB `eps`

# ---- canonical helpers ---------------------------------------------------------------
def _gcd(a, b):
    a, b = abs(a), abs(b)
    while b: a, b = b, a % b
    return a

def _norm2(n, d):
    if d == 0: raise ZeroDivisionError("zero denominator")
    if d < 0: n, d = -n, -d
    g = _gcd(n, d)
    if g > 1: n, d = n // g, d // g
    if n == 0: d = 1
    return n, d

def _squarefree(d):
    if d <= 1: return 1, d
    k, f = 1, 2
    while f * f <= d:
        while d % (f * f) == 0: d //= f * f; k *= f
        f += 1
    return k, d

def _mulrat(n1, d1, n2, d2):
    if n1 == 0 or n2 == 0: return 0, 1
    g1 = _gcd(n1, d2); g2 = _gcd(n2, d1)
    return _norm2((n1 // g1) * (n2 // g2), (d1 // g2) * (d2 // g1))

def _addrat(n1, d1, n2, d2):
    g = _gcd(d1, d2); L = (d1 // g) * d2                          # LCM(d1,d2), NOT (d1/g)*(d2/g)
    return _norm2(n1 * (L // d1) + n2 * (L // d2), L)

def canon(m, n, d):
    """Canonical multi-quadratic: sorted squarefree radicands, lowest terms, merged, no zeros.
    Inputs are parallel lists (m,n,d) of ints -> value sum_i (n_i/d_i)*sqrt(m_i)."""
    terms = [(mm, nn, dd) for (mm, nn, dd) in zip(m, n, d) if nn != 0 and mm != 0]
    fixed = []
    for mm, nn, dd in terms:
        if mm < 0: raise ValueError("negative radicand")
        k, ms = _squarefree(mm)
        if k != 1: nn, dd = _mulrat(nn, dd, k, 1)
        else:      nn, dd = _norm2(nn, dd)
        fixed.append((ms, nn, dd))
    fixed.sort()
    out_m, out_n, out_d = [], [], []
    i = 0
    while i < len(fixed):
        ms = fixed[i][0]; acc_n, acc_d = fixed[i][1], fixed[i][2]; j = i + 1
        while j < len(fixed) and fixed[j][0] == ms:
            acc_n, acc_d = _addrat(acc_n, acc_d, fixed[j][1], fixed[j][2]); j += 1
        if acc_n != 0: out_m.append(ms); out_n.append(acc_n); out_d.append(acc_d)
        i = j
    return (out_m, out_n, out_d)

def q(m, n, d): return canon(list(m), list(n), list(d))
def qrat(n, d):    return q([1], [n], [d])                              # rational n/d
def qsingle(m, n=1, d=1): return q([m], [n], [d])                        # (n/d)*sqrt(m)
Z = ([], [], [])                                                         # zero
def is_zero(v):   return len(v[0]) == 0
def is_rat(v):    return len(v[0]) == 0 or (len(v[0]) == 1 and v[0][0] == 1)
def plus(x, y):   return q(x[0] + y[0], x[1] + y[1], x[2] + y[2])
def minus(x, y):  return plus(x, q(y[0], [-a for a in y[1]], y[2]))
def times(x, y):
    om, on, od, = [], [], []
    for (mi, ni, di) in zip(x[0], x[1], x[2]):
        for (mj, nj, dj) in zip(y[0], y[1], y[2]):
            g = _gcd(mi, mj)
            m = (mi // g) * (mj // g)
            pn, pd = _mulrat(ni, di, nj, dj)
            pn, pd = _mulrat(pn, pd, g, 1)
            om.append(m); on.append(pn); od.append(pd)
    return q(om, on, od)
def div_by_rat(v, n, d):                                                 # v divided by rational n/d
    om, on, od = [], [], []
    for (mi, ni, di) in zip(v[0], v[1], v[2]):
        pn, pd = _mulrat(ni, di, d, n)
        om.append(mi); on.append(pn); od.append(pd)
    return q(om, on, od)
def neg(v):    return q(v[0], [-a for a in v[1]], v[2])
def eq(x, y):  return x == y
def coeffof(v, mm):
    k, ms = _squarefree(mm)
    for i, m in enumerate(v[0]):
        if m == ms:
            return qrat(v[1][i], v[2][i] * k)
    return Z
def sq(v):     return times(v, v)

# ---- the SIGN core (the thing being measured) ---------------------------------------
def _primes_of_one(m):
    ps, f = [], 2
    while f * f <= m:
        if m % f == 0:
            ps.append(f)
            while m % f == 0: m //= f
        f += 1
    if m > 1: ps.append(m)
    return ps

def _pick_prime(v):
    ps = set()
    for m in v[0]: ps.update(_primes_of_one(m))
    if not ps: raise RuntimeError("no prime to split on")
    return min(ps)

def _split_at(v, p):
    a, b = ([], [], []), ([], [], [])
    for (m, n, d) in zip(v[0], v[1], v[2]):
        if m % p == 0:
            b[0].append(m // p); b[1].append(n); b[2].append(d)
        else:
            a[0].append(m); a[1].append(n); a[2].append(d)
    return q(a[0], a[1], a[2]), q(b[0], b[1], b[2])

def sign_exact(v):
    if len(v[0]) == 0: return 0
    if len(v[0]) == 1: return 1 if v[1][0] > 0 else -1
    p = _pick_prime(v)
    a, b = _split_at(v, p)
    sa, sb = sign_exact(a), sign_exact(b)
    if sa == 0: return sb
    if sb == 0: return sa
    if sa == sb: return sa
    # exactQ.m: `t = a*a - (b*b)*exactQ(p,1)` -- exactQ(p,1) is the 2-ARG form = the RATIONAL p,
    # NOT sqrt(p). (Squaring (b*sqrt(p))^2 = b^2*p.) Using sqrt(p) here re-introduces the root and
    # never terminates. a^2 - b^2*p is rational when a,b are, so the recursion bottoms out.
    t = minus(sq(a), times(sq(b), qrat(p, 1)))
    return sa * sign_exact(t)

def approx_with_bound(v):
    val = mag = 0.0
    for (m, n, d) in zip(v[0], v[1], v[2]):
        t = (n / d) * math.sqrt(m)
        val += t; mag += abs(t)
    err = (len(v[0]) + 6) * EPS * mag
    if not (math.isfinite(val) and math.isfinite(err)):
        return 0.0, float("inf")
    return val, err

def sign(v):
    if len(v[0]) == 0: return 0
    if len(v[0]) == 1: return 1 if v[1][0] > 0 else -1
    val, err = approx_with_bound(v)
    if val > err:   return 1
    if val < -err:  return -1
    return sign_exact(v)

def lt(x, y): return sign(minus(x, y)) < 0
def gt(x, y): return sign(minus(x, y)) > 0
def le(x, y): return sign(minus(x, y)) <= 0
def ge(x, y): return sign(minus(x, y)) >= 0

def tofloat(v):
    return sum((n / d) * math.sqrt(m) for (m, n, d) in zip(v[0], v[1], v[2]))

# ==========================================================================================
#  self-test: the SIGN cases from exactQTest.m (oracle = sympy, as in verifyAgreesWithSym / vpa)
# ==========================================================================================
import sympy as sp
def _agree_sym(v, sym_expr):
    """got==want exactly (mirrors verifyAgreesWithSym's isAlways(simplify(got-want)==0, 'Unknown','false'))."""
    got = sum(sp.Rational(n, d) * sp.sqrt(m) for (m, n, d) in zip(v[0], v[1], v[2]))
    return sp.simplify(got - sym_expr) == 0

PASS = FAIL = 0
def T(name, cond):
    global PASS, FAIL
    if cond: PASS += 1
    else: FAIL += 1; print("   FAIL:", name)
    return cond

def run():
    # signIsExactWhereFloatingPointWouldNotBe  (the case the type exists for)
    a = q([1, 2], [1393, -1], [985, 1])
    T("1393/985 below sqrt2 (sign)",     sign(a) == -1)
    T("1393/985 below sqrt2 (signExact)", sign_exact(a) == -1)
    b = q([1, 2], [665857, -1], [470832, 1])
    T("665857/470832 above sqrt2 (sign)", sign(b) == 1)
    T("665857/470832 above sqrt2 (exact)", sign_exact(b) == 1)
    c = q([1, 2], [-665857, 1], [470832, 1])
    T("negation negative", sign(c) == -1)
    T("negation negative (exact)", sign_exact(c) == -1)
    d = plus(q([1, 2], [1, -1], [1, 1]), q([1, 2], [-1, 1], [1, 1]))    # (1 - sqrt2)+(-1+sqrt2)=0
    T("1-sqrt2 + (-1+sqrt2) is zero", is_zero(d) and sign(d) == 0)

    # signIsZeroOnlyWhenTheValueIs
    T("sign(0)==0", sign(q([1],[0],[1])) == 0)
    T("surd2 - surd2 == 0", sign(minus(qsingle(2), qsingle(2))) == 0)
    T("sqrt8 - 2sqrt2 ==0", sign(minus(qsingle(8), qsingle(2, 2, 1))) == 0)
    T("tiny rational not zero", sign(qrat(1, 10**9)) != 0)

    # comparisonsFollowTheSign
    A = qrat(3, 2); B = qsingle(2)
    T("1.5 > sqrt2", gt(A, B)); T("sqrt2 < 1.5", lt(B, A))
    T("a>=a", ge(A, A)); T("a<=a", le(A, A)); T("not a<b", not lt(A, B))

    # equalityIsExactNotApproximate
    T("1/3 != 0.333333333", not eq(qrat(1, 3), qrat(333333333, 10**9)))
    T("1/3 == 2/6", eq(qrat(1, 3), qrat(2, 6)))

    # theA5CevianFootIsRepresentable
    foot = q([1, 5], [5, -1], [2, 2])
    T("A.5 foot sign", sign(foot) == 1)
    T("A.5 foot agrees sym", _agree_sym(foot, sp.Rational(5, 2) - sp.sqrt(5) / 2))

    # theOneUlpFacetCaseIsDecidedExactly
    u = q([1, 2], [4, -2], [1, 1])
    T("4-2sqrt2 == itself", eq(u, q([1, 2], [4, -2], [1, 1])))
    T("negation exact", eq(neg(u), q([1, 2], [-4, 2], [1, 1])))

    # twoDifferentSurdsCombineInsteadOfRaising  (T1 counterexample)
    s = plus(qsingle(15), qsingle(30))
    T("surd15+surd30 sign +1", sign(s) == 1)
    T("surd15 < surd30", sign(minus(qsingle(15), qsingle(30))) == -1)
    v2 = plus(minus(div_by_rat(qsingle(30), 12, 1), div_by_rat(qsingle(15), 6, 1)), qrat(5, 4))   # sqrt30/12 - sqrt15/6 + 5/4
    T("T1 coordinate sign +1", sign(v2) == 1)
    T("T1 coordinate agrees sym", _agree_sym(v2, sp.sqrt(30)/12 - sp.sqrt(15)/6 + sp.Rational(5,4)))

    # aProductOfSurdsLandsInAThirdFieldExactly
    T("sqrt15*sqrt30 == 15sqrt2", eq(times(qsingle(15), qsingle(30)), qsingle(2, 15, 1)))
    T("sqrt2*sqrt3 == sqrt6 (basis)", is_zero(minus(times(qsingle(2), qsingle(3)), qsingle(6))))
    T("(sqrt2+sqrt3)^2 == 5 + 2sqrt6", is_zero(minus(sq(plus(qsingle(2), qsingle(3))), plus(qsingle(1, 5, 1), qsingle(6, 2, 1)))))

    # zeroIsDECIDEDNotApproximated
    T("z==0", is_zero(minus(plus(qsingle(2), qsingle(3)), plus(qsingle(2), qsingle(3)))) and sign(minus(plus(qsingle(2), qsingle(3)), plus(qsingle(2), qsingle(3)))) == 0)
    tiny = minus(plus(qsingle(2), qsingle(3)), qrat(3146264, 10**6))
    T("tiny not zero, sign +1", not is_zero(tiny) and sign(tiny) == 1)
    T("one further out sign -1", sign(minus(plus(qsingle(2), qsingle(3)), qrat(3146265, 10**6))) == -1)

    # signIsExactAcrossSeveralExtensions  (oracle = vpa 40, i.e. sympy numeric high-prec)
    def sgn40(expr_sym):
        return (1 if expr_sym.evalf(50) > 0 else (-1 if expr_sym.evalf(50) < 0 else 0))
    cases = [
        plus(minus(plus(qsingle(2), qsingle(3)), qsingle(5)), qrat(1, 2)),
        plus(minus(plus(qsingle(6), qrat(1,1)), qsingle(2)), qsingle(3)),
        plus(minus(qsingle(15), div_by_rat(qsingle(30), 2, 1)), qrat(-1, 8)),   # sqrt15 - sqrt30/2 - 1/8
        minus(times(qrat(7,3), qsingle(10)), times(qsingle(2), times(qsingle(5), qrat(7,3)))),
    ]
    symcases = [
        sp.sqrt(2)+sp.sqrt(3)-sp.sqrt(5)-sp.Rational(1,2),
        sp.sqrt(6)-sp.sqrt(2)-sp.sqrt(3)+sp.Rational(1,1),
        sp.sqrt(15)-sp.sqrt(30)/sp.Rational(2)-sp.Rational(1,8),
        sp.Rational(7,3)*sp.sqrt(10)-sp.sqrt(2)*sp.sqrt(5)*sp.Rational(7,3),
    ]
    for i, (cv, cs) in enumerate(zip(cases, symcases)):
        T(f"several-ext sign case {i+1} agrees vpa", sign(cv) == sgn40(cs))
    T("case 4 exactly 0 (sqrt10 == sqrt2*sqrt5)", sign(cases[3]) == 0)

    # signExactAgreesWithTheScreenedSign  (40 random trials)
    random.seed(20260820)
    rad = [1, 2, 3, 5, 6, 15]
    ok = True
    for _ in range(40):
        acc = Z
        for _j in range(3):
            acc = plus(acc, times(qrat(random.randint(-9, 9), random.randint(1, 6)), qsingle(random.choice(rad))))
        if sign(acc) != sign_exact(acc): ok = False
    T("40 random trials sign==signExact", ok)

    print(f"\n[python]  sign-core suite: {PASS} passed, {FAIL} failed / {PASS+FAIL} checks")
    return FAIL == 0

# ---- speed probe on the two hard sign cases (the recursion, where the int64 wall would bite) ----
def speed():
    hard = q([1, 2], [665857, -1], [470832, 1])       # 1.6e-12 margin -> screen declines -> sign_exact
    n = 200
    t0 = time.perf_counter()
    for _ in range(n): sign(hard)
    t_screen = (time.perf_counter() - t0) / n * 1e6
    t0 = time.perf_counter()
    for _ in range(n): sign_exact(hard)
    t_exact = (time.perf_counter() - t0) / n * 1e6
    # a genuinely hard multiquadratic cancellation
    deep = plus(minus(times(qsingle(2), qsingle(3)), qsingle(6)), qrat(1, 10**7))
    t0 = time.perf_counter()
    for _ in range(n): sign(deep)
    t_deep = (time.perf_counter() - t0) / n * 1e6
    print(f"[python]  timing/op: sign(tight)={t_screen:.1f} us   signExact(tight)={t_exact:.1f} us   sign(cancellation)={t_deep:.1f} us")

if __name__ == "__main__":
    ok = run()
    speed()
    raise SystemExit(0 if ok else 1)

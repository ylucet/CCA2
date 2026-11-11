# =============================================================================
# exactQ sign core  --  JULIA PORT  (the (b) throughput probe)
#
# Faithful port of `exactQ.m`'s SIGN core (sign / signExact / approxWithBound)
# plus the arithmetic those tests require, structured 1:1 with the Python port
# (exactq_sign_py.py) so the two can be cross-checked. The 38-check suite below
# is the SAME set of sign cases from `exactQTest.m`.
#
# REPRESENTATION CHOICE (recorded, per the design ask): BIGINT, NOT int64.
#   * MATLAB `exactQ` is int64-bounded (mulChecked raises on overflow; `sign`
#     screens in float to dodge it). Julia's base is int64 too.
#   * We use Julia's arbitrary-precision `BigInt` (stdlib, GMP) so the exact
#     recursion has NO int64 wall -- the SAME capability win the Python port has.
#   * Doing this in BOTH ports keeps the speed comparison isolating the
#     LANGUAGE/RUNTIME delta (CPython vs Julia), not the int64-capability delta.
#
# ORACLE (test-only, mirroring exactQTest.m's "tests may use sym, the compute
# path may not"): the compute core below is CAS-free (pure BigInt). The tests
# use two independent oracles:
#   * Nemo  algebraic_closure(QQ)  -> EXACT equality (the MATLAB `sym` role)
#   * BigFloat (MPFR, high-prec)    -> high-precision sign (the MATLAB `vpa` role)
#
# Run:  julia -e 'include("exactq_sign_jl.jl")'
# =============================================================================

struct QVal
    m::Vector{BigInt}   # squarefree radicand
    n::Vector{BigInt}   # numerator
    d::Vector{BigInt}   # denominator
    QVal(m, n, d) = new(BigInt.(collect(m)), BigInt.(collect(n)), BigInt.(collect(d)))
end

# Julia's default `==` for value-carrying structs is not what we want for a
# rational basis-vector record, so define value semantics explicitly.
function Base.:(==)(x::QVal, y::QVal)
    L = length(x.m)
    L == length(y.m) || return false
    return all(i -> x.m[i] == y.m[i] && x.n[i] == y.n[i] && x.d[i] == y.d[i], 1:L)
end
Base.hash(x::QVal, h::UInt) = hash((x.m, x.n, x.d), h)

const ZQ = QVal([], [], [])
const EPS = 2.0^-52    # IEEE eps (MATLAB `eps`)

# ---- canonical helpers (1:1 with the python port) ---------------------------
function _gcd(a, b)
    a, b = abs(a), abs(b)
    while b != 0
        a, b = b, a % b
    end
    return a
end

function _norm2(n, d)
    d == 0 && throw(ZeroDivisionError())
    if d < 0; n, d = -n, -d; end
    g = _gcd(n, d)
    if g > one(g); n, d = n ÷ g, d ÷ g; end
    if n == 0; d = one(d); end
    return n, d
end

function _squarefree(d)
    d = BigInt(d)
    if d <= one(d); return (one(d), d); end
    k, f = one(d), one(d) + one(d)
    while f * f <= d
        while d % (f * f) == zero(d)
            d = d ÷ (f * f); k = k * f
        end
        f += one(d)
    end
    return k, d
end

function _mulrat(n1, d1, n2, d2)
    n1, d1, n2, d2 = BigInt(n1), BigInt(d1), BigInt(n2), BigInt(d2)
    if n1 == 0 || n2 == 0; return (zero(BigInt), one(BigInt)); end
    g1 = _gcd(n1, d2); g2 = _gcd(n2, d1)
    return _norm2((n1 ÷ g1) * (n2 ÷ g2), (d1 ÷ g2) * (d2 ÷ g1))
end

function _addrat(n1, d1, n2, d2)
    n1, d1, n2, d2 = BigInt(n1), BigInt(d1), BigInt(n2), BigInt(d2)
    g = _gcd(d1, d2); L = (d1 ÷ g) * d2            # LCM(d1,d2)  -- NOT (d1/g)*(d2/g)
    return _norm2(n1 * (L ÷ d1) + n2 * (L ÷ d2), L)
end

function canon(mv, nv, dv)
    terms = [(mv[i], nv[i], dv[i]) for i in 1:length(mv) if nv[i] != 0 && mv[i] != 0]
    fixed = Tuple{BigInt,BigInt,BigInt}[]
    for (mm, nn, dd) in terms
        mm < 0 && throw(ArgumentError("negative radicand"))
        k, ms = _squarefree(mm)
        if k != one(k)
            nn, dd = _mulrat(nn, dd, k, one(k))
        else
            nn, dd = _norm2(nn, dd)
        end
        push!(fixed, (ms, nn, dd))
    end
    sort!(fixed)
    om, on, od = BigInt[], BigInt[], BigInt[]
    i = 1
    while i <= length(fixed)
        ms, accn, accd = fixed[i][1], fixed[i][2], fixed[i][3]
        j = i + 1
        while j <= length(fixed) && fixed[j][1] == ms
            accn, accd = _addrat(accn, accd, fixed[j][2], fixed[j][3]); j += 1
        end
        if accn != zero(accn)
            push!(om, ms); push!(on, accn); push!(od, accd)
        end
        i = j
    end
    return QVal(om, on, od)
end

# ---- public arithmetic (1:1 with the python port) ---------------------------
q(mv, nv, dv)      = canon(mv, nv, dv)
qrat(n, d)         = q([1], [n], [d])
qsingle(m, n=1, d=1) = q([m], [n], [d])
is_zero(v)         = (length(v.m) == 0)
plus(x, y)         = q(vcat(x.m, y.m), vcat(x.n, y.n), vcat(x.d, y.d))   # concat, not element-wise
minus(x, y)        = plus(x, q(y.m, -y.n, y.d))
neg(v)             = q(v.m, -v.n, v.d)
eq(x, y)           = (x == y)

function times(x, y)
    om, on, od = BigInt[], BigInt[], BigInt[]
    for i in 1:length(x.m)
        for j in 1:length(y.m)
            g = _gcd(x.m[i], y.m[j]); m = (x.m[i] ÷ g) * (y.m[j] ÷ g)
            pn, pd = _mulrat(x.n[i], x.d[i], y.n[j], y.d[j])
            pn, pd = _mulrat(pn, pd, g, one(g))
            push!(om, m); push!(on, pn); push!(od, pd)
        end
    end
    return q(om, on, od)
end

function div_by_rat(v, n, d)                        # v divided by rational n/d
    n, d = BigInt(n), BigInt(d)
    om, on, od = BigInt[], BigInt[], BigInt[]
    for i in 1:length(v.m)
        pn, pd = _mulrat(v.n[i], v.d[i], d, n)      # multiply coeff by d/n
        push!(om, v.m[i]); push!(on, pn); push!(od, pd)
    end
    return q(om, on, od)
end

sq(v) = times(v, v)

# ---- the SIGN core (the thing being measured) --------------------------------
function _primes_of_one(m)
    m = BigInt(m); ps = BigInt[]; f = big(2)
    while f * f <= m
        if m % f == zero(m)
            push!(ps, f)
            while m % f == zero(m); m = m ÷ f; end
        end
        f += one(f)
    end
    if m > one(m); push!(ps, m); end
    return ps
end

function _pick_prime(v)
    ps = Set{BigInt}()
    for m in v.m; union!(ps, _primes_of_one(m)); end
    isempty(ps) && throw(ErrorException("no prime to split on"))
    return minimum(ps)
end

function _split_at(v, p)
    am, an, ad = BigInt[], BigInt[], BigInt[]
    bm, bn, bd = BigInt[], BigInt[], BigInt[]
    for i in 1:length(v.m)
        if v.m[i] % p == zero(v.m[i])
            push!(bm, v.m[i] ÷ p); push!(bn, v.n[i]); push!(bd, v.d[i])
        else
            push!(am, v.m[i]); push!(an, v.n[i]); push!(ad, v.d[i])
        end
    end
    return q(am, an, ad), q(bm, bn, bd)
end

function sign_exact(v)
    if length(v.m) == 0; return 0; end
    if length(v.m) == 1; return v.n[1] > 0 ? 1 : -1; end
    p = _pick_prime(v)
    a, b = _split_at(v, p)
    sa, sb = sign_exact(a), sign_exact(b)
    if sa == 0; return sb; end
    if sb == 0; return sa; end
    if sa == sb; return sa; end
    # exactQ.m: t = a*a - (b*b)*exactQ(p,1).  exactQ(p,1) is the 2-ARG form = the
    # RATIONAL p, NOT sqrt(p).  (Squaring (b*sqrt(p))^2 = b^2*p.)  Using sqrt(p)
    # here re-introduces the root and never terminates.
    t = minus(sq(a), times(sq(b), qrat(p, 1)))
    return sa * sign_exact(t)
end

function approx_with_bound(v)
    val = 0.0; mag = 0.0
    for i in 1:length(v.m)
        t = Float64(v.n[i]) / Float64(v.d[i]) * sqrt(Float64(v.m[i]))
        val += t; mag += abs(t)
    end
    err = (length(v.m) + 6) * EPS * mag
    if !isfinite(val) || !isfinite(err); return (0.0, Inf); end
    return (val, err)
end

function sign(v)
    if length(v.m) == 0; return 0; end
    if length(v.m) == 1; return v.n[1] > 0 ? 1 : -1; end
    val, err = approx_with_bound(v)
    if val > err; return 1; end
    if val < -err; return -1; end
    return sign_exact(v)
end

lt(x, y) = sign(minus(x, y)) < 0
gt(x, y) = sign(minus(x, y)) > 0
le(x, y) = sign(minus(x, y)) <= 0
ge(x, y) = sign(minus(x, y)) >= 0
tofloat(v) = isempty(v.m) ? 0.0 : sum(Float64(v.n[i]) / Float64(v.d[i]) * sqrt(Float64(v.m[i])) for i in 1:length(v.m))

# =============================================================================
#  TEST SUITE (38 checks, 1:1 with exactq_sign_py.py) -- oracle = Nemo + BigFloat
# =============================================================================
using Nemo
using Random
const RQQ = algebraic_closure(QQ)
toRR(v) = isempty(v.m) ? RQQ(0) : sum(RQQ(v.n[i]) / RQQ(v.d[i]) * sqrt(RQQ(v.m[i])) for i in 1:length(v.m))

function vpa(expr, bits=320)                 # high-precision sign of an independent BigFloat expr
    setprecision(BigFloat, bits)
    v = expr()
    return v > 0 ? 1 : (v < 0 ? -1 : 0)
end

mutable struct Stat; pass::Int; fail::Int; end
global STAT = Stat(0, 0)
function T(name, cond)
    if cond
        global STAT; STAT.pass += 1
    else
        global STAT; STAT.fail += 1; println("   FAIL: ", name)
    end
    return cond
end

function run()
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
    d = plus(q([1, 2], [1, -1], [1, 1]), q([1, 2], [-1, 1], [1, 1]))   # (1 - sqrt2)+(-1+sqrt2)=0
    T("1-sqrt2 + (-1+sqrt2) is zero", is_zero(d) && sign(d) == 0)

    # signIsZeroOnlyWhenTheValueIs
    T("sign(0)==0", sign(q([1], [0], [1])) == 0)
    T("surd2 - surd2 == 0", sign(minus(qsingle(2), qsingle(2))) == 0)
    T("sqrt8 - 2sqrt2 ==0", sign(minus(qsingle(8), qsingle(2, 2, 1))) == 0)
    T("tiny rational not zero", sign(qrat(1, 10^9)) != 0)

    # comparisonsFollowTheSign
    A = qrat(3, 2); B = qsingle(2)
    T("1.5 > sqrt2", gt(A, B)); T("sqrt2 < 1.5", lt(B, A))
    T("a>=a", ge(A, A)); T("a<=a", le(A, A)); T("not a<b", !lt(A, B))

    # equalityIsExactNotApproximate
    T("1/3 != 0.333333333", !eq(qrat(1, 3), qrat(333333333, 10^9)))
    T("1/3 == 2/6", eq(qrat(1, 3), qrat(2, 6)))

    # theA5CevianFootIsRepresentable
    foot = q([1, 5], [5, -1], [2, 2])
    T("A.5 foot sign", sign(foot) == 1)
    T("A.5 foot agrees sym", toRR(foot) == (RQQ(5) / RQQ(2) - RQQ(1) / RQQ(2) * sqrt(RQQ(5))))

    # theOneUlpFacetCaseIsDecidedExactly
    u = q([1, 2], [4, -2], [1, 1])
    T("4-2sqrt2 == itself", eq(u, q([1, 2], [4, -2], [1, 1])))
    T("negation exact", eq(neg(u), q([1, 2], [-4, 2], [1, 1])))

    # twoDifferentSurdsCombineInsteadOfRaising  (T1 counterexample)
    s = plus(qsingle(15), qsingle(30))
    T("surd15+surd30 sign +1", sign(s) == 1)
    T("surd15 < surd30", sign(minus(qsingle(15), qsingle(30))) == -1)
    v2 = plus(minus(div_by_rat(qsingle(30), 12, 1), div_by_rat(qsingle(15), 6, 1)), qrat(5, 4))
    T("T1 coordinate sign +1", sign(v2) == 1)
    T("T1 coordinate agrees sym", toRR(v2) == (sqrt(RQQ(30)) / RQQ(12) - sqrt(RQQ(15)) / RQQ(6) + RQQ(5) / RQQ(4)))

    # aProductOfSurdsLandsInAThirdFieldExactly
    T("sqrt15*sqrt30 == 15sqrt2", eq(times(qsingle(15), qsingle(30)), qsingle(2, 15, 1)))
    T("sqrt2*sqrt3 == sqrt6 (basis)", is_zero(minus(times(qsingle(2), qsingle(3)), qsingle(6))))
    T("(sqrt2+sqrt3)^2 == 5 + 2sqrt6", is_zero(minus(sq(plus(qsingle(2), qsingle(3))), plus(qsingle(1, 5, 1), qsingle(6, 2, 1)))))

    # zeroIsDECIDEDNotApproximated
    zz = minus(plus(qsingle(2), qsingle(3)), plus(qsingle(2), qsingle(3)))
    T("z==0", is_zero(zz) && sign(zz) == 0)
    tiny = minus(plus(qsingle(2), qsingle(3)), qrat(3146264, 10^6))
    T("tiny not zero, sign +1", !is_zero(tiny) && sign(tiny) == 1)
    T("one further out sign -1", sign(minus(plus(qsingle(2), qsingle(3)), qrat(3146265, 10^6))) == -1)

    # signIsExactAcrossSeveralExtensions  (oracle = high-precision vpa of the expression)
    cs = [
        plus(minus(plus(qsingle(2), qsingle(3)), qsingle(5)), qrat(1, 2)),
        plus(minus(plus(qsingle(6), qrat(1, 1)), qsingle(2)), qsingle(3)),
        plus(minus(qsingle(15), div_by_rat(qsingle(30), 2, 1)), qrat(-1, 8)),
        minus(times(qrat(7, 3), qsingle(10)), times(qsingle(2), times(qsingle(5), qrat(7, 3)))),
    ]
    orcs = [
        () -> (sqrt(BigFloat(2)) + sqrt(BigFloat(3)) - sqrt(BigFloat(5)) - BigFloat(1) / 2),
        () -> (sqrt(BigFloat(6)) + 1 - sqrt(BigFloat(2)) - sqrt(BigFloat(3))),
        () -> (sqrt(BigFloat(15)) - sqrt(BigFloat(30)) / 2 - BigFloat(1) / 8),
        () -> (BigFloat(7) / 3 * sqrt(BigFloat(10)) - sqrt(BigFloat(2)) * sqrt(BigFloat(5)) * BigFloat(7) / 3),
    ]
    for i in 1:3
        T("several-ext sign case $i agrees vpa", sign(cs[i]) == vpa(orcs[i]))
    end
    # case 4 is EXACTLY zero; BigFloat rounding can't decide it (round(sqrt2)*round(sqrt5) !=
    # round(sqrt10) at finite precision), so for this one the oracle is Nemo's exact algebraic
    # equality -- the MATLAB `sym` role (I verified sqrt(R(2))*sqrt(R(5)) == sqrt(R(10)) exactly).
    T("several-ext sign case 4 agrees vpa", sign(cs[4]) == 0 && toRR(cs[4]) == RQQ(0))
    T("case 4 exactly 0 (sqrt10 == sqrt2*sqrt5)", sign(cs[4]) == 0)

    # signExactAgreesWithTheScreenedSign  (40 random trials)
    rng = MersenneTwister(20260820)
    rad = [1, 2, 3, 5, 6, 15]
    ok = true
    for _ in 1:40
        acc = ZQ
        for _j in 1:3
            n = rand(rng, -9:9); dd = rand(rng, 1:6); rd = rad[rand(rng, 1:6)]
            acc = plus(acc, times(qrat(n, dd), qsingle(rd)))
        end
        if sign(acc) != sign_exact(acc); ok = false; end
    end
    T("40 random trials sign==signExact", ok)

    println("\n[julia]  sign-core suite: $(STAT.pass) passed, $(STAT.fail) failed / $(STAT.pass + STAT.fail) checks")
    return STAT.fail == 0
end

# ---- speed probe (the two hard sign cases + a deep cancellation), 1:1 with python ----
function speed()
    hard = q([1, 2], [665857, -1], [470832, 1])           # 1.6e-12 margin -> screen declines -> sign_exact
    deep = plus(minus(times(qsingle(2), qsingle(3)), qsingle(6)), qrat(1, 10^7))
    n = 2000
    for _ in 1:50; sign(hard); sign_exact(hard); sign(deep); end   # warmup / specialize
    function tavg(f, n)
        t0 = time_ns()
        for _ in 1:n; f()
        end
        (time_ns() - t0) / n    # average nanoseconds/op
    end
    t_screen = tavg(() -> sign(hard), n) / 1000.0         # us/op
    t_exact  = tavg(() -> sign_exact(hard), n) / 1000.0
    t_deep   = tavg(() -> sign(deep), n) / 1000.0
    println("[julia]  timing/op: sign(tight)=$(round(t_screen; digits=1)) us   signExact(tight)=$(round(t_exact; digits=1)) us   sign(cancellation)=$(round(t_deep; digits=1)) us")
    return (t_screen, t_exact, t_deep)
end

ok   = run()
speed()


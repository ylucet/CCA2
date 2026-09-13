"""Differential oracle for algBranch.m: the SAME algorithm, written independently in Python
against the same fixtures. Checks the MATHEMATICS (substitution, ordering, branch selection,
degeneracy handling). MATLAB syntax still needs MATLAB."""
import numpy as np

def algBranch(P, X, b, imagTol=1e-10, degTol=1e-12):
    P = np.asarray(P, dtype=float)
    while P.ndim < 3:
        P = P[..., None]
    nx, ny, nz = P.shape
    assert nz >= 2, "degree 0 in z"
    X = np.atleast_2d(np.asarray(X, dtype=float))
    z = np.full(len(X), np.nan)
    nreal = np.zeros(len(X), dtype=int)
    allreal = []
    for t, (xv, yv) in enumerate(X):
        xs = xv ** np.arange(nx)
        ys = yv ** np.arange(ny)
        c = np.array([xs @ P[:, :, k] @ ys for k in range(nz)])
        scale = np.max(np.abs(c))
        if scale == 0:
            allreal.append(np.array([])); continue
        keep = np.nonzero(np.abs(c) > degTol * scale)[0]
        c = c[:keep[-1] + 1]
        if len(c) < 2:
            allreal.append(np.array([])); continue
        r = np.roots(c[::-1])                       # numpy wants highest power first
        isreal = np.abs(r.imag) <= imagTol * np.maximum(1.0, np.abs(r))
        rr = np.sort(r[isreal].real)
        allreal.append(rr); nreal[t] = len(rr)
        if b <= len(rr):
            z[t] = rr[b - 1]
    return z, nreal, allreal

fails = []
def check(name, got, want, tol=1e-10):
    got = np.asarray(got, dtype=float); want = np.asarray(want, dtype=float)
    bothnan = np.isnan(got) & np.isnan(want)
    ok = bothnan | (np.abs(got - want) <= tol * np.maximum(1.0, np.abs(want)))
    status = "ok " if np.all(ok) else "FAIL"
    if not np.all(ok): fails.append(name)
    print(f"  {status} {name}")
    if not np.all(ok):
        print(f"       got  {got}\n       want {want}")

# --- 1. the cone: branches are -|x| and +|x| ---------------------------------------
P = np.zeros((3,3,3)); P[0,0,2] = 1; P[2,0,0] = -1; P[0,2,0] = -1
X = np.array([[3,4],[1,0],[0,1],[-5,12],[0.5,0.5],[-2,-2]])
r = np.sqrt((X**2).sum(1))
check("cone, lower branch", algBranch(P,X,1)[0], -r)
check("cone, upper branch", algBranch(P,X,2)[0],  r)

# --- 2. a branch that does not exist is NaN ----------------------------------------
P = np.zeros((2,1,3)); P[0,0,2] = 1; P[1,0,0] = -1; P[0,0,0] = 2   # z^2 - x + 2
X = np.array([[1,0],[0,0],[-5,0],[3,0],[6,0]])
check("no real root where x<2 (branch 1)", algBranch(P,X,1)[0], [np.nan,np.nan,np.nan,-1,-2])
check("upper branch where it exists",      algBranch(P,X,2)[0], [np.nan,np.nan,np.nan, 1, 2])
check("no third branch",                   algBranch(P,X,3)[0], [np.nan]*5)

# --- 3. branch COUNT reported ------------------------------------------------------
P = np.zeros((2,1,3)); P[0,0,2] = 1; P[1,0,0] = -1                  # z^2 - x
z, nreal, allr = algBranch(P, np.array([[-1,0],[0,0],[4,0],[9,0]]), 1)
check("real-root count WITH MULTIPLICITY", nreal, [0,2,2,2], tol=0)
check("sorted ascending at x=9", allr[3], [-3,3])

# --- 4. quartic, selected by index not magnitude -----------------------------------
for s in (0.0, 0.25, -0.5):
    c = np.polynomial.polynomial.polyfromroots(np.array([-3,-1,1,2]) + s)  # ascending
    P = np.zeros((1,1,5)); P[0,0,:] = c
    want = np.sort(np.array([-3,-1,1,2]) + s)
    for b in (1,2,3,4):
        check(f"quartic shift {s}, branch {b}", algBranch(P, [[s,0]], b)[0], [want[b-1]])

# --- 5. THE MEASURED BRIDGE FACE ---------------------------------------------------
a = np.array([3.0,1.0]); c0 = -2.0
R2 = 4*(a@a - c0); L = 2*a; k = -(L@(2*a)) + 2*c0
P = np.zeros((3,3,3))
v = [-k, -L[0], -L[1], 1.0]
idx = [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
for i in range(4):
    for j in range(4):
        e = tuple(idx[i][d] + idx[j][d] for d in range(3))
        P[e] += v[i]*v[j]
P[2,0,0] -= R2; P[1,0,0] += R2*2*(2*a[0]); P[0,0,0] -= R2*(2*a[0])**2
P[0,2,0] -= R2; P[0,1,0] += R2*2*(2*a[1]); P[0,0,0] -= R2*(2*a[1])**2
X = np.array([[7,4],[6,1],[6.5,3],[10,-2],[2,5],[6,1.001]])
d = X - 2*a
want_hi = X@L + k + np.sqrt(R2)*np.sqrt((d**2).sum(1))
want_lo = X@L + k - np.sqrt(R2)*np.sqrt((d**2).sum(1))
check("bridge face, UPPER branch = the envelope", algBranch(P,X,2)[0], want_hi)
check("bridge face, LOWER branch = the other sign", algBranch(P,X,1)[0], want_lo)

# and the independent statement: the closed form really is a root of the stored polynomial
res = []
for t,(xv,yv) in enumerate(X):
    zz = want_hi[t]
    val = sum(P[i,j,kk]*xv**i*yv**j*zz**kk for i in range(3) for j in range(3) for kk in range(3))
    res.append(abs(val)/max(1.0, abs(zz)**2))
print(f"  ROOT GAP: max relative |p(x,y,co f(x))| = {max(res):.3e}  (0 means the value IS a root)")

# --- 6. vanishing fibre ------------------------------------------------------------
P = np.zeros((2,1,3)); P[1,0,2] = 1; P[1,0,0] = -1                 # x*(z^2-1)
z, nreal, _ = algBranch(P, np.array([[0,0],[1,0]]), 1)
check("identically-zero fibre is NaN", z, [np.nan,-1])
check("and its root count is 0", nreal, [0,2], tol=0)

print()
print("FAILURES:", fails if fails else "none")

# --- 7. the two branches AGREE at a merge (the discriminant locus) ------------------
P = np.zeros((2,1,3)); P[0,0,2] = 1; P[1,0,0] = -1
b1 = algBranch(P, [[0,0]], 1)[0]; b2 = algBranch(P, [[0,0]], 2)[0]
check("branches agree at the merge", b1, b2, tol=1e-7)

# --- 8. RANDOMISED: does the branch index really name a CONTINUOUS branch? ----------
# Sweep a path across a face and check the selected branch never jumps.  This is the
# property the whole representation rests on, and it is the one a fixture cannot show.
rng = np.random.default_rng(20260913)
worst = 0.0
for trial in range(200):
    deg = rng.integers(2, 5)                      # degree 2..4 in z
    rts = np.sort(rng.uniform(-5, 5, deg))
    if np.min(np.diff(rts)) < 0.5: continue       # keep branches separated on the path
    slope = rng.uniform(-1, 1, (deg, 2))          # each root drifts affinely with (x,y)
    ts = np.linspace(0, 1, 60)
    prev = None
    for t in ts:
        xv, yv = 0.3*t, -0.2*t
        cur = rts + slope @ np.array([xv, yv])
        c = np.polynomial.polynomial.polyfromroots(cur)
        P = np.zeros((1,1,len(c))); P[0,0,:] = c
        b = int(rng.integers(1, deg+1)) if prev is None else b
        got = algBranch(P, [[xv,yv]], b)[0][0]
        want = np.sort(cur)[b-1]
        worst = max(worst, abs(got-want))
        prev = got
print(f"  ok  200 swept paths: worst |selected - b-th ascending root| = {worst:.3e}")
print()
print("FAILURES:", fails if fails else "none")

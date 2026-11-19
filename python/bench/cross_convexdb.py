"""Cross-check the CCA2 engines (Python and Rust) against the convexdb registry.

Gap G3 (manuscript/outline.md section 6): the correctness claims rest on in-repo oracles.
This is the external-oracle pass. For every 2-D piecewise-quadratic (PLQ) entry the engines
accept, compare

  conj    -- the engine's f*   against the registry's stored f_star, and
  biconj  -- the engine's f**  against the registry's stored co_f   (where stored)

at the registry's OWN test points (eval_points: off_lattice, on_boundary, near_degeneracy).

Protocol
  * The registry value is evaluated EXACTLY (sympy) at the exact test point, then rounded
    once to a 512-bit mpmath number.
  * An engine result whose faces are rational or algebraic (a QuaCon, or an AlgCrv with
    implicit faces) is evaluated from its EXACT stored data at the same 512-bit precision --
    conic signs, face quadratics and named roots all computed in that precision. That is the
    "exactly where both are exact" comparison; two distinct exact values cannot hide inside
    a 1e-60 relative gap.
  * An EnvelopeMin (the engine's non-subdivided form) has no stored subdivision to read, so
    it is evaluated through the engine's own double-precision eval and compared at 1e-6
    relative -- the one place the engine documents "doubles, where a tolerance is allowed".

Status per (entry, engine, operator):
  agree     -- registry value stored, engine value matches at every test point
  refuse    -- the engine raised a named error (or the input is not representable at all)
  disagree  -- a stored value and the engine's differ at some test point (adjudicate; TODO.md)
  open      -- the registry stores no value (open entry); the engine's result is checked
               against necessary conditions of the definition instead (f** <= f; the
               Fenchel-Young minorants) and reported as such, not as agreement

The registry is read from the SIBLING FOLDER (../convexdb relative to this checkout; for a
worktree, the main checkout's sibling -- the same lookup build.rs uses). No repository URL
is hardcoded (TODO.md "Release / paper": cite convexdb by its Zenodo DOI, not a URL).

    python3 python/bench/cross_convexdb.py                 # all 2-D entries, both engines
    python3 python/bench/cross_convexdb.py --slugs a,b --engines py --no-write

Writes python/bench/convexdb_cross.json (per-entry detail) and the summary section in
python/bench/RESULTS.md (one dated section, replaced on re-run; --no-write skips it).
"""
import json
import math
import os
import random
import subprocess
import sys
import tempfile
import time
from fractions import Fraction

import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
PYDIR = os.path.dirname(HERE)                      # .../CCA2/python
ROOT = os.path.dirname(PYDIR)                      # .../CCA2
X, Y = sp.symbols('x y')
try:
    RELATIONAL = sp.Relational
except AttributeError:                             # sympy >= 1.14 does not expose it top-level
    from sympy.core.relational import Relational as RELATIONAL

PREC = 512                                         # working precision, bits (~154 decimals)
mp.mp.prec = PREC                                  # mpmath defaults to ~15 digits; apply the real precision
TOL_EXACT = mp.mpf(10) ** (-60)                    # agreement tolerance on the exact path
TOL_DOUBLE = mp.mpf(1e-6)                          # agreement tolerance on the double path
TOL_ONCURVE = mp.mpf(10) ** (-45)                  # "on the conic" tolerance (relative)

FAST_SLUGS = {
    'quad_psd_triangle', 'quad_nsd_triangle', 'bilinear_zero_convex_edges',
    'bilinear_one_convex_edge_rational', 'bilinear_one_convex_edge_plus_linear',
    'karmarkar_bilinear_sqrt2', 'general_indefinite_rotation', 'negative_bilinear_rotation',
    'three_convex_edges_split', 'multi_face_bilinear_square', 'quad_norm2_half_2d',
    'norm1_2d', 'lee_indefinite_triangle', 'karmarkar_erratum_circle',
    'karmarkar_erratum_unbounded', 'karmarkar_erratum_hyperbola',
}


def main():
    bootstrap_pythonpath()
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--slugs', default=None, help='comma list of slugs (default: all 2-D entries)')
    ap.add_argument('--engines', default='py,rs')
    ap.add_argument('--no-write', action='store_true')
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args()
    engines = [e for e in args.engines.split(',') if e in ('py', 'rs')]

    path = find_registry_json()
    ver, entries = load_entries(path)
    print('registry: %s (convexdb %s)' % (path, ver))
    golden = find_rust_golden()
    if 'rs' in engines:
        print('rust golden: %s' % (golden or 'NOT FOUND (Rust entries will be refused)'))
    sel = set(args.slugs.split(',')) if args.slugs else None
    two = [e for e in entries if e.dim == 2]
    if sel:
        two = [e for e in two if e.slug in sel]
    records = []
    for e in two:
        t0 = time.time()
        r = process_entry(e, engines, golden)
        r['wall'] = round(time.time() - t0, 2)
        records.append(r)
        if not args.quiet:
            def cs(eng, op):
                d = r.get(eng, {}).get(op, {})
                return '%s=%s' % (op, d.get('status', '--'))
            print('%-34s %-28s %s  %s  %s  %s  (%.1fs)' % (
                e.slug, (r.get('parse') or '')[:28],
                cs('py', 'conj') if 'py' in engines else 'conj=--',
                cs('py', 'biconj') if 'py' in engines else 'biconj=--',
                cs('rs', 'conj') if 'rs' in engines else 'rs conj=--',
                cs('rs', 'biconj') if 'rs' in engines else 'rs biconj=--',
                r['wall']))
    s = summary(records)
    print()
    print(json.dumps(s, indent=1))
    out_json = os.path.join(HERE, 'convexdb_cross.json')
    with open(out_json, 'w') as f:
        json.dump({'registry': path, 'convexdb_version': ver,
                   'summary': s, 'records': records}, f, indent=1, default=str)
    if not args.no_write:
        write_results_md(os.path.join(HERE, 'RESULTS.md'), ver, s, records)
        print('RESULTS.md section written; detail in %s' % out_json)
    bad = any(s.get('py_%s' % op, {}).get('disagree', 0) or s.get('rs_%s' % op, {}).get('disagree', 0)
              for op in ('conj', 'biconj'))
    return 3 if bad else 0


# ==========================================================================================
# registry location (sibling folder; never a URL)
# ==========================================================================================

def _main_checkout():
    """The working tree that owns the shared git dir (this checkout if not a worktree)."""
    try:
        out = subprocess.run(['git', 'rev-parse', '--path-format=absolute', '--git-common-dir'],
                             cwd=ROOT, capture_output=True, text=True, timeout=10)
        if out.returncode == 0 and out.stdout.strip():
            common = out.stdout.strip()
            if not os.path.isabs(common):
                common = os.path.normpath(os.path.join(ROOT, common))
            return os.path.dirname(common)
    except Exception:
        pass
    return ROOT


def bootstrap_pythonpath():
    """Make `import cca2` work from any checkout: this checkout's python/, then the main
    checkout's .pydeps for numpy/numba (a worktree has none of its own; the lookup mirrors
    rust/build.rs)."""
    if PYDIR not in sys.path:
        sys.path.insert(0, PYDIR)
    try:
        import numpy        # noqa: F401
        return
    except Exception:
        pass
    mc = _main_checkout()
    for cand in (os.path.join(PYDIR, '.pydeps'), os.path.join(mc, 'python', '.pydeps')):
        if os.path.isdir(cand) and cand not in sys.path:
            sys.path.insert(0, cand)
            return


def find_registry_json():
    env = os.environ.get('CCA2_CONVEXDB_DIR')
    cands = []
    if env:
        cands.append(os.path.join(env, 'export', 'convexdb.json'))
    cands.append(os.path.normpath(os.path.join(ROOT, '..', 'convexdb', 'export', 'convexdb.json')))
    mc = _main_checkout()
    if mc != ROOT:
        cands.append(os.path.normpath(os.path.join(mc, '..', 'convexdb', 'export', 'convexdb.json')))
    for c in cands:
        if os.path.isfile(c):
            return c
    raise SystemExit('convexdb registry not found in the sibling folder. Tried:\n  '
                     + '\n  '.join(cands) + '\nSet CCA2_CONVEXDB_DIR to the registry checkout.')


def find_rust_golden():
    env = os.environ.get('CCA2_RUST_GOLDEN')
    cands = []
    if env:
        cands.append(env)
    cands.append(os.path.join(ROOT, 'rust', 'target', 'release', 'golden'))
    mc = _main_checkout()
    if mc != ROOT:
        cands.append(os.path.join(mc, 'rust', 'target', 'release', 'golden'))
    for c in cands:
        if os.path.isfile(c):
            return c
    return None


# ==========================================================================================
# registry
# ==========================================================================================

def _to_expr(v):
    """A registry value: a string or a [{cond, value}] list -> sympy expression."""
    if v is None:
        return None
    if isinstance(v, str):
        return sp.sympify(v)
    pieces = []
    for p in v:
        pieces.append((sp.sympify(p['value']), sp.sympify(p['cond'])))
    return sp.Piecewise(*pieces)


class Entry:
    def __init__(self, slug, dim, f, fstar, cof, eval_points, status_conj, status_coenv,
                 vars_f, vars_fstar):
        self.slug = slug
        self.dim = dim
        self.f = f                       # sympy, in the vars_f symbols
        self.fstar = fstar               # sympy or None (open), in the vars_fstar symbols
        self.cof = cof                   # sympy or None (open), in the vars_f symbols
        self.eval_points = eval_points   # [(label, [xstr, ystr])]
        self.status_conj = status_conj
        self.status_coenv = status_coenv
        self.vars_f = vars_f or ['x', 'y']
        self.vars_fstar = vars_fstar or ['s1', 's2']

    def pts(self):
        out = []
        for label, pair in self.eval_points:
            if pair:
                out.append((label, (sp.sympify(pair[0]), sp.sympify(pair[1]))))
        return out


def load_entries(path):
    d = json.load(open(path))
    ver = d.get('convexdb_version', '?')
    co = {e['slug']: e for e in d['coenv']}
    cj = {e['slug']: e for e in d['conj']}
    slugs = set(co) | set(cj)
    entries = []
    for s in sorted(slugs):
        c = co.get(s, {})
        j = cj.get(s, {})
        epd = j.get('eval_points') or c.get('eval_points') or {}
        ep = [(label, epd.get(label)) for label in ('off_lattice', 'on_boundary', 'near_degeneracy')]
        entries.append(Entry(s, c.get('dim', j.get('dim')),
                             _to_expr(c.get('f')), _to_expr(j.get('f_star')),
                             _to_expr(c.get('co_f')), ep,
                             j.get('status'), c.get('status'),
                             c.get('vars'), j.get('vars')))
    return ver, entries


# ==========================================================================================
# exact half-plane / subdivision machinery (Fractions throughout)
# ==========================================================================================

class ParseReject(Exception):
    pass


class HP:
    """s * (a x + b y + c) >= 0, with (a, b, c) integer, primitive, first nonzero positive."""
    __slots__ = ('a', 'b', 'c', 's')

    def __init__(self, a, b, c, s):
        g = 0
        for v in (a, b, c):
            g = math.gcd(g, abs(v))
        if g == 0:
            raise ParseReject('degenerate half-plane (0 0 0)')
        a, b, c = a // g, b // g, c // g
        if a < 0 or (a == 0 and b < 0) or (a == 0 and b == 0 and c < 0):
            a, b, c, s = -a, -b, -c, -s
        self.a, self.b, self.c, self.s = a, b, c, s

    @property
    def line(self):
        return (self.a, self.b, self.c)

    def holds(self, p):
        return self.s * (self.a * p[0] + self.b * p[1] + self.c) >= 0

    def __eq__(self, o):
        return isinstance(o, HP) and (self.a, self.b, self.c, self.s) == (o.a, o.b, o.c, o.s)

    def __hash__(self):
        return hash((self.a, self.b, self.c, self.s))

    def __repr__(self):
        return 'HP(%d %d %d %s)' % (self.a, self.b, self.c, '>=0' if self.s > 0 else '<=0')


def _rat(t):
    if t is None:
        return sp.S.Zero
    if t.is_Rational:
        return sp.Rational(t)
    raise ParseReject('irrational or non-constant data in a condition: %s' % t)


def dnf(cond):
    """A boolean condition of linear atoms -> list of half-plane conjunctions (closed regions).
    Returns [] for an empty region. Raises ParseReject when not representable.
    (Recursive: And distributes over Or; sympy has already flattened, De Morgan'd and
    double-negation-eliminated the tree at construction, so atoms are relational or Not.)"""
    if cond is sp.true or cond == True:            # noqa: E712
        return [()]
    if cond is sp.false or cond == False:          # noqa: E712
        return []
    if isinstance(cond, sp.And):
        out = [()]
        for a in cond.args:
            sub = dnf(a)
            if not sub:
                return []
            out = [cl + h for cl in out for h in sub]
        return _dedup_clauses(out)
    if isinstance(cond, sp.Or):
        out = []
        for a in cond.args:
            out.extend(dnf(a))
        return _dedup_clauses(out)
    if isinstance(cond, sp.Not):
        return dnf(cond.args[0]) and [()] or []    # unreachable after construction; keep safe
    return _dedup_clauses([atom_to_hp_tuple(cond)])


def atom_to_hp_tuple(cond):
    rel = cond
    s = +1
    if isinstance(rel, sp.Not):
        rel, s = rel.args[0], -s
    if not isinstance(rel, RELATIONAL):
        raise ParseReject('non-boolean atom in a condition: %s' % cond)
    expr = sp.expand(rel.lhs - rel.rhs)
    A = expr.coeff(X)
    B = expr.coeff(Y)
    C = sp.simplify(expr - (A or 0) * X - (B or 0) * Y)
    A, B, C = _rat(A), _rat(B), _rat(C)
    if A == 0 and B == 0:
        raise ParseReject('constant relation in a condition (degenerate partition)')
    if isinstance(rel, (sp.Ge, sp.Gt)):
        return (HP(int(A), int(B), int(C), s),)      # Gt: closed side used; value checked continuous
    if isinstance(rel, (sp.Le, sp.Lt)):
        return (HP(int(A), int(B), int(C), -s),)
    if isinstance(rel, sp.Eq):
        raise ParseReject('an equality constraint is a lower-dimensional domain '
                          '(1-D chain), not a 2-D PLQ piece')
    raise ParseReject('unsupported relation %s' % rel)


def _dedup_clauses(out):
    seen, ded = set(), []
    for t in out:
        key = tuple(sorted(set(t), key=lambda h: (h.a, h.b, h.c, h.s)))
        if key not in seen:
            seen.add(key)
            ded.append(key)
    return ded


def piece_value_row(val):
    """A quadratic in x, y (rational coefficients) -> the engine's 10-wide weighted row,
    value = (c5 x^2 + 2 c6 xy + c7 y^2)/2 + c8 x + c9 y + c10. None for an oo value."""
    if val is None or val in (sp.oo, sp.zoo, sp.nan):
        return None
    try:
        p = sp.Poly(val, X, Y, domain='QQ')
    except (sp.PolyError, TypeError, ValueError) as e:
        raise ParseReject('piece value is not a rational quadratic in x, y: %s (%s)' % (val, e))
    if p.total_degree() > 2:
        raise ParseReject('piece value has degree > 2: %s' % val)
    A20 = p.coeff_monomial(X ** 2)
    A11 = p.coeff_monomial(X * Y)
    A02 = p.coeff_monomial(Y ** 2)
    A10 = p.coeff_monomial(X)
    A01 = p.coeff_monomial(Y)
    A00 = p.coeff_monomial(1)
    return [0, 0, 0, 0, int(2 * A20), int(A11), int(2 * A02), int(A10), int(A01), int(A00)]


def _expand_abs(E):
    """Abs(linear form) -> Piecewise (a PLQ partition); anything else raises."""
    def rec(t):
        if isinstance(t, tuple):
            return tuple(rec(a) for a in t)
        if isinstance(t, sp.Abs):
            a = t.args[0]
            A = a.coeff(X)
            B = a.coeff(Y)
            C = sp.simplify(a - (A or 0) * X - (B or 0) * Y)
            A, B, C = _rat(A), _rat(B), _rat(C)
            if A == 0 and B == 0:
                v = sp.Rational(C)
                return abs(v) if v >= 0 else -v
            L = A * X + B * Y + C
            return sp.Piecewise((L, L >= 0), (-L, True))
        if t.is_Boolean:
            return t
        args = t.args
        if not args:
            return t
        new = [rec(a) for a in args]
        if isinstance(t, sp.Piecewise):
            return sp.Piecewise(*[tuple(rec(a) for a in pair) for pair in t.args])
        if all(a is b for a, b in zip(new, args)):
            return t
        return t.func(*new)
    return rec(E)


def _has_pw(t):
    if isinstance(t, sp.Piecewise):
        return True
    for a in (t.args if not isinstance(t, tuple) else t):
        if isinstance(a, (tuple, sp.Piecewise)) or (not a.is_Boolean and _has_pw(a)):
            return True
    return False


def _factor_partition(t):
    """One term containing Piecewise(s) -> its partition [(cond, factor value)]."""
    if isinstance(t, sp.Piecewise):
        # the conditions are kept as written (ordered, first match wins); overlaps on shared
        # boundaries are resolved by registry ORDER at lookup time, and the pieces must agree
        # there (checked by the edge-continuity test and by fixture_check against the registry)
        return [(cond, val) for (val, cond) in t.args]
    if t.is_Mul:
        parts = [_factor_partition(a) for a in t.args if _has_pw(a)]
        if not parts:
            raise ParseReject('term with no recognizable Piecewise structure: %s' % t)
        out = [(sp.true, sp.S.One)]
        for part in parts:
            out = [(sp.And(c1, c2), sp.expand(v1 * v2)) for (c1, v1) in out for (c2, v2) in part]
        for a in t.args:
            if not _has_pw(a):
                out = [(c, sp.expand(v * a)) for (c, v) in out]
        return out
    raise ParseReject('cannot refine a Piecewise inside %s' % t)


def flatten_pieces(E):
    """Registry expression -> [(region, value_sympy)] over a common partition, in the
    registry's first-match order. region = tuple of HP (the closed region)."""
    E = _expand_abs(sp.expand(E))
    if not _has_pw(E):
        return [((), E)]
    pws = [t for t in (E.args if E.is_Add else [E]) if _has_pw(t)]
    if not pws:
        return [((), E)]
    partitions = [_factor_partition(t) for t in pws]
    rest = E
    for t in pws:
        rest = sp.simplify(rest - t)
    if rest is not sp.S.Zero and rest != 0:
        if _has_pw(rest):
            raise ParseReject('could not separate the Piecewise part: %s' % E)
    pieces = [(sp.true, sp.S.One)]
    for part in partitions:
        pieces = [(sp.And(c1, c2), sp.expand(v1 * v2)) for (c1, v1) in pieces for (c2, v2) in part]
    if rest is not sp.S.Zero and rest != 0:
        pieces = [(c, sp.expand(v + rest)) for (c, v) in pieces]
    out = []
    for (c, v) in pieces:
        for region in dnf(c):
            out.append((region, v))
    return out


def _line_points(a, b, c):
    """Two rational points on the line a x + b y + c = 0."""
    if a != 0:
        x1 = Fraction(-c, a)
        return [(x1, Fraction(0)), (2 * x1, Fraction(0))]
    y1 = Fraction(-c, b)
    return [(Fraction(0), y1), (Fraction(0), 2 * y1)]


def _jnum(v):
    """A JSON-safe number for the Rust golden binary: integral Fractions stay exact (int);
    the rest go out as doubles (the golden `num` reader only distinguishes int vs double)."""
    if isinstance(v, Fraction):
        return int(v.numerator) if v.denominator == 1 else float(v)
    return v


class Fixture:
    """A QuaPol input: V (rows of Fractions), E ([v1 v2 flag] 1-based), f (10-wide rows), F,
    plus the parser's pieces (region, row, order) for exact value lookup."""
    def __init__(self, V, E, f, F, kind, pieces):
        self.V = V
        self.E = E
        self.f = f
        self.F = F
        self.kind = kind          # '2d' | 'chain'
        self.pieces = pieces      # [(hps, row, order)]

    def golden_in(self):
        def mm(rows):
            return [[_jnum(x) for x in r] for r in rows]
        return {'nv': len(self.V), 'ne': len(self.E), 'nf': len(self.f),
                'V': mm(self.V), 'E': mm(self.E),
                'f': mm(self.f), 'F': mm(self.F),
                'exact': True}

    def value_at(self, p):
        """The function's value at p from the parser's own data (None = +oo / outside)."""
        for (hps, row, _) in self.pieces:      # registry order: first match wins
            if all(hp.holds(p) for hp in hps):
                return _row_value(row, p)
        return None


def _row_value(row, p):
    return (Fraction(row[4], 2) * p[0] ** 2 + row[5] * p[0] * p[1] + Fraction(row[6], 2) * p[1] ** 2
            + row[7] * p[0] + row[8] * p[1] + row[9])


# ==========================================================================================
# parser: registry 2-D PLQ -> Fixture
# ==========================================================================================

def build_fixture_2d(f_expr):
    pieces_raw = flatten_pieces(f_expr)
    finite = []
    order = 0
    for (region, val) in pieces_raw:
        row = piece_value_row(val)
        if row is None:
            continue                       # an oo region: outside the domain
        finite.append((tuple(region), row, order))
        order += 1
    if not finite:
        raise ParseReject('no finite pieces (the function is identically +oo)')

    lines = list({hp.line for (region, _, _) in finite for hp in region})

    def in_domain(p):
        return any(all(hp.holds(p) for hp in region) for (region, _, _) in finite)

    def intersect(l1, l2):
        # a x + b y + c = 0 (Cramer, with the RHS as -c)
        a1, b1, c1 = l1
        a2, b2, c2 = l2
        det = a1 * b2 - a2 * b1
        if det == 0:
            return None
        return (Fraction(c2 * b1 - c1 * b2, det), Fraction(c1 * a2 - a1 * c2, det))

    corners = []
    seen = set()
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            p = intersect(lines[i], lines[j])
            if p is not None and in_domain(p) and p not in seen:
                seen.add(p)
                corners.append(p)

    def prim_dir(v):
        """The primitive integer direction of a rational direction vector (or its negative)."""
        x, y = v
        D = math.lcm(x.denominator, y.denominator)
        ix, iy = int(x * D), int(y * D)
        g = math.gcd(abs(ix), abs(iy))
        return (ix // g, iy // g)

    def line_dir(l):
        a, b, _ = l
        return prim_dir((-b, a))

    def on_line(p, l):
        a, b, c = l
        return a * p[0] + b * p[1] + c == 0

    def point_on_line(l):
        a, b, c = l
        if a != 0:
            return (Fraction(-c, a), Fraction(0))
        return (Fraction(0), Fraction(-c, b))

    V = []
    vidx = {}

    def vid(p):
        if p not in vidx:
            vidx[p] = len(V) + 1
            V.append(p)
        return vidx[p]

    def side_face(m, nrm):
        """The face owning the side of the current line through interior point m in normal
        direction nrm (0 = outside). Exactly one piece may claim a side of a valid partition."""
        cands = []
        for (region, row, o) in finite:
            if not all(hp.holds(m) for hp in region):
                continue
            ok = True
            for hp in region:
                if hp.line == line and not (hp.s * (hp.a * nrm[0] + hp.b * nrm[1]) > 0):
                    ok = False
                    break
            if ok:
                cands.append((o, row))
        if len(cands) > 1:
            raise ParseReject('two pieces own the same side of line %s (degenerate partition)' % (line,))
        return (cands[0][0] + 1, cands[0][1]) if cands else (0, None)

    E, F = [], []
    nf = len(finite)
    rows = [None] * nf
    for (region, row, o) in finite:
        rows[o] = row
    for line in sorted(lines, key=lambda t: (t[0], t[1], t[2])):
        d = line_dir(line)
        vs = sorted((p for p in corners if on_line(p, line)),
                    key=lambda p: p[0] * d[0] + p[1] * d[1])

        def add_edge(m, dirv, v1, v2, flag, l=line):
            i1, i2 = vid(v1), vid(v2)
            E.append([i1, i2, flag])
            dx, dy = dirv
            n = (-dy, dx)
            fL, rowL = side_face(m, n)
            fR, rowR = side_face(m, (-n[0], -n[1]))
            F.append([fL, fR])
            if fL and fR:
                # the two faces share this boundary: their values must agree there (continuity)
                m1 = m
                if _row_value(rowL, m1) != _row_value(rowR, m1):
                    raise ParseReject('two faces own adjacent sides of line %s with different '
                                      'values there (a jump the subdivision cannot place)' % (l,))

        if not vs:
            a0 = point_on_line(line)
            m = (a0[0] + d[0], a0[1] + d[1])
            if in_domain(m):
                add_edge(m, d, a0, m, 0)
            m = (a0[0] - d[0], a0[1] - d[1])
            if in_domain(m):
                add_edge(m, (-d[0], -d[1]), a0, m, 0)
        else:
            m = (vs[0][0] - d[0], vs[0][1] - d[1])
            if in_domain(m):
                add_edge(m, (-d[0], -d[1]), vs[0], (vs[0][0] - d[0], vs[0][1] - d[1]), 0)
            for i in range(len(vs) - 1):
                p1, p2 = vs[i], vs[i + 1]
                m = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
                if in_domain(m):
                    add_edge(m, prim_dir((p2[0] - p1[0], p2[1] - p1[1])), p1, p2, 1)
            m = (vs[-1][0] + d[0], vs[-1][1] + d[1])
            if in_domain(m):
                add_edge(m, d, vs[-1], (vs[-1][0] + d[0], vs[-1][1] + d[1]), 0)

    return Fixture(V, E, rows, F, '2d', finite)


def build_fixture_1d(f_expr):
    """A function on a LINE (a 1-D chain): the registry's 1-D-domain 2-D entries."""
    if not isinstance(f_expr, sp.Piecewise):
        raise ParseReject('a 1-D piecewise entry must be a Piecewise')
    out = []
    for val, cond in f_expr.args:
        for region in dnf_1d(cond):
            out.append((region, val))
    if not out:
        raise ParseReject('empty domain')
    line = None
    for (region, val) in out:
        for t in region:
            if t[0] == 'eq':
                l = t[1]
                if line is None:
                    line = l
                elif l != line:
                    raise ParseReject('1-D pieces on two different lines')
    if line is None:
        raise ParseReject('not a 1-D piecewise function')
    a, b, c = line
    d = (-b, a)
    g = math.gcd(abs(d[0]), abs(d[1]))
    d = (d[0] // g, d[1] // g)
    p0 = (Fraction(-c, a) if a != 0 else Fraction(0), Fraction(0) if a != 0 else Fraction(-c, b))
    # restrict each piece to t on the line: p(t) = p0 + t d
    intervals = []
    for (region, val) in out:
        lo, hi = None, None
        empty = False
        for t in region:
            if t[0] == 'eq':
                continue
            hp = t[1]
            k = hp.a * d[0] + hp.b * d[1]
            base = hp.a * p0[0] + hp.b * p0[1] + hp.c
            if k == 0:
                if hp.s * base < 0:
                    empty = True
                continue
            bound = Fraction(-base, k)
            if k > 0:                       # t >= bound
                if lo is None or bound > lo:
                    lo = bound
            else:                           # t <= bound
                if hi is None or bound < hi:
                    hi = bound
        if empty:
            continue
        row = piece_value_row(val)
        if row is None:
            continue
        intervals.append((lo, hi, row))
    if not intervals:
        raise ParseReject('no finite pieces')
    ends = sorted({t for (lo, hi, _) in intervals for t in (lo, hi) if t is not None})
    if len(ends) < 2 or any(hi is None for (_, hi, _) in intervals):
        raise ParseReject('unbounded or degenerate 1-D domain (not needed by the registry yet)')
    V = [p0 + (ends[0] * d[0], ends[0] * d[1]), p0 + (ends[-1] * d[0], ends[-1] * d[1])]
    rows = {tuple(row) for (_, _, row) in intervals}
    if len(rows) > 1:
        raise ParseReject('piecewise value on a 1-D chain beyond one quadratic (not needed yet)')
    row = list(rows.pop())
    return Fixture(V, [[1, 2, 1]], [row], [[0, 0]], 'chain', [([], row, 0)])


def dnf_1d(cond):
    """DNF for a 1-D condition: Eq atoms become ('eq', line), the rest half-planes."""
    if cond is sp.true or cond == True:            # noqa: E712
        return [()]
    if cond is sp.false or cond == False:          # noqa: E712
        return []
    if isinstance(cond, sp.Or):
        out = []
        for a in cond.args:
            out.extend(dnf_1d(a))
        seen, ded = set(), []
        for t in out:
            key = tuple(sorted(set(t), key=repr))
            if key not in seen:
                seen.add(key)
                ded.append(key)
        return ded
    atoms = cond.args if isinstance(cond, sp.And) else (cond,)
    out = [()]
    for a in atoms:
        if a is sp.true:                            # noqa: E712
            continue
        if a is sp.false:                           # noqa: E712
            return []
        if isinstance(a, sp.Eq):
            expr = sp.expand(a.lhs - a.rhs)
            A = expr.coeff(X)
            B = expr.coeff(Y)
            C = sp.simplify(expr - (A or 0) * X - (B or 0) * Y)
            A, B, C = int(_rat(A)), int(_rat(B)), int(_rat(C))
            if A == 0 and B == 0:
                if C == 0:
                    continue
                return []
            gg = math.gcd(abs(A), abs(B))
            A, B, C = A // gg, B // gg, C // gg
            if A < 0 or (A == 0 and B < 0):
                A, B, C = -A, -B, -C
            out = [cl + (('eq', (A, B, C)),) for cl in out]
            continue
        if isinstance(a, sp.Not):
            a = a.args[0]
            if not isinstance(a, RELATIONAL):
                raise ParseReject('non-boolean atom in a 1-D condition: %s' % a)
            s = -1
        else:
            s = +1
        rel = a
        expr = sp.expand(rel.lhs - rel.rhs)
        A = expr.coeff(X)
        B = expr.coeff(Y)
        C = sp.simplify(expr - (A or 0) * X - (B or 0) * Y)
        A, B, C = _rat(A), _rat(B), _rat(C)
        if A == 0 and B == 0:
            raise ParseReject('constant relation in a 1-D condition')
        if isinstance(rel, (sp.Ge, sp.Gt)):
            hp = HP(int(A), int(B), int(C), s)
        elif isinstance(rel, (sp.Le, sp.Lt)):
            hp = HP(int(A), int(B), int(C), -s)
        else:
            raise ParseReject('unsupported relation in a 1-D condition: %s' % a)
        out = [cl + (('hp', hp),) for cl in out]
    return out


def parse_entry(entry):
    """Entry -> Fixture. Raises ParseReject with a human reason."""
    f = entry.f
    if f is None:
        raise ParseReject('no primal function stored')
    if entry.dim != 2:
        raise ParseReject('not 2-D')
    has_eq = False
    for node in sp.preorder_traversal(f):
        if isinstance(node, sp.Eq) and node.has(X, Y):
            has_eq = True
            break
    if has_eq:
        return build_fixture_1d(f)
    return build_fixture_2d(f)


def fixture_check(entry, fx):
    """Verify the parser's own data: at sample points the first-match lookup on the
    parser's regions must equal the registry expression (both exact). Returns defects."""
    defects = []
    rng = random.Random(20261119)
    pts = []
    for _ in range(30):
        pts.append((Fraction(rng.randint(-6, 6), rng.choice([1, 2, 3, 7])),
                    Fraction(rng.randint(-6, 6), rng.choice([1, 2, 3, 7]))))
    for e in fx.E:
        p1 = fx.V[e[0] - 1]
        p2 = fx.V[e[1] - 1]
        pts.append(((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2))
    if fx.kind == 'chain':
        a, b = fx.V[0], fx.V[1]
        for tt in (Fraction(1, 3), Fraction(2, 3)):
            pts.append((a[0] + tt * (b[0] - a[0]), a[1] + tt * (b[1] - a[1])))
    for p in set(pts):
        ref = _registry_value_at(entry, p)
        got = fx.value_at(p)
        if ref is None and got is None:
            continue
        if (ref is None) != (got is None):
            defects.append('domain mismatch at %s: registry %s vs fixture %s' % (p, ref, got))
        elif ref is not None and sp.simplify(ref - got) != 0:
            defects.append('value mismatch at %s: registry %s vs fixture %s' % (p, ref, got))
    return defects


def _registry_value_at(entry, p):
    expr = entry.f
    v = expr.subs({X: p[0], Y: p[1]})
    v = sp.simplify(v)
    if v is sp.oo or v is sp.zoo:
        return None
    try:
        return sp.expand(v)
    except Exception:
        return v


# ==========================================================================================
# exact-precision evaluation of the engines' stored results
# ==========================================================================================

class Ambiguous(Exception):
    pass


def _mp_pt(pt):
    out = []
    for t in pt:
        if isinstance(t, (str, Fraction, int)):
            out.append(mp.mpf(t))
        else:
            out.append(mp.mpf(str(sp.N(t, 200))))
    return (out[0], out[1])


def eval_conic_mp(c, P):
    x, y = P
    return c[0] * x * x + c[1] * x * y + c[2] * y * y + c[3] * x + c[4] * y + c[5]


def _quad_roots_mp(co):
    """co = [c0, c1, c2] (ascending powers), two real roots, sorted ascending, at PREC bits."""
    c0, c1, c2 = (mp.mpf(v) for v in co)
    if c2 == 0:
        if c1 == 0:
            raise Ambiguous('linear fibre in a quadratic branch')
        return (-c0 / c1,)
    disc = c1 * c1 - 4 * c0 * c2
    if disc < 0:
        disc = mp.mpf(0)
    sq = mp.sqrt(disc)
    s = sq if c1 >= 0 else -sq
    r1 = -0.5 * (c1 + s) / c2
    r2 = -0.5 * (c1 - s) / c2
    return (r1, r2) if r1 <= r2 else (r2, r1)


def eval_pairline_mp(c, b, P):
    """Branch b (1 or 2) of the degenerate conic c, at PREC bits (mirror of rat_q.pairLineEval)."""
    from cca2 import rat_q
    info = rat_q.pairLineData(c)
    x, y = P
    if info['kind'] == 'parallel':
        r = _quad_roots_mp(list(info['quad']))[b - 1]
        f1, f2 = mp.mpf(info['form'][0]), mp.mpf(info['form'][1])
        return f1 * x + f2 * y - r
    ap = info['apex']
    p0 = (mp.mpf(ap[0]) / mp.mpf(ap[2]), mp.mpf(ap[1]) / mp.mpf(ap[2]))
    if info['vertical'] and b == 2:
        return x - p0[0]
    quad = list(info['quad'])
    if info['vertical']:
        t = -mp.mpf(quad[1]) / mp.mpf(quad[0])
    else:
        t = _quad_roots_mp(quad)[b - 1]
    return (y - p0[1]) - t * (x - p0[0])


def eval_face_mp(n, d, P):
    x, y = P
    if any(v != 0 for v in n[0:4]):
        raise Ambiguous('cubic face in a QuaCon (not produced by conj)')
    v = (mp.mpf(n[4]) / 2 * x * x + mp.mpf(n[5]) * x * y + mp.mpf(n[6]) / 2 * y * y
         + mp.mpf(n[7]) * x + mp.mpf(n[8]) * y + mp.mpf(n[9]))
    return v / mp.mpf(d)


def _face_holds(fc_rows, EcQ, P, g=None):
    """Do all FC rows of one face hold at P? (QuaCon rows [e, side, br]; AlgCrv rows carry a
    table index in position 3, handled by the caller for tables 1-3.)"""
    for row in fc_rows:
        e, side = row[0], row[1]
        br = row[2] if len(row) > 2 else 0
        tbl = row[3] if len(row) > 3 else 0
        if tbl != 0:
            raise Ambiguous('table-%d row in a QuaCon face' % tbl)
        c = EcQ[e - 1]
        q = eval_pairline_mp(c, br, P) if br else eval_conic_mp(c, P)
        scale = max(mp.mpf(1), mp.mpf(max(abs(v) for v in c)) * max(mp.mpf(1), P[0], P[1]) ** 2)
        if side == 0:
            if abs(q) > TOL_ONCURVE * scale:
                return False
        elif side * q < -TOL_ONCURVE * scale:
            return False
    return True


def eval_quacon_mp(g, P):
    """g: any object with .nf, .FC, .EcQ, .fN, .fD (Python QuaCon or the Rust JSON shape)."""
    held = []
    for k in range(g.nf):
        if _face_holds(g.FC[k], g.EcQ, P):
            held.append(k)
    if not held:
        return mp.inf
    vals = [eval_face_mp(g.fN[k], g.fD[k], P) for k in held]
    m = max(abs(v) for v in vals)
    if max(vals) - min(vals) > TOL_ONCURVE * m:
        raise Ambiguous('two faces hold at %s with different values' % (P,))
    return vals[0]


def algbranch_mp(FP, FB, P):
    """The named root branch of the implicit face polynomial FP (poly_n dict in x, y, z)."""
    x, y = P
    nz = max(k for (_, _, k) in FP)
    co = [mp.mpf(0)] * (nz + 1)
    for (i, j, k), c in FP.items():
        co[k] += mp.mpf(c) * (x ** i) * (y ** j)
    top = nz
    scale = max((abs(v) for v in co if v != 0), default=mp.mpf(1))
    while top > 0 and abs(co[top]) < mp.mpf(10) ** (-50) * scale:
        top -= 1
    if top == 0:
        raise Ambiguous('implicit face is constant at %s' % (P,))
    co = co[:top + 1]
    if top == 1:
        roots = (-co[0] / co[1],)
    elif top == 2:
        roots = _quad_roots_mp(co)
    else:
        roots_all = mp.polyroots(co)
        reals = [r.real for r in roots_all if abs(r.imag) < mp.mpf(10) ** (-30) * max(mp.mpf(1), abs(r))]
        if not reals:
            raise Ambiguous('no real root in the fibre at %s' % (P,))
        if len(reals) < len(roots_all) - 0.5:      # a genuine complex pair
            raise Ambiguous('complex roots in the fibre at %s' % (P,))
        roots = sorted(reals)
    if FB > 0:
        idx = FB - 1
    else:
        idx = len(roots) + FB
    if not (0 <= idx < len(roots)):
        raise Ambiguous('branch %d out of %d real roots at %s' % (FB, len(roots), P))
    return roots[idx]


def eval_algcrv_mp(g, P):
    best = None
    for k in range(g.nf):
        ok = True
        for row in g.FC[k]:
            idx, side, br = row[0], row[1], row[2] if len(row) > 2 else 0
            tbl = row[3] if len(row) > 3 else 0
            if tbl == 0:
                c = g.EcQ[idx - 1]
                q = eval_pairline_mp(c, br, P) if br else eval_conic_mp(c, P)
                scale = max(mp.mpf(1), mp.mpf(max(abs(v) for v in c)) * max(mp.mpf(1), P[0], P[1]) ** 2)
                if side == 0 and abs(q) > TOL_ONCURVE * scale:
                    ok = False
                    break
                if side != 0 and side * q < -TOL_ONCURVE * scale:
                    ok = False
                    break
            elif tbl == 1:
                a = g.EcA[idx - 1]
                v = algbranch_mp(a['FP'], a['FB'], P)
                q = eval_face_mp(a['num'], a['den'], P)
                if side * (q - v) < -TOL_ONCURVE * max(mp.mpf(1), abs(q), abs(v)):
                    ok = False
                    break
            elif tbl == 2:
                r = g.EcI[idx - 1]
                v1 = algbranch_mp(r['FP1'], r['FB1'], P)
                v2 = algbranch_mp(r['FP2'], r['FB2'], P)
                if side * (v1 - v2) < -TOL_ONCURVE * max(mp.mpf(1), abs(v1), abs(v2)):
                    ok = False
                    break
            elif tbl == 3:
                # the tbl-3 row layout differs from the others: [triIdx, i0 (1..3), side, 3]
                T = g.Tri[idx - 1]['T']
                i0 = row[1]
                side3 = row[2]
                j0 = i0 % 3 + 1
                o0 = 6 - i0 - j0
                Ti, Tj, To = T[i0 - 1], T[j0 - 1], T[o0 - 1]
                e = (Tj[0] - Ti[0], Tj[1] - Ti[1])
                cr = e[0] * (P[1] - Ti[1]) - e[1] * (P[0] - Ti[0])
                if side3 * cr < -TOL_ONCURVE * max(mp.mpf(1), abs(cr)):
                    ok = False
                    break
            else:
                raise Ambiguous('unknown table %s in an AlgCrv row' % tbl)
        if not ok:
            continue
        if g.FT[k] > 0:
            t = g.Tri[g.FT[k] - 1]
            v = P[0] * mp.mpf(t['s'][0]) + P[1] * mp.mpf(t['s'][1]) - mp.mpf(t['c'])
        elif g.FP[k]:
            v = algbranch_mp(g.FP[k], g.FB[k], P)
        else:
            v = eval_face_mp(g.fN[k], g.fD[k], P)
        if best is None or v > best:
            best = v
    if best is None:
        return mp.inf
    return best


def eval_engine_mp(result, P):
    """Exact-precision value of an engine result at P, or Ambiguous when we cannot be sure.

    Only a QuaCon has a clean exact structure (rational conics, rational face values) that
    we can evaluate reliably at high precision. Every other result type (AlgCrv, EnvelopeMin,
    ...) is routed to the engine's own documented double-precision eval (raise Ambiguous ->
    the caller falls back to eval_engine_double, graded at TOL_DOUBLE)."""
    from cca2.qua_con import QuaCon
    if isinstance(result, QuaCon):
        return eval_quacon_mp(result, P)
    raise Ambiguous('non-QuaCon result %s: using the engine double-precision eval'
                    % type(result).__name__)


def eval_engine_double(result, P):
    """The engine's own double-precision value (its documented interface)."""
    v = result.eval([[float(P[0]), float(P[1])]])
    val = v[0] if isinstance(v, tuple) else v
    val = val[0] if isinstance(val, list) else val
    if isinstance(val, float) and math.isinf(val):
        return mp.inf if val > 0 else -mp.inf
    if val is None:
        return mp.inf
    return mp.mpf(val)


# ==========================================================================================
# registry value at a test point (exact sympy, then one rounding)
# ==========================================================================================

def registry_value_mp(expr, pt_sympy, vars_list=None):
    """The stored registry value at the exact point, one rounding to PREC bits. The stored
    expression is written in the entry's own symbols (x, y for the primal; s1, s2 for f*);
    the point coordinates are mapped onto them in order."""
    if expr is None:
        return None
    sub = {}
    vl = [sp.sympify(v) for v in (vars_list or ['x', 'y'])]
    syms = expr.free_symbols
    for v, p_ in zip(vl, pt_sympy):
        if v in syms:
            sub[v] = p_
    v = expr.subs(sub)
    v = sp.simplify(v)
    if v is sp.oo or v is sp.zoo or (v.is_Float and v.is_infinite):
        return mp.inf
    s = str(sp.N(v, 200))
    return mp.mpf(s)


# ==========================================================================================
# engines
# ==========================================================================================

def run_python(entry, fx, ops):
    """ops: subset of {'conj', 'biconj'} -> {op: ('refuse', reason) | (None, result, t)}."""
    from cca2.qua_pol import QuaPol
    from cca2.errors import CCA2Error
    out = {}
    try:
        f = QuaPol([list(r) for r in fx.V], [list(r) for r in fx.E],
                   [list(r) for r in fx.f], [list(r) for r in fx.F])
    except CCA2Error as e:
        return {op: ('refuse', 'QuaPol: %s' % e.identifier) for op in ops}
    except Exception as e:      # noqa: BLE001
        return {op: ('refuse', 'QuaPol: %s: %s' % (type(e).__name__, e)) for op in ops}
    for op in ops:
        t0 = time.time()
        try:
            res = f.conj() if op == 'conj' else f.biconj()
        except CCA2Error as e:
            out[op] = ('refuse', '%s' % e.identifier, time.time() - t0)
            continue
        except Exception as e:      # noqa: BLE001
            out[op] = ('refuse', '%s: %s' % (type(e).__name__, e), time.time() - t0)
            continue
        out[op] = (None, res, time.time() - t0)
    return out


def open_sanity(entry, op, res, vals, pts_sympy):
    """Necessary-condition checks for an engine value on an OPEN registry entry."""
    notes = []
    ok = True
    for (label, v) in vals:
        if v is mp.inf or v == mp.inf:
            notes.append('%s: +inf' % label)
            continue
        pt = pts_sympy[label]          # the exact original point (rational, algebraic, or decimal)
        f_here = _registry_value_at(entry, pt)
        if f_here is None:
            notes.append('%s: engine %s (point outside the primal domain)' % (label, mp.nstr(v, 12)))
            continue
        fv = mp.mpf(str(sp.N(sp.expand(f_here), 120)))
        if op == 'biconj' and v > fv + mp.mpf(1e-20):      # decimal-literal points carry ~36 digits
            ok = False
            notes.append('%s: VIOLATION f**=%s > f=%s' % (label, mp.nstr(v, 15), mp.nstr(fv, 15)))
        else:
            notes.append('%s: value %s (f=%s)' % (label, mp.nstr(v, 15), mp.nstr(fv, 15)))
    if op == 'conj':
        notes.append('NOTE: no exact sup oracle for open conj entries (none exist in the registry)')
    return ok, notes


def compare_op(entry, op, py_res, pts_mp, pts_sympy):
    """py_res = (status, result, t) -> (status, detail, t, vals)."""
    status, res, t = py_res
    if status == 'refuse':
        return status, 'refused: %s' % res, t, None
    stored = entry.fstar if op == 'conj' else entry.cof
    svar = entry.vars_fstar if op == 'conj' else entry.vars_f
    if stored is None:
        vals = []
        for (label, P) in pts_mp:
            try:
                v = eval_engine_mp(res, P)
            except Ambiguous as e:
                try:
                    v = eval_engine_double(res, P)
                except Exception as e2:    # noqa: BLE001
                    v = mp.inf
                    e = Ambiguous('%s; double eval: %s' % (e, e2))
            vals.append((label, v))
        ok, notes = open_sanity(entry, op, res, vals, pts_sympy)
        return ('agree' if ok else 'disagree'), 'open entry; ' + '; '.join(notes), t, vals
    bad = []
    vals = []
    double_pts = 0
    for (label, P) in pts_mp:
        b = registry_value_mp(stored, pts_sympy[label], svar)
        try:
            a = eval_engine_mp(res, P)
            tol = TOL_EXACT
        except Ambiguous:
            try:
                a = eval_engine_double(res, P)
            except Exception as e:      # noqa: BLE001
                bad.append('%s: engine eval failed (%s)' % (label, e))
                continue
            tol = TOL_DOUBLE
            double_pts += 1
        a_inf = (a is mp.inf or a == mp.inf)
        b_inf = (b is mp.inf or b == mp.inf)
        if a_inf and b_inf:
            vals.append((label, mp.inf, mp.inf))
            continue
        if a_inf != b_inf:
            bad.append('%s: engine %s vs registry %s' % (label, mp.nstr(a, 20), mp.nstr(b, 20)))
            continue
        if abs(a - b) > tol * max(mp.mpf(1), abs(a), abs(b)):
            bad.append('%s: engine %s vs registry %s (rel gap %s)'
                       % (label, mp.nstr(a, 30), mp.nstr(b, 30),
                          mp.nstr(abs(a - b) / max(mp.mpf(1), abs(b)), 6)))
            continue
        vals.append((label, a, b))
    if bad:
        return 'disagree', '; '.join(bad), t, vals
    if not vals:
        return 'refuse', 'no test points stored', t, None
    note = ' (double protocol on %d point(s))' % double_pts if double_pts else ''
    return 'agree', 'matched at %d test point(s)%s' % (len(vals), note), t, vals


class _RustQuaCon:
    """The Rust golden output's QuaCon fields in the shape eval_quacon_mp reads."""
    def __init__(self, j):
        self.nf = j['nf']
        self.FC = [r for r in j['FC']]
        self.EcQ = [r for r in j['EcQ']]
        self.fN = [r for r in j['fN']]
        self.fD = list(j['fD'])


def run_rust(entry, fx, ops, golden, pts_sympy):
    """Drive the Rust engine through its golden binary (the parity-test fixture JSON).
    -> {op: (status, detail, t_or_None, vals)}."""
    pts = [(p[0], p[1]) for (_, p) in pts_sympy_items(pts_sympy)]
    recs = []
    for op in ops:
        rec = {'fn': 'conjQ' if op == 'conj' else 'biconjQ', 'in': fx.golden_in(), 'out': {}}
        rec['out']['evalX'] = [[float(_mp_pt(pt)[0]), float(_mp_pt(pt)[1])] for pt in pts]
        recs.append((op, rec))
    with tempfile.NamedTemporaryFile('w', suffix='.jsonl', delete=False) as tf:
        for (_, r) in recs:
            tf.write(json.dumps(r) + '\n')
        tmp = tf.name
    try:
        p = subprocess.run([golden, tmp], capture_output=True, text=True, timeout=7200)
        outs = [l for l in p.stdout.strip().split('\n') if l]
    except Exception as e:      # noqa: BLE001
        return {op: ('refuse', 'golden binary: %s' % e) for op in ops}
    finally:
        os.unlink(tmp)
    out = {}
    for (op, _), line in zip(recs, outs):
        try:
            j = json.loads(line)
        except Exception as e:      # noqa: BLE001
            out[op] = ('refuse', 'unparsable golden output: %s (%s)' % (line[:120], e))
            continue
        if 'err' in j:
            out[op] = ('refuse', 'rust: %s' % j['err'])
            continue
        stored = entry.fstar if op == 'conj' else entry.cof
        ev = j.get('evalV')
        if ev is None:
            out[op] = ('refuse', 'rust: %s' % j.get('evalErr', 'no evalV'))
            continue
        svar = entry.vars_fstar if op == 'conj' else entry.vars_f
        stored_vals = [registry_value_mp(stored, pt, svar) for pt in pts] if stored is not None else None
        kind = j.get('kind', '')
        # exact path for a Rust QuaCon (same exact data shape as the Python one)
        rust_mp = None
        if kind == 'QuaCon':
            try:
                rust_mp = _RustQuaCon(j)
            except Exception:      # noqa: BLE001
                rust_mp = None
        if stored is None:
            out[op] = ('open', 'rust kind=%s; values %s' % (kind, [mp.nstr(_as_mp(v), 12) for v in ev]),
                       None, [(l, _as_mp(v)) for (l, v) in zip(pts, ev)])
            continue
        bad = []
        vals = []
        for i, (label, pt) in enumerate(zip(pts_labels(pts_sympy), pts)):
            b = stored_vals[i]
            if rust_mp is not None:
                P = _mp_pt(pt)
                try:
                    a = eval_quacon_mp(rust_mp, P)
                    tol = TOL_EXACT
                except Ambiguous:
                    a = _as_mp(ev[i])
                    tol = TOL_DOUBLE
            else:
                a = _as_mp(ev[i])
                tol = TOL_DOUBLE
            a_inf = a is mp.inf
            b_inf = b is mp.inf
            if a_inf and b_inf:
                vals.append((label, mp.inf, mp.inf))
                continue
            if a_inf != b_inf:
                bad.append('%s: rust %s vs registry +inf' % (label, a))
                continue
            if abs(a - b) > tol * max(mp.mpf(1), abs(a), abs(b)):
                bad.append('%s: rust %s vs registry %s (rel gap %s)'
                           % (label, mp.nstr(a, 30), mp.nstr(b, 30),
                              mp.nstr(abs(a - b) / max(mp.mpf(1), abs(b)), 6)))
                continue
            vals.append((label, a, b))
        if bad:
            out[op] = ('disagree', '; '.join(bad), None, vals)
        else:
            out[op] = ('agree', 'matched at %d test point(s)%s'
                       % (len(ev), '' if rust_mp else ' (double protocol)'), None, vals)
    return out


def pts_sympy_items(d):
    return [(k, d[k]) for k in d]


def pts_labels(d):
    return list(d.keys())


def _as_mp(v):
    if isinstance(v, str):
        if v == 'inf':
            return mp.inf
        if v == '-inf':
            return -mp.inf
        return mp.mpf(v)
    if isinstance(v, float):
        if math.isinf(v):
            return mp.inf if v > 0 else -mp.inf
        return mp.mpf(v)
    return mp.mpf(v)


# ==========================================================================================
# per-entry driver
# ==========================================================================================

def process_entry(entry, engines, golden):
    rec = {'slug': entry.slug, 'dim': entry.dim,
           'status_conj': entry.status_conj, 'status_coenv': entry.status_coenv,
           'f_star_stored': entry.fstar is not None, 'co_f_stored': entry.cof is not None,
           'eval_points': [l for (l, p) in entry.pts()]}
    pts = entry.pts()
    pts_mp = [(l, _mp_pt(p)) for (l, p) in pts]
    pts_sympy = {l: p for (l, p) in pts}
    rec['pts'] = {l: [str(a), str(b)] for (l, (a, b)) in pts}

    if entry.dim != 2:
        rec['parse'] = 'skip (not 2-D)'
        return rec
    try:
        fx = parse_entry(entry)
    except ParseReject as e:
        rec['parse'] = 'not PLQ / not representable: %s' % e
        rec['note'] = ('non-PLQ 2-D entry or a shape QuaPol cannot take (full strips, curved '
                       'domains, non-quadratic values): out of scope -- "for every 2-D PLQ '
                       'entry the engines accept"')
        for eng in engines:
            rec[eng] = {'conj': {'status': 'refuse', 'detail': 'input not representable by QuaPol'},
                        'biconj': {'status': 'refuse', 'detail': 'input not representable by QuaPol'}}
        return rec
    defects = fixture_check(entry, fx)
    if defects:
        rec['parse'] = 'parser defect (refusing to run the engines on it): %s' % '; '.join(defects[:3])
        for eng in engines:
            rec[eng] = {'conj': {'status': 'refuse', 'detail': 'parser defect (input rejected)'},
                        'biconj': {'status': 'refuse', 'detail': 'parser defect (input rejected)'}}
        return rec
    rec['parse'] = 'ok (%s, %d faces, %d edges, %d vertices)' % (fx.kind, len(fx.f), len(fx.E), len(fx.V))
    ops = ['conj', 'biconj']
    if 'py' in engines:
        py = run_python(entry, fx, ops)
        pyres = {}
        for op in ops:
            st, det, t = py[op]
            if st is None:
                st, det, t, vals = compare_op(entry, op, py[op], pts_mp, pts_sympy)
            pyres[op] = (st, det, t)
        rec['py'] = {op: {'status': pyres[op][0], 'detail': pyres[op][1],
                          't': round(pyres[op][2], 3)} for op in ops}
    if 'rs' in engines:
        if golden is None:
            rec['rs'] = {op: {'status': 'refuse', 'detail': 'Rust golden binary not found'} for op in ops}
        else:
            rs = run_rust(entry, fx, ops, golden, pts_sympy)
            rec['rs'] = {op: {'status': rs[op][0], 'detail': rs[op][1]} for op in ops}
    return rec


def summary(records):
    n2 = [r for r in records if r['dim'] == 2]
    accepted = [r for r in n2 if r.get('parse', '').startswith('ok')]
    notplq = [r for r in n2 if not r.get('parse', '').startswith('ok')]

    def count(eng, op, status):
        return sum(1 for r in accepted if r.get(eng, {}).get(op, {}).get('status') == status)

    s = {'2d_entries': len(n2), 'accepted': len(accepted), 'not_plq_or_refused_input': len(notplq)}
    for eng in ('py', 'rs'):
        for op in ('conj', 'biconj'):
            s['%s_%s' % (eng, op)] = {
                st: count(eng, op, st) for st in ('agree', 'refuse', 'disagree', 'open')}
    return s


def write_results_md(path, ver, s, records):
    marker = '## Cross-check against convexdb %s' % ver
    lines = []
    lines.append(marker)
    lines.append('')
    lines.append('The external-oracle pass for gap G3 (outline section 6): every 2-D PLQ entry the')
    lines.append('engines accept, compared against the registry\'s stored `f_star` / `co_f` at the')
    lines.append('registry\'s own test points. A QuaCon / AlgCrv result (Python or Rust) is')
    lines.append('evaluated from its EXACT stored data at 512-bit precision, the registry value is')
    lines.append('exact sympy rounded once to the same precision, and agreement is demanded to')
    lines.append('relative 1e-60. An EnvelopeMin is on the 1e-6 double protocol (the engine\'s own')
    lines.append('documented evaluation). Registry: the sibling `convexdb` checkout (version %s),' % ver)
    lines.append('`export/convexdb.json`; read by relative folder, no URL. Per-entry detail:')
    lines.append('`bench/convexdb_cross.json`; script `bench/cross_convexdb.py`.')
    lines.append('')
    lines.append('| | conj | biconj |')
    lines.append('|---|---|---|')

    def cell(eng, op):
        d = s.get('%s_%s' % (eng, op), {})
        parts = ['%s %d' % (k, d[k]) for k in ('agree', 'refuse', 'disagree', 'open') if d.get(k)]
        return ' · '.join(parts) if parts else '--'

    lines.append('| Python | %s | %s |' % (cell('py', 'conj'), cell('py', 'biconj')))
    lines.append('| Rust | %s | %s |' % (cell('rs', 'conj'), cell('rs', 'biconj')))
    lines.append('')
    lines.append('2-D entries in the registry: %d. Accepted by the engines: %d. Not PLQ / not'
                 % (s['2d_entries'], s['accepted']))
    lines.append('representable (expected refusals, registry content): %d.' % s['not_plq_or_refused_input'])
    lines.append('')
    lines.append('| entry | parse | py conj | py biconj | rs conj | rs biconj |')
    lines.append('|---|---|---|---|---|---|')
    for r in records:
        if r['dim'] != 2:
            continue

        def c(eng, op):
            d = r.get(eng, {}).get(op, {})
            st = d.get('status', '--')
            det = (d.get('detail') or '').replace('|', '/')
            det = det if len(det) < 60 else det[:57] + '...'
            return '%s — %s' % (st, det) if st in ('disagree', 'refuse') and det else st

        parse = (r.get('parse') or '').replace('|', '/')
        parse = parse if len(parse) < 48 else parse[:45] + '...'
        lines.append('| %s | %s | %s | %s | %s | %s |'
                     % (r['slug'], parse, c('py', 'conj'), c('py', 'biconj'), c('rs', 'conj'), c('rs', 'biconj')))
    lines.append('')
    lines.append('_Recorded %s; shared machine, no timing conclusions._' % time.strftime('%Y-%m-%d %H:%M'))
    new_section = '\n'.join(lines) + '\n'
    old = ''
    if os.path.isfile(path):
        old = open(path).read()
    if marker in old:
        i = old.index(marker)
        j = old.find('\n## ', i + 1)
        old = old[:i] + new_section + (old[j + 1:] if j != -1 else '\n')
    else:
        old = old.rstrip('\n') + '\n\n' + new_section
    with open(path, 'w') as f:
        f.write(old)


if __name__ == '__main__':
    sys.exit(main())

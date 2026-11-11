"""io_matlab -- read the golden traces MATLAB wrote (python/tools/export_golden.sh) and compare
a Python result against one.

A trace is JSON lines; each line is one call of conjQ or biconjQ: the input QuaPol (the caller's
doubles as %.17g, plus MATLAB's exact reading of them), and either the exact output or the
identifier of the named refusal. Index values are 1-based in both worlds, so nothing is shifted.
"""
import json
import math

from .qua_pol import QuaPol


def _num(x):
    if isinstance(x, str):
        return {'inf': math.inf, '-inf': -math.inf, 'nan': math.nan}[x]
    return x


def _mat(A):
    if A and not isinstance(A[0], list):
        A = [A]
    return [[_num(x) for x in r] for r in A]


def load(path):
    """All records of a trace, deduplicated on (function, input)."""
    out = []
    seen = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            key = (rec['fn'], json.dumps(rec['in'], sort_keys=True))
            if key in seen:
                continue
            seen.add(key)
            out.append(rec)
    return out


def quaPolFrom(j):
    """Build the Python QuaPol from a record's input exactly as the MATLAB caller did."""
    f = _mat(j['f'])
    if j['nv'] == 0 and not j['V']:
        P = QuaPol(f[0] if len(f) == 1 else f)
    else:
        P = QuaPol(_mat(j['V']), _mat(j['E']) if j['E'] else [], f, _mat(j['F']) if j['F'] else [])
    if not j['exact']:
        # MATLAB's object was INEXACT -- built by the float pipeline, or stripped by a test --
        # which the trace cannot spell otherwise: the numbers themselves may look rational.
        P.fN = P.fD = P.VN = P.VD = None
    if 'fIsConvex' in j:
        P.fIsConvex = bool(j['fIsConvex'])
    return P


def inputMismatch(P, j):
    """'' when the Python QuaPol read the input as MATLAB did, else what differs."""
    if P.isExact() != bool(j['exact']):
        return 'isExact %s vs MATLAB %s' % (P.isExact(), j['exact'])
    if not j['exact']:
        return ''
    if (P.nv, P.ne, P.nf) != (j['nv'], j['ne'], j['nf']):
        return 'counts %s vs MATLAB %s' % ((P.nv, P.ne, P.nf), (j['nv'], j['ne'], j['nf']))
    if P.fN != _mat(j['fN']) or P.fD != _mat(j['fD'])[0]:
        return 'fN/fD differ'
    if j['VN'] and (P.VN != _mat(j['VN']) or P.VD != _mat(j['VD'])[0]):
        return 'VN/VD differ'
    return ''


def quaConFields(g):
    return {'nv': g.nv, 'ne': g.ne, 'nf': g.nf, 'Vname': g.Vname, 'E': g.E, 'EcQ': g.EcQ,
            'fN': g.fN, 'fD': g.fD, 'F': g.F, 'FC': g.FC}


def quaConMismatch(g, j):
    """'' when the QuaCon g equals MATLAB's field for field, else the first field that differs."""
    mine = quaConFields(g)
    theirs = {'nv': j['nv'], 'ne': j['ne'], 'nf': j['nf'],
              'Vname': _mat(j['Vname']) if j['Vname'] else [],
              'E': _mat(j['E']) if j['E'] else [],
              'EcQ': _mat(j['EcQ']) if j['EcQ'] else [],
              'fN': _mat(j['fN']) if j['fN'] else [],
              'fD': _mat(j['fD'])[0] if j['fD'] else [],
              'F': _mat(j['F']) if j['F'] else [],
              'FC': [(_mat(c) if c else []) for c in j['FC']]}
    for k in ('nf', 'ne', 'nv', 'fN', 'fD', 'EcQ', 'FC', 'Vname', 'E', 'F'):
        if mine[k] != theirs[k]:
            return '%s: python %r vs MATLAB %r' % (k, _short(mine[k]), _short(theirs[k]))
    return ''


def _short(x, n=400):
    s = repr(x)
    return s if len(s) <= n else s[:n] + '...'


def evalMismatch(g, j, tol=1e-9):
    """'' when g.eval reproduces MATLAB's sampled values (and faces, when recorded)."""
    if 'evalX' not in j:
        return ''
    X = _mat(j['evalX'])
    V = _mat(j['evalV'])[0]
    res = g.eval(X)
    val = res[0] if isinstance(res, tuple) else res
    for i, (a, b) in enumerate(zip(val, V)):
        if math.isinf(a) or math.isinf(b) or math.isnan(a) or math.isnan(b):
            if not (a == b or (math.isnan(a) and math.isnan(b))):
                return 'eval at %r: python %r vs MATLAB %r' % (X[i], a, b)
        elif abs(a - b) > tol * max(1.0, abs(b)):
            return 'eval at %r: python %r vs MATLAB %r' % (X[i], a, b)
    return ''


# ------------------------------------------------------------------------------------------------
# generic traces (python/tools/export_trace.sh): arguments and results of intermediate functions
# ------------------------------------------------------------------------------------------------
def dec(x):
    """Decode one serialised MATLAB value, STRICTLY by shape:
      numeric m x n      -> list of m rows (ints where MATLAB's value was integral); 0 x k -> []
      numeric N-D        -> {'nd': size, 'v': column-major values}
      cell               -> flat list (vector) or list of rows
      struct 1x1         -> dict;  struct array -> list of dicts
      object             -> dict of its properties plus 'cls'
      bigZ               -> like numeric, with Python ints;  BigInteger -> int
    Use sc / vec below to squeeze a 1x1 or a vector."""
    if isinstance(x, list):
        return [dec(v) for v in x]
    if not isinstance(x, dict):
        return x
    t = x.get('t')
    if t in ('n', 'z'):
        sz = x['sz']
        v = [int(a) if t == 'z' else _num(a) for a in x['v']]
        if len(sz) > 2:
            return {'nd': sz, 'v': v}
        m, n = sz
        if m == 0 or n == 0:
            return []
        return [[v[j * m + i] for j in range(n)] for i in range(m)]
    if t == 'j':
        return int(x['v'])
    if t == 'c':
        sz = x['sz']
        v = [dec(a) for a in x['v']]
        if len(sz) == 2 and sz[0] > 1 and sz[1] > 1:
            m, n = sz
            return [[v[j * m + i] for j in range(n)] for i in range(m)]
        return v
    if t == 's':
        v = [{k: dec(a) for k, a in e.items()} for e in x['v']]
        if x['sz'] == [1, 1]:
            return v[0]
        return v
    if t == 'o':
        d = {k: dec(a) for k, a in x['p'].items()}
        d['cls'] = x['cls']
        return d
    if t == 'f':
        return x['v']
    return x


def sc(M):
    """The scalar of a 1x1 matrix (or the value itself)."""
    while isinstance(M, list):
        if not M:
            return None
        M = M[0]
    return M


def vec(M):
    """A row or column matrix as a flat list."""
    if not M:
        return []
    if not isinstance(M[0], list):
        return list(M)
    if len(M) == 1:
        return list(M[0])
    return [r[0] for r in M]


def loadTrace(path, fn=None):
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if fn is not None and rec.get('fn') != fn:
                continue
            if 'serr' in rec:
                continue
            out.append(rec)
    return out


def quaConFromTrace(d):
    """A QuaCon from its traced property dict."""
    from .qua_con import QuaCon
    FC = [[list(r) for r in (c or [])] for c in d['FC']]
    return QuaCon(d['Vname'] or [], d['E'] or [], d['EcQ'] or [], d['fN'] or [], vec(d['fD']),
                  d['F'] or [], FC)


def quaPolFromTrace(d):
    f = d['f']
    if d['nv'] == [[0]] or not d['V']:
        P = QuaPol(f[0] if len(f) == 1 else f)
    else:
        P = QuaPol(d['V'], d['E'] or [], f, d['F'] or [])
    if not d.get('fN'):
        P.fN = P.fD = P.VN = P.VD = None
    if d.get('fIsConvex'):
        P.fIsConvex = bool(sc(d['fIsConvex']))
    return P


def polyFromTrace(x, nv):
    """A traced MATLAB polynomial (N-D array, 2-D array, bigZ, or ratQ's {sz, c} struct) as a
    poly_n dict in nv variables."""
    from . import poly_n
    if isinstance(x, dict) and 'nd' in x:
        return poly_n.fromND(x['nd'], x['v'], nv)
    if isinstance(x, dict) and 'sz' in x and 'c' in x:          # ratQ.isBigPoly
        sz = vec(x['sz'])
        v = [0 if (c is None or c == []) else int(c) for c in x['c']]
        return poly_n.fromND(sz, v, nv)
    if isinstance(x, dict) and x.get('cls') == 'bigZ':
        return polyFromTrace(x['c'], nv)
    if isinstance(x, list):
        if not x:
            return {}
        if x and isinstance(x[0], list):
            m, n = len(x), len(x[0])
            return poly_n.fromND([m, n], [x[i][j] for j in range(n) for i in range(m)], nv)
    raise ValueError('unrecognised polynomial in trace: %r' % (x,))

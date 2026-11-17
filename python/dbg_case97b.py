"""Dump caseD cells for conj[97] before and after assembleQuaConCells to find the dropped cell."""
import json
from cca2 import io_matlab
from cca2 import conj_q as _q
from cca2 import rat_q
from cca2.assemble_qua_con_cells import assembleQuaConCells

raw = [json.loads(l) for l in open('/home/ylucet/AI/CCA2/python/golden/conj.jsonl') if l.strip()]
rec = raw[97]
P = io_matlab.quaPolFrom(rec['in'])

# Reproduce caseDBoundaryMax's cell construction manually (mirror conj_q.py caseDBoundaryMax)
fN, fD = P.faceQ(1)
Hn = [[fN[4], fN[5]], [fN[5], fN[6]]]
Ln = [fN[7], fN[8]]
kn = fN[9]
Vi, vd = rat_q.combineDen(P.VN, P.VD)
sh = _q.pieceShape(P, 1, Vi, vd)
m = len(sh['vs'])

from cca2 import cell_interior as _ci

denV = 2 * fD * vd ** 2
qv = []
numV = []
for vi in sh['vs']:
    v = Vi[vi - 1]
    q = rat_q.quad(v, Hn, v) + 2 * vd * rat_q.dot(Ln, v) + 2 * vd ** 2 * kn
    qv.append(q)
    numV.append([0, 0, 0, 0, 0, 0, 0, 2 * fD * vd * v[0], 2 * fD * vd * v[1], -q])

# Edges with alpha > 0
Ed = []
for e in sh['ed']:
    d = e['d']
    alpha = rat_q.quad(d, Hn, d)
    if alpha <= 0:
        continue
    a = Vi[e['a'] - 1]
    ga = [x + vd * y for x, y in zip(rat_q.mv(Hn, a), Ln)]
    qa = rat_q.quad(a, Hn, a) + 2 * vd * rat_q.dot(Ln, a) + 2 * vd ** 2 * kn
    Td = [fD * vd * d[0], fD * vd * d[1], -rat_q.dot(ga, d)]
    lin = [2 * alpha * fD * vd * a[0], 2 * alpha * fD * vd * a[1]]
    Ed.append({'Td': Td, 'alpha': alpha,
               'num': [0, 0, 0, 0, 2 * Td[0] ** 2, 2 * Td[0] * Td[1], 2 * Td[1] ** 2,
                       lin[0] + 2 * Td[0] * Td[2], lin[1] + 2 * Td[1] * Td[2],
                       Td[2] ** 2 - alpha * qa],
               'den': 2 * alpha * fD * vd ** 2, 'isRay': e['isRay']})
print("m=%d Ed=%d" % (m, len(Ed)))
for ei, ed in enumerate(Ed):
    print("  edge[%d] num=%s den=%s alpha=%s" % (ei, ed['num'], ed['den'], ed['alpha']))

# cells0: vertex-leading cells
cells0 = []
for i in range(m):
    rows = []
    vI = Vi[sh['vs'][i] - 1]
    for j in range(m):
        if j == i:
            continue
        vJ = Vi[sh['vs'][j] - 1]
        r = [2 * fD * vd * (vI[0] - vJ[0]), 2 * fD * vd * (vI[1] - vJ[1]), -(qv[i] - qv[j])]
        rows.append(r)
    if rat_q.feasible2(rows, True):
        cells0.append({'lin': rows, 'vi': i, 'ei': []})
    else:
        print("  vertex cell i=%d INFEASIBLE (dropped)" % i)
print("cells0 (vertex-leading feasible) = %d" % len(cells0))
for c in cells0:
    print("   vi=%d ei=%s" % (c['vi'], c['ei']))

# Now apply edge regimes (like conj_q.py)
for e, ed in enumerate(Ed):
    Td = ed['Td']
    alpha = ed['alpha']
    neg = [-Td[0], -Td[1], -Td[2]]
    if ed['isRay']:
        regimes = [[neg], [Td]]
        avail = [False, True]
    else:
        regimes = [[neg], [Td, [-Td[0], -Td[1], alpha - Td[2]]], [[Td[0], Td[1], Td[2] - alpha]]]
        avail = [False, True, False]
    nxt = []
    for c in cells0:
        for r in range(len(regimes)):
            lin = c['lin'] + regimes[r]
            if not rat_q.feasible2(lin, True):
                continue
            ei = c['ei'] + [e] if avail[r] else c['ei']
            nxt.append({'lin': lin, 'vi': c['vi'], 'ei': ei})
    cells0 = nxt
print("after edge regimes: cells0 = %d" % len(cells0))
for c in cells0:
    print("   vi=%d ei=%s" % (c['vi'], c['ei']))

# Build the out cells
out = []
for c in cells0:
    base, ok = _q.toCon(c['lin'])
    if not ok:
        continue
    cand = [{'num': numV[c['vi']], 'den': denV}]
    for e in c['ei']:
        cand.append({'num': Ed[e]['num'], 'den': Ed[e]['den']})
    for i in range(len(cand)):
        con = [list(r) for r in base]
        for j in range(len(cand)):
            if j == i:
                continue
            dn, _ = rat_q.sub(cand[i]['num'], cand[i]['den'], cand[j]['num'], cand[j]['den'])
            if all(x == 0 for x in dn):
                continue
            Rr = [dn[4], 2 * dn[5], dn[6], 2 * dn[7], 2 * dn[8], 2 * dn[9]]
            con.append(rat_q.conic(Rr) + [_q.sgnOf(Rr, +1)])
        out.append({'num': cand[i]['num'], 'den': cand[i]['den'], 'con': con})
print("out cells (pre-assemble) = %d" % len(out))
from collections import Counter
cn = Counter(tuple(c['num']) for c in out)
for v, c in cn.items():
    print("   %s x%d" % (v, c))

# Now run assembleQuaConCells and compare
res = assembleQuaConCells(out)
print("after assemble: nf=%d" % res.nf)
cn2 = Counter(tuple(r) for r in res.fN)
print("MATLAB extra cell (0, 0, 0, 0, 392, 224, 128, 4548, 4640, 40305):")
print("   in out cells:", (0, 0, 0, 0, 392, 224, 128, 4548, 4640, 40305) in cn)
print("   in result:", (0, 0, 0, 0, 392, 224, 128, 4548, 4640, 40305) in cn2)

"""Investigate conj[97] (raw line): pentagon, indefinite H -> caseD. Dump the cell pipeline."""
import json
from cca2 import io_matlab
from cca2 import conj_q as _q
from cca2.conj_q import conjQ, conjPieceQ, caseDBoundaryMax, pieceShape
from cca2 import rat_q
from cca2.assemble_qua_con_cells import assembleQuaConCells

raw = [json.loads(l) for l in open('/home/ylucet/AI/CCA2/python/golden/conj.jsonl') if l.strip()]
rec = raw[97]
P = io_matlab.quaPolFrom(rec['in'])

fN, fD = P.faceQ(1)
Hn = [[fN[4], fN[5]], [fN[5], fN[6]]]
D = rat_q.detExact(Hn)
Vi, vd = rat_q.combineDen(P.VN, P.VD)
sh = pieceShape(P, 1, Vi, vd)
print("convex=%s bounded=%s nv=%d" % (sh['convex'], sh['bounded'], len(sh['vs'])))
print("Hn=%s det=%s  H[0][0]>0=%s det>0=%s" % (Hn, D, Hn[0][0] > 0, D > 0))
print("vs=%s" % sh['vs'])
print("ed=%s" % [(e['a'], e['b'], e['d']) for e in sh['ed']])
print("hp=%s" % sh['hp'])
print("rays=%s" % sh['rays'])

# Is the pentagon convex or not? If convex -> caseD directly. If not -> earClip split.
if not sh['convex']:
    print("NON-CONVEX: earClip")
    cyc = _q.boundaryCycle(sh, Vi)
    tris = _q.earClip(Vi, cyc)
    print("boundaryCycle=%s" % cyc)
    print("triangles=%s" % tris)
    for ti, tri in enumerate(tris):
        sht = _q.triangleShape(Vi, vd, tri)
        g = caseDBoundaryMax(fN, fD, Vi, vd, sht, [])
        print("  tri[%d]=%s -> nf=%d" % (ti, tri, g.nf))
        for k in range(g.nf):
            print("      fN=%s" % g.fN[k])
else:
    print("CONVEX: caseD directly")
    g = caseDBoundaryMax(fN, fD, Vi, vd, sh, [])
    print("  caseD nf=%d" % g.nf)

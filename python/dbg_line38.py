"""Trace conj[38] (1D line domain) construction in Python."""
import json
from cca2 import io_matlab
from cca2 import conj_q as _q
from cca2 import rat_q
from cca2.conj_q import conjQ, degenerateShape, recessionConditions, caseDBoundaryMax

raw = [json.loads(l) for l in open('/home/ylucet/AI/CCA2/python/golden/conj.jsonl') if l.strip()]
rec = raw[38]
P = io_matlab.quaPolFrom(rec['in'])
print("dom dim =", P.dom['dim'], " nv=%d ne=%d nf=%d" % (P.nv, P.ne, P.nf))
fN, fD = P.faceQ(1)
Hn = [[fN[4], fN[5]], [fN[5], fN[6]]]
Ln = [fN[7], fN[8]]
Vi, vd = rat_q.combineDen(P.VN, P.VD)
print("Vi=", Vi, " vd=", vd, " Hn=", Hn, " Ln=", Ln, " fD=", fD)
sh = degenerateShape(P, Vi)
print("sh.vs=", sh['vs'], " sh.ed=", [(e['a'], e['b'], e['d'], e['isRay']) for e in sh['ed']])
print("sh.rays=", sh['rays'], " sh.hp=", sh['hp'], " bounded=", sh['bounded'])
ok, dom, why = recessionConditions(sh, Hn, Ln, fD, Vi, vd)
print("recession ok=%s dom=%s why=%r" % (ok, dom, why))
g = caseDBoundaryMax(fN, fD, Vi, vd, sh, dom)
print("caseD nf=%d" % g.nf)
for k in range(g.nf):
    print("  fN=%s fD=%s FC=%s" % (g.fN[k], g.fD[k], g.FC[k]))
print("  EcQ=", g.EcQ)

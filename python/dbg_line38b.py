"""Trace conj[38] via the real conjPieceQ path (pieceShape + linealityConj)."""
import json
from cca2 import io_matlab
from cca2 import conj_q as _q
from cca2 import rat_q

raw = [json.loads(l) for l in open('/home/ylucet/AI/CCA2/python/golden/conj.jsonl') if l.strip()]
rec = raw[38]
P = io_matlab.quaPolFrom(rec['in'])
fN, fD = P.faceQ(1)
Hn = [[fN[4], fN[5]], [fN[5], fN[6]]]
Ln = [fN[7], fN[8]]
D = rat_q.detExact(Hn)
Vi, vd = rat_q.combineDen(P.VN, P.VD)
print("Hn=", Hn, " D=", D, " Ln=", Ln, " fD=", fD)
sh = _q.pieceShape(P, 1, Vi, vd)
print("pieceShape: vs=%s ed=%s rays=%s hp=%s bounded=%s convex=%s" % (
    sh['vs'], [(e['a'],e['b'],e['d'],e['isRay']) for e in sh['ed']], sh['rays'], sh['hp'], sh['bounded'], sh['convex']))
ok, dom, why = _q.recessionConditions(sh, Hn, Ln, fD, Vi, vd)
print("recession ok=%s dom=%s why=%r" % (ok, dom, why))
gT, isThin = _q.linealityConj(P, 1, sh, Hn, fN, fD, Vi, vd)
print("linealityConj isThin=%s" % isThin)
if isThin and gT is not None:
    print("  gT nf=%d" % gT.nf)
    for k in range(gT.nf):
        print("   fN=%s fD=%s FC=%s" % (gT.fN[k], gT.fD[k], gT.FC[k]))
    print("  EcQ=", gT.EcQ)
print("strictly convex? Hn[0][0]>0 and D>0:", Hn[0][0] > 0 and D > 0)

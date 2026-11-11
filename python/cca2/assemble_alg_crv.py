"""assemble_alg_crv -- turn envelopeSubdivisionAlg's face list into an AlgCrv (port of
assembleAlgCrv.m): tabulate the conics, the algebraic rows and the value ties, and drop faces
that repeat another's rows and value.
"""
from .alg_crv import AlgCrv
from .errors import error


def assembleAlgCrv(subdivision):
    if not isinstance(subdivision, list):
        error('assembleAlgCrv:shape',
              'subdivision must be envelopeSubdivisionAlg\'s list of faces.')
    EcQ = []
    EcA = []
    EcI = []
    FC = []
    fN = []
    fD = []
    FP = []
    FB = []
    FT = []
    Tri = []
    if subdivision and subdivision[0].get('tfs'):
        Tri = list(subdivision[0]['tfs'])
    for f in subdivision:
        rows = []
        for r in f['con']:
            c = list(r[0:6])
            if all(x == 0 for x in c[0:5]):
                continue                           # a constant is not a curve
            try:
                hit = EcQ.index(c) + 1
            except ValueError:
                EcQ.append(c)
                hit = len(EcQ)
            br = r[7] if len(r) >= 8 else 0
            rows.append([hit, r[6], br, 0])
        for a in f.get('alg') or []:
            hit = 0
            for j, e in enumerate(EcA):
                if e['FB'] == a['FB'] and e['FP'] == a['FP'] and e['num'] == a['num'] \
                        and e['den'] == a['den']:
                    hit = j + 1
                    break
            if hit == 0:
                EcA.append({'FP': a['FP'], 'FB': a['FB'], 'num': a['num'], 'den': a['den']})
                hit = len(EcA)
            rows.append([hit, a['side'], 0, 1])
        for r in f.get('ii') or []:                # a minimum over pairs: table 2
            hit = 0
            for j, e in enumerate(EcI):
                if e['FB1'] == r['FB1'] and e['FB2'] == r['FB2'] and e['FP1'] == r['FP1'] \
                        and e['FP2'] == r['FP2']:
                    hit = j + 1
                    break
            if hit == 0:
                EcI.append({'FP1': r['FP1'], 'FB1': r['FB1'], 'FP2': r['FP2'], 'FB2': r['FB2']})
                hit = len(EcI)
            rows.append([hit, 1, 0, 2])
        for t in f.get('tri') or []:
            rows.append([t[0], t[1], t[2], 3])
        FC.append(rows)
        if f['kind'] == 'triangle':
            FP.append(None)
            FB.append(0)
            fN.append([0] * 10)
            fD.append(1)
            FT.append(next(j + 1 for j, u in enumerate(Tri)
                           if u['T'] == f['tf']['T'] and u['s'] == f['tf']['s']))
        elif f['kind'] == 'bridge':
            FP.append(f['FP'])
            FB.append(f['FB'])
            fN.append([0] * 10)
            fD.append(1)
            FT.append(0)
        else:
            FP.append(None)
            FB.append(0)
            fN.append(list(f['num']))
            fD.append(f['den'])
            FT.append(0)

    # A face that repeats another's rows AND value is kept once.
    keep = [True] * len(FC)
    keys = [tuple(sorted(set(tuple(r) for r in R))) for R in FC]
    for k in range(1, len(FC)):
        for j in range(k):
            if keep[j] and keys[j] == keys[k] and FT[j] == FT[k] and FB[j] == FB[k] \
                    and FP[j] == FP[k] and fN[j] == fN[k] and fD[j] == fD[k]:
                keep[k] = False
                break

    def sel(L):
        return [x for x, kp in zip(L, keep) if kp]
    return AlgCrv([], [], EcQ, sel(fN), sel(fD), [], sel(FC), sel(FP), sel(FB), EcA, EcI, Tri,
                  sel(FT))

"""Classify every biconjQ golden record: input shape, golden outcome, Python outcome."""
import json
from cca2 import io_matlab, rat_q
from cca2.biconj_q import biconjQ, isChain, isChainSubdivision
from cca2.is_convex_q import isConvexQ

raw = [json.loads(l) for l in open('/home/ylucet/AI/CCA2/python/golden/biconj.jsonl') if l.strip()]

def Hn_of(fN):
    r = fN
    return [[r[4], r[5]], [r[5], r[6]]]

rows = []
for i, r in enumerate(raw):
    if r.get('fn') != 'biconjQ':
        continue
    inp = r['in']
    exact = inp.get('exact', True)
    nf = len(inp.get('fN', []))
    allPSD = all(rat_q.isPSD2(Hn_of(f)) for f in inp.get('fN', [])) if inp.get('fN') else '-'
    try:
        P = io_matlab.quaPolFrom(inp)
        dim = P.dom['dim']
        chain = isChain(P)
        chainsub = isChainSubdivision(P)
        conv, why = (isConvexQ(P) if dim >= 2 else (False, 'dim<2'))
    except Exception as e:
        dim, chain, chainsub, conv = 'ERR', False, False, False
    golden = r.get('out') or {}
    gkind = golden.get('kind', '-')
    gnf = golden.get('nf', '-')
    gerr = r.get('err') or ''
    # python
    try:
        g = biconjQ(P)
        pkind, pnf = type(g).__name__, getattr(g, 'nf', '?')
        # does python match the golden field-for-field?
        ok_fields = ''
        if gkind == 'QuaCon':
            ok_fields = 'fields=' + ('OK' if io_matlab.quaConMismatch(g, golden) == '' else 'DIFF')
        ok_eval = ''
        try:
            ok_eval = 'eval=' + ('OK' if io_matlab.evalMismatch(g, golden) == '' else 'DIFF')
        except Exception:
            ok_eval = 'eval=ERR'
    except Exception as e:
        pkind, pnf = type(e).__name__, str(e)[:40]
        ok_fields, ok_eval = '', ''
    rows.append((i, dim, nf, inp['nv'], inp['ne'], allPSD, chain, chainsub, conv, gkind, gnf, gerr, pkind, pnf, ok_fields, ok_eval))

print("%-5s %-4s %-3s %-3s %-3s %-5s %-5s %-5s %-5s %-6s %-4s %-12s %-16s %s" % (
    'raw', 'dim', 'nf', 'nv', 'ne', 'PSD', 'chain', 'chsub', 'conv', 'gkind', 'gnf', 'gerr', 'pyresult', 'fields/eval'))
for t in rows:
    i, dim, nf, nv, ne, psd, ch, chs, cv, gk, gn, ge, pk, pn, of, oe = t
    print("%-5s %-4s %-3s %-3s %-3s %-5s %-5s %-5s %-5s %-6s %-4s %-12s %-16s %s" % (
        i, dim, nf, nv, ne, psd, ch, chs, cv, gk, gn, ge[:12], (str(pk)[:10]+' '+str(pn))[:16], (of+' '+oe).strip() or '-'))

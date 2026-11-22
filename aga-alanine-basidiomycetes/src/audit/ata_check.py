"""Decoy follow-up: are ATA codons placed at Met columns start codons (refpos 1) or internal?
Reuses repro_star.place(); compares with ATT/ATC (Ile) placement at M columns. argv: gdir pdir"""
import sys, re, subprocess
from collections import Counter
gdir, pdir = sys.argv[1], sys.argv[2]
sys.argv = ['x', gdir, pdir, '/dev/null/none']
src = open('repro_star.py').read().split('CL = {')[0]   # definitions only (Q, REF, place, core_faa)
exec(src)
for q, label, role, cls, rels in Q:
    if q not in ('MZ439919.1', 'BK074411.1', 'PP995850.1', 'NC_054273.1', 'BK074202.1'): continue
    faa, n = core_faa(REF[cls]); ref = {}; name = None
    for l in open(faa):
        if l[0] == '>': name = l[1:].strip(); ref[name] = ''
        else: ref[name] += l.strip()
    qp = place(faa, q); rp = [place(faa, r) for r in rels if r != q]
    out = {'ATA': Counter(), 'ATH': Counter()}; starts = 0; ex = []
    for (prot, pos), (cod, aa) in qp.items():
        res = [ref[prot][pos - 1]] + [r[(prot, pos)][1] for r in rp if (prot, pos) in r]
        res = [x for x in res if x not in 'X*-']
        if not (len(res) >= 3 and len(set(res)) == 1): continue
        k = 'ATA' if cod == 'ATA' else 'ATH' if cod in ('ATT', 'ATC') else None
        if not k: continue
        out[k][res[0]] += 1
        if k == 'ATA' and res[0] == 'M':
            starts += pos <= 2; ex.append(f'{prot}:{pos}')
    print(q, label, '| ATA unanimous', dict(out['ATA']), '| M at refpos<=2:', starts, '| ATT/ATC', dict(out['ATH'].most_common(5)), '|', ' '.join(ex[:12]))

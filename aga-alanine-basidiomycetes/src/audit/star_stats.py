"""Exact counts from the reference-protein (tblastn) reproduction: per query, codon class x unanimous residue
(A / R / other), pooled Fisher tests. Reuses repro_star.py definitions (Q, REF, place, core_faa). argv: gdir pdir"""
import sys
from collections import Counter, defaultdict
from scipy.stats import fisher_exact
gdir, pdir = sys.argv[1], sys.argv[2]
sys.argv = ['x', gdir, pdir, '/dev/null/none']
exec(open('repro_star.py').read().split('CL = {')[0])
CL = {'AGA': {'AGA'}, 'CGN': {'CGT', 'CGC', 'CGA', 'CGG'}, 'GCN': {'GCT', 'GCC', 'GCA', 'GCG'}}
pool = defaultdict(Counter)
print('query\tlabel\trole\t' + '\t'.join(f'{c}_unan\t{c}_A\t{c}_R' for c in CL))
for q, label, role, cls, rels in Q:
    faa, n = core_faa(REF[cls]); ref = {}; name = None
    for l in open(faa):
        if l[0] == '>': name = l[1:].strip(); ref[name] = ''
        else: ref[name] += l.strip()
    qp = place(faa, q); rp = [place(faa, r) for r in rels]
    t = defaultdict(Counter)
    for (prot, pos), (cod, aa) in qp.items():
        res = [ref[prot][pos - 1]] + [r[(prot, pos)][1] for r in rp if (prot, pos) in r]
        res = [x for x in res if x not in 'X*-']
        if len(res) >= 3 and len(set(res)) == 1:
            for c, cs in CL.items():
                if cod in cs: t[c][res[0]] += 1; t[c]['n'] += 1
    for c in CL: pool[(role, c)].update(t[c])
    print(f'{q}\t{label}\t{role}\t' + '\t'.join(f"{t[c]['n']}\t{t[c]['A']}\t{t[c]['R']}" for c in CL))
for role in ('in 32', 'control', 'extra'):
    a, g, cg = pool[(role, 'AGA')], pool[(role, 'GCN')], pool[(role, 'CGN')]
    print(f'POOLED {role}: AGA n{a["n"]} A{a["A"]} R{a["R"]} | CGN n{cg["n"]} A{cg["A"]} R{cg["R"]} | GCN n{g["n"]} A{g["A"]} R{g["R"]}')
    if a['n']:
        print('  AGA vs CGN at A vs R columns, one-sided p =', f"{fisher_exact([[a['A'], a['R']], [cg['A'], cg['R']]], alternative='greater')[1]:.2g}")
        print('  Ala share AGA vs GCN:', f"{a['A'] / a['n']:.2f} vs {g['A'] / g['n']:.2f}, two-sided p =",
              f"{fisher_exact([[a['A'], a['n'] - a['A']], [g['A'], g['n'] - g['A']]])[1]:.2g}")
        print('  AGA non-A/R residues:', dict((k, v) for k, v in a.items() if k not in 'nAR'))

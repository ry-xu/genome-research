import numpy as np, csv
from collections import Counter, defaultdict
from scipy.stats import fisher_exact
S = 'scratch/'
C = np.load(S + 'mycache.npz')
AA = 'ACDEFGHIKLMNPQRSTVWY'; B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
sg, sc, col, strict, cons = C['sg'], C['sc'], C['col'], C['strict'], C['cons']
acc, tax = C['acc'].tolist(), C['tax'].tolist(); K = int(C['K']); models = C['models'].tolist()
G = len(acc); ai = {a: i for i, a in enumerate(acc)}
ala = np.load(S + 'ala32.npy'); arg = np.load(S + 'arg.npy')
st = strict[col]; cc = cons[col]
A, R = AA.index('A'), AA.index('R'); AGA = COD.index('AGA'); AGG = COD.index('AGG')
SYN = [COD.index(c) for c in ('CGT','CGC','CGA','CGG','AGG')]
# --- tRNA: re-derive G3:U70 and A73 from the stems/sequence
comp = {('G','T'),('T','G'),('G','C'),('C','G'),('A','T'),('T','A')}
tr = defaultdict(list); fungal_acc = set(acc[g] for g in range(G) if 'Fungi' in tax[g].split(';'))
n_rows = 0; mism = 0
for r in csv.DictReader(open('lab/work/new-stops-euk/results/aga_ala/fungi.aragorn.tsv'), delimiter='\t'):
    n_rows += 1
    if r['anticodon'] != 'TCT': continue
    s5, s3, seq = r['stem5'], r['stem3'], r['seq']
    p = f'{s5[2]}:{s3[4]}'
    if p != r['pair3_70']: mism += 1
    # check stems are at seq start and just before N73
    i3 = seq.rfind(s3)
    n73 = seq[i3 + 7] if i3 >= 0 and i3 + 7 < len(seq) else '?'
    tr[r['acc']].append((p, n73, seq.startswith(s5), r['N73']))
print('aragorn rows', n_rows, 'genomes with any row', len(set(r for r in tr)), 'pair mismatches', mism)
accs_in_aragorn = set(r['acc'] for r in csv.DictReader(open('lab/work/new-stops-euk/results/aga_ala/fungi.aragorn.tsv'), delimiter='\t'))
print('fungal genomes with any aragorn output', len(accs_in_aragorn & fungal_acc), 'of', len(fungal_acc))
def g3u70(a): return [t for t in tr.get(a, []) if t[0] == 'G:T']
c28 = [acc[g] for g in ala if g3u70(acc[g])]
print('32 with G3:U70 UCU:', len(c28), '; all A73 (own seq parse)?', all(any(t[1] == 'A' for t in g3u70(a)) for a in c28),
      '; without:', [acc[g] for g in ala if not g3u70(acc[g])])
argacc = [acc[g] for g in arg]
withucu = [a for a in argacc if tr.get(a)]
print('arg class', len(argacc), 'with UCU tRNA', len(withucu), 'with G3:U70', sum(1 for a in withucu if g3u70(a)),
      'no aragorn output', sum(1 for a in argacc if a not in accs_in_aragorn))
# all fungal genomes carrying a G3:U70 UCU tRNA
carr = [a for a in fungal_acc if g3u70(a)]
print('all fungal G3:U70 UCU carriers:', len(carr), sorted((tax[ai[a]].split(';')[-1], a) for a in carr if ai[a] not in set(ala.tolist())))
# --- AGG at strict columns among 32 and among the three classes
CL = ('Tremellomycetes', 'Exobasidiomycetes', 'Cystobasidiomycetes', 'Microbotryomycetes')
cls = [g for g in range(G) if any(k in tax[g].split(';') for k in CL)]
x = np.isin(sg, cls) & (sc == AGG) & st
print('AGG at strict cols in classes:', [(acc[s], tax[s].split(';')[-1], AA[c], models[k // K]) for s, c, k in zip(sg[x], cc[x], col[x])])
# --- selection-free per-order distinct-column contrast: all genomes in each order, regardless of AGA use
for l in ('Trichosporonales', 'Tremellales', 'Tilletiales', 'Cystobasidiales', 'Sporidiobolales'):
    gs = [g for g in range(G) if l in tax[g].split(';')]
    mm = np.isin(sg, gs) & st
    def dc(cm, r): return len(set(col[mm & cm & (cc == r)].tolist()))
    t = [dc(sc == AGA, A), dc(sc == AGA, R), dc(np.isin(sc, SYN), A), dc(np.isin(sc, SYN), R)]
    s = [int((mm & (sc == AGA) & (cc == A)).sum()), int((mm & (sc == AGA) & (cc == R)).sum()), int((mm & np.isin(sc, SYN) & (cc == A)).sum()), int((mm & np.isin(sc, SYN) & (cc == R)).sum())]
    agaR = Counter(acc[g] for g in sg[mm & (sc == AGA) & (cc == R)])
    print(f'{l}: all {len(gs)} genomes; distinct {t} p={fisher_exact([[t[0],t[1]],[t[2],t[3]]], "greater")[1]:.2g}; sites {s}; AGA@R from {dict(agaR)}')
# --- BK074393.1's AGA at Arg: which column, and what codons do the other Tilletia have there
g = ai['BK074393.1']
x = (sg == g) & (sc == AGA) & st & (cc == R)
for k in set(col[x].tolist()):
    print('BK074393.1 AGA@R column', models[k // K], k % K)
    y = (col == k) & np.isin(sg, [ai[a] for a in ('NC_009880.1','NC_010651.1','NC_073546.1','BK074390.1','BK074391.1','BK074392.1','BK074121.1')])
    print('   other Tilletiales codons there:', [(acc[s], COD[c]) for s, c in zip(sg[y], sc[y])])
# --- shared Ala columns across lineages (9 shared?)
LIN = ('Trichosporonales', 'Tremellales', 'Tilletiales', 'Cystobasidiales', 'Sporidiobolales')
per = {l: set(col[np.isin(sg, [g for g in ala if l in tax[g].split(';')]) & (sc == AGA) & st & (cc == A)].tolist()) for l in LIN}
sub = {'Agaricomycotina': per['Trichosporonales'] | per['Tremellales'], 'Ustilaginomycotina': per['Tilletiales'], 'Pucciniomycotina': per['Cystobasidiales'] | per['Sporidiobolales']}
allc = set().union(*sub.values())
shared = [k for k in allc if sum(k in v for v in sub.values()) >= 2]
print('sizes by subphylum', {k: len(v) for k, v in sub.items()}, 'columns in >=2 subphyla', len(shared))

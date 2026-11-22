"""Independent re-derivation of the AGA=Ala census numbers (own loader; same strict definition as the post:
one residue >= 90% of family-weighted weight, residues read with each genome's declared table, sense-to-stop erased)."""
import sys, numpy as np, csv
from collections import Counter, defaultdict
from Bio.Data import CodonTable
from scipy.stats import fisher_exact, hypergeom
D = 'lab/work/new-stops-euk/data/codes/'
AA = 'ACDEFGHIKLMNPQRSTVWY'; B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
STD = CodonTable.unambiguous_dna_by_id[1].forward_table
def tab(t):
    try: f = CodonTable.unambiguous_dna_by_id[t].forward_table
    except KeyError: f = STD
    return np.array([AA.index(f[c]) if c in f else (AA.index(STD[c]) if c in STD else -1) for c in COD], np.int8)
acc, tax, table, SG, SC, SM, SK = [], [], [], [], [], [], []
for p in ('mito', 'gbmito', 'gbmitoorf'):
    d = np.load(D + p + '.census.npz')
    off = len(acc)
    acc += d['acc'].tolist(); tax += d['taxonomy'].tolist(); table += [int(t) for t in d['table']]
    SG.append((d['site_genome'] + off).astype(np.int32)); SC.append((d['site_codon'].astype(np.int64) % 64).astype(np.int8))
    SM.append(d['site_model'].astype(np.int32)); SK.append(d['site_col'].astype(np.int32))
    models = d['models'].tolist()
    print(p, len(d['acc']), 'genomes', len(d['site_genome']), 'sites', file=sys.stderr)
sg = np.concatenate(SG); sc = np.concatenate(SC); sm = np.concatenate(SM); sk = np.concatenate(SK); del SG, SC, SM, SK
K = int(sk.max()) + 1; M = len(models)
col = sm * K + sk; del sm, sk
G = len(acc)
fam = Counter(t.rsplit(';', 1)[0] for t in tax)
w = np.array([1 / fam[t.rsplit(';', 1)[0]] for t in tax])
gt = np.array(table)
res = np.full(len(sg), -1, np.int8)
for t in set(table):
    gm = np.flatnonzero(gt == t)
    m = np.isin(sg, gm)
    res[m] = tab(t)[sc[m]]
ok = res >= 0
freq = np.zeros((M * K, 20))
np.add.at(freq, (col[ok], res[ok]), w[sg[ok]])
tot = freq.sum(1)
strict = (tot > 0) & (freq.max(1) >= 0.9 * np.maximum(tot, 1e-12))
cons = freq.argmax(1).astype(np.int8)
print('genomes', G, 'unique acc', len(set(acc)), 'strict cols', strict.sum(), file=sys.stderr)
np.savez('scratch/mycache.npz',
         sg=sg, sc=sc, col=col, strict=strict, cons=cons, acc=np.array(acc), tax=np.array(tax), table=gt, K=K, models=np.array(models), freq=freq)

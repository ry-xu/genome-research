"""Audit re-count of the census numbers in the post, with own code (not sitescan_org.load), streaming one census file
at a time to keep memory low.

Definition re-implemented from the README/sitescan_org docstring: column = (model, column); each genome's residue read
with its declared table (declared stops that are standard sense -> standard residue; stops not counted); each genome
weighted 1/(genomes sharing its taxonomy minus the last rank); strict = top residue >= 0.9 of the column weight.

usage: python -I census_check.py mito.census.npz gbmito.census.npz gbmitoorf.census.npz
"""
import sys
from collections import Counter, defaultdict
import numpy as np
from Bio.Data import CodonTable

AA = 'ACDEFGHIKLMNPQRSTVWY'
CODONS = [a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG']
STD = CodonTable.unambiguous_dna_by_id[1].forward_table
ARGC = ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT']


def resvec(t):
    try:
        f = CodonTable.unambiguous_dna_by_id[int(t)].forward_table
    except KeyError:
        f = STD
    v = []
    for c in CODONS:
        a = f.get(c) or STD.get(c)
        v.append(AA.index(a) if a else -1)
    return np.array(v, dtype=np.int8)


paths = sys.argv[1:]
meta = []  # (acc, tax, table)
for p in paths:
    d = np.load(p)
    meta += list(zip(d['acc'].tolist(), d['taxonomy'].tolist(), [int(x) for x in d['table']]))
    K = d['emax'].shape[1]
fam = Counter(t.rsplit(';', 1)[0] for _, t, _ in meta)
W = np.array([1 / fam[t.rsplit(';', 1)[0]] for _, t, _ in meta])
TAB = np.array([t for _, _, t in meta])
RV = {t: resvec(t) for t in set(TAB.tolist())}
RED = np.array(['Rhodophyta' in t for _, t, _ in meta])
ACC = [a for a, _, _ in meta]
print('genomes', len(meta), 'red algal', RED.sum(), file=sys.stderr)

freq = np.zeros((120 * (K + 1), 20))
red_sites = []  # (genome, codon, col)
off = 0
for p in paths:
    d = np.load(p)
    sg = d['site_genome'].astype(np.int64) + off
    sc = d['site_codon'].astype(np.int64) % 64
    col = d['site_model'].astype(np.int64) * (K + 1) + d['site_col'].astype(np.int64)
    res = np.full(len(sg), -1, dtype=np.int16)
    tg = TAB[sg]
    for t, v in RV.items():
        m = tg == t
        res[m] = v[sc[m]]
    ok = res >= 0
    np.add.at(freq, (col[ok], res[ok].astype(np.int64)), W[sg[ok]])
    r = RED[sg]
    red_sites.append(np.stack([sg[r], sc[r], col[r]], 1))
    off += len(d['acc'])
    del sg, sc, col, res, tg, ok, r, d
tot = freq.sum(1)
strict = (tot > 0) & (freq.max(1) >= 0.9 * np.maximum(tot, 1e-12))
cons = freq.argmax(1)
R = np.concatenate(red_sites)
print(f'red algal genomes with sites: {len(set(R[:, 0].tolist()))}; red sites {len(R)}')

M, Ri = AA.index('M'), AA.index('R')
st = strict[R[:, 2]]
cc = cons[R[:, 2]]
cod = np.array(CODONS)[R[:, 1]]
g_acc = np.array(ACC)[R[:, 0]]
red_tables = Counter(TAB[RED].tolist())
print('declared tables of red algal genomes:', dict(red_tables))
for c in ('AGG', 'AGA', 'ATG'):
    m = cod == c
    gm = Counter(g_acc[m & st & (cc == M)].tolist())
    print(f'{c}: red sites {m.sum()}, at strict Met {np.sum(m & st & (cc == M))} '
          f'(in {len(gm)} genomes; {dict(gm.most_common(3))}...), at strict Arg {np.sum(m & st & (cc == Ri))}')

for acc in ('NC_062302.1', 'MT483996.1', 'MT483995.1', 'NC_062401.1', 'NC_062400.1', 'NC_062303.1', 'NC_062397.1'):
    if acc not in ACC:
        print(acc, 'not in census'); continue
    m = g_acc == acc
    print(f'\n== {acc}: {m.sum()} sites')
    for res_name, ri in (('Met', M), ('Arg', Ri)):
        k = m & st & (cc == ri)
        print(f'  strict {res_name} columns by codon: {dict(Counter(cod[k].tolist()).most_common())}')
    k = m & (cod == 'AGG')
    print(f'  AGG sites {k.sum()}; consensus at AGG sites (all): {dict(Counter(AA[x] for x in cc[k]).most_common())}')
    print(f'  AGG at strict columns by consensus: {dict(Counter(AA[x] for x in cc[k & st]).most_common())}')
    k = m & (cod == 'AGA')
    print(f'  AGA at strict columns by consensus: {dict(Counter(AA[x] for x in cc[k & st]).most_common())}')
    # per-genome Arg codon usage at strict Arg
    print(f'  Arg codons at strict Arg: {dict(Counter(cod[m & st & (cc == Ri) & np.isin(cod, ARGC)].tolist()))}')

# per-genome table for subset checks
with open(sys.argv[0].rsplit('/', 2)[0] + '/red_per_genome.tsv', 'w') as fh:
    fh.write('acc\tsource\ttable\ttaxonomy\tsites\tAGG_sites\tAGG_strictMet\tAGG_strictArg\n')
    src = []
    for i, p in enumerate(paths):
        src += [p.split('/')[-1]] * len(np.load(p)['acc'])
    for g in sorted(set(R[:, 0].tolist())):
        m = R[:, 0] == g
        a = m & (cod == 'AGG')
        fh.write(f"{ACC[g]}\t{src[g]}\t{TAB[g]}\t{meta[g][1]}\t{m.sum()}\t{a.sum()}\t{np.sum(a & st & (cc == M))}\t{np.sum(a & st & (cc == Ri))}\n")

# five relatives' codons at SAG's AGG columns
rel = ['MT483996.1', 'NC_062401.1', 'NC_062400.1', 'NC_062303.1', 'NC_062397.1']
sag = (g_acc == 'NC_062302.1') & (cod == 'AGG')
cols = set(R[sag, 2].tolist())
n = atg = arg = 0
per = defaultdict(Counter)
for a in rel:
    m = (g_acc == a) & np.isin(R[:, 2], list(cols))
    cs = cod[m]
    # a genome can have >1 codon per column only if duplicated genes; count distinct (col) entries
    n += len(cs); atg += np.sum(cs == 'ATG'); arg += np.sum(np.isin(cs, ARGC))
    per[a].update(cs.tolist())
print(f'\nSAG AGG columns: {len(cols)}; relatives: {n} codons at them, ATG {atg}, Arg codons {arg}')
for a in rel:
    print(' ', a, sum(per[a].values()), 'ATG', per[a]['ATG'], 'Arg', sum(per[a][c] for c in ARGC))

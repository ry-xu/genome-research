"""Check the acceptor-stem parse in fungi.aragorn.tsv.
(1) sensitivity: share of tRNA-Ala(TGC) with G:T at 3:70 (should be ~all).
(2) independent parse: for each TCT tRNA, all (o,t) readings with >=6 of 7 acceptor pairs + U8 or T-loop TTC; does any
    reading give G3:U70 when the table says otherwise (missed carriers), or vice versa?
(3) anticodon loop: 7-nt loop around the anticodon and anticodon-stem pairing, per G3:U70 UCU tRNA."""
import csv, numpy as np
from collections import Counter, defaultdict
A = 'lab/work/new-stops-euk/results/aga_ala/'
PAIRS = {('G', 'C'), ('C', 'G'), ('A', 'T'), ('T', 'A'), ('G', 'T'), ('T', 'G')}
C = np.load(A + 'audit/cache.npz', allow_pickle=True)
ala = set(C['acc'][np.load(A + 'audit/ala_class_idx.npy')]); arg = set(C['acc'][np.load(A + 'audit/arg_class_idx.npy')])
ala_t = Counter(); ucu = []
for r in csv.DictReader(open(A + 'fungi.aragorn.tsv'), delimiter='\t'):
    if r['aa'] == 'Ala':
        ala_t[r['pair3_70']] += 1
    if r['anticodon'] == 'TCT':
        ucu.append(r)
print('tRNA-Ala 3:70 pairs:', ala_t.most_common(6))
def readings(s):
    out = []
    for o in range(0, 4):
        for t in range(0, 6):
            e = len(s) - t
            if e - o < 55: continue
            st5, st3 = s[o:o + 7], s[e - 8:e - 1]
            wc = sum((st5[i], st3[6 - i]) in PAIRS for i in range(7))
            if wc >= 6:
                out.append((o, t, st5[2] + ':' + st3[4], s[e - 1], wc, s[o + 7] == 'T', s[e - 20:e - 17] == 'TTC'))
    return out
miss, cls = Counter(), Counter()
for r in ucu:
    rd = readings(r['seq'])
    alt_gt = any(x[2] == 'G:T' for x in rd)
    tab_gt = r['pair3_70'] == 'G:T'
    k = 'ala' if r['acc'] in ala else 'arg' if r['acc'] in arg else 'other'
    cls[(k, tab_gt, alt_gt)] += 1
    if alt_gt and not tab_gt:
        miss[k] += 1
        if miss[k] <= 5:
            print('alt G:T reading not chosen:', k, r['acc'], r['seq'], [x for x in rd if x[2] == 'G:T'], r['stem5'], r['stem3'])
print('(class, table G:T, any >=6-pair reading G:T):', dict(cls))
# anticodon-loop check for the G3:U70 UCU tRNAs
print('\nanticodon arm of G3:U70 TCT tRNAs (5stem loop 3stem, pairs):')
seen = set()
for r in ucu:
    if r['pair3_70'] != 'G:T': continue
    s = r['seq']; i = s.find('TCT')
    hits = []
    for j in range(len(s)):
        if s[j:j + 3] != 'TCT' or j < 25 or j > len(s) - 30: continue
        loop = s[j - 2:j + 5]; st5 = s[j - 7:j - 2]; st3 = s[j + 5:j + 10]
        p = sum((st5[q], st3[4 - q]) in PAIRS for q in range(5))
        hits.append((p, st5, loop, st3))
    best = max(hits) if hits else None
    key = (r['seq'])
    if key in seen: continue
    seen.add(key)
    print(r['acc'], best, 'U33' if best and best[2][1] == 'T' else '', 'purine37' if best and best[2][5] in 'AG' else '')

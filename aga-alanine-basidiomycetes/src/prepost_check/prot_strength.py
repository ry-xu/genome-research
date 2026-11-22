import csv, sys
from collections import Counter
RES = sys.argv[1]
thr = {'A14 (scan 3271)': 1.218192, 'A13 (scan 3307)': 0.9537924}
dec_scores, tgt_scores, a_dec = [], [], []
for r in csv.DictReader(open(RES), delimiter='\t'):
    s = float(r['sage_discriminant_score'])
    if r['label'] == '-1':
        dec_scores.append(s)
        m = [p for p in r['proteins'].split(';') if p.replace('rev_', '').startswith('mito|')]
        if m and {p.rsplit('|', 1)[1] for p in m} == {'A'}: a_dec.append((s, r['peptide'], r['spectrum_q']))
    else:
        tgt_scores.append(s)
nd = len(dec_scores)
for k, t in thr.items():
    d = sum(s >= t for s in dec_scores); tg = sum(s >= t for s in tgt_scores)
    print(f'{k}: score {t}; decoys >= {d}, targets >= {tg}; expected false A-only at this score ~ {d * 45 / nd:.4f}')
print('best A-only decoy PSMs:', sorted(a_dec, reverse=True)[:3])
